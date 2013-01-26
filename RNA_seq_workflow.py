import dxpy
import subprocess
import logging
import os
import multiprocessing

from dxpy.dxlog import DXLogHandler

def run_shell(command):
    print "Running "+command
    subprocess.check_call(command, shell=True)

def make_indexed_reference( ref_ID ):

    run_shell("dx-contigset-to-fasta %s reference.fasta" % ref_ID)
    ref_details = dxpy.DXRecord(ref_ID).get_details()
    ref_name = dxpy.DXRecord(ref_ID).describe()['name']

    # call bowtie2-build
    run_shell("bowtie2-build reference.fasta indexed_ref")
    # package it into an archive for uploading
    run_shell("XZ_OPT=-0 tar -cJf reference.tar.xz indexed_ref*")

    indexed_ref_dxfile = dxpy.upload_local_file("reference.tar.xz", hidden=True, wait_on_close=True)
    
    indexed_ref_record = dxpy.new_dxrecord(name=ref_name + " (indexed for Bowtie2)",
                                           types=["BowtieLetterContigSetV2"],
                                           details={'index_archive': dxpy.dxlink(indexed_ref_dxfile.get_id()),
                                                    'original_contigset': dxpy.dxlink(ref_ID)})
    indexed_ref_record.close()
    
    '''
    # TODO: dxpy project workspace convenience functions
    if "projectWorkspace" in job:
        indexed_ref_record.clone(job["projectWorkspace"])
    '''

    return indexed_ref_record.get_id()

def upload_transcripts_file( trans_file ):
    with open(trans_file, 'r') as fh:
        # eat column header line
        line = fh.readline().rstrip('\n')

        line = line.split('\t')

        trans_schema = [("chr", "string"),
                        ("lo", "int32"),
                        ("hi", "int32"),
                        ("tracking_id", "string"),
                        ("class_code", "string"),
                        ("nearest_ref_id", "string"),
                        ("gene_id", "string"),
                        ("gene_short_name", "string"),
                        ("tss_id", "string"),
                        ("length", "int32"),
                        ("coverage", "float"),
                        ("FPKM", "float"),
                        ("FPKM_lo", "float"),
                        ("FPKM_hi", "float"),
                        ("status", "string")]

        column_descriptors = [dxpy.DXGTable.make_column_desc(name, type) for name, type in trans_schema]

        gri_index = dxpy.DXGTable.genomic_range_index("chr", "lo", "hi")
        transcripts = dxpy.new_dxgtable(column_descriptors, indices=[gri_index])
        transcripts.rename("FPKM_per_gene")

        while True:
            line = fh.readline()
            line = line.rstrip('\n')
            if line == '':
                break

            line = line.split('\t')

            try:
                chrom = line[6].split(":")[0]
                lo = int(line[6].split(":")[1].split("-")[0]) - 1
                hi = int(line[6].split(":")[1].split("-")[1])
                # no length set, set to 0
                if line[7] == '-':
                    line[7] = 0
                if line[8] == '-':
                    line[8] = -1

                trans_row = [chrom, lo, hi, 
                             line[0], 
                             line[1], 
                             line[2], 
                             line[3], 
                             line[4], 
                             line[5], 
                             int(line[7]), 
                             float(line[8]), 
                             float(line[9]), 
                             float(line[10]), 
                             float(line[11]), 
                             line[12]]

                transcripts.add_row(trans_row)
            except IndexError:
                raise dxpy.AppError("Error parsing transcript file from cufflinks.  Line: "+line)

    transcripts.close(block = True)

    return transcripts

def check_reads( reads_tables ):
    # validate that tables contain data that can be used together (all paired or all unpaired, etc)

    if len(reads_tables) == 0:
        raise AppError("Please enter at least one Reads table as input")

    single = 0
    paired = 0

    for table in reads_tables:
        if 'sequence2' in dxpy.DXGTable(table).get_col_names():
            paired = paired + 1
        else:
            single = single + 1

    if single > 0 and paired > 0:
        raise AppError("Found both single and paired-end reads.  Please only input one type.")

    return

def dump_fastqa( reads_ID, output_base ):
    if 'sequence2' in dxpy.DXGTable(reads_ID).get_col_names():
        paired = True
    else:
        paired = False

    if paired:
        run_shell(" ".join(["dx-reads-to-fastq", reads_ID, "--output "+output_base+"_1", "--output2 "+output_base+"_2"]))
    else:
        run_shell(" ".join(["dx-reads-to-fastq", reads_ID, "--output "+output_base+"_1"]))

    run_shell("head "+output_base+"_1")

    if paired:
        return output_base+"_1", output_base+"_2"
    else:
        return output_base+"_1", None

@dxpy.entry_point('main')
def main(**job_inputs):
    print "Beginning processing of RNA data"

    output = {}

    check_reads( job_inputs['reads'] )

    # Convert reads tables to FASTQ/FASTA files
    left_reads = []
    right_reads = []

    current_reads = 0
    for reads in job_inputs['reads']:
        print "Converting reads table "+str(reads['$dnanexus_link'])
        left, right = dump_fastqa( reads['$dnanexus_link'], "reads_"+str(current_reads) )

        left_reads.append( left )
        if right != None:
            right_reads.append( right )

        current_reads += 1
    
    # Convert Genes Object to GFF file 

    run_shell("dx-genes-to-gtf --output genes.gtf "+job_inputs['gene_model']['$dnanexus_link'])

    # Create or download indexed genome
    genome = dxpy.DXRecord(job_inputs['reference'])

    if not 'indexed_reference' in job_inputs:
        output['indexed_reference'] = dxpy.dxlink(make_indexed_reference(genome.get_id()))
    else:
        output['indexed_reference'] = job_inputs['indexed_reference']
        indexed_genome = dxpy.DXRecord(job_inputs['indexed_reference'])
        dxpy.download_dxfile(indexed_genome.get_details()['index_archive'], "reference.tar.xz")
        run_shell("tar -xJf reference.tar.xz")

    # call tophat
    num_cpus = multiprocessing.cpu_count()

    cmd = " ".join(['tophat', "-p", str(num_cpus), job_inputs['tophat_options'], "-G genes.gtf", "--transcriptome-index=./genes", "-T", "indexed_ref", " ", ",".join(left_reads)])

    if len(right_reads) != 0:
        cmd += " " + ",".join(right_reads)

    # Invoke tophat2 with FASTQ/A file(s) and indexed reference    
    try:
        run_shell(cmd)
    except:
        raise dxpy.AppError("Error while running Tophat.  This could be caused by an incompatible gene model and reference or incorrect optional parameters.  Please check that these are all correct")

    # upload and import the BAM as a Mappings table
    accepted_hits_file = dxpy.upload_local_file('tophat_out/accepted_hits.bam', wait_on_close=True)
    name = job_inputs.get('output_name', "RNA-seq mappings")
    sam_importer = dxpy.DXApp(name="sam_importer")
    print "Importing BAM output of Tophat"
    import_job = sam_importer.run(app_input={"file":dxpy.dxlink(accepted_hits_file.get_id()), 
                                             "reference_genome":dxpy.dxlink(genome.get_id()),
                                             "name":name})

    cuff_cmd = " ".join(['cufflinks', '-p', str(num_cpus), '-G genes.gtf', '-o cuff'])

    if 'cufflinks_options' in job_inputs:
        cuff_cmd += " "+job_inputs['cufflinks_options']

    cuff_cmd += " tophat_out/accepted_hits.bam"

    # now with mapped reads in hand we can run cufflinks
    try:
        run_shell(cuff_cmd)
    except:
        raise dxpy.AppError("Error while running Cufflinks.  Please check that your parameters are valid")

    print "Packing, uploading, and parsing cufflinks output"
    # package cufflinks output
    run_shell("tar -czf cufflinks_output.tar.gz cuff/")
    orig_trans_file = dxpy.upload_local_file('cufflinks_output.tar.gz')
    transcripts_table = upload_transcripts_file('cuff/genes.fpkm_tracking')

    output['mappings'] = {"job":import_job.get_id(), "field": "mappings"}
    output['transcripts'] = dxpy.dxlink(transcripts_table.get_id())
    output['cufflinks_output'] = dxpy.dxlink(orig_trans_file.get_id())

    print "DONE!"

    return output

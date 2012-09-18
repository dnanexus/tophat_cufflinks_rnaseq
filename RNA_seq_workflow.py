
import dxpy
import subprocess
import logging
import os
import multiprocessing

from dxpy.dxlog import DXLogHandler

def run_shell(command):
    logging.debug("Running "+command)
    subprocess.check_call(command, shell=True)

def make_indexed_reference( ref_ID ):
    
    run_shell("contigset2fasta %s reference.fasta" % ref_ID)
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


def parse_tophat_options( options ):
    # take job input and translate into command line options for tophat
    # also do some sanity checking?

    opt_string = []

    if not 'tophap_options' in options:
        seg_len = 100
        multi_hit = 1 

        opt_list = " ".join(["-g", str(multi_hit), "--b2-very-fast", "--no-coverage-search", "--segment-length", str(seg_len)])

    else:
        opt_list = options['tophat_options']

    return opt_list

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
                        ("q0_FPKM", "float"),
                        ("q0_FPKM_lo", "float"),
                        ("q0_FPKM_hi", "float"),
                        ("q0_status", "string")]

        column_descriptors = [dxpy.DXGTable.make_column_desc(name, type) for name, type in trans_schema]

        gri_index = dxpy.DXGTable.genomic_range_index("chr", "lo", "hi")
        transcripts = dxpy.new_dxgtable(column_descriptors, indices=[gri_index])

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
    pass

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
    
    logging.getLogger().setLevel(logging.DEBUG)

    logging.debug("Beginning processing of RNA data")

    output = {}

    options = parse_tophat_options(job_inputs)

    check_reads( job_inputs['reads'] )

    resources_ID = os.environ.get("DX_RESOURCES_ID")

    # Convert reads tables to FASTQ/FASTA files
    left_reads = []
    right_reads = []

    current_reads = 0
    for reads in job_inputs['reads']:
        logging.debug("Converting reads table "+str(reads['$dnanexus_link']))
        left, right = dump_fastqa( reads['$dnanexus_link'], "reads_"+str(current_reads) )

        left_reads.append( left )
        if right != None:
            right_reads.append( right )

        current_reads += 1
    
    # hard code hg19 and genes tracks into analysis

    resources_id = os.environ['DX_RESOURCES_ID']
    resource_bundle_id = dxpy.find_one_data_object(classname="file", name="tophat_resources.tar.gz", project=resources_id, return_handler = False)['id']
    genome_id = dxpy.find_one_data_object(classname="record", name="hg19", project=resources_id, return_handler = False)['id']

    #resource_bundle_id = job_inputs['resources']

    logging.debug("Downloading hg19 and transcript information")
    dxpy.download_dxfile(dxpy.dxlink(resource_bundle_id, project_id=resources_id), "tophat_resources.tar.gz")

    logging.debug("Unpacking resource bundle")
    run_shell("tar -xzf tophat_resources.tar.gz")


    # if we're taking in a reference from the user then 
    '''
    # download reference (and make index if necessary)
    if "ContigSet" in dxpy.DXRecord(job_inputs['reference']).describe()['types']:
        output['indexed_reference'] = dxpy.dxlink(make_indexed_reference(job_inputs['reference']['$dnanexus_link']))

    elif "BowtieLetterContigSetV2" in dxpy.DXRecord(job_inputs['reference']).describe()['types']:
        dxpy.download_dxfile(dxpy.get_details(job_inputs["reference"])['index_archive'], "indexed_ref.tar.xz")
        run_shell("tar -xJf indexed_ref.tar.xz")
        output['indexed_reference'] = job_inputs['reference']

    run_shell("ls -l")

    '''

    num_cpus = multiprocessing.cpu_count()

    tophat_options = parse_tophat_options( job_inputs )

    cmd = " ".join(['tophat', "-p", str(num_cpus), tophat_options, "--transcriptome-index=./genes", "--no-novel-juncs", "-T", "genome", " ", ",".join(left_reads)])

    if len(right_reads) != 0:
        cmd += " " + ",".join(right_reads)

    # Invoke tophat2 with FASTQ/A file(s) and indexed reference    
    run_shell(cmd)

    #ref = dxpy.DXRecord(output['indexed_reference']).get_details()['original_contigset']['$dnanexus_link']
    #ref_proj = dxpy.DXRecord(ref).describe()['project']

    # upload and import the BAM as a Mappings table
    accepted_hits_file = dxpy.upload_local_file('tophat_out/accepted_hits.bam', wait_on_close=True)
    name = job_inputs.get('output name', "RNA-seq mappings")
    sam_importer = dxpy.DXApp(name="sam_bam_importer")
    logging.debug("Importing BAM output of Tophat")
    import_job = sam_importer.run(app_input={"file":dxpy.dxlink(accepted_hits_file.get_id()), 
                                             "reference_genome":dxpy.dxlink(genome_id, project_id=resources_id),
                                             "name":name})

    cuff_cmd = " ".join(['cufflinks', '-p', str(num_cpus), '-G genes.gff', '-o cuff', 'tophat_out/accepted_hits.bam'])    

    # now with mapped reads in hand we can run cufflinks
    run_shell(cuff_cmd)

    logging.debug("Packing, uploading, and parsing cufflinks output")
    # package cufflinks output
    run_shell("tar -czf cufflinks_output.tar.gz cuff/")
    orig_trans_file = dxpy.upload_local_file('cufflinks_output.tar.gz')
    transcripts_table = upload_transcripts_file('cuff/genes.fpkm_tracking')

    output['mappings'] = {"job":import_job.get_id(), "field": "mappings"}
    output['transcripts'] = dxpy.dxlink(transcripts_table.get_id())
    output['cufflinks_output'] = dxpy.dxlink(orig_trans_file.get_id())

    logging.debug("DONE!")

    return output

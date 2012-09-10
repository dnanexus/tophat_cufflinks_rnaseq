
import dxpy
import subprocess
import logging


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
    
    indexed_ref_record = dxpy.new_dxrecord(name=ref_name + " (indexed for Bowtie)",
                                           types=["BowtieLetterContigSetV2"],
                                           details={'index_archive': dxpy.dxlink(indexed_ref_dxfile.get_id()),
                                                    'original_contigset': job['input']['reference']})
    indexed_ref_record.close()
    
    '''
    # TODO: dxpy project workspace convenience functions
    if "projectWorkspace" in job:
        indexed_ref_record.clone(job["projectWorkspace"])
    '''

    return indexed_ref_record


def parse_tophat_options( ):
    # take job input and translate into command line options for tophat
    # also do some sanity checking?
    pass

def check_reads( reads_tables ):
    # validate that tables contain data that can be used together (all paired or all unpaired, etc)
    pass

def dump_fastqa( reads_ID, output_base ):

    if 'sequence2' in dxpy.DXGTable("reads_ID").get_col_names():
        paired = True
    else:
        paired = False

    if paired:
        run_shell(" ".join(["dx-reads-to-fastq", "--output "+output_base+"_1", "--output2 "+output_base+"_2"]))
    else:
        run_shell(" ".join(["dx-reads-to-fastq", "--output "+output_base+"_1"]))

@dxpy.entry_point
def main(**job_inputs):
    
    output = {}

    options = parse_tophat_options( job_inputs )

    check_reads( job_inputs['reads'] )


    # Convert reads tables to FASTQ/FASTA files
    left_reads = []
    right_reads = []

    current_reads = 0
    for reads in job_inputs['reads']:
        left, right = dump_fastqa( reads['$dnanexus_link'], "reads_"+str(current_reads) )

        left_reads.append( left )
        if right_reads != None:
            right_reads.append( right )


    # download reference (and make index if necessary)
    if "ContigSet" in dxpy.DXRecord(job['input']['reference']).describe()['types']:
        job['output']['indexed_reference'] = make_indexed_reference()

    elif "BowtieLetterContigSetV2" in dxpy.DXRecord(job['input']['reference']).describe()['types']:
        dxpy.download_dxfile(dxpy.get_details(job["input"]["reference"])['index_archive'], "indexed_ref.tar.xz")
        run_shell("tar -xJf indexed_ref.tar.xz")

    output['indexed_reference'] = job['input']['reference']

    # Invoke tophat2 with FASTQ/A file(s) and indexed reference
    run_shell(" ".join(['tophat',"indexed_ref"," ".join(left_reads)," ".join(right_reads)]))


    return outputs

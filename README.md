RNA-Seq Analysis Pipeline
=========================

This app bundles Bowtie2, Tophat2 and Cufflinks to map RNA-Seq reads and quantitate expression.  These result in a Mappings table containing all mapped reads and a table containing per-gene expression level represented in FPKM values (Fragments Per Kilobase of transcript per Million mapped reads).  Full output of the cufflinks program is also output as a tar file which also contains expression on the per isoform level.

This pipeline does not support novel gene or isoform discovery.  Reads will only be mapped to transcripts found in the input Genes object.  This corresponds to running Tophat with the "-G" and "--transcriptome-index" options.  While these options will always be passed to Tophat, further options modifying both Tophat and Cufflinks steps are accepted by the app.

Links to the source and information about the underlying programs can be found here:
[Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
[Tophat](http://tophat.cbcb.umd.edu/)
[Cufflinks](http://cufflinks.cbcb.umd.edu/)

---------

**Inputs:**

*reads*: An array of gtables of type "Reads" which contain the RNA-seq reads.  These can be generated  from FASTQ or FASTA files with the Reads Importer app.  The pipeline support both paired and unpaired reads but all inputs must be of the same type:  all paired or all unpaired.

*Tophat Options*: A string containing all additional options to be passed to Tophat during execution.  Tophat uses the Bowtie program to align RNA-Seq reads to the given transcriptome.  A guide to Tophat options can be found [here](http://tophat.cbcb.umd.edu/manual.html#toph).  As mentioned above the "-G" and "--transcriptome-index" options will always be passed to Tophat.  The string input here will be passed directly to the Tophat program and therefore must be formatted as it would be on the command line.  The default value is "-g 1 --b2-very-fast --no-coverage-search --segment-length 100 --no-novel-juncs".  Contradictory or invalid parameters will not be caught and will cause the app to fail.

*Cufflinks Options*: (OPTIONAL) A string containing all additional options to be passed to Cufflinks during execution.  Cufflinks takes the mappings generate by Tophat and calculates the level of expression of each gene and transcript.  A guide to Cufflinks options can be found [here](http://cufflinks.cbcb.umd.edu/manual.html#cufflinks) The pipeline uses the "-G", "-p", and "-o" options to input the gene model, use all processors, and capture the output.  Do not include these as further parameters.  No additional options are set by default.

*Reference*: Reference genome for mapping.  This must be in the form of a [ContigSet](http://wiki.dnanexus.com/Types/ContigSet) object.  These can be created by importing a FASTA file using the "Genome Importer" app.  The genome must be compatible with the gene model supplied.

*Gene Model*: A [Genes](http://wiki.dnanexus.com/create/Types/Genes) object describing the transcripts to map reads to.  A GTF, GFF, or BED file can be imported to a Genes object using an importer app.

*Indexed Reference*: (OPTIONAL) An Bowtiev2 indexed copy of the reference genome.  If not supplied an index will be generated, a process that may take up to several hours.  It will then be included as part of the app output.  Including this indexed genome in later runs of the app will speed processing by skipping this indexing step.

*output name*: (OPTIONAL) Name of Mappings table output.  If not supplied the Mappings output will be "RNA-seq mappings".

**Outputs:**

*Mappings*: A gtable of type [Mappings](http://wiki.dnanexus.com/Types/Mappings) that contains the reads after alignment to the transcriptome.

*Transcripts*: A gtable containing the expression levels found for each gene as calculated by cufflinks.  This table contains the information output in the file "genes.fpkm_tracking".  The format of this file is described [here](http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format).

*Cufflinks Output*: File object, called "cufflinks_output.tar.gz", which is a compressed archive of all files output by cufflinks.  Useful for downstream analysis of expression values.

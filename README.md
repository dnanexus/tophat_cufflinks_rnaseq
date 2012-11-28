RNA-Seq Analysis Pipeline
=========================

This app bundles Bowtie2, Tophat2 and Cufflinks to map RNA-Seq reads to the human genome and quantitate expression.  Expression is represented in FPKM values (Fragments Per Kilobase of transcript per Million mapped reads).  This is calculated both on the per gene and per isoform level.

Also includes the human genome and transcriptome (from the [iGenome](http://cufflinks.cbcb.umd.edu/igenomes.html) project - version UCSC hg19).

Links to the source and information about the underlying programs can be found here:
[Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
[Tophat](http://tophat.cbcb.umd.edu/)
[Cufflinks](http://cufflinks.cbcb.umd.edu/)

---------

**Inputs:**

*reads*: An array of gtables of type "Reads" which contain the RNA-seq reads.  These can be generated  from FASTQ or FASTA files with the Reads Importer app.  The pipeline support both paired and unpaired reads but all inputs must be of the same type:  all paired or all unpaired.

*output name*: Name of Mappings table output.  This is optional.  If not supplied the Mappings output will be "RNA-seq mappings".

**Outputs:**

*mappings*: A gtable of type [Mappings](http://wiki.dnanexus.com/Types/Mappings) that contains the reads after alignment to the human transcriptome.

*transcripts*: A gtable containing the expression levels found for each gene as calculated by cufflinks.  This table contains the information output in the file "genes.fpkm_tracking".  The format of this file is described [here](http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format).

*cufflinks_output*: File object, called "cufflinks_output.tar.gz", which is a compressed archive of all files output by cufflinks.  Useful for downstream analysis of expression values.

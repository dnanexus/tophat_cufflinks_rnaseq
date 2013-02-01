all: bowtie bowtie2 cufflinks tophat samtools

bowtie:
	make -C bowtie
	cp -a bowtie/bowtie bowtie/bowtie-build bowtie/bowtie-inspect resources/usr/local/bin

bowtie2:
	make -C bowtie2
	cp -a bowtie2/bowtie2 bowtie2/bowtie2-align bowtie2/bowtie2-build bowtie2/bowtie2-inspect resources/usr/local/bin

cufflinks: libbam libeigen
	cd cufflinks2; ./configure --with-bam-lib=`pwd`/../libbam --with-bam=`pwd`/../libbam --with-eigen=`pwd`/../libeigen; make
	cp -a cufflinks2/src/cuffcompare cufflinks2/src/cuffdiff cufflinks2/src/cufflinks cufflinks2/src/cuffmerge cufflinks2/src/gffread cufflinks2/src/gtf_to_sam resources/usr/local/bin

tophat: libbam
	cd tophat2; ./configure --with-bam-lib=`pwd`/../libbam --with-bam=`pwd`/../libbam; make
	cp -a tophat2/src/bam2fastx tophat2/src/bam_merge tophat2/src/bed_to_juncs tophat2/src/closure_juncs tophat2/src/contig_to_chr_coords tophat2/src/fix_map_ordering tophat2/src/gtf_juncs tophat2/src/gtf_to_fasta tophat2/src/juncs_db tophat2/src/long_spanning_reads tophat2/src/map2gtf tophat2/src/prep_reads tophat2/src/sam_juncs tophat2/src/segment_juncs tophat2/src/sra_to_solid tophat2/src/tophat tophat2/src/tophat2 tophat2/src/tophat-fusion-post tophat2/src/tophat_reports resources/usr/local/bin
	chmod +x resources/usr/local/bin/tophat2

libeigen:
	mkdir -p libeigen/include
	cp -a eigen/* libeigen/include

libbam: samtools
	mkdir -p libbam/bin
	mkdir -p libbam/include/bam
	cp -a samtools/*.h libbam/include/bam
	cp -a samtools/libbam.a libbam/bin

samtools:
	make -C samtools
	cp -a samtools/samtools resources/usr/local/bin

clean:
	rm -f resources/usr/local/bin/*
	rm -rf libbam
	rm -rf libeigen

.PHONY: all bowtie bowtie2 cufflinks tophat libbam samtools clean

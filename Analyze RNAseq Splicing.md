# Analyze tie-1 for splicing defects

SNP analysis has suggested that there is a splice-acceptor mutation in Solyc02g080120.  Here we examine possible splicing defects  splicing in the *tie-1* GA2 Oxidase gene Solyc02g080120.

fq files for lanes 1,2,3,7 from /iplant/home/shared/ucd.plantbio/maloof.lab/members/amanda/tiernaseqfinal

fq files for lanes 4-6 from

icd /iplant/home/shared/ucd.plantbio/maloof.lab/members/upendra/kazu_palmer_libraries/processed_data/amanda_stuff/JKNE3_4re
icd /iplant/home/shared/ucd.plantbio/maloof.lab/members/upendra/kazu_palmer_libraries/processed_data/amanda_stuff/JKNE3_5re
icd /iplant/home/shared/ucd.plantbio/maloof.lab/members/upendra/kazu_palmer_libraries/processed_data/amanda_stuff/JKNE3_6re


### combine fq files

First rename to prevent clobbering, then move to common directory

in lane 4 directory:

	for file in *.fq
	    do
	        base=`basename $file .fq`
	        newName="$base.4.fq"
	        mv $file $newName
	    done
	    
repeat for all lanes, move files to a common directory and concatenate

### next run tophat
	    
    bowtie2-build S_lycopersicum_chromosomes.2.40.fa S_lycopersicum_chromosomes.2.40.fa
    
Working in directory: /home/jnmaloof/ebs5/tie1RNAseqSNPs/newTophat020114

    /usr/local/src/tophat-2.0.8b.Linux_x86_64/tophat2 --GTF ../ITAG2.3_gene_models.gff3 --output-dir tophat_M82 ../S_lycopersicum_chromosomes.2.40 ../fqfiles/M82_all.fq
    
     /usr/local/src/tophat-2.0.8b.Linux_x86_64/tophat2 --GTF ../ITAG2.3_gene_models.gff3 --output-dir tophat_tie1 ../S_lycopersicum_chromosomes.2.40 ../fqfiles/Tie_all.fq

### View in IGV.  

It appears that the first intron of Solyc02g080120 is not spliced in tie-1.  If true this would lead to a premature stop codon.

## Quantify splicing defect

Working in `newTophat020114` directory

### determine the number of reads spanning the splice junction or spliced in M82 and tie-1

For the first of these I am just going to hand count in IGV

| location     | M82 splice | M82 span | tie-1 splice | tie-1 span |
|---------------|----------------|---------------|----------------|---------------|
| 5' junction | 25              | 1               | 0                | 33             | 
| 3' junction | 25              | 2               | 0                | 15             |


### determine coverage before, in, and after splice junction in M82 and tie-1

Use samtools to calculate coverage:

    cd "/Users/jmaloof/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_tie1"
    samtools depth -r SL2.40ch02:39009924-39009943 accepted_hits.bam > tie1_5prime.tsv
    samtools depth -r SL2.40ch02:39009944-39010179 accepted_hits.bam > tie1_intron.tsv
    samtools depth -r SL2.40ch02:39010180-39010199 accepted_hits.bam > tie1_3prime.tsv
    samtools idxstats accepted_hits.bam > tie1_hits.tsv
    
    cd ../tophat_M82
    samtools depth -r SL2.40ch02:39009924-39009943 accepted_hits.bam > M82_5prime.tsv
    samtools depth -r SL2.40ch02:39009944-39010179 accepted_hits.bam > M82_intron.tsv
    samtools depth -r SL2.40ch02:39010180-39010199 accepted_hits.bam > M82_3prime.tsv
    samtools idxstats accepted_hits.bam > M82_hits.tsv
    
For subsequent analysis see `Analyze RNAseq Splicing.R`


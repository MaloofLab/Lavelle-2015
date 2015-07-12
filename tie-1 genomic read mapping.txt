# Mapping tie-1 mutation

## Background
Amanda and Natalie found tie-1 in the Zamir collection.  This mutant has elongated internodes.

tie-1 is in M82

Amanda and Natalie crossed to Heinz, selfed F1s, and grew out F2s

~ 250 mutant F2s were pooled and sequenced

Now I will try to map the mutation

## Sequence files

working in directory ~/ebs5/tie1 on an atmosphere instance

get files from irods with:

    iget /iplant/home/shared/ucd.tomato/processed.seq.data/filtered_fq.gz_files/163.20130501.JMAS003.PE1.fltr.fq
    iget /iplant/home/shared/ucd.tomato/processed.seq.data/filtered_fq.gz_files/163.20130501.JMAS003.PE2.fltr.fq
    
used my runBWA.pl script to map, in SE mode

Remove duplicates with picard
    
    java -Xmx3g -jar /usr/local/bin/MarkDuplicates.jar INPUT=163.20130501.JMAS003.PE1.fltr.bam \
    	 OUTPUT=PE1_nodup.bam \
    	 ASSUME_SORTED=true \
    	 REMOVE_DUPLICATES=true \
    	 METRICS_FILE=PE1_duplicates.txt \
    	 CREATE_INDEX=true \
    	 VALIDATION_STRINGENCY=LENIENT \
    	 TMP_DIR=./
    
        java -Xmx3g -jar /usr/local/bin/MarkDuplicates.jar INPUT=163.20130501.JMAS003.PE2.fltr.bam \
    	 OUTPUT=PE2_nodup.bam \
    	 ASSUME_SORTED=true \
    	 REMOVE_DUPLICATES=true \
    	 METRICS_FILE=PE2_duplicates.txt \
    	 CREATE_INDEX=true \
    	 VALIDATION_STRINGENCY=LENIENT \
    	 TMP_DIR=./
    	 
Use samtools to do a mpileup:

    samtools mpileup -B -f S_lycopersicum_chromosomes.2.40.fa PE2_nodup.bam PE1_nodup.bam > PE1_PE2.mpileup
    
And then Varscan to call SNPs (with default values)

    java -jar ~/bin/VarScan.v2.3.5.jar mpileup2snp PE1_PE2.mpileup > PE1_PE2.varscan
    
Plotting with the script [plot_tie1.R](<file:///Users/jmaloof/Documents/Lab%20Notebook%20support/2013/tie1_phyE/plot_tie1.R>) clearly shows that tie-1 lies on chromosome 2.

![](file:///Users/jmaloof/Documents/Lab%20Notebook%20support/2013/tie1_phyE/tie1AllChrom.png)

It is very near phyE (red line in the plot below).

![](file:///Users/jmaloof/Documents/Lab%20Notebook%20support/2013/tie1_phyE/tie1closeup.png)

But strangely I find no mutations in phyE itself in spite of adequate coverage.  I am missing coverage upstream and in the 5' UTR.  I am now going to try mapping in paired-end mode to see if I can get more information about what is going on upstream.  An alternative explanation is that there is a deletion upstream…

### Paired end mapping

    bwa sampe S_lycopersicum_chromosomes.2.40.fa 163.20130501.JMAS003.PE1.fltr.sai 163.20130501.JMAS003.PE2.fltr.sai \
        163.20130501.JMAS003.PE1.fltr.fq 163.20130501.JMAS003.PE2.fltr.fq | samtools view -bS - | samtools sort -m 300000000 - PE1PE2
        
    samtools index PE1PE2.bam PE1PE2.bai
    
    samtools idxstats PE1PE2.bam > PE1PE2.counts
    
        java -Xmx3g -jar /usr/local/bin/MarkDuplicates.jar INPUT=PE1PE2.bam \
    	 OUTPUT=PE1PE2_nodup.bam \
    	 ASSUME_SORTED=true \
    	 REMOVE_DUPLICATES=true \
    	 METRICS_FILE=PE2_duplicates.txt \
    	 CREATE_INDEX=true \
    	 VALIDATION_STRINGENCY=LENIENT \
    	 TMP_DIR=./
    	 
    	  samtools mpileup -B -f S_lycopersicum_chromosomes.2.40.fa PE1PE2.bam > PE1PE2_paired.mpileup
    	  
    	  java -jar ~/bin/VarScan.v2.3.5.jar mpileup2snp PE1PE2_paired.mpileup > PE1PE2_paired.varscan
    	  
### I also want to compare this to M82.  Use reads at /Volumes/OuterColomaData/Mike/map.ILs_to_Slyc.bwa_tophat/M82_n05/merged/bwa_tophat_M82_n05-Slyc.sorted.dupl_rm.bam

assign readgroups:

    java -Xmx3g -jar /usr/local/bin/AddOrReplaceReadGroups.jar \
    	INPUT=PE1PE2_nodup.bam \
    	OUTPUT=PE1PE2_nodup.RG.bam \
    	RGID=tie1 \
    	RGSM=tie1 \
    	RGLB=1 \
    	RGPU=NA \
    	RGPL=Illumina \
    	CREATE_INDEX=true \
      	VALIDATION_STRINGENCY=LENIENT \
    	TMP_DIR=./ 
    	
    java -Xmx3g -jar /usr/local/bin/AddOrReplaceReadGroups.jar \
    	INPUT=bwa_tophat_M82_n05-Slyc.sorted.dupl_rm.bam \
    	OUTPUT=M82_nodup.RG.bam \
    	RGID=M82 \
    	RGSM=M82 \
    	RGLB=1 \
    	RGPU=NA \
    	RGPL=Illumina \
    	CREATE_INDEX=true \
      	VALIDATION_STRINGENCY=LENIENT \
    	TMP_DIR=./ 

    	 

Combined pileup on chrom 2 only

Make bedfile to specify region
	
	cat > chr2.bed
	SL2.40ch02      1       49918234

    samtools mpileup -B -l chr2.bed -f S_lycopersicum_chromosomes.2.40.fa PE1PE2_nodup.RG.bam  M82_nodup.RG.bam > PE1PE2_paired_M82.mpileup
    
    java -Xmx3g -jar ~/bin/VarScan.v2.3.5.jar mpileup2snp PE1PE2_paired_M82.mpileup > PE1PE2_paired_M82.varscan

Also generate chr2 bcf file

    samtools mpileup -BgD -l chr2.bed -f S_lycopersicum_chromosomes.2.40.fa PE1PE2_nodup.RG.bam  M82_nodup.RG.bam > PE1PE2_paired_M82.bcf

first make filter for variant sites only and convert to VCF

    bcftools view -Nv PE1PE2_paired_M82.bcf > PE1PE2_paired_M82.vcf

filter: not needed

	#sed s/,*X/N/g PE1PE2_paired_M82.vcf > PE1PE2_paired_M82_filtered.vcf
	
    java -Xmx3g -jar ../snpEff/snpEff.jar eff -c ../snpEff/snpEff.config Solyc2.40 PE1PE2_paired_M82.vcf > annotated_PE1PE2_M82.vcf
    
MEanwhile, scanning through BAM files on IGV, looking for mutations in tie1 that are not in M82.  Going from 35,000,000 to 35,500,000 on Chrom 2

	* Solyc02g071100 has possible deletion in tie1 not in M82 (?)
	
	* Solyc02g069570 has a CTA -> TTA mutation in  most tie-1 reads that is not in Heinz or M82.  This is a silent mutation.
    	  
    	 * Solyc02g069670.2.1 Alpha-glucosidase has a mutation in tie1, not in M82, but not quite homozygous in F2s

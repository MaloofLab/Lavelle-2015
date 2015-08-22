## Background

Trying to ID *tie-1* mutation.

Will map both _tie-1_ BSA F2 genomic reads and M82 Genomic reads to Heinz, and then filter SNPs.

### _tie-1_ Paired end mapping

Note reads already filtered

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

### M82 Concatenate reads

    cat /mnt/ebs2/uploads/EAS517_0034_FC619A6AAXX/*1_sequence.txt.gz /mnt/ebs2/uploads/EAS67_0101_PEFC619A7AAXX/*1_sequence.txt.gz > M82G_PE1.txt.gz &
    cat /mnt/ebs2/uploads/EAS517_0034_FC619A6AAXX/*2_sequence.txt.gz /mnt/ebs2/uploads/EAS67_0101_PEFC619A7AAXX/*2_sequence.txt.gz > M82G_PE2.txt.gz &
    
### M82 filter reads

Looks like M82 reads are Illumina 1.5 format (Phred+64)

    zcat M82G_PE1.txt.gz | fastq_quality_trimmer -t 20 -l 50  | fastq_quality_filter -q 30 -p 90 -z > M82G_PE1_filtered.fq.gz &
    zcat M82G_PE2.txt.gz | fastq_quality_trimmer -t 20 -l 50  | fastq_quality_filter -q 30 -p 90 -z > M82G_PE2_filtered.fq.gz &

### M82 map reads

     bwa aln -t 4 -I ../S_lycopersicum_chromosomes.2.40.fa  M82G_PE1_filtered.fq.gz > M82G_PE1.sai
     bwa aln -t 4 -I ../S_lycopersicum_chromosomes.2.40.fa  M82G_PE2_filtered.fq.gz > M82G_PE2.sai
     
     bwa sampe -P -r '@RG\tID:M82G\tSM:M82G' ../S_lycopersicum_chromosomes.2.40.fa M82G_PE1.sai M82G_PE2.sai M82G_PE1_filtered.fq.gz M82G_PE2_filtered.fq.gz > M82G_PE1PE2.sam
    
This seems not to be working; BWA is estimating insert sizes in the tens of thousands but Tony says it should be a couple of hundred.

I maybe should come back to this, but for now will use samse instead

    bwa samse -r '@RG\tID:M82G\tSM:M82G' ../S_lycopersicum_chromosomes.2.40.fa M82G_PE1.sai M82G_PE1_filtered.fq.gz  > M82G_PE1.sam &
    
    bwa samse -r '@RG\tID:M82G\tSM:M82G' ../S_lycopersicum_chromosomes.2.40.fa M82G_PE2.sai M82G_PE2_filtered.fq.gz  > M82G_PE2.sam &
     
### M82 remove duplicates

        samtools view -buS M82G_PE1.sam | samtools sort -o -m 1000000000 - M82GPE1tmpsort | samtools rmdup - M82G_PE1_rmdup.bam
        samtools index  M82G_PE1_rmdup.bam 
             
        samtools view -buS M82G_PE2.sam | samtools sort -o -m 1000000000 - M82GPE2tmpsort | samtools rmdup - M82G_PE2_rmdup.bam
        samtools index  M82G_PE2_rmdup.bam
     
### tie-1 and M82 realign

For TIE1:

     java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T RealignerTargetCreator \
       -R S_lycopersicum_chromosomes.2.40.fa \
       -I PE1PE2_nodup.RG.bam \
       -o PE1PE2_nodup.RG.forIndelRealigner.intervals 
       
    java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T IndelRealigner \
       -R S_lycopersicum_chromosomes.2.40.fa \
       -I PE1PE2_nodup.RG.bam \
       -targetIntervals PE1PE2_nodup.RG.forIndelRealigner.intervals \
       -o PE1PE2_nodup.RG.realigned.bam 
       
    samtools index PE1PE2_nodup.RG.realigned.bam PE1PE2_nodup.RG.realigned.bai
       
For M82PE1:

     java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T RealignerTargetCreator \
       -R ../S_lycopersicum_chromosomes.2.40.fa \
       -I M82G_PE1_rmdup.bam \
       -o M82G_PE1_rmdup.bam.forIndelRealigner.intervals 
       
    java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T IndelRealigner \
       -R ../S_lycopersicum_chromosomes.2.40.fa \
       -I M82G_PE1_rmdup.bam \
       -targetIntervals M82G_PE1_rmdup.bam.forIndelRealigner.intervals  \
       -o M82G_PE1_rmdup.realigned.bam
       
    samtools index M82G_PE1_rmdup.realigned.bam M82G_PE1_rmdup.realigned.bai
       
For M82PE2:

     java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T RealignerTargetCreator \
       -R ../S_lycopersicum_chromosomes.2.40.fa \
       -I M82G_PE2_rmdup.bam \
       -o M82G_PE2_rmdup.bam.forIndelRealigner.intervals 
       
    java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T IndelRealigner \
       -R ../S_lycopersicum_chromosomes.2.40.fa \
       -I M82G_PE2_rmdup.bam \
       -targetIntervals M82G_PE2_rmdup.bam.forIndelRealigner.intervals  \
       -o M82G_PE2_rmdup.realigned.bam
       
    samtools index M82G_PE2_rmdup.realigned.bam M82G_PE2_rmdup.realigned.bai
       
#### merge the two M82G files

    samtools merge -h M82G_PE1_rmdup.realigned.bam M82G_PE1PE2_rmdup.realigned.bam M82G_PE1_rmdup.realigned.bam M82G_PE2_rmdup.realigned.bam
    
    samtools index M82G_PE1PE2_rmdup.realigned.bam M82G_PE1PE2_rmdup.realigned.bai
        
#### tie-1 and M82 pileup

    samtools mpileup -B -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.realigned.bam  M82G_PE1PE2_rmdup.realigned.bam > PE1PE2_paired_M82G.realigned.all.mpileup
    
    java -Xmx3g -jar /usr/local/bin/VarScan.v2.3.5.jar mpileup2snp PE1PE2_paired_M82G.realigned.all.mpileup > PE1PE2_paired_M82G.realigned.all.varscan
    

### Generate bcf for SnpEff to use

    samtools mpileup -BgD -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.realigned.bam  M82G_PE1PE2_rmdup.realigned.bam > PE1PE2_paired_M82G.realigned.bcf  
 
    bcftools view -Nv PE1PE2_paired_M82G.realigned.bcf   > PE1PE2_paired_M82G.realigned.vcf  

But if we use the pre-built database in SnpEff we need to adjust the chromnames:

    sed 's/SL2.40ch0*//' PE1PE2_paired_M82G.realigned.vcf  > PE1PE2_paired_M82G.realigned.renamed.vcf  
  
    java -Xmx3g -jar ../snpEff/snpEff.jar eff  SL2.40.26  PE1PE2_paired_M82G.realigned.renamed.vcf   >  annotated_PE1PE2_paired_M82G.realigned.renamed,vcf  
  
  ### Generate tie-1 only bcf for SnpEff to use (probably there is a better way, but...)
  
      samtools mpileup -BgD -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.realigned.bam   >  tie1_only_PE1PE2_paired.realigned.bcf  
      
       bcftools view -Nv tie1_only_PE1PE2_paired.realigned.bcf    > tie1_only_PE1PE2_paired.realigned.vcf   

       sed 's/SL2.40ch0*//' tie1_only_PE1PE2_paired.realigned.vcf > tie1_only_PE1PE2_paired.realigned.renamed.vcf
  
      java -Xmx3g -jar ../snpEff/snpEff.jar eff  -no-downstream -no-intergenic -t SL2.40.26 tie1_only_PE1PE2_paired.realigned.renamed.vcf   >  annotated_tie1_only_PE1PE2_paired.realigned.renamed.vcf

## Redo with FreeBayes

    bamaddrg -b M82G_PE1PE2_rmdup.realigned.bam -r M82 > M82G_PE1PE2_rmdup.realigned.rg.bam

    freebayes -f --use-best-n-alleles 4 -C 5 ../ M82G_PE1PE2_rmdup.realigned.rg.bam ../ 

### Further Analysis

See `plot_tie1M82G_clean.R`

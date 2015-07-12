## Background

Trying to ID *tie-1* mutation.
Mapping against Tony Bolger's M82 reference seemed really noisy.  I now have M82 genomic raw reads from him.  Will try mapping these along with tie-1 to the Heinz reference.  

### Concatenate reads

    cat /mnt/ebs2/uploads/EAS517_0034_FC619A6AAXX/*1_sequence.txt.gz /mnt/ebs2/uploads/EAS67_0101_PEFC619A7AAXX/*1_sequence.txt.gz > M82G_PE1.txt.gz &
    cat /mnt/ebs2/uploads/EAS517_0034_FC619A6AAXX/*2_sequence.txt.gz /mnt/ebs2/uploads/EAS67_0101_PEFC619A7AAXX/*2_sequence.txt.gz > M82G_PE2.txt.gz &
    
### filter reads

Looks like M82 reads are Illumina 1.5 format (Phred+64)

    zcat M82G_PE1.txt.gz | fastq_quality_trimmer -t 20 -l 50  | fastq_quality_filter -q 30 -p 90 -z > M82G_PE1_filtered.fq.gz &
    zcat M82G_PE2.txt.gz | fastq_quality_trimmer -t 20 -l 50  | fastq_quality_filter -q 30 -p 90 -z > M82G_PE2_filtered.fq.gz &

### map reads

     bwa aln -t 4 -I ../S_lycopersicum_chromosomes.2.40.fa  M82G_PE1_filtered.fq.gz > M82G_PE1.sai
     bwa aln -t 4 -I ../S_lycopersicum_chromosomes.2.40.fa  M82G_PE2_filtered.fq.gz > M82G_PE2.sai
     
     bwa sampe -P -r '@RG\tID:M82G\tSM:M82G' ../S_lycopersicum_chromosomes.2.40.fa M82G_PE1.sai M82G_PE2.sai M82G_PE1_filtered.fq.gz M82G_PE2_filtered.fq.gz > M82G_PE1PE2.sam
    
This seems not to be working; BWA is estimating insert sizes in the tens of thousands but Tony says it should be a couple of hundred.

I maybe should come back to this, but for now will use samse instead

    bwa samse -r '@RG\tID:M82G\tSM:M82G' ../S_lycopersicum_chromosomes.2.40.fa M82G_PE1.sai M82G_PE1_filtered.fq.gz  > M82G_PE1.sam &
    
    bwa samse -r '@RG\tID:M82G\tSM:M82G' ../S_lycopersicum_chromosomes.2.40.fa M82G_PE2.sai M82G_PE2_filtered.fq.gz  > M82G_PE2.sam &
     
### remove duplicates

        samtools view -buS M82G_PE1.sam | samtools sort -o -m 1000000000 - M82GPE1tmpsort | samtools rmdup - M82G_PE1_rmdup.bam
        samtools index  M82G_PE1_rmdup.bam 
             
        samtools view -buS M82G_PE2.sam | samtools sort -o -m 1000000000 - M82GPE2tmpsort | samtools rmdup - M82G_PE2_rmdup.bam
        samtools index  M82G_PE2_rmdup.bam
     
### pileup

    samtools mpileup -B -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.bam  M82G_PE1_rmdup.bam M82G_PE2_rmdup.bam > PE1PE2_paired_M82G.mpileup
    
    java -Xmx3g -jar /usr/local/bin/VarScan.v2.3.5.jar mpileup2snp PE1PE2_paired_M82G.mpileup > PE1PE2_paired_M82G.varscan
    
Also generate chr2 bcf file

    samtools mpileup -BgD -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.bam  M82G_PE1_rmdup.bam M82G_PE2_rmdup.bam > PE1PE2_paired_M82G.mpileup.bcf
    
    bcftools view -Nv PE1PE2_paired_M82G.mpileup.bcf > PE1PE2_paired_M82G.mpileup.vcf
    	
    java -Xmx3g -jar /home/jnmaloof/ebs5/snpEff/snpEff.jar eff -c /home/jnmaloof/ebs5/snpEff/snpEff.config PE1PE2_paired_M82G.mpileup.vcf > annotated_PE1PE2_paired_M82G.mpileup.vcf

### realign

Going to go back and do a realignment on the bam files

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
        
#### redo the pileup

    samtools mpileup -B -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.realigned.bam  M82G_PE1PE2_rmdup.realigned.bam > PE1PE2_paired_M82G.realigned.mpileup
    
    java -Xmx3g -jar /usr/local/bin/VarScan.v2.3.5.jar mpileup2snp PE1PE2_paired_M82G.realigned.mpileup > PE1PE2_paired_M82G.realigned.varscan
    
### look for indels

        java -Xmx3g -jar /usr/local/bin/VarScan.v2.3.5.jar mpileup2indel PE1PE2_paired_M82G.realigned.mpileup > PE1PE2_paired_M82G.realigned.varscan_indel

### Generate bcf for SnpEff to use

    samtools mpileup -BgD -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.realigned.bam  M82G_PE1PE2_rmdup.realigned.bam > PE1PE2_paired_M82G.realigned.bcf  
 
  bcftools view -Nv PE1PE2_paired_M82G.realigned.bcf   > PE1PE2_paired_M82G.realigned.vcf  
  
  java -Xmx3g -jar ~/ebs5//snpEff/snpEff.jar eff -c ~/ebs5/snpEff/snpEff.config Solyc2.40  PE1PE2_paired_M82G.realigned.vcf   >  annotated_PE1PE2_paired_M82G.realigned.vcf  
  
  ### Generate tie-1 only bcf for SnpEff to use (probably there is a better way, but...)
  
      samtools mpileup -BgD -l ../chr2.bed -f ../S_lycopersicum_chromosomes.2.40.fa ../PE1PE2_nodup.RG.realigned.bam   >  tie1_only_PE1PE2_paired.realigned.bcf  
      
       bcftools view -Nv tie1_only_PE1PE2_paired.realigned.bcf    > tie1_only_PE1PE2_paired.realigned.vcf   
  
  java -Xmx3g -jar ~/ebs5//snpEff/snpEff.jar eff -c ~/ebs5/snpEff/snpEff.config -no-downstream -no-intergenic -t Solyc2.40  tie1_only_PE1PE2_paired.realigned.vcf   >  annotated_tie1_only_PE1PE2_paired.realigned.vcf
  
### Next

* what column is SnpEff acting on?
	* maybe I should just feed it the tie1 info
* Consider places where M82 is het.
* 


     

    


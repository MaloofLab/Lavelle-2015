# Analyze RNAseq for SNPs
Here we use RNAseq data from _tie-1_ and M82 to look for mutations in the _tie-1_ region.  For most genes, coverage should be higher than what we had in genomic.

### sort and merge

    samtools sort M82-1.genome.bam M82-1.sorted &
    samtools sort M82-2.genome.bam M82-2.sorted &
    samtools sort M82-3.genome.bam M82-3.sorted &
    
    samtools sort Tie-1.genome.bam Tie-1.sorted &
    samtools sort Tie-2.genome.bam Tie-2.sorted &
    samtools sort Tie-3.genome.bam Tie-3.sorted &

    samtools merge M82All.bam M82-*sorted.bam &
    samtools merge tie1All.bam Tie-*.sorted.bam &
    
### remove duplicate reads

    samtools rmdupM82All.bam M82All.rmdup.bam &
    samtools rmdup tie1All.bam tie1All.rmdup.bam &
    
## Assign Readgroups

    java -Xmx3g -jar ~/bin/picard-tools-1.103/AddOrReplaceReadGroups.jar \
    	INPUT=M82All.rmdup.bam \
    	OUTPUT=M82All.rmdup.RG.bam \
    	RGID=M82 \
    	RGSM=M82 \
    	RGLB=1 \
    	RGPU=NA \
    	RGPL=Illumina \
    	CREATE_INDEX=true \
      	VALIDATION_STRINGENCY=LENIENT \
    	TMP_DIR=./ 
    	
    java -Xmx3g -jar ~/bin/picard-tools-1.103/AddOrReplaceReadGroups.jar \
    	INPUT=tie1All.rmdup.bam \
    	OUTPUT=tie1All.rmdup.RG.bam \
    	RGID=tie1 \
    	RGSM=tie1 \
    	RGLB=1 \
    	RGPU=NA \
    	RGPL=Illumina \
    	CREATE_INDEX=true \
      	VALIDATION_STRINGENCY=LENIENT \
    	TMP_DIR=./ 
    	
index seemed not to work, so

    samtools index tie1All.rmdup.RG.bam
    samtools index M82All.rmdup.RG.bam

### Realign

I get an error because of bad reads.

Writing short python script to discard bad reads (in this case reads without adequate quality info) `filterBam.py`

Now...

     java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T RealignerTargetCreator \
       -R ../tie1/S_lycopersicum_chromosomes.2.40.fa \
       -I M82All.rmdup.RG.filtered.bam \
       -o M82.forIndelRealigner.intervals \
       -DBQ 1
       
    java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T IndelRealigner \
       -R ../tie1/S_lycopersicum_chromosomes.2.40.fa \
       -I M82All.rmdup.RG.filtered.bam \
       -targetIntervals M82.forIndelRealigner.intervals  \
       -o M82.rmdup.realigned.bam \
       -DBQ 1
       
    samtools index M82.rmdup.realigned.bam M82.rmdup.realigned.bai
    
    java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T RealignerTargetCreator \
       -R ../tie1/S_lycopersicum_chromosomes.2.40.fa \
       -I tie1All.rmdup.RG.filter.bam \
       -o tie1.forIndelRealigner.intervals \
       -DBQ 1
   
    java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar \
       -T IndelRealigner \
       -R ../tie1/S_lycopersicum_chromosomes.2.40.fa \
       -I tie1All.rmdup.RG.filter.bam \
       -targetIntervals tie1.forIndelRealigner.intervals  \
       -o tie1.rmdup.realigned.bam \
       -DBQ 1
       
    samtools index tie1.rmdup.realigned.bam tie1.rmdup.realigned.bai
    
### mpileup and snpEFF

    samtools mpileup -B -D  -f ../tie1/S_lycopersicum_chromosomes.2.40.fa -l ../tie1/chr2.bed M82All.rmdup.RG.bam tie1All.rmdup.RG.bam > M82_tie1_mpileupNEW.mpileup
        
    java -Xmx3g -jar /usr/local/bin/VarScan.v2.3.5.jar mpileup2cns M82_tie1_mpileupNEW.mpileup \
         --variants 1 --min-coverage 4 --min-var-freq .2 --output-vcf 1 --vcf-sample-list > M82_tie1_mpileupNEW.varscan
         
    java -Xmx3g -jar /home/jnmaloof/ebs5/snpEff/snpEff.jar -c /home/jnmaloof/ebs5/snpEff/snpEff.config -no-downstream -no-intergenic -no-upstream -no-intron -no-utr  -v -t Solyc2.40 M82_tie1_mpileupNEW.varscan > anotated_M82_tie1.varscanNEW.vcf
    
### Further Analysis

See `Analyze RNAseq for SNPs.R`

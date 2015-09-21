# Take a fresh start on genomic mapping of tie-1 F2 BSA reads

## Tie1 Reads:

    gzip 163.20130501.JMAS003.PE1.fltr.fq
    gzip 163.20130501.JMAS003.PE2.fltr.fq
    cat 163.20130501.JMAS003.PE1.fltr.fq.gz 163.20130501.JMAS003.PE2.fltr.fq.gz > tie1_PE1PE2.fltr.fq.gz
    bwa mem -R @RG\tID:tie1\tSM:tie1 ../../S.lyc/S_lycopersicum_chromosomes.2.40.fa tie1_PE1PE2.fltr.fq.gz | samtools view -Sbu - | samtools rmdup - - | samtools sort - tie1_PE1PE2_rmdup
    samtools index tie1_PE1PE2_rmdup.bam


files are `163.20130501.JMAS003.PE1.fltr.fq.gz` and `163.20130501.JMAS003.PE2.fltr.fq.gz`

## M82 Reads

files are `M82G_PE1_filtered.fq.gz` and `M82G_PE2_filtered.fq.gz`

I am going to treat these as SE reads because they aren't matched.

    cat * > M82G_PE1PE2_filtered.fq.gz

    bwa mem -R '@RG\tID:M82\tSM:M82' ../../S.lyc/S_lycopersicum_chromosomes.2.40.fa M82G_PE1PE2_filtered.fq.gz | samtools view -Sbu - | samtools rmdup - - | samtools sort - M82G_PE1PE2_rmdup
    samtools index M82G_PE1PE2_rmdup.bam

## heinz

    fastq-dump --gzip SRR404081

    bwa mem -t 8 -R '@RG\tID:Heinz\tSM:Heinz' ../../S.lyc/S_lycopersicum_chromosomes.2.40.fa SRR404081.fastq.gz | samtools view -Sbu - | samtools rmdup - - | samtools sort -m 4000000000 - HeinzG_PE1PE2_rmdup

## Call SNPs

    freebayes -f ../S.lyc/S_lycopersicum_chromosomes.2.40.fa tie1F2G/tie1_PE1PE2_rmdup.bam M82G/M82G_PE1PE2_rmdup.bam heinz/HeinzG_PE1PE2_rmdup.bam > tie1_M82_Heinz_fb.vcf

Do some filtering

    vcfqualfilter -c 20 < tie1_M82_Heinz_fb.vcf > tie1_M82_Heinz_fb.filter20.vcf

## Annotation

Making my own genome database file for SnpEff seems not to be working well.  Using the pre-built one is annoying because the chromosome names are different, but that is what I am going to try.  This requires changing the chromosome names in the vcf file to match.

Download pre-built database: 
    
    java -jar /usr/local/bin/snpEff.jar download -v SL2.40.26

Take a look at nomenclature

    java -jar /usr/local/bin/snpEff.jar dump SL2.40.26 | more

Chromosomes are labelled from 1 to 12.  Chromosome 0 is ''

    gunzip -c tie1_M82_Heinz_fb.filter20.vcf.gz |  awk '{ sub(/^SL2.40ch0{0,2}/,""); print }' | gzip -c > tie1_M82_Heinz_fb.filter20.vcf.sub_chrom.gz

    java -jar /usr/local/bin/snpEff.jar ann -t -fi tie1_region.bed  -noStats -no-intergenic -no-upstream -no-downstream -v SL2.40.26 tie1_M82_Heinz_fb.filter20.vcf.sub_chrom.gz > annotated.tie1_M82_Heinz_fb.filter20.vcf

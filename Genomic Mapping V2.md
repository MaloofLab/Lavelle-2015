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
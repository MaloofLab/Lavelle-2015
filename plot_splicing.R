# Script to make plot of altered GA2OX7/8 transcript in tie-1
# Julin Maloof
# Aug 19, 2015

setwd("~/git/Lavelle-tie1--2015/")

library(ggbio)
library(Biostrings)
library(GenomicFeatures)
library(genomeIntervals)
library(biovizBase)


Slchroms <- readDNAStringSet("~/Sequences/ref_genomes/tomato/ITAG2.3_Chromo2.4/S_lycopersicum_chromosomes.2.40.fa")
  # Can be downloaded from SGN ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/wgs/assembly/build_2.40/S_lycopersicum_chromosomes.2.40.fa.gz

chrom.info <- data.frame(
  chrom=names(Slchroms),
  length=nchar(Slchroms),
  is_circular=FALSE)

chrom.info

txdb <- makeTxDbFromGFF(file = "~/Sequences/ref_genomes/tomato/ITAG2.3_Chromo2.4/ITAG2.3_gene_models.gff3",
                        chrominfo = chrom.info,
                        dataSource="ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG2.3_release/ITAG2.3_gene_models.gff3")
        
#test out select function and see what is in txdb
select(txdb,keys = "Solyc02g080120.1",keytype="GENEID",columns = columns(txdb))

#Get a Grange object for the TIE1 gene
TIE1.Grange <- crunch(txdb,which=list(gene_id="Solyc02g080120.1"))

#limit this to exons
TIE1.Grange <- TIE1.Grange[TIE1.Grange$type=="exon",]

#now make a Grange object for the tie1 mutant (first intron skipped)
tie1.Grange <- TIE1.Grange[-2]
ranges(tie1.Grange[1]) <- IRanges(start=start(ranges(TIE1.Grange[1])),
                                 end=end(ranges(TIE1.Grange[2])))

#Combine the two Grange objects
TIE1.Grange$gt <- "wild type"
tie1.Grange$gt <- "tie1"
TIE1.tie1.gr <- GRangesList(tie1=tie1.Grange,WT=TIE1.Grange)

#plot it
pl <- ggplot(TIE1.tie1.gr) + geom_alignment(aes(y=gt,fill=gt),
                                            gap.geom="chevron",
                                            label=FALSE)
pl

#or for separate plots

TIE1.splice <- ggplot(TIE1.Grange) + geom_alignment(fill="skyblue")

tie1.splice <- ggplot(tie1.Grange) + geom_alignment(fill="red")
#can we add bam files?

# not working...
plotSpliceSum(data="~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_tie1/accepted_hits.bam",
                   model=TIE1.tie1.gr)

#try another way

#create a "which" object that encompasses the TIE1 gene
wh <- range(TIE1.Grange,ignore.strand=TRUE)

ranges(wh) <- IRanges(start=start(ranges(wh))-100, #expand the range a bit
                      end=end(ranges(wh))+100)

tie1.bam <- autoplot("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_tie1/accepted_hits.bam",
         which=wh,
         geom="gapped.pair") + ylim(0,40)

TIE1.bam <-  autoplot("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_M82/accepted_hits.bam",
                 which=wh,
                 geom="gapped.pair") + ylim(0,40)

tie1.bam.coverage <- autoplot("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_tie1/accepted_hits.bam",
                     which=wh)

TIE1.bam.coverage <-  autoplot("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_M82/accepted_hits.bam",
                      which=wh) 

tracks(WT=TIE1.splice, "WT coverage"=TIE1.bam.coverage, tie1=tie1.splice, "tie1 reads"=tie1.bam.coverage, heights=c(1,5,1,5))

tracks(WT=TIE1.splice, "WT reads"=TIE1.bam, tie1=tie1.splice, "tie1 reads"=tie1.bam, heights=c(1,5,1,5))

ggsave("tie1.splicing.pdf",height=6.5,width=6.5)


# try blow-up of intron one region to make it easier to see

wh.small <- range(TIE1.Grange[1:2],ignore.strand=TRUE)

ranges(wh.small) <- IRanges(start=start(ranges(wh.small)), #tweak the range a bit...
                      end=end(ranges(wh.small)-50))

tie1.bam.small <- autoplot("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_tie1/accepted_hits.bam",
                     which=wh.small,
                     geom="gapped.pair")

TIE1.bam.small <-  autoplot("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_M82/accepted_hits.bam",
                      which=wh.small,
                      geom="gapped.pair")

tracks(TIE1.bam.small,tie1.bam.small) # could add this as a blow-up to the paper
ggsave("tie1_intron1.pdf",height=6.5,width=6.5)




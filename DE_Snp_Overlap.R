setwd("~/git/Lavelle-tie-1--2015/")
library(rtracklayer)
library(plyr)

DE.data <- read.csv("tie1_DEgenes.csv",row.names=1,stringsAsFactors = FALSE)

head(DE.data)

names(DE.data)[1] <- "ITAG"

gff <- import.gff3("ITAG2.3_genes_only.gff3")
#create this file with 
# grep "\tgene\t" ITAG2.3_Chromo2.4/ITAG2.3_gene_models.gff3 > ITAG2.3_genes_only.gff3
# this avoids records with errors and all we care about here is gene space

gff.DE.genes <- gff[unlist(gff$Alias) %in% DE.data$ITAG]

#exapnd window around each gene

start(gff.DE.genes) <- start(gff.DE.genes) - 2000
end(gff.DE.genes) <- end(gff.DE.genes) + 2000

snps <- read.delim("tie1_M82_Heinz_fb.filter20.vcf.gz",stringsAsFactors = FALSE,comment.char = "#",header=FALSE)

colnames(snps) <- unlist(
  strsplit(
    sub("#","",system("zgrep '^#[A-Z,a-z]' tie1_M82_Heinz_fb.filter20.vcf.gz",intern = TRUE)),
    split = "\t")
)

snps <- snps[snps$CHROM=="SL2.40ch02",]

# the data we are interested in is in colon-delimited fields.  Here we create a function to split it up.
split.gt <- function(gt.data,format,prefix) {
  # data is a character vector with the data to be split,
  # format is a string indicating the labels for each column, 
  # and prefix is a character string to be pre-pended to each new label
  
  gt.data <- sub("^\\.$","NA:NA:NA:NA:NA:NA:NA",gt.data)
  
  tmp <- ldply(strsplit(gt.data,split=":"))
  
  colnames(tmp) <- paste(prefix,unlist(strsplit(unique(format)[1],split=":")),sep=".")
  tmp[,2:6] <- apply(tmp[,2:6],2,as.numeric)
  
  tmp <- within(tmp,assign(paste(prefix,"percent.alt",sep="."),tmp[,5]/tmp[,2]*100))
  
  tmp
}

snps <- cbind(snps,
              split.gt(snps$tie1,snps$FORMAT,"tie1"),
              split.gt(snps$M82,snps$FORMAT,"M82"),
              split.gt(snps$Heinz,snps$FORMAT,"Heinz")
)

head(snps[,-8])

#filter snps to focus on relevant SNPs

#filter for  coverage.  

snps.small <- snps[snps$M82.DP > 4 & snps$M82.DP < 200 ,] 

snps.small <- snps.small[snps.small$Heinz.DP > 4 & snps.small$Heinz.DP < 200 ,] 

snps.small <- snps.small[snps.small$tie1.DP > 4 & snps.small$tie1.DP < 150 ,]

# additional filters

#get rid of weird SNPs
snps.small <- snps.small[grepl("0/0",snps.small$Heinz.GT),]

#focus on SNPs where M82 is fixed different from Hz
snps.small <- snps.small[snps.small$M82.percent.alt> 90,]

#only look at SNPs where we have information for all genotypes
snps.small <- snps.small[!is.na(snps.small$tie1.DP),] 
snps.small <- snps.small[!is.na(snps.small$M82.DP),] 
snps.small <- snps.small[!is.na(snps.small$Heinz.DP),] 

snps.range <- RangedData(ranges=IRanges(start = snps.small$POS,width = 1),snps.small[,c(-1,-2)],space = snps.small$CHROM)

snps.range

#for some reason subsetByOverlaps() isn't working and I can't figure our why not
snpOverlaps <- as.data.frame(findOverlaps(snps.range,gff.DE.genes))
snpOverlaps$ITAG <- unlist(gff.DE.genes$Alias)[snpOverlaps$subjectHits]
head(snpOverlaps)

snpOverlaps <- cbind(snpOverlaps,snps.range[snpOverlaps$queryHits,])
head(snpOverlaps)

tie1.percent.gene <- tapply(snpOverlaps$tie1.percent.alt,snpOverlaps$ITAG,mean)

DE.data$tie.percent.alt <- NA

DE.data$tie.percent.alt[DE.data$ITAG %in% names(tie1.percent.gene)] <- tie1.percent.gene

DE.data[1:10,]

write.csv(DE.data,file="tie1_DEgenes_plus_snp_percent.csv",row.names=FALSE)



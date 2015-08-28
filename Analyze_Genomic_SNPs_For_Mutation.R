
#get nearby genes

library(rtracklayer)
library(IRanges)
library(stringr)

#unfortunately there are errors in the ITAG2.3 file so I can't use import.gff3()
genes <- read.delim("~/Documents/Lab Notebook support/2013//ITAG2.3_gene_models.gff3",sep="\t",comment.char="#",as.is=T,header=F) #Can be downloaded from SGN
#ftp://ftp.sgn.cornell.edu/genomes/Solanum_lycopersicum/annotation/ITAG2.3_release/ITAG2.3_gene_models.gff3
#for finding candidate causal mutations underlying tie1
#include M82 and Heinz for clean-up

library(ggplot2)
library(Biostrings)
library(plyr)
library(reshape2)

setwd("~/git/Lavelle-tie1--2015/")

data <- read.delim("annotated.tie1_M82_Heinz_fb.filter20.vcf.gz",stringsAsFactors = FALSE,comment.char = "#",header=FALSE)
colnames(data) <- unlist(
  strsplit(
    sub("#","",system("zgrep '^#[A-Z,a-z]' annotated.tie1_M82_Heinz_fb.filter20.vcf.gz",intern = TRUE)),
    split = "\t")
  )

head(data)
tail(data)
names(data)

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

data <- cbind(data,
              split.gt(data$tie1,data$FORMAT,"tie1"),
              split.gt(data$M82,data$FORMAT,"M82"),
              split.gt(data$Heinz,data$FORMAT,"Heinz")
)

head(data[,-8])

#filter data to reduce to relevant SNPs
#filters
#M82 coverage
qplot(data$M82.DP,geom="histogram") + scale_x_log10()
qplot(data$M82.DP,geom="histogram") + xlim(0,100)
qplot(data$M82.DP,geom="histogram") + xlim(0,200)

data.small <- data[data$M82.DP > 4 & data$M82.DP < 200 ,] 

#Heinz coverage
qplot(data$Heinz.DP,geom="histogram") + scale_x_log10()
qplot(data$Heinz.DP,geom="histogram") + xlim(0,200)

data.small <- data.small[data.small$Heinz.DP > 4 & data.small$Heinz.DP < 200 ,] 

data.small <- data.small[grepl("0/0",data.small$Heinz.GT),] #get rid of weird SNPs

data.small <- data.small[data.small$M82.percent.alt> 80 | data.small$M82.percent.alt < 20,] #focus on SNPs where M82 is fixed
#tie1 coverage

qplot(data.small$tie1.DP,geom="histogram") + scale_x_log10()
qplot(data.small$tie1.DP,geom="histogram") + xlim(0,200)


data.small <- data.small[data.small$tie1.DP > 4 & data.small$tie1.DP < 150 ,]

data.small <- data.small[!is.na(data.small$tie1.DP),] 
data.small <- data.small[!is.na(data.small$M82.DP),] 
data.small <- data.small[data.small$tie1.percent.alt > 90,]
data.small$annotation <- ifelse(grepl("ANN=",data.small$INFO),
                                sub("^.+(ANN=.+$)","\\1",data.small$INFO,),
                                NA)

data.small <- data.small[!is.na(data.small$annotation),]

data.small <- data.small[!grepl("intron_variant|synonymous_variant|UTR_variant",data.small$annotation),]

data.small <- data.small[data.small$tie1.GT!=data.small$M82.GT,]


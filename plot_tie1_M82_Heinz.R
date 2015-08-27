#for plotting segregating SNPs of tie1 F2 BSA plants
#include M82 and Heinz for clean-up

library(ggplot2)
library(caTools) #Could also try using IRanges or zoo packages
library(Biostrings)
library(plyr)
library(reshape2)

setwd("~/git/Lavelle-tie1--2015/")

data <- read.delim("tie1_M82_Heinz_fb.filter20.vcf.gz",stringsAsFactors = FALSE,comment.char = "#",header=FALSE)
colnames(data) <- unlist(
  strsplit(
    sub("#","",system("zgrep '^#[A-Z,a-z]' tie1_M82_Heinz_fb.filter20.vcf.gz",intern = TRUE)),
    split = "\t")
  )

head(data)
tail(data)
names(data)

#limit to biallelic events

data <- data[!grepl(",",data$ALT),]

split.gt <- function(gt.data,format,prefix) {
  # data is a character vector with the data to be split,
  # format is a string indicating the labels for each column, 
  # and prefix is a character string to be pre-pended to each new label
  
  gt.data <- sub("^\\.$","NA:NA:NA:NA:NA:NA:NA",gt.data)
  
  tmp <- ldply(strsplit(gt.data,split=":"))
                  
  colnames(tmp) <- paste(prefix,unlist(strsplit(unique(format)[1],split=":")),sep=".")
  tmp[,2:6] <- apply(tmp[,2:6],2,as.numeric)
  
  tmp <- within(tmp,assign(paste(prefix,"percent.alt",sep="."),tmp[,5]/tmp[,2]*100))
  
  tmp[,7] <- sub("NA","NA,NA,NA",tmp[,7])
  
  tmp.gl <- ldply(strsplit(tmp[,7],split=","))
  
  colnames(tmp.gl) <- paste(prefix,"gl",c("ref","het","alt"),sep=".")
  
  tmp.gl <- apply(tmp.gl,2,as.numeric)

  cbind(tmp,tmp.gl)
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

data.small <- data[data$M82.DP > 9 & data$M82.DP < 200 ,] 

#Heinz coverage
qplot(data$Heinz.DP,geom="histogram") + scale_x_log10()
qplot(data$Heinz.DP,geom="histogram") + xlim(0,200)

data.small <- data.small[data.small$Heinz.DP > 9 & data.small$Heinz.DP < 200 ,] 

data.small <- data.small[grepl("1/1",data.small$M82.GT),] #Focus on places where M82 is homozygous different from Heinz
data.small <- data.small[grepl("0/0",data.small$Heinz.GT),] #get rid of weird SNPs

qplot(data.small$M82.percent,geom="histogram")
qplot(data.small$M82.percent,geom="histogram",binwidth=1) + xlim(80,101)


data.small <- data.small[data.small$M82.percent>95,] #focus on SNPs where we are confident that M82 is different from Heinz

#tie1 coverage

qplot(data.small$tie1.DP,geom="histogram") + scale_x_log10()
qplot(data.small$tie1.DP,geom="histogram") + xlim(0,200)


data.small <- data.small[data.small$tie1.DP > 4 & data.small$tie1.DP < 150 ,]

data.small <- data.small[!is.na(data.small$tie1.DP),] 
data.small <- data.small[!is.na(data.small$M82.DP),] 

#checkout the quality

qplot(data.small$QUAL,geom="histogram") 
qplot(data.small$QUAL,geom="histogram") + xlim(0,200)


#instead of using a smoothed line, compute running average.  The reason is that 
#smoothed lines are overly influenced by the SNP density

data.small$tie1.percent.run10 <- caTools::runmean(data.small$tie1.percent.alt,10)
data.small$tie1.percent.run20 <- caTools::runmean(data.small$tie1.percent.alt,20)

data.small <- data.small[data.small$CHROM!="SL2.40ch00",]

##plot some of this 

plotSnp <- function(start=NULL,finish=NULL,data=data.small,title=NULL,alpha=0.5,chrom=NULL) {
  if(! is.null(chrom)) {
    data <- data[data$CHROM==chrom,]
  }
  pl <- ggplot(data=data,mapping=aes(x=POS))
  pl <- pl + geom_point(aes(y=tie1.percent.alt),alpha=alpha)
  if(!is.null(start) & ! is.null(finish)) {
    pl <- pl + xlim(c(start,finish))
    #pl <- pl + scale_size_discrete(range=c(2,6))
  }
  if(!is.null(title)) {
    pl <- pl + ggtitle(title)
  }
  if(is.null(chrom)) {
    pl <- pl + facet_wrap(~ CHROM, ncol = 2)
    #pl <- pl + geom_smooth(aes(y=tie1.percent.alt),color="skyblue",lwd=1) + ylim(0,100)
}# else {
  pl <- pl + geom_line(aes(y=tie1.percent.run20),color="skyblue",lwd=1)
#}
  print(pl)
}

plotSnp(title="tie1 F2 BSA Analysis")
ggsave("tie1_F2_BSA_All_Chroms.pdf",height=11,width=8.5)

plotSnp(title="tie1 F2 BSA Analysis Chrom 2",chrom="SL2.40ch02") #essentially two regions, around 3.4- 3.6 and 3.9
ggsave("tie1_F2_BSA_Ch02.pdf",width=8.5,height=2)

plotSnp(3e07,4.5e07,title="Chrom 2 Enlargement",chrom="SL2.40ch02")
ggsave("tie1_F2_BSA_Ch02_enlargement.pdf",width=8.5,height=2)

plotSnp(3.3e07,3.6e07,title="Chrom 2 Enlargement Region 1",chrom="SL2.40ch02")

plotSnp(3.8e07,4.0e07,title="Chrom 2 Enlargement Region 2",chrom="SL2.40ch02")

#get nearby genes

library(rtracklayer)
library(IRanges)
library(stringr)

#unfortunately there are errors in the ITAG2.3 file so I can't use import.gff3()
genes <- read.delim("~/Documents/Lab Notebook support/2013//ITAG2.3_gene_models.gff3",sep="\t",comment.char="#",as.is=T,header=F) #Can be downloaded from SGN
#ftp://ftp.sgn.cornell.edu/genomes/Solanum_lycopersicum/annotation/ITAG2.3_release/ITAG2.3_gene_models.gff3

head(genes)
colnames(genes) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
head(genes)

genes <- genes[genes$type=="gene",]

genes$locus <- regmatches(genes$attributes,regexpr("Solyc[01][0-9]g[0-9]{6}",genes$attributes))

head(genes)

genes.rd <- RangedData(IRanges(start=genes$start, end=genes$end, names=genes$locus),space=genes$seqid)

head(genes.rd)

#make it a ranged data object

region1.rd <- RangedData(IRanges(start=3.45e07,end=3.55e07,names="tie1Region1"),space="SL2.40ch02")
region1.rd

region2.rd <- RangedData(IRanges(start=3.85e07,end=3.95e07,names="tie1Region2"),space="SL2.40ch02")

overlaps1 <- as.data.frame(as.matrix(findOverlaps(ranges(region1.rd),ranges(genes.rd))))

head(overlaps1)

overlaps2 <- as.data.frame(as.matrix(findOverlaps(ranges(region2.rd),ranges(genes.rd))))

overlaps1$locus <- genes$locus[overlaps1$subjectHits]

overlaps2$locus <- genes$locus[overlaps2$subjectHits]

write.csv(overlaps1,"tie1Overlaps1.csv")
write.csv(overlaps2,"tie1Overlaps2.csv")

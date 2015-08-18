#for plotting cleaned up SNPs of tie1 chrom2 region

library(ggplot2)
library(caTools) #Could also try using IRanges or zoo packages
library(Biostrings)

setwd("~/git/Lavelle-tie1--2015/")

data <- read.delim("PE1PE2_paired_M82G.realigned.all.varscan.gz",sep="\t",as.is=T)

head(data)

names(data)

data <- data[,c(1:4,11)]

split <- matrix(unlist(strsplit(data$Cons.Cov.Reads1.Reads2.Freq.P.value.1,split=" ")),ncol=2,byrow=T)

colnames(split) <- c("tie1F2","M82")

tie1F2 <- matrix(unlist(strsplit(split[,"tie1F2"],split=":")),ncol=6,byrow=T)
head(tie1F2)
colnames(tie1F2) <- c("tie1.cons","tie1.cov","tie1.ref","tie1.alt","tie1.percent","tie1.P")

M82 <- matrix(unlist(strsplit(split[,"M82"],split=":")),ncol=6,byrow=T)
head(M82)
colnames(M82) <- c("M82.cons","M82.cov","M82.ref","M82.alt","M82.percent","M82.P")

data <- cbind(data,tie1F2,M82)
head(data)
data[data=="-"] <- NA
head(data)
summary(data)

data$tie1.cov <- as.numeric(as.character(data$tie1.cov))
data$tie1.ref <- as.numeric(as.character(data$tie1.ref))
data$tie1.alt <- as.numeric(as.character(data$tie1.alt))
data$tie1.P <- as.numeric(as.character(data$tie1.P))
data$tie1.percent <- as.numeric(sub("%","",data$tie1.percent))
data$tie1.percent2 <- data$tie1.alt/data$tie1.cov*100
data$tie1.cons.bases <- IUPAC_CODE_MAP[as.character(data$tie1.cons)]
data$tie1.cons.bases[is.na(data$tie1.cons.bases)] <- data$tie1.cons[is.na(data$tie1.cons.bases)] #if there is no IUPAC code (ie deletion) put in the call
data$tie1.het.Pvalue <- pbinom(data$tie1.alt,data$tie1.ref + data$tie1.alt,.5,lower=F)
data$tie1.het.Pvalue.adj <- p.adjust(data$tie1.het.Pvalue,method="fdr")

data$M82.cov <- as.numeric(as.character(data$M82.cov))
data$M82.ref <- as.numeric(as.character(data$M82.ref))
data$M82.alt <- as.numeric(as.character(data$M82.alt))
data$M82.percent <- as.numeric(sub("%","",data$M82.percent))
data$M82.percent2 <- data$M82.alt/(data$M82.cov)*100
data$M82.cons.bases <- IUPAC_CODE_MAP[as.character(data$M82.cons)]
data$M82.cons.bases[is.na(data$M82.cons.bases)] <- data$M82.cons[is.na(data$M82.cons.bases)] #if there is no IUPAC code (ie deletion) put in the call

head(data)
summary(data)

#filter data to reduce to relevant SNPs
#filters
#M82 coverage
qplot(data$M82.cov,geom="histogram") + scale_x_log10()
qplot(data$M82.cov,geom="histogram") + xlim(0,100)
qplot(data$M82.cov,geom="histogram") + xlim(0,200)

data.small <- data[data$M82.cov > 9 & data$M82.cov < 200 ,] 
qplot(data.small$M82.cov,geom="histogram") + xlim(0,200)


data.small <- data.small[grepl("[GATC]",data.small$M82.cons),] #get rid of M82 hets

qplot(data.small$M82.percent2,geom="histogram")
qplot(data.small$M82.percent2,geom="histogram",binwidth=1) + xlim(80,101)


data.small <- data.small[data.small$M82.percent2>98,] #focus on SNPs where we are confident that M82 is different from Heinz

#tie1 coverage

qplot(data.small$tie1.cov,geom="histogram") + scale_x_log10()
qplot(data.small$tie1.cov,geom="histogram",binwidth=1) + xlim(0,25)
qplot(data.small$tie1.cov,geom="histogram") + xlim(0,100)


data.small <- data.small[data.small$tie1.cov > 4 & data.small$tie1.cov < 50 ,]

data.small <- data.small[!is.na(data.small$tie1.ref),] 
data.small <- data.small[!is.na(data.small$M82.ref),] 

#instead of using a smoothed line, compute running average.  The reason is that 
#smoothed lines are overly infleunced by the SNP density

data.small$tie1.percent.run <- caTools::runmean(data.small$tie1.percent,20)
#data.small$tie1.het.Pvalue.run <- caTools::runmean(data.small$tie1.het.Pvalue,5)

data.small <- data.small[data.small$Chrom!="SL2.40ch00",]
##plot some of this 

plotSnp <- function(start=NULL,finish=NULL,data=data.small,title=NULL,alpha=0.5,chrom=NULL) {
  if(! is.null(chrom)) {
    data <- data[data$Chrom==chrom,]
  }
  pl <- ggplot(data=data,mapping=aes(x=Position))
  pl <- pl + geom_point(aes(y=tie1.percent2),alpha=alpha)
  if(!is.null(start) & ! is.null(finish)) {
    pl <- pl + xlim(c(start,finish))
    #pl <- pl + scale_size_discrete(range=c(2,6))
  }
  if(!is.null(title)) {
    pl <- pl + ggtitle(title)
  }
  if(is.null(chrom)) {
    pl <- pl + facet_wrap(~ Chrom, ncol = 2)
    #pl <- pl + geom_smooth(aes(y=tie1.percent2),color="skyblue",lwd=1) + ylim(0,100)
}# else {
  pl <- pl + geom_line(aes(y=tie1.percent.run),color="skyblue",lwd=1)
#}
  print(pl)
}
  
plotSnp(title="tie1 F2 BSA Analysis")
ggsave("tie1_F2_BSA_All_Chroms.pdf",height=11,width=8.5)

plotSnp(title="tie1 F2 BSA Analysis Chrom 2",chrom="SL2.40ch02") #essentially two regions, around 3.4- 3.6 and 3.9
ggsave("tie1_F2_BSA_Ch02.pdf",width=4,height=4)

plotSnp(2.8e07,4.2e07,title="Chrom 2 Enlargement",chrom="SL2.40ch02")
ggsave("tie1_F2_BSA_Ch02.pdf",width=4,height=4)

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

region1.rd <- RangedData(IRanges(start=3.3e07,end=3.6e07,names="tie1Region1"),space="SL2.40ch02")
region1.rd

region2.rd <- RangedData(IRanges(start=3.9e07,end=4.1e07,names="tie1Region2"),space="SL2.40ch02")

overlaps1 <- as.data.frame(as.matrix(findOverlaps(ranges(region1.rd),ranges(genes.rd))))

head(overlaps1)

overlaps2 <- as.data.frame(as.matrix(findOverlaps(ranges(region2.rd),ranges(genes.rd))))

overlaps1$locus <- genes$locus[overlaps1$subjectHits]

overlaps2$locus <- genes$locus[overlaps2$subjectHits]

write.csv(overlaps1,"tie1Overlaps1.csv")
write.csv(overlaps2,"tie1Overlaps2.csv")

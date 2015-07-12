setwd("~/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/")
library(plyr)
library(reshape2)
library(ggplot2)

files <- dir(pattern="(prime)|(intron)",recursive = TRUE)

files

for (f in files) {
  shortf <- sub("\\.tsv","",basename(f)) 
  assign(shortf,cbind(read.delim(f,header=F,col.names=c("chr","pos","depth")),file=shortf))
}

counts <- ldply(lapply(ls(pattern = "(prime)|(intron)"),function(x) get(x)))

counts$gt <- unlist(strsplit(as.character(counts$file),split="_"))[c(T,F)]
counts$position <- unlist(strsplit(as.character(counts$file),split="_"))[c(F,T)]

counts.melt <- melt(tapply(counts$depth,list(counts$gt,counts$position),mean),id.vars = 0)
names(counts.melt) <- c("genotype","position","depth")
counts.melt$position <- factor(counts.melt$position,
                               levels = c("5prime","intron","3prime"),
                               labels=c("5-prime","intron","3-prime"))

pl <- ggplot(counts.melt,aes(x=position,y=depth,fill=genotype))
pl <- pl + geom_bar(stat="Identity",position="dodge")
pl + ggtitle("Read depth in M82 and tie-1")

ggsave("tie-1_M82_depth.pdf",height=6,width=6)

## another figure based on junction spanning

span <- data.frame(
  genotype=rep(c("M82","tie-1"),each=4),
  location=rep(c("5-prime junction","3-prime junction"),each = 2,length.out = 8),
  type=rep(c("spliced","unspliced"),each=1,length.out=8),
  counts=c(25,1,25,2,0,33,0,15))

span

span$location <- relevel(span$location,ref="5-prime junction")

pl2 <- ggplot(span,aes(x=genotype,y=counts,fill=type))
pl2 <- pl2 + facet_wrap( ~ location,ncol=2)
pl2 <- pl2 + geom_bar(stat="Identity")
pl2 + ggtitle("Number of reads spliced or unspliced at splice site")

ggsave("tie-1_M82_junctions.pdf",height=6,width=6)

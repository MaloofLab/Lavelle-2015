# Script to filter annotated SNPs to find candidates for causing tie-1
# Julin Maloof

library(ggplot2)
library(plyr)
library(reshape2)

setwd("~/git/Lavelle-tie1--2015/")

data <- read.delim("annotated.tie1_M82_Heinz_fb.filter20.vcf",stringsAsFactors = FALSE,comment.char = "#",header=FALSE)

#extract the column headers from the VCF file
colnames(data) <- unlist(
  strsplit(
    sub("#","",system("zgrep '^#[A-Z,a-z]' annotated.tie1_M82_Heinz_fb.filter20.vcf",intern = TRUE)),
    split = "\t")
  )

head(data)
tail(data)
names(data)

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

data <- cbind(data,
              split.gt(data$tie1,data$FORMAT,"tie1"),
              split.gt(data$M82,data$FORMAT,"M82"),
              split.gt(data$Heinz,data$FORMAT,"Heinz")
)

head(data[,-8])

#filter data to focus on relevant SNPs

#filter for M82 coverage.  First plot some histograms to see what is reasonable.
qplot(data$M82.DP,geom="histogram") + scale_x_log10()
qplot(data$M82.DP,geom="histogram") + xlim(0,100)
qplot(data$M82.DP,geom="histogram") + xlim(0,200)

data.small <- data[data$M82.DP > 4 & data$M82.DP < 200 ,] 

#filter for Heinz coverage
qplot(data$Heinz.DP,geom="histogram") + scale_x_log10()
qplot(data$Heinz.DP,geom="histogram") + xlim(0,200)

data.small <- data.small[data.small$Heinz.DP > 4 & data.small$Heinz.DP < 200 ,] 



#filter for tie1 coverage

qplot(data.small$tie1.DP,geom="histogram") + scale_x_log10()
qplot(data.small$tie1.DP,geom="histogram") + xlim(0,200)

data.small <- data.small[data.small$tie1.DP > 4 & data.small$tie1.DP < 150 ,]

# additional filters

#get rid of weird SNPs
data.small <- data.small[grepl("0/0",data.small$Heinz.GT),]

#focus on SNPs where M82 is fixed
data.small <- data.small[data.small$M82.percent.alt> 80 | data.small$M82.percent.alt < 20,]

#only look at SNPs where we have information for all tie1 and M82 genotypes
data.small <- data.small[!is.na(data.small$tie1.DP),] 
data.small <- data.small[!is.na(data.small$M82.DP),] 

#reduce to near homozygous tie-1 SNPs
data.small <- data.small[data.small$tie1.percent.alt > 90,]

#create a separate column with the annotation information
data.small$annotation <- ifelse(grepl("ANN=",data.small$INFO),
                                sub("^.+(ANN=.+$)","\\1",data.small$INFO,),
                                NA)

#remove non-annotated SNPs
data.small <- data.small[!is.na(data.small$annotation),]

#remove SNPs that don't change coding sequence
data.small <- data.small[!grepl("MODIFIER",data.small$annotation),]

#only keep SNPs where M82 and tie-1 have a large difference in the percentage of the alternate allele
data.small2 <- data.small[abs(data.small$M82.percent.alt-data.small$tie1.percent.alt) > 60,] 

write.csv(data.small2,"tie1_candidate_SNPs_from_gDNA.csv",row.names = FALSE)
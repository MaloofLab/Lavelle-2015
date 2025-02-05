#Looking for interesting SNPs in RNAseq data
#Julin Maloof
#December 2, 2013

data <- read.delim("anotated_M82_tie1.varscanNEW.vcf",comment.char="#",sep="\t",header=F,as.is=T)
head(data)
colnames(data) <- unlist(strsplit("CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  M82 tie1",split=" +"))
head(data,20)

data$tie1 <- sub("(\\./\\.:\\.:[0-9]+$)","\\1::::::::::::",data$tie1)

tie1 <- as.data.frame(matrix(unlist(strsplit(data[,"tie1"],split=":")),ncol=14,byrow=T))
head(tie1)
colnames(tie1) <- paste("tie1",c("genotype","genotypeQuality","raw","cov",
                                 "ref","alt","percent","pval","r.qual","alt.qual",
                                 "ref.for","ref.rev","alt.for","alt.rev"),sep=".")
head(tie1)

tie1$tie1.percent <- sub("%","",tie1$tie1.percent)

tie1[,-1] <- apply(tie1[,-1],2,as.numeric)

data$M82 <- sub("(\\./\\.:\\.:[0-9]+$)","\\1::::::::::::",data$M82)

M82 <- as.data.frame(matrix(unlist(strsplit(data[,"M82"],split=":")),ncol=14,byrow=T))
head(M82)
colnames(M82) <- paste("M82",c("genotype","genotypeQuality","raw","cov",
                                 "ref","alt","percent","pval","r.qual","alt.qual",
                                 "ref.for","ref.rev","alt.for","alt.rev"),sep=".")
head(M82)

M82$M82.percent <- sub("%","",M82$M82.percent)

M82[,-1] <- apply(M82[,-1],2,as.numeric)

summary(M82)

data.large <- cbind(data,tie1[,1:8],M82[,1:8])

data.small <- data.large[data.large$M82.genotype!="./." & data.large$tie1.genotype!="./.",]
head(data.small,30)
summary(data.small)
data.small <- data.small[as.character(data.small$M82.genotype) != as.character(data.small$tie1.genotype),]
head(data.small,30)
data.small <- data.small[data.small$tie1.genotype=="1/1",]
data.small <- data.small[data.small$tie1.percent>90,]

data.small <- data.small[grep("EFF",data.small$INFO),]
data.small <- data.small[!grepl("=SYNONYMOUS_CODING",data.small$INFO),]

data.small <- data.small[data.small$M82.percent < 20 | data.small$M82.percent > 80,] # limit to nonsegregating M82 sites

head(data.small)

data.small$ITAG <- regmatches(data.small$INFO,regexpr("Solyc[0-9]{2}g[0-9]{6}",data.small$INFO))

head(data.small)

HRD <- read.delim("/Users/jmaloof/Documents/Lab Notebook support/2011/Tomato Annotation and GO/ITAG2.3 annotation/ITAG2.3_HRD.tsv",
                  sep="\t")

head(HRD)

HRD$ITAG <- substr(HRD$ITAG,1,14)

data.small.annotated <- merge(data.small,HRD,by="ITAG",all.x=T)

write.csv(data.small.annotated,"tie1Candidates_RNAseq.csv")


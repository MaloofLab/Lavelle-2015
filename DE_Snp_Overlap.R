DE.data <- read.csv("~/git/Lavelle-tie1--2015/tie1_DEgenes.csv")

library(rtracklayer)
library(genomeIntervals)

readGff3("~/Sequences/ref_genomes/tomato/ITAG2.3_Chromo2.4/ITAG2.3_gene_models.gff3")

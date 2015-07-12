Working in `newTophat020114` directory

Mike suggested

* number of reads spanning the splice junction or spliced in M82 and tie-1
* coverage before, in, and after splice junction in M82 and tie-1

For the first of these I am just going to hand count in IGV

| location     | M82 splice | M82 span | tie-1 splice | tie-1 span |
|---------------|----------------|---------------|----------------|---------------|
| 5' junction | 25              | 1               | 0                | 33             | 
| 3' junction | 25              | 2               | 0                | 15             |

For the second use samtools to calculate coverage:

    cd "/Users/jmaloof/Documents/Lab Notebook support/2013/tie1_phyE/RNAseq/newTophat020114/tophat_tie1"
    samtools depth -r SL2.40ch02:39009924-39009943 accepted_hits.bam > tie1_5prime.tsv
    samtools depth -r SL2.40ch02:39009944-39010179 accepted_hits.bam > tie1_intron.tsv
    samtools depth -r SL2.40ch02:39010180-39010199 accepted_hits.bam > tie1_3prime.tsv
    samtools idxstats accepted_hits.bam > tie1_hits.tsv
    
    cd ../tophat_M82
    samtools depth -r SL2.40ch02:39009924-39009943 accepted_hits.bam > M82_5prime.tsv
    samtools depth -r SL2.40ch02:39009944-39010179 accepted_hits.bam > M82_intron.tsv
    samtools depth -r SL2.40ch02:39010180-39010199 accepted_hits.bam > M82_3prime.tsv
    samtools idxstats accepted_hits.bam > M82_hits.tsv
    
And then normalize by library size?

In R: <file:///Users/jmaloof/Documents/Lab%20Notebook%20support/2013/tie1_phyE/RNAseq/newTophat020114/Figures_for_paper.R>

![](<file:///Users/jmaloof/Documents/Lab%20Notebook%20support/2013/tie1_phyE/RNAseq/newTophat020114/tie-1_M82_depth.pdf>)
![](<file:///Users/jmaloof/Documents/Lab%20Notebook%20support/2013/tie1_phyE/RNAseq/newTophat020114/tie-1_M82_junctions.pdf>)






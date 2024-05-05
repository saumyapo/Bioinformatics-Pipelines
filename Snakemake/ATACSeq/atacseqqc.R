#!/usr/bin/env Rscript

library(ATACseqQC)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(BSgenome.Hsapiens.UCSC.hg38)


## Getting the BAM file name and sample ID
args<- commandArgs(trailingOnly=TRUE)

bamfile <- args[1]
output_path <- args[2]

package = "ATACseqQC"
bamFileLabels <- gsub("_sorted_shifted.bam", "", basename(bamfile))

## Plotting size distribution of fragments (Figure 1G)
jpeg(file.path(output_path, paste0(bamFileLabels, ".fragment.size.distribution.jpeg"))) 

fragSize <- fragSizeDist(bamfile, bamFileLabels)

dev.off()

#!/usr/bin/env Rscript

library(RColorBrewer)
myPalette <- brewer.pal(8, "Set2")

args<- commandArgs(trailingOnly=TRUE)

annotated <- args[1]
output_path <- args[2]

anno = read.delim(annotated, sep="\t", header=T, quote="")
FileLabel <- gsub(".txt", "", basename(annotated))


pietable = table(unlist(lapply(strsplit(as.character(anno$Annotation), " \\("),"[[",1))) 

newtable = table(c("3' UTR","5' UTR","exon","Intergenic","intron","promoter-TSS","TTS"))
newtable[names(newtable)] = 0 #reset everything to 0
newtable[names(pietable)] = pietable

names(newtable) = paste(names(newtable), "(", round(newtable/sum(newtable)*100), "%, ", newtable, ")", sep="")

jpeg(file.path(output_path, paste0(FileLabel, "_piechart.jpeg")))
pie(newtable, main="Distribution of HOMER annotation", col=myPalette)

dev.off()

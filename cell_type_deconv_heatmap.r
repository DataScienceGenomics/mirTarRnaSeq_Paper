#Paper correlationHeatmapEBV miRNA
library(pheatmap)
library(dplyr)
library(mirTarRnaSeq)
library(readxl)

##CellTypeDeconv

CellTypeDecon<-read.table("~/Desktop/CellTypeDeconv_GEDIT.txt", as.is = TRUE, header = T, row.names = 1)
summary(CellTypeDecon)#What is max-what mean and what is min
breaks<-seq(0,1,length.out=2001)#One element longer than the color vector
col_heat<-colorRampPalette(c("#FFFBF3","red","purple","blue"))(2000)
pheatmap(t(log2(CellTypeDecon+1)),breaks=breaks,col=col_heat,
         fontsize_col=5,fontsize_row=5,fontsize = 6,
         cellwidth=10,cellheight=10)

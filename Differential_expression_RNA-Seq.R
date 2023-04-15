## Author - John Philip George
## Date - 09/27/2023
## Differential gene expression analysis


setwd("PATH")
library (DESeq2)
library(ggplot2)

#1. preparing count data
file <- read.csv("COVID19_exp_count_092222.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = TRUE)
head(file)
countdata <- as.matrix(file)
countdata
condition <- factor(c(rep("Normal",77), rep("COVID-19",277)))
condition
# assignement of samples to Covid and normal
coldata <- data.frame(row.names = colnames(countdata), condition)

# Mapping with matrix
ddsFull <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design =~condition)
ddsFull

dds = DESeq(ddsFull)
res <- results(dds)
summary (res)
head(res)
# Visualization
res_ordered <- res[order(res$padj),]
res_x <- as.data.frame(res_ordered)
res_x
res_x$diffexpressed <- "NO"

# if log2FoldChange > 0.6 & Pvalue <0.05, set as 'UP'
res_x$diffexpressed[res_x$log2FoldChange > 0.6 & res_x$pvalue < 0.05] <-"Upregulated"
res_x$diffexpressed[res_x$log2FoldChange < -0.6 & res_x$pvalue < 0.05] <-"Downregulated"


res_x <- cbind(rownames(res_x), data.frame(res_x, row.names = NULL))
res_x
colnames(res_x)[1] <- "genes"
res_x

#Volcano plot``
plot1 <- ggplot(res_x, aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(colour = diffexpressed), size = 1.5, alpha = 2)+ 
  scale_color_manual(values = c("darkblue", "grey", "red"))+
  theme_update(legend.position='bottom')+
  geom_vline(xintercept = c(-1,1), lty=4, col="black", lwd=0.8)+
  geom_hline(yintercept = 1.301, lty=4, col="black", lwd=0.8)

plot1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))

# Save outputs
a = res_x[which(res_x$diffexpressed=="Upregulated"),]
b = res_x[which(res_x$diffexpressed=="Downregulated"),]
write.table (a, file = "Upregulated_covid19_092222.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)
write.table (b, file = "Downregulated_covid19_092222.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)

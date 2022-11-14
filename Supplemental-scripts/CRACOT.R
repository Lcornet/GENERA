library(dplyr)
library(ggplot2)
library(ggpubr)
#Conta part
contam <- read.delim("conta.txt")
contam$level <- factor(contam$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
conta_median <- contam %>% group_by(level) %>% summarize(Avg = median(percent))
#CheckM
results <- read.delim("resultsCheckM.txt")
results$level <- factor(results$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
subresults <- results[ which(results$category=='CK-conta'), ]
results_median <- subresults %>% group_by(level) %>% summarize(Avg = median(percent))
Lresults <- subresults[ which(subresults$level=='phylum'), ]
Lcontam <- contam[ which(contam$level=='phylum'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Pcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='class'), ]
Lcontam <- contam[ which(contam$level=='class'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ccor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='order'), ]
Lcontam <- contam[ which(contam$level=='order'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ocor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='family'), ]
Lcontam <- contam[ which(contam$level=='family'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Fcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='genus'), ]
Lcontam <- contam[ which(contam$level=='genus'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Gcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='species'), ]
Lcontam <- contam[ which(contam$level=='species'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Scor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lcor <- c(Pcor, Ccor, Ocor, Fcor, Gcor, Scor)
CK1 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -1, size = 2, colour = "red") + theme(legend.position="none") + ylim(-1,NA)
#CK1 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue")
#BUSCO
results <- read.delim("resultsBUSCO.txt")
results$level <- factor(results$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
subresults <- results[ which(results$category=='BUSCO-dup'), ]
results_median <- subresults %>% group_by(level) %>% summarize(Avg = median(percent))
Lresults <- subresults[ which(subresults$level=='phylum'), ]
Lcontam <- contam[ which(contam$level=='phylum'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Pcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='class'), ]
Lcontam <- contam[ which(contam$level=='class'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ccor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='order'), ]
Lcontam <- contam[ which(contam$level=='order'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ocor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='family'), ]
Lcontam <- contam[ which(contam$level=='family'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Fcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='genus'), ]
Lcontam <- contam[ which(contam$level=='genus'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Gcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='species'), ]
Lcontam <- contam[ which(contam$level=='species'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Scor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lcor <- c(Pcor, Ccor, Ocor, Fcor, Gcor, Scor)
BUSCO <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -1, size = 2, colour = "red") +  annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -1, size = 2, colour = "red") + theme(legend.position="none") + ylim(-1,NA)
#BUSCO <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue")
#GUNC
results <- read.delim("resultsGUNC.txt")
results$level <- factor(results$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
subresults <- results[ which(results$category=='GUNC-conta'), ]
results_median <- subresults %>% group_by(level) %>% summarize(Avg = median(percent))
Lresults <- subresults[ which(subresults$level=='phylum'), ]
Lcontam <- contam[ which(contam$level=='phylum'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Pcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='class'), ]
Lcontam <- contam[ which(contam$level=='class'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ccor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='order'), ]
Lcontam <- contam[ which(contam$level=='order'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ocor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='family'), ]
Lcontam <- contam[ which(contam$level=='family'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Fcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='genus'), ]
Lcontam <- contam[ which(contam$level=='genus'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Gcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='species'), ]
Lcontam <- contam[ which(contam$level=='species'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Scor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lcor <- c(Pcor, Ccor, Ocor, Fcor, Gcor, Scor)
GUNC <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -1, size = 2, colour = "red") +  annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -1, size = 2, colour = "red") + theme(legend.position="none") + ylim(-1,NA)
#GUNC <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue")
#Physeter
results <- read.delim("resultsPhyseter.txt")
results$level <- factor(results$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
subresults <- results[ which(results$category=='Physeter-conta'), ]
results_median <- subresults %>% group_by(level) %>% summarize(Avg = median(percent))
Lresults <- subresults[ which(subresults$level=='phylum'), ]
Lcontam <- contam[ which(contam$level=='phylum'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Pcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='class'), ]
Lcontam <- contam[ which(contam$level=='class'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ccor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='order'), ]
Lcontam <- contam[ which(contam$level=='order'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ocor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='family'), ]
Lcontam <- contam[ which(contam$level=='family'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Fcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='genus'), ]
Lcontam <- contam[ which(contam$level=='genus'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Gcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='species'), ]
Lcontam <- contam[ which(contam$level=='species'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Scor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lcor <- c(Pcor, Ccor, Ocor, Fcor, Gcor, Scor)
Physeter <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -1, size = 2, colour = "red") +  annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -1, size = 2, colour = "red") + theme(legend.position="none") + ylim(-1,NA)
#Physeter <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue")
#Kraken2
results <- read.delim("resultsKraken2.txt")
results$level <- factor(results$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
subresults <- results[ which(results$category=='Kraken2-conta'), ]
results_median <- subresults %>% group_by(level) %>% summarize(Avg = median(percent))
Lresults <- subresults[ which(subresults$level=='phylum'), ]
Lcontam <- contam[ which(contam$level=='phylum'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Pcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='class'), ]
Lcontam <- contam[ which(contam$level=='class'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ccor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='order'), ]
Lcontam <- contam[ which(contam$level=='order'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ocor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='family'), ]
Lcontam <- contam[ which(contam$level=='family'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Fcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='genus'), ]
Lcontam <- contam[ which(contam$level=='genus'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Gcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='species'), ]
Lcontam <- contam[ which(contam$level=='species'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Scor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lcor <- c(Pcor, Ccor, Ocor, Fcor, Gcor, Scor)
Kraken2 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -1, size = 2, colour = "red") +  annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -1, size = 2, colour = "red") + theme(legend.position="none") + ylim(-1,NA)
#Kraken2 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue")
#CK2
results <- read.delim("resultsCheckM2.txt")
results$level <- factor(results$level, levels=c("phylum", "class", "order", "family", "genus", "species"))
subresults <- results[ which(results$category=='CheckM2-conta'), ]
results_median <- subresults %>% group_by(level) %>% summarize(Avg = median(percent))
Lresults <- subresults[ which(subresults$level=='phylum'), ]
Lcontam <- contam[ which(contam$level=='phylum'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Pcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='class'), ]
Lcontam <- contam[ which(contam$level=='class'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ccor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='order'), ]
Lcontam <- contam[ which(contam$level=='order'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Ocor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='family'), ]
Lcontam <- contam[ which(contam$level=='family'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Fcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='genus'), ]
Lcontam <- contam[ which(contam$level=='genus'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Gcor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lresults <- subresults[ which(subresults$level=='species'), ]
Lcontam <- contam[ which(contam$level=='species'), ]
df_merge <- merge(Lresults, Lcontam,by="ID")
Scor <- cor(df_merge$percent.x, df_merge$percent.y, method = "spearman")
Lcor <- c(Pcor, Ccor, Ocor, Fcor, Gcor, Scor)
CK2 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -1, size = 2, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -1, size = 2, colour = "red") +  annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -1, size = 2, colour = "red") + theme(legend.position="none") + ylim(-1,NA)
#CK2 labl font size
# CK2 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + theme(text = element_text(size = 22)) + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue", size = 0.2) + annotate("text", label = sprintf("%0.3f", round(Lcor[1], digits = 3)), x =1, y = -2, size = 4, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[2], digits = 3)), x =2, y = -2, size = 4, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[3], digits = 3)), x =3, y = -2, size = 4, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[4], digits = 3)), x =4, y = -2, size = 4, colour = "red") + annotate("text", label = sprintf("%0.3f", round(Lcor[5], digits = 3)), x =5, y = -2, size = 4, colour = "red") +  annotate("text", label = sprintf("%0.3f", round(Lcor[6], digits = 3)), x =6, y = -2, size = 4, colour = "red")
#CK2 <- ggplot(results, aes(x=level, y=percent, fill = category)) + geom_violin() + geom_line(data = conta_median, mapping = aes(x = level, y = Avg, fill=NULL), group = 1, color="blue")
#Combine plot
A <- ggarrange(CK1, BUSCO, GUNC, Physeter, Kraken2, CK2, labels = c("CheckM", "BUSCO", "GUNC", "Physeter", "Kraken2", "CheckM2"), font.label = list(size = 10, color = "black"), align = "hv", heights=20, hjust=-0.5, vjust = 0.1, ncol = 1, nrow = 6) + theme(plot.margin = margin(1,1,1,1, "cm"))
#Print pdf
#pdf("results.pdf",  width = 4, height = 8)
pdf("results.pdf")
A
dev.off()

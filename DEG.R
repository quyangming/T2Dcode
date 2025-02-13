
data <- read.table("data.csv",header=T, sep=",",row.names = 1)
zero_counts <- apply(data == 0, 1, sum)   
data <- data[-which(zero_counts>ncol(data)*0.25),]
coldata <- read.table("coldata.csv",header=T, sep=",",row.names = 1)
library(DESeq2)
data<- round(as.matrix(data)) 
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition)
head(dds)
dds <- DESeq(dds)	
res <- results(dds)


result <- as.data.frame(results(dds)) 
result[which(result$log2FoldChange>0),'up_down']<-'up'
result[which(result$log2FoldChange<0),'up_down']<-'down'
head(result)  

DGE <-subset(result,pvalue < 0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))
DGE <- DGE[order(DGE$log2FoldChange),]
DGE_matrix <- data[rownames(DGE),]

write.csv(result,'DESeq2_DEGs_result.csv',row.names = TRUE)
write.csv(DGE,'DESeq2_DEGs.csv',row.names = TRUE)
write.csv(DGE_matrix,'DESeq2_DEGs_matrix.csv',row.names = TRUE)

library(ggplot2)
library(ggrepel)

DEG = result
head(DEG)
DEG$logFC=DEG$log2FoldChange
DEG$PValue=DEG$pvalue
logFC_cutoff= 0.5 
p_cutoff= 0.05 
DEG$change = as.factor(ifelse(DEG$PValue < p_cutoff & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT') )

p=ggplot(data = DEG, aes(x = logFC, y = -log10(PValue),color=change )) +    
geom_hline(yintercept=-log10(p_cutoff),lty=4,lwd=0.6,alpha=0.8) + 
geom_vline(xintercept = c(logFC_cutoff,-logFC_cutoff),lty=4,lwd=0.6,alpha=0.8)+ 
geom_point( alpha=0.8)+      
scale_color_manual(name = c("change"), values = c("#FE180A", "#263394", "#AAAAAA"), limits = c("UP", "DOWN", "NOT"))+
  ylim(c(0,8))+
  guides(size = F)+
  labs(x="log2 (Fold Change)",
       y="-log10 (P Value)",
       title = "DESeq2"          
  )+
  theme_bw(base_size = 12) + 
  theme(plot.title = element_text(size=15,hjust = 0.5),  
        legend.position = c(0.9,0.9),     
        legend.title = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))

ggsave("DESeq2_Volcano_plot.png", p, width = 4, height = 5)



group_list <-read.table("coldata.csv", header=T,sep=",",row.names = 1)
colnames(group_list)[1] <- "group"
need_DEG =result
head(need_DEG)
need_DEG$logFC=need_DEG$log2FoldChange
need_DEG$PValue=need_DEG$pvalue
logFC_cutoff= 0.5 
p_cutoff= 0.05 
need_DEG$change = as.factor(ifelse(need_DEG$PValue < p_cutoff & abs(need_DEG$logFC) > logFC_cutoff,
                                   ifelse(need_DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT') )

down_gene=need_DEG[need_DEG$change=="DOWN",]
down_gene=down_gene[order(down_gene$PValue,decreasing = F),]
up_gene=need_DEG[need_DEG$change=="UP",]
up_gene=up_gene[order(up_gene$PValue,decreasing = F),]


library(pheatmap)
choose_gene<- c(row.names(down_gene[order(down_gene$PValue)[1:200],]),row.names(up_gene[order(up_gene$PValue)[1:200],]) )
choose_matrix=expr[choose_gene,]
choose_matrix = na.omit(choose_matrix)
choose1=choose_matrix[,group_list=="control"]  
choose2=choose_matrix[,group_list=="T2D"]      
choose_matrix =cbind(choose1,choose2)
choose_matrix[1:4,1:4]
n=t(scale(t(log2(choose_matrix+1))))  
n[n>2]=2 
n[n< -2]= -2 
n[1:4,1:4]

p=pheatmap(n,
           show_colnames = F,  
           show_rownames = F, 
           annotation_col = group_list,
           cluster_cols = T 
)

ggsave("DESeq2_pheatmap.png", p, width = 5, height = 4)






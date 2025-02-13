
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
samples<- list.files("seurat1/")
samples
seurat_list<- list()
for(sample in samples) {
  data.path<- paste0("seurat1/", sample)
  seurat_data<- Read10X(data.dir = data.path)
  seurat_obj<- CreateSeuratObject(counts = seurat_data,
                                  project = sample, min.features= 200, min.cells= 3)
  seurat_list<- append(seurat_list, seurat_obj)
}

seurat_combined<- merge(seurat_list[[1]], y= seurat_list[-1], add.cell.ids= samples)

pbmc= JoinLayers(seurat_combined)

abc789=pbmc@meta.data  
head(abc789)
meta=abc789

a<- c("GSM6846488_MS17001","GSM6846490_MS17003","GSM6846491_MS17004","GSM6846492_MS17005","GSM6846493_MS17006",
        "GSM6846497_MS18004","GSM6846499_MS18006","GSM6846501_MS19002","GSM6846504_MS19005","GSM6846508_MS19009",
        "GSM6846518_MS19019","GSM6846526_MS19035","GSM6846530_MS19043","GSM6846532_MS20003","GSM6846539_MS21001",
        "GSM6846537_MS20087","GSM6846541_MS21009")
meta$group<-ifelse(meta$orig.ident %in% a,'T2D','contrl')
head(meta)

pbmc<- AddMetaData(object = pbmc, metadata = meta, col.name = "group")
table(meta$orig.ident,meta$group)

dim(pbmc)  

table(str_split(colnames(pbmc),'_',simplify = T)[,2])
pbmc<- AddMetaData(object = pbmc,
                   metadata = str_split(colnames(pbmc),'_',simplify = T)[,2],
                   col.name = "orig.ident")

saveRDS(pbmc,"0pbmc.rds")

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
pbmc= readRDS("0pbmc.rds")
pbmc[["percent.mt"]]<- PercentageFeatureSet(pbmc,pattern = "^MT-")
HB.genes<- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes<- CaseMatch(HB.genes, row.names(pbmc))
pbmc[["percent.HB"]]<- PercentageFeatureSet(pbmc, features= HB.genes)
FeatureScatter(pbmc, "nCount_RNA", "percent.mt", group.by = "orig.ident") 
FeatureScatter(pbmc,"nCount_RNA","nFeature_RNA", group.by = "orig.ident")  
theme.set2= theme(axis.title.x= element_blank())
plot.featrures= c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")
group= "orig.ident"
plots= list()
for(i in c(1:length(plot.featrures))) {
  plots[[i]]= VlnPlot(pbmc, group.by = group, pt.size = 0,
                     features = plot.featrures[i]) + theme.set2 + NoLegend() }
violin <- wrap_plots(plots=plots, nrow = 2)
violin
ggsave("lvlnplot_before_qc.pdf", plot = violin, width = 14, height = 8)



quantile(pbmc$nFeature_RNA, seq(0.01,0.1,0.01))
quantile(pbmc$nFeature_RNA, seq(0.9,1,0.01))
quantile(pbmc$nCount_RNA, seq(0.01,0.1,0.01))
quantile(pbmc$nCount_RNA, seq(0.9,1,0.01))
quantile(pbmc$percent.mt, seq(0.9,1,0.01))
quantile(pbmc$percent.HB, seq(0.9,1,0.01))


minGene=200
maxGene=10000
minUMI=600
pctMT=10
pctHB=1


pbmc<- subset(pbmc, subset = nFeature_RNA> minGene & nFeature_RNA< maxGene &
                nCount_RNA> minUMI & percent.mt < pctMT & percent.HB< pctHB)
plot= list()
for(i in seq_along(plot.featrures)){
  plots[[i]]=VlnPlot(pbmc,group.by= group, pt.size=0,
                    features= plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots=plots, nrow = 2)   
violin
ggsave("2lvlnplot_after_qc.pdf", plot = violin, width = 14, height = 8)
dim(pbmc)  
saveRDS(pbmc,"1pbmc_qc.rds")


library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
pbmc= readRDS("1pbmc_qc.rds")


pbmc= pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

pbmc= pbmc %>%
  RunPCA() %>%
  RunUMAP(dims=1:30) %>%
  RunTSNE(dims=1:30) %>%
  FindNeighbors(dims=1:30) %>%
  FindClusters(resolution=0.1)

sweep.res.list<- paramSweep(pbmc, PCs=1:30, sct=FALSE)
sweep.stats<- summarizeSweep(sweep.res.list,GT= FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best= bcmvn %>%
  dplyr::arrange(desc(BCmetric))%>%
  dplyr::pull(pK) %>%
  .[1] %>% as.character() %>% as.numeric()

annotations<- pbmc$seurat_clusters
homotypic.prop<- modelHomotypic(annotations)
print(homotypic.prop)


nExp_poi <- round(0.07*nrow(pbmc@meta.data))
nExp_poi.adj<- round(nExp_poi*(1-homotypic.prop))


pbmc<- doubletFinder(pbmc,PCs=1:30,
                     pN=0.25,pK= pk_best, nExp = nExp_poi.adj,
                     reuse.pANN = FALSE, sct=FALSE)
colnames(pbmc@meta.data)

colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))-1] <- "Double_score"
colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))] <- "Is_Double"


head(pbmc@meta.data[,c("Double_score","Is_Double")])


DimPlot(pbmc, reduction = "tsne", group.by = "Is_Double")


VlnPlot(pbmc, group.by = "Is_Double",
        features = c("nCount_RNA", "nFeature_RNA"),
        pt.size = 0, ncol = 2)

saveRDS(pbmc, "2pbmc_double.rds")



library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
setwd("/home/data/t090502/Analyse/3000GSE221156")
pbmc= readRDS("2pbmc_double.rds")


pbmc<- NormalizeData(pbmc)


g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(pbmc))


s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(pbmc))


pbmc<- CellCycleScoring(pbmc, g2m.features = g2m_genes, s.features = s_genes)

colnames(pbmc@meta.data)
table(pbmc$Phase)


DimPlot(pbmc, group.by = "Phase", reduction = "tsne")
FeaturePlot(pbmc, features = "FAM114A2", ncol=3)
saveRDS(pbmc,"3pbmc_Cellcycle.rds")



library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
setwd("/home/data/t090502/Analyse/3000GSE221156")
pbmc= readRDS("3pbmc_Cellcycle.rds")


pbmc<- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc<- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
pbmc<- ScaleData(pbmc, vars.to.regress = c("S.Score","G2M.Score")) #用高变基因进行数据放缩
pbmc <- RunPCA(pbmc)
saveRDS(pbmc, "4pbmc_Normalize.rds")



library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(Rogue)
library(clustree)
library(harmony)
pbmc= readRDS("4pbmc_Normalize.rds")


DefaultAssay(pbmc)
table(Idents(pbmc))
Idents(pbmc)="orig.ident"
table(Idents(pbmc))
top10<- head(VariableFeatures(pbmc),10)


pdf(file="1.pdf",width = 7, height = 6)
VariableFeaturePlot(object = pbmc)
dev.off()
pdf(file="2.pdf",width = 7, height = 6)
LabelPoints(plot = VariableFeaturePlot(object = pbmc),points=top10, repel = TRUE)
dev.off()
pdf(file = "3.pdf",width = 7, height = 6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="4.pdf",width = 10, height = 9)
VizDimLoadings(object = pbmc, dims=1:4, reduction = "pca", nfeatures = 20)
dev.off()
pdf(file = "5.pdf", width = 10, height = 9)
DimHeatmap(object = pbmc, dims=1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()
pdf(file = "6.pdf", width = 7, height = 6)
ElbowPlot(pbmc, ndims = 50)
dev.off()
pct<- pbmc[["pca"]]@stdev/sum(pbmc[["pca"]]@stdev) *100
pct
cumu<- cumsum(pct)
cumu
pcs=1:41
pbmc<- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="RNA", max.iter.hamony=20)
table(pbmc@meta.data$orig.ident)


seq= seq(0.1, 2, by=0.1)
pbmc<- FindNeighbors(pbmc, dims = pcs)
for(res in seq) {
  pbmc= FindClusters(pbmc, resolution = res,graph.name = "RNA_snn") 
}
head(pbmc[[]])
p1= clustree(pbmc, prefix= "RNA_snn_res.")+coord_flip()  
p= p1+plot_layout(widths = c(3,1))
ggsave("RNA_sun_res.png", p, width = 30, height = 14)

pbmc<- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.5)
pbmc<- RunUMAP(pbmc, reduction = "harmony", dims = pcs) %>% RunTSNE(dims= pcs, resolution = "harmony")
colnames(pbmc@meta.data)

pdf(file = "7.pdf",width = 7, height = 6)
DimPlot(pbmc, reduction = "umap", label = T)
dev.off()

pdf(file = "8.pdf",width = 7, height = 6)
DimPlot(pbmc, reduction = "umap", label = F,group.by = "orig.ident")
dev.off()

pdf(file = "9.pdf",width = 7, height = 6)
DimPlot(pbmc, reduction = "umap", label = F,group.by = "Is_Double")
dev.off()

pdf(file = "10.pdf",width = 7, height = 6)
DimPlot(pbmc, reduction = "tsne", label = T)
dev.off()

pdf(file = "11.pdf",width = 7, height = 6)
DimPlot(pbmc, reduction = "tsne", label = F,group.by = "orig.ident")
dev.off()

pdf(file = "12.pdf",width = 7, height = 6)
DimPlot(pbmc, reduction = "tsne", label = F,group.by = "Is_Double")
dev.off()

saveRDS(pbmc,"5pbmc_UMPA.TSNE.rds")



library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(Rogue)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
pbmc= readRDS("5pbmc_UMPA.TSNE.rds")

DefaultAssay(pbmc)
table(Idents(pbmc))
pbmc.markers1<- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 1,)  
write.csv(pbmc.markers1, file = "markers.1.SCT.csv")  

DefaultAssay(pbmc)="RNA"
pbmc.markers2<- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 1)
write.csv(pbmc.markers2, file = "markers.2.SCT.csv")  
DefaultAssay(pbmc)
DefaultAssay(pbmc)="RNA"

top5_markers <- pbmc.markers2 %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5_markers
DotPlot(pbmc, features = unique(top5_markers$gene)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) + NoLegend()
DoHeatmap(pbmc, features = unique(top5_markers$gene)) +
  theme(axis.text.y = ggplot2::element_text(size = 8)) + NoLegend()
save(top5_markers, file = "top5_markers.csv")


markers<- c( "ATP2A3", "BMP5", "NKX6-1", "ADCYAP1",  #beta cells
             "ARRDC4", "ARX","DPP4", "PLCE1",   #alpha cells
             "MS4A8","BCHE", "HHEX", #delta cells
             "CD86","SDS", #macrophages
             "PROM1",#ducts
             "FGB","MEIS1","ENTPD2","THSD7A","CARTPT", #PPY/Gamma
             "VWF", #Endothelial cell
             "COL6A3","SFRP2","COL1A2","TNFAIP6", #stellate cells
             "ALB", "ALDOB" ,"PNLIP") #Acinar cell



p<- FeaturePlot(pbmc, features = markers, ncol = 3)
ggsave("48.pdf",p, width = 15, height = 32)
p<-DotPlot(pbmc, features = markers)+ RotatedAxis()
ggsave("49.pdf", p, width = 14, height=6)
p<- VlnPlot(pbmc, features = markers, stack = T, flip = T)+ NoLegend()
ggsave("50.pdf", p, width = 14, height=6)
p= DimPlot(pbmc, reduction = "umap", label = T, group.by="celltype.1")
ggsave("51.pdf",p , width = 7, height= 6)

p1= DimPlot(pbmc, reduction = "umap", label = T, group.by="celltype.1")
p2= DimPlot(pbmc, reduction = "tsne", label = T, group.by="celltype.1")
ggsave("52.pdf",p1|p2 , width = 15, height= 6)
table(Idents(pbmc))
DefaultAssay(pbmc)

pbmc.markers3<- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.3) #min.pct必须在30%的细胞中表达
write.csv(pbmc.markers3, file = "markers.celltype.RNA.csv")
top10<- pbmc.markers3 %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
markers= as.data.frame(top10[,"gene"])
pbmc<- ScaleData(pbmc, features = as.character(unique(markers$gene)))

p= DoHeatmap(pbmc, features = as.character(unique(markers$gene)),
             group.by = "celltype.1")
ggsave("53.pdf", p, width = 10, height=9)


allCells= names(Idents(pbmc))
allType= levels(Idents(pbmc))
Choose_Cells = unlist(lapply(allType, function(x){
  cgCells= allCells[Idents(pbmc)==x]
  cg=sample(cgCells, min(table(pbmc@meta.data$celltype.1)))
  cg
}))


cg_sce= pbmc[, allCells %in% Choose_Cells]
table(Idents(cg_sce))

p= DoHeatmap(cg_sce, features = as.character(unique(markers$gene)),
               group.by = "celltype.1")
ggsave("54.pdf", p, width = 10, height=10)

p=DimPlot(pbmc, reduction = "tsne", label = F,group.by = "group")
ggsave("56.pdf", p, width = 7, height=6)

p= DimPlot(pbmc, reduction = "umap", label = T,split.by = 'group')
ggsave("57.pdf",p , width = 15, height= 6)
saveRDS(pbmc, "7pbmc_celltype.rds")



library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(Rogue)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
pbmc= readRDS("7pbmc_celltype.rds")

DimPlot(pbmc)
names(pbmc@meta.data)
unique(pbmc$group)
DimPlot(pbmc,split.by = 'group')


table(Idents(pbmc))  
Idents(pbmc) <- "celltype.1"
pbmc$celltype.group <- paste(pbmc$celltype.1, pbmc$group, sep = "_")  
unique(pbmc$celltype.group) 
Idents(pbmc) <- "celltype.group"  

mydeg <- FindMarkers(pbmc,ident.1 = 'beta_cells_T2D',ident.2 = 'beta_cells_contrl', verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg, file = "mydeg_beta_cells.csv")

cellfordeg<-levels(pbmc$celltype.1)  
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(pbmc, ident.1 = paste0(cellfordeg[i],"_T2D"), ident.2 = paste0(cellfordeg[i],"_contrl"), verbose = T)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}
list.files()
write.csv(CELLDEG, file = "CELLDEG.csv")


library(dplyr)
top10 <- CELLDEG  %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10
pbmc <- ScaleData(pbmc, features =  rownames(pbmc))
DoHeatmap(pbmc,features = top10,size=3)  

Idents(pbmc) <- "celltype.1"
P<-VlnPlot(pbmc,features = top10,split.by = 'group') 


t2D_causal_genes11<-c("ABCC5", "APOL4", "B3GAT1", "ESPN", "NAP1L4", "PDK2",  "PROCR", "SAMD15", "SH2D2A", "SORBS1","TENT5C")
t2D_causal_genes11 <- c("B3GAT1", "NAP1L4", "NTAN1", "PROCR","SORBS1","TMCC2")
library(viridis)
DoHeatmap(pbmc,features = t2D_causal_genes11)  
p<- DotPlot(pbmc,features = t2D_causal_genes11,split.by = 'group')
ggsave("t2D_causal_genes1.pdf", p, width = 14, height=6)
p<- FeaturePlot(pbmc, features = "TENT5C",split.by = 'group')
ggsave("TENT5C.pdf",p, width = 6, height = 3)


library(readxl)
diff_genes <- read_excel("diff_genes.xlsx")
library(scRNAtoolVis)
head(diff_genes)
p<- jjVolcano(diffData = diff_genes)
ggsave("diff_genes.pdf",p, width = 10, height = 8)

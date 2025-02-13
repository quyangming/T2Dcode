
library(readxl)
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) 
library(data.table) 


enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 

if (T) {
  fpkm00 <- read_excel("Expressed_mRNA_genes_fpkm.xlsx") 
  ### 合并具有相同基因名的行
  table(duplicated(fpkm00$gene_name)) #统计重复基因名的基因
  gene <- fpkm00$gene_name ; fpkm00 <- fpkm00[,-1]
  fpkm0 <- aggregate(fpkm00, by=list(gene), FUN=sum)
  fpkm <- column_to_rownames(fpkm0,"Group.1")
}
data <- log2(fpkm+1)
zero_counts <- apply(data == 0, 1, sum)   
keep_data <- data[-which(zero_counts>ncol(data)*0.25),] #移除那些在大于30%的样本中表达量为0的miRNA
datTraits <- data.frame(row.names = colnames(data),group=colnames(data))
fix(datTraits) 
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
datTraits
datExpr0 <- as.data.frame(t(keep_data))

gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK


if(T){
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 4,height = 3)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  dev.off()
}


datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")



R.sq_cutoff = 0.85
if(T){
  powers <- c(seq(1,10,by = 1), seq(10,20,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  pdf("step2_power-value.pdf",width = 8,height = 6)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste( "Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=R.sq_cutoff ,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste( "Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}

sft$powerEstimate  


power = sft$powerEstimate
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average");



load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", 
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 50,    
    mergeCutHeight = 0.25, 
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
}

  
if(T){
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 12,height = 9)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")



rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")


if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group) 
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap.pdf",
      width = 3*length(colnames(design)), 
      height = 0.4*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) 
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}


if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") 
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              legend = "") +stat_compare_means()
  }
  
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width = 7.5,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) 
  dev.off()
}



levels(datTraits$group)
choose_group <- "T2D"  

if(T){
  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  

  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  selectModule <- modNames  
  pdf("step4_gene-Module-trait-significance.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2)) 
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}

 


load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
  if(T){
    geneTree = net$dendrograms[[1]]
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png("step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
 
  if(F){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf("step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}



if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)

  if(T){
    design
    primed = as.data.frame(design[,2])   
    names(primed) = "T2D"
    MET = orderMEs(cbind(MEs, primed))
  }
  pdf("step5_module_cor_Eigengene-dendrogram.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="", 
                        marDendro = c(0,4,1,4),  
                        marHeatmap = c(5,5,1,2), 
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}



rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step3_genes_modules.Rdata")
table(moduleColors)

module = "red"

if(T){
  dat=datExpr[,moduleColors==module]
  library(pheatmap)
  n=t(scale(dat)) 
  
  group_list=datTraits$group
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  
  pheatmap::pheatmap(n,
                     fontsize = 8,
                     show_colnames =T,
                     show_rownames = F,
                     cluster_cols = T,
                     annotation_col =ac,
                     width = 8,
                     height = 6,
                     angle_col=45,
                     main = paste0("module_",module,"-gene heatmap"),
                     filename = paste0("step7_module_",module,"_Gene-heatmap.pdf"))
  
}




rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
DOWN <- c("darkmagenta","darkred","skyblue","violet")
UP <- c("darkgreen","lightyellow","salmon","sienna3","yellow","turquoise")
module = "red"
 
gene <- colnames(datExpr) 
gene <- colnames(datExpr) 
inModule <- moduleColors=="turquoise"
modgene <- gene[inModule]
modgene<- data.frame(modgene)
write.csv(modgene,'turquoise.csv',row.names = T)



library(caret)
train_dataset<-datExpr
m<-c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
train_dataset<-cbind(m,train_dataset)
colnames(train_dataset)[1]<-"group"
inTrain<-createDataPartition(y=train_dataset$group,p=0.7,list=FALSE)
train<-train_dataset[inTrain,-1]
test<-train_dataset[-inTrain,-1]
setLabels = c("Train", "Test")
multiExpr = list(Train = list(data = train), Test = list(data = test))
multiColor = list(Train = moduleColors)
nSets = 2

mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1,
                        nPermutations = 200,
                        randomSeed = 1,
                        quickCor = 0,
                        verbose = 3)
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]
row.names(statsZ[statsZ$Zsummary.pres<10,])

plotMods = !(modColors %in% row.names(statsZ[statsZ$Zsummary.pres<0,]))  

text = modColors[plotMods]
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])


mains = c("Preservation Median rank", "Preservation Zsummary")

sizeGrWindow(10, 5)
pdf("step9_preservation.pdf",width = 20, height = 10)

par(mfrow = c(1,2))

par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2){
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  if (p==2){
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 3.0,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.8, cex.axis = 1.8, cex.main =2.0)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1.4, offs = 0.08);
  if (p==2){
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()

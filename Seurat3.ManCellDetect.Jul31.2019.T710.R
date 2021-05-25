library(Seurat)
library(Matrix)
library(dplyr)
library(cowplot)
#library(harmony)
library(MAST)


options(bitmapType='cairo') #to be able to run png() on severs without X11

setwd("/home/lolab/")
##args to supply = dirList, pctMitoCutoff, prefix for the Seurat combined object file
dirList = "scRNAseq.List2.cleaned.txt"
pctMitoCutoff = 0.1
prefix = "Prins.allSample.scRNAseq.Jul31.2019"

# Read in the data
dirs <- read.delim(dirList, header=T, sep="\t")
ID = rownames(dirs)
#samples = dirs$Sample
ob.list = list()
objectlist = c()


# Update Jun 2019: A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#convert human to mouse genes
require("biomaRt")
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

AnnotFile<-read.delim("HOM_MouseHumanSequence.rpt")
convertHumanGeneList <- function(x){
  human = AnnotFile[AnnotFile$Common.Organism.Name=="human",]
  mouse = AnnotFile[AnnotFile$Common.Organism.Name=="mouse, laboratory",]
  gene_human = x
  
  human_subset <- AnnotFile[AnnotFile$Symbol %in% x,]$HomoloGene.ID
  gene_mouse = as.character(droplevels(mouse[mouse$HomoloGene.ID %in% human_subset,]$Symbol))
  
  # Print the first 6 genes found to the screen
  print(head(gene_mouse))
  return(gene_mouse)
}

s.genes <- convertHumanGeneList(s.genes)
g2m.genes <- convertHumanGeneList(g2m.genes)

## read in each single sample for pre-filtering ################
##note Jul23 2018: We download the raw_gene_bc_matrices from T710 and 
#rename it as filtered_gene_bc_matrices, so that cellranger can read it
#i=1
for (i in 1:length(ID)){
	dirID <- as.character(dirs[i,3])
  sampleID <- as.character(dirs[i,2])
  ID <- as.character(dirs[i,1])
  cellranger_pipestance_path <- paste("/media/lolab/Gabri/scRNA_cellranger_out/",ID,"/",dirID,sep="")
  
	gbm <- Read10X(data.dir =cellranger_pipestance_path) #load it for the first time. will save a new gbm file after cleaning. Everytime afterwards will load the new gbm
	
	gbm = gbm[rowSums(gbm != 0) >= 20,]
	gbm = gbm[,colSums(gbm != 0) >= 200]
	dim(gbm)
	
	nUMI <- Matrix::colSums(gbm)
	nGene <- Matrix::colSums(gbm != 0)
	
	mito.genes <- grep(pattern = "^MT-", x = rownames(gbm), value = TRUE)
	percent.mito <- Matrix::colSums(gbm[which(rownames(gbm) %in% mito.genes),])/Matrix::colSums(gbm)
	
	#png(paste(sampleID,".UMI.pctMito.png",sep=""))
	#plot(nUMI,percent.mito,xlab="nUMI",ylab="pctMito")
	#dev.off()
	
	ribosome.genes <- grep(pattern = "^RP[LS]", x = rownames(gbm), value = TRUE)
	percent.ribosome <- Matrix::colSums(gbm[which(rownames(gbm) %in% ribosome.genes),])/Matrix::colSums(gbm)
	
	#png(paste(sampleID,".UMI.pctRibo.png",sep=""))
	#plot(nUMI,percent.ribosome,xlab="nUMI",ylab="pctRibo")
	#dev.off()
	
	#gbm = gbm[,percent.mito<=pctMitoCutoff] #do not run this since we already regress out mitochondrial later
	gbm.seurat <- CreateSeuratObject(counts = gbm, project = sampleID, min.cells = 5)
	gbm.seurat$ID <- sampleID
	gbm.seurat$condition <- dirs[i,5]
	gbm.seurat$pembro_response <- dirs[i,6]
	
	gbm.seurat <- subset(x = gbm.seurat, subset = nFeature_RNA > 200)
	#update Jun 2019: we use the new SCTransform function to do normalization and variance(e.g.,mitochondrial mapping percentage) stablization, 
	gbm.seurat <- NormalizeData(object = gbm.seurat, verbose = FALSE)
	gbm.seurat <- FindVariableFeatures(object = gbm.seurat, selection.method = "vst", nfeatures = 2000)
	gbm.seurat <- PercentageFeatureSet(gbm.seurat, pattern = "^MT-", col.name = "percent.mt")
	gbm.seurat <- CellCycleScoring(gbm.seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #We assign scores in the CellCycleScoring function, which stores S and G2/M scores in object meta data, along with the predicted classification of each cell in either G2M, S or G1 phase.
	gbm.seurat$CC.Difference <- gbm.seurat$S.Score-gbm.seurat$G2M.Score
	#gbm.seurat <- SCTransform(object = gbm.seurat, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = FALSE)
	
	ob.list = append(ob.list,gbm.seurat)
	objectlist = c(objectlist,sampleID)
	gc()
}

rm(gbm)
rm(gbm.seurat)
gc()

k.filterMax <- min(200, min(sapply(ob.list, ncol)))
#dirs$cellNumbers = sapply(ob.list, ncol)
#write.table(dirs,file="scRNAseq.List3.txt",sep="\t")

#anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:20,k.filter = 150)
anchors <- FindIntegrationAnchors(object.list = ob.list)
gbm.combined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(object = gbm.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#gbm.combined <- ScaleData(object = gbm.combined, verbose = FALSE, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
gbm.combined <- ScaleData(object = gbm.combined, verbose = FALSE, vars.to.regress = c("percent.mt", "CC.Difference"))

gbm.combined <- RunPCA(object = gbm.combined, npcs = 30, verbose = FALSE)

gbm.combined <- RunHarmony(gbm.combined, group.by.vars="old.ident")
harmony_embeddings <- Embeddings(gbm.combined, 'harmony')

saveRDS(gbm.combined,file=paste(prefix, "gbm.combined.rds", sep="_"))
#gbm.combined=readRDS(file=paste(prefix, "gbm.combined.rds", sep="_"))



# t-SNE and Clustering
gbm.combined <- RunUMAP(object = gbm.combined, reduction = "harmony", dims = 1:30)
#gbm.combined <- RunUMAP(object = gbm.combined, reduction = "pca", dims = 1:30, min_dist = 0.3)
gbm.combined <- FindNeighbors(object = gbm.combined, reduction = "harmony", dims = 1:30)
#gbm.combined <- FindClusters(gbm.combined, resolution = 0.3)
gbm.combined <- FindClusters(gbm.combined, resolution = 0.2)

saveRDS(gbm.combined,file=paste(prefix, "gbm.combined.rds", sep="_"))
gbm.combined=readRDS(file=paste(prefix, "gbm.combined.rds", sep="_"))


#####Plot clustering pattern

#level<-unique(gbm.combined@meta.data$old.ident)

#gbm.combined$ID <- factor(gbm.combined$ID, 
#                          levels = level) #so that each sample will have the same number clusters tablized (even no cells from that cluster)


outfile = paste(prefix,"_tSNEplot.byBatches.png",sep="")
png(filename = outfile,width = 800, height = 1200)
DimPlot(object = gbm.combined, reduction = "umap", group.by = "ID", pt.size = 1)
dev.off()

outfile = paste(prefix,"_tSNEplot.Res0.2.png",sep="")
png(filename = outfile,width = 800, height = 800)
DimPlot(object = gbm.combined, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE, label.size = 8, pt.size = 1)
dev.off()

outfile = paste(prefix,"_tSNEplot.nUMI.png",sep="")
png(filename = outfile,width = 800, height = 800)
FeaturePlot(object = gbm.combined, features = c("nCount_RNA"), min.cutoff = "q9", pt.size = 2) #umi plot
dev.off()

outfile = paste(prefix,"_tSNEplot.nGene.png",sep="")
png(filename = outfile,width = 800, height = 800)
FeaturePlot(object = gbm.combined, features = c("nFeature_RNA"), min.cutoff = "q9", pt.size = 2) #umi plot
dev.off()

#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
DefaultAssay(object = gbm.combined) <- "RNA" #Need to assign RNA instead of integrated as default assay, coz integrated only contains the top 2000 features

outfile = paste(prefix,"_tSNEplot.Markers.png",sep="")
png(filename = outfile,width = 1200, height = 800)

FeaturePlot(object = gbm.combined, features = c("Sox9","Mitf","Axl","Ptprc","Thy1","Col1a1","Luciferase"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
FeaturePlot(object = gbm.combined, features = c("Ccnd1","Mki67","Cdk1","Pcna"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
FeaturePlot(object = gbm.combined, features = c("Cd3d","Cd4","Cd8a","Il2ra","Foxp3","Il7r"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
FeaturePlot(object = gbm.combined, features = c("Pdcd1","Ctla4","Cd274","Pdcd1lg2","Havcr2","Icos"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
FeaturePlot(object = gbm.combined, features = c("Sell","Lag3","Tnfrsf9","Tnfrsf4","Tnfrsf18"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
FeaturePlot(object = gbm.combined, features = c("Cd69","Ncr1","Klrg1","Fcgr3","Il2rb","Itga2"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red")) #NKp46 (also called NCR1 or CD335) has been put forward as the most specific marker for NK cells in mammalian cells (Westgaard et al. 2004; Walzer et al. 2007).
#FeaturePlot(object = gbm.combined, features = c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2","Tgfbr3"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
#S100A9/8 a new marker for monocytic human myeloid-derived suppressor cells, Cxcl2 is a macrophage marker
FeaturePlot(object = gbm.combined, features = c("Cd14","Itgam","Cd68","Fcgr3","Itgax","H2-Aa","S100a8","S100a9","Cxcl2"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))

#Pecam1, Vegfa (endothelial marker) ; Cd19, Cd79a (B cell marker); Batf3(DC marker); 
FeaturePlot(object = gbm.combined, features = c("Pecam1","Cd19","Cd79a","Batf3","Apoe","Vegfa"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
###denritic cell
FeaturePlot(object = gbm.combined, features = c("Ccr7","Cd7","Cd24a","Cd209a","Irf8","Batf"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))

#Pdpn(lymphatic vessel marker),Cd248(tumor associated pericytes) 
FeaturePlot(object = gbm.combined, features = c("Pdpn","Cd248","Acta2","Pdgfrb","Fap","Mmp2","S100a4","Gas6","Ogn"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))
FeaturePlot(object = gbm.combined, features = c("mt-Nd4","mt-Nd5","mt-Co1","mt-Cytb"), min.cutoff = "q9",  pt.size = 0.5, ncol=3, cols = c("lightgrey","red"))

dev.off()

levels(gbm.combined)

# find markers for every cluster compared to all remaining cells, report
# only the positive ones

DefaultAssay(gbm.combined) <- "RNA"
Idents(object = gbm.combined) <- "integrated_snn_res.0.2"
cluster.markers<- FindAllMarkers(object = gbm.combined, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.585,test.use="MAST")

outfile = paste(prefix,"_ClusterDiffGene.Res0.2.csv",sep="")
write.csv(cluster.markers, file=outfile)

#generates an expression heatmap for given cells and genes. In this case, we are plotting 
#the top 5 markers (or all markers if less than 10) for each cluster.
top <- cluster.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

#assign each cluster cell type
current.cluster.ids <- c(0:17)
new.cluster.ids <- c("Tumor","Tumor","Tcell","Tumor","Tumor","Tumor",
                     "Myeloid","TumorDoublets","NK","Myeloid","Tcell","Myeloid",
                     "TumorDoublets","Fibroblast","Myeloid","Endothelial","Myeloid","Tumor")
Idents(object = gbm.combined) <- plyr::mapvalues(x = Idents(object = gbm.combined), from = current.cluster.ids, to = new.cluster.ids)

outfile = paste(prefix,"_tSNEplot.Celltype.png",sep="")
png(filename = outfile,width = 800, height = 800)
DimPlot(gbm.combined, label = TRUE, pt.size = 1,label.size = 5)
dev.off()

gbm.combined@meta.data$celltype <- Idents(object = gbm.combined)

saveRDS(gbm.combined,file=paste(prefix, "gbm.combined.rds", sep="_"))
gbm.combined=readRDS(file=paste(prefix, "gbm.combined.rds", sep="_"))

#Cluster 6 (NK) was not well defined. We can manualy select it
# p <- DimPlot(object=gbm.combined, reduction = "umap")
# select.cells <- CellSelector(plot = p,object=NULL)
# Idents(gbm.combined, cells = select.cells) <- "Myeloid"
# 
# saveRDS(gbm.combined,file=paste(prefix, "gbm.combined.rds", sep="_"))
#gbm.combined=readRDS(file=paste(prefix, "gbm.combined.rds", sep="_"))


###### Percentage of each celltype in each sample ###########################################
gbm.combined@meta.data$celltype <- Idents(object = gbm.combined)

gbm.combined@meta.data$celltype <- factor(gbm.combined@meta.data$celltype, 
                                                 levels = c("Tumor","Tcell","NK","Myeloid","Fibroblast","Endothelial","TumorDoublets")) #so that each sample will have the same number clusters tablized (even no cells from that cluster)


##Tabulate each cell type's numbers
totalNum=c(table(gbm.combined@meta.data$ID)[[1]])
Group = names(table(gbm.combined@meta.data$ID))
table_sub = data.frame()
table_sub = subset(gbm.combined@meta.data, ID==Group[1])
percentages=as.data.frame((table(table_sub$celltype)))$Freq
table1 = percentages

for(i in 2:length(Group)){
  id = Group[i]
  n_total = table(gbm.combined@meta.data$ID)[[i]]
  totalNum=c(totalNum,n_total)
  
  table_sub = data.frame()
  table_sub = subset(gbm.combined@meta.data, ID==id)
  
  percentages=as.data.frame((table(table_sub$celltype)))$Freq
  #percentages=countBatchInCluster(table_sub$celltype,levels(table_sub$celltype))
  
  #print(length(percentages))
  #print(percentages)
  table1=as.data.frame(rbind(table1,percentages))
}

table2 = as.data.frame(cbind(table1,Total=totalNum))
rownames(table2)=Group
colnames(table2)=c(levels(gbm.combined@meta.data$celltype),"Total")
outfile = paste(prefix,"_Celltype.txt",sep="")
write.table(table2,file=outfile,sep="\t")


################################################################
# extract subset cells (perhaps, to load into another package)
###############################################################
ids = rownames(subset(gbm.combined@meta.data,celltype %in% c("Tcell","NK")))
fileConn<-file(paste(prefix,".TcellandNK.ID.txt",sep=""))
writeLines(ids, fileConn)

ids = rownames(subset(gbm.combined@meta.data,celltype %in% c("Myeloid")))
fileConn<-file(paste(prefix,".Myeloid.ID.txt",sep=""))
writeLines(ids, fileConn)

ids = rownames(subset(gbm.combined@meta.data,celltype %in% c("Tumor")))
fileConn<-file(paste(prefix,".Tumor.ID.txt",sep=""))
writeLines(ids, fileConn)




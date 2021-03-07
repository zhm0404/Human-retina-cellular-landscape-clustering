library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

dir <- "/Users/zhm/R/785project/GSE137537_counts/"

r <- readMM(paste0(dir,"matrix.mtx"))
barcodes <- read.csv(paste0(dir,"barcodes.tsv"), header = F, sep = "\t", skip = 1)
genenames <- read.table(paste0(dir,"features.txt.gz"),header = F,sep = "\t")
gene_name <- genenames[,1]
gene_name_u <- gene_name[which(!duplicated(gene_name))]
r_u <- r[which(!duplicated(gene_name)),]
dim(r_u)

rownames(r_u) <- gene_name_u
colnames(r_u) <- barcodes$V1

#QC and selecting cells for fuether analysis
retina <- CreateSeuratObject(counts = r_u, project = "retina", min.cells = 3, min.features = 200)
r.violin <-  VlnPlot(retina,features = c("nFeature_RNA","nCount_RNA"),ncol = 3)
r.scatter <- FeatureScatter(retina,feature1 = "nCount_RNA","nFeature_RNA")
r.violin + r.scatter
retina <- subset(retina,nFeature_RNA<5000&nCount_RNA<20000)

r_filterd <- GetAssayData(retina,slot="data",assay="RNA")

#Normalization
retina <- NormalizeData(retina, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
retina <- FindVariableFeatures(retina, selection.method = "vst", nfeatures = 2000)
retina_top10 <- head(VariableFeatures(retina),10)
retina_plot1 <- VariableFeaturePlot(retina)
retina_plot2 <- LabelPoints(plot = retina_plot1, points = retina_top10, repel = T)
retina_plot2

#sclaing the data
all.genes <- rownames(retina)
retina <- ScaleData(retina, features = all.genes)

#Perform linear dimensional reduction
retina <- RunPCA(retina,features = VariableFeatures(object = retina))

VizDimLoadings(retina, dims = 1:2, reduction = "pca")
DimPlot(retina, reduction = "pca")
DimHeatmap(retina, dims = 1:20, cells = 500, balanced = TRUE)

ElbowPlot (retina, ndims = 50, reduction = "pca")

#Cluster the cells
retina <- FindNeighbors(retina, dims = 1:40)
retina <- FindClusters(retina, resolution = 0.5)

head(Idents(retina),5)

retina <- RunUMAP(retina, dims = 1:40)
DimPlot(retina, reduction = "umap")
FeaturePlot(retina, features = c("PDE6A","GNAT2","NEFL","CAMK2B","GAD1","ONECUT1","GLUL"
                                  ,"C1QA","TM4SF1"),pt.size = 1)
DoHeatmap(retina, features = c("PDE6A","GNAT2","NEFL","CAMK2B","GAD1","ONECUT1","GLUL"
                             ,"C1QA","TM4SF1")) + NoLegend()
gene.viz <- c("PDE6A","GNAT2","NEFL","CAMK2B","GAD1","ONECUT1","GLUL","C1QA","TM4SF1")
#Find and Visualize Cell-Type Specific Gene Expression
DotPlot(retina, features = gene.viz) + RotatedAxis()
#VlnPlot(retina, features = gene.viz, combine = TRUE)
new.cluster.ids<-c("Rods","macroglia","retinal ganglion cells","Rods","Rods",
                   "bipolar cells","bipolar cells","macroglia","bipolar cells","macroglia",
                   "amacrine cells","bipolar cells","bipolar cells","macroglia","macroglia",
                   "macroglia","Cones","horizontal cells","microglia","vascular cells","retinal ganglion cells")
names(new.cluster.ids) <- levels(retina)
retina <- RenameIdents(retina, new.cluster.ids)
DimPlot(retina, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# install.packages("devtools")
# install.packages("igraph")
# devtools::install_github("VCCRI/CIDR")
#devtools::install_github("yycunc/SAMEclustering")
library("SAMEclustering")

cluster.result1 <- individual_clustering(r_u1, mt_filter = TRUE,
                                        percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = TRUE,
                                        nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, SEED = 123)
#cluster.result4 <- individual_clustering(inputTags = r_u, mt_filter = TRUE,
     #                                    percent_dropout = 10, SC3 = FALSE, CIDR = FALSE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
      #                                   nPC.seurat = NULL, resolution = 0.7, tSNE = FALSE, dimensions = 2, perplexity = 30, SIMLR = FALSE, diverse = TRUE, SEED = 123)
dim(r_u1)
rarray<- as.array(r_u1)
indi_cluster_tsne1<- tSNE_kmeans_SAME(rarray,percent_dropout = 10,dimensions = 2,perplexity = 30,k.min = 2,k.max = 10, var_genes = 2000,SEED=123)

cluster.results1 <- readRDS(file = "/Users/zhm/R/785Project/GSE137846_counts/indi_cluster_same.rds")

cluster.results1[,1:40]

unique(cluster.result1[4,])

indi_cluster_tsne1$cluster[1:40]
sc3_annotation<- cluster.results1[1,]
cidr_annotation <- cluster.results1[2,]
seurat_annotation<- cluster.results1[3,]
cluster.results1[3,] <- seurat0
simlr_annotation <- cluster.results1[4,]
tsne_annotation <- indi_cluster_tsne1$cluster
ensemble

#gpmetis and shmetis are stored in program.dir "/Users/zhm/R/785project"
#cluster.ensemble <- SAFE(cluster_results = cluster.results1, program.dir = "/Users/zhm/R/785project", 
#                         MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = 123)
#cluster.results1[1:4, 1:10]

cluster.ensemble <- SAMEclustering(Y = t(cluster.results1), rep = 3, SEED = 123)
cluster.ensemble$BIC_result_summary
library(cidr)
data("data_SAME")
data_SAME$Zheng.celltype

adjustedRandIndex(annotation,cluster.ensemble$BICcluster)





r1 <- readMM("/Users/zhm/R/785project/GSE137846_counts/GSE137846_Seq-Well_counts.mtx.gz")
dim(r1)
barcodes1 <- read.csv("/Users/zhm/R/785project/GSE137846_counts/GSE137846_Seq-Well_sample_annotations.txt.gz", header = F, sep = "\t",skip=1)
genenames1 <- read.table("/Users/zhm/R/785project/GSE137846_counts/GSE137846_Seq-Well_gene_names.txt.gz",header = F,sep = "\t")
gene_name1 <- genenames1[,1]
gene_name_u1 <- gene_name1[which(!duplicated(gene_name1))]
r_u1 <- r1[which(!duplicated(gene_name1)),]
dim(r_u1)

rownames(r_u1) <- gene_name_u1
colnames(r_u1) <- barcodes1$V1
annotation<- barcodes1$V49

retina1@meta.data[,"annotation"] <- annotation
retina1@meta.data[,"sc3annotation"] <- sc3_annotation
retina1@meta.data[,"cidrannotation"] <- cidr_annotation
retina1@meta.data[,"tsneannotation"] <- tsne_annotation
retina1@meta.data[,"simlrannotation"] <- simlr_annotation
retina1@meta.data[,"sameannotation"] <- cluster.ensemble$BICcluster
#QC and selecting cells for fuether analysis
retina1 <- CreateSeuratObject(counts = r_u1, project = "retina1", min.cells = 3, min.features = 200)
r.violin1 <-  VlnPlot(retina1,features = c("nFeature_RNA","nCount_RNA"),ncol = 3)
r.scatter1 <- FeatureScatter(retina1,feature1 = "nCount_RNA","nFeature_RNA")
r.violin1 + r.scatter1
retina1 <- subset(retina1,nFeature_RNA<4000&nCount_RNA<15000)

#Normalization
retina1 <- NormalizeData(retina1, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
retina1 <- FindVariableFeatures(retina1, selection.method = "vst", nfeatures = 2000)
retina1_top10 <- head(VariableFeatures(retina1),10)
retina1_plot1 <- VariableFeaturePlot(retina1)
retina1_plot2 <- LabelPoints(plot = retina1_plot1, points = retina1_top10, repel = T)
retina1_plot2

#sclaing the data
all.genes <- rownames(retina1)
retina1 <- ScaleData(retina1, features = all.genes)

#Perform linear dimensional reduction
retina1 <- RunPCA(retina1,features = VariableFeatures(object = retina1))

VizDimLoadings(retina1, dims = 1:2, reduction = "pca")
DimPlot(retina1, reduction = "pca")
DimHeatmap(retina1, dims = 1:20, cells = 500, balanced = TRUE)

ElbowPlot (retina1, ndims = 50, reduction = "pca")

#Cluster the cells
retina1 <- FindNeighbors(retina1, dims = 1:40)
retina1 <- FindClusters(retina1, resolution = 0.5)

head(Idents(retina1),5)

retina1 <- RunUMAP(retina1, dims = 1:40)
DimPlot(retina1, reduction = "umap")
FeaturePlot(retina1, features = c("PDE6A","GNAT2","NEFL","CAMK2B","GAD1","ONECUT1","GLUL"
                                 ,"C1QA","TM4SF1"),pt.size = 1)
DoHeatmap(retina1, features = c("PDE6A","GNAT2","NEFL","CAMK2B","GAD1","ONECUT1","GLUL"
                               ,"C1QA","TM4SF1")) + NoLegend()
gene.viz <- c("PDE6A","GNAT2","NEFL","CAMK2B","GAD1","ONECUT1","GLUL","C1QA","TM4SF1")
#Find and Visualize Cell-Type Specific Gene Expression
DotPlot(retina1, features = gene.viz) + RotatedAxis()
VlnPlot(retina1, features = gene.viz, combine = TRUE)

plotmenon<- DimPlot(retina1, reduction = "umap",group.by = "annotation", pt.size = 0.1)
plotseurat <- DimPlot(retina1, reduction = "umap", pt.size = 0.1)
plottsne <- DimPlot(retina1, reduction = "umap",group.by = "tsneannotation", pt.size = 0.1)
plotsame <- DimPlot(retina1, reduction = "umap",group.by = "sameannotation", pt.size = 0.1)
plotmenon
plotseurat+plotmenon
plottsne+plotmenon
plotsame+plotmenon
seurat0<- retina1@active.ident


##################### scRNA-seq
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(cowplot)
library(data.table)
library(limma)
library(harmony)
#install.packages("Matrix")
#install.packages("Matrix.utils")
library(org.Hs.eg.db)
library(clusterProfiler)
library(Matrix)
library(Matrix.utils)
library(dplyr)
library(hdf5r)
library(cowplot)
library(data.table)
library(DoubletFinder)

setwd("./20.kidney_Capsule_sc/01.scRNA_satandard_pipline") 
# sample list
sample_ids <- c("YKRRD230303001", "YKRRD230303002",
                "YKRRD230309001", "YKRRD230309002",
                "YKRRD230315003", "YKRRD230315004",
                "YKRRD230320001", "YKRRD230320002")

# raw data directory
base_dir <- "./18.kidney_Capsule/01.Data/01.scRNA"
read_one_sample <- function(sample_id, base_dir) {
  mtx_file  <- file.path(base_dir, sample_id, "matrix.mtx")
  gene_file <- file.path(base_dir, sample_id, "genes.tsv")
  bc_file   <- file.path(base_dir, sample_id, "barcodes.tsv")
  # read matrix
  mat <- Matrix::readMM(file = mtx_file)
  feature.names <- read.delim(gene_file, header = FALSE, stringsAsFactors = FALSE)
  barcode.names <- read.delim(bc_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V1
  # ENSEMBL -> SYMBOL
  gsym.id <- bitr(rownames(mat), fromType = "ENSEMBL", 
                  toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
  mat <- as.matrix(mat)
  mat <- cbind(ENSEMBL = rownames(mat), mat)
  gsym.fc.id <- merge(gsym.id, mat, by = "ENSEMBL", all = FALSE)
  gsym.fc.id <- gsym.fc.id[, 2:ncol(gsym.fc.id)]
  
  rt1 <- as.matrix(gsym.fc.id)
  rownames(rt1) <- rt1[, 1]
  exp1 <- rt1[, 2:ncol(rt1)]
  dimnames1 <- list(rownames(exp1), colnames(exp1))
  data1 <- matrix(as.numeric(as.matrix(exp1)), nrow = nrow(exp1), dimnames = dimnames1)
  mat <- avereps(data1)  
  # create seurat object
  seu <- CreateSeuratObject(
    count = mat, 
    names.delim = "_",
    min.cells = 1,
    min.features = 1,
    project = sample_id
  )
  seu[['percent.mt']] <- PercentageFeatureSet(seu, pattern = "^MT-")
  return(seu)
}

# deal with all samples
scRNAlist <- list()
for (i in seq_along(sample_ids)) {
    scRNAlist[[i]] <- read_one_sample(sample_ids[i], base_dir)
    print(paste0("Finished: ", sample_ids[i]))
}

saveRDS(scRNAlist, "kidney_Capsule_Integrated.rds")
######################################################################################
# scRNAlist <- readRDS("./18.kidney_Capsule/01.Data/01.scRNA/kidney_Capsule_Integrated.rds")
scRNA_harmony <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]],scRNAlist[[7]],scRNAlist[[8]]))
dim(scRNA_harmony)
saveRDS(scRNA_harmony, "kidney_Capsule_merged.rds")

######################################################################################
# scRNA_harmony <- readRDS("kidney_Capsule_merged.rds")
dim(scRNA_harmony)
scRNA_harmony <- subset(scRNA_harmony, subset = nCount_RNA > 200 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
dim(scRNA_harmony)
table(scRNA_harmony@meta.data$orig.ident)
scRNA_harmony@meta.data$Region<-"unknown"
scRNA_harmony@meta.data$Sample<-"unknown"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230303001","Region"]="Capsule"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230303002","Region"]="Cancer"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230309001","Region"]="Cancer"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230309002","Region"]="Capsule"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230315003","Region"]="Capsule"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230315004","Region"]="Cancer"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230320001","Region"]="Cancer"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230320002","Region"]="Capsule"

scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230303001","Sample"]="Sample1"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230303002","Sample"]="Sample1"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230309001","Sample"]="Sample2"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230309002","Sample"]="Sample2"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230315003","Sample"]="Sample3"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230315004","Sample"]="Sample3"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230320001","Sample"]="Sample4"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$orig.ident=="YKRRD230320002","Sample"]="Sample4"
table(scRNA_harmony@meta.data$Sample,scRNA_harmony@meta.data$Region)
## PCA reduction
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scRNA_harmony <- CellCycleScoring(scRNA_harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scRNA_harmony@meta.data[1:5,]
scRNA_harmony<-NormalizeData(scRNA_harmony,verbose = T) 
scRNA_harmony<-FindVariableFeatures(scRNA_harmony,selection.method = "vst", nfeatures = 2000)
scRNA_harmony<-ScaleData(scRNA_harmony,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = FALSE)
scRNA_harmony<-RunPCA(scRNA_harmony,verbose = T,npcs = 50)
pdf(file="1.ElbowPlot.pdf",width=6,height=6)
ElbowPlot(scRNA_harmony,ndims = 50)
dev.off()
saveRDS(scRNA_harmony, "kidney_PreProcess.rds")
######################################################################################
# scRNA_harmony <- readRDS("kidney_PreProcess.rds")
scRNA_harmony<-RunHarmony(scRNA_harmony,group.by.vars = c("orig.ident","Sample","Region"), plot_convergence = TRUE)

scRNA_harmony <- scRNA_harmony %>% 
                RunUMAP(reduction = "harmony", dims = 1:50) %>% 
                RunTSNE(reduction = "harmony", dims = 1:50) %>%
                FindNeighbors(reduction = "harmony", dims = 1:50) %>%
                FindClusters(resolution = 0.4)

scRNA_harmony@meta.data[1:4,]
saveRDS(scRNA_harmony, "kidney_Harmony2.rds")
######################################################################################
kidney <- readRDS("kidney_Harmony2.rds")
mycols <- c("0"='#E5D2DD', 
			"1"='#53A85F',
			"2"='#F1BB72',
			"3"='#F3B1A0',
			"4"='#D6E7A3',
			"5"='#57C3F3',
			"6"='#476D87',
			"7"='#E95C59',
			"8"='#E59CC4',
			"9"='#AB3282',
			"10"='#23452F',
			"11"='#BD956A',
			"12"='#8C549C',
			"13"='#585658',
			"14"='#9FA3A8',
			"15"='#E0D4CA',
			"16"='#5F3D69',
			"17"='#C5DEBA',
			"18"='#58A4C3',
			"19"='#E4C755',
			"20"='#F7F398',
			"21"='#AA9A59',
			"22"='blue')
cols =  mycols
	
pdf(file="1.UMAP_newCol2.pdf",width=6.5,height=5)
DimPlot(kidney, 
        reduction = "umap",
        pt.size = 0.1,
		label.size = 4,
		cols =  mycols,
        label = T)
dev.off()

kidney_val_db <- kidney
sweep.res.list_kidney_val <- paramSweep_v3(kidney_val_db, PCs = 1:50)
sweep.stats_kidney_val <- summarizeSweep(sweep.res.list_kidney_val, GT = FALSE)
bcmvn_kidney_val<-find.pK(sweep.stats_kidney_val)
pK_value <- as.numeric(as.character(bcmvn_kidney_val$pK[bcmvn_kidney_val$BCmetric == max(bcmvn_kidney_val$BCmetric)]))
annotations <- kidney_val_db@meta.data$RNA_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)  # get the estimated number
nExp_poi <- round(homotypic.prop*length(kidney_val_db@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
kidney_val_db <- doubletFinder_v3(kidney_val_db, PCs = 1:50, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
kidney_val_db <- doubletFinder(kidney_val_db, PCs = 1:50, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
kidney_val_db@meta.data[1:5,]
saveRDS(kidney_val_db, "kidney_Capsule_0.4_doublets_UMAP2.rds")
######################################################################################
kidney<-readRDS("kidney_Capsule_0.4_doublets_UMAP2.rds")
kidney$Doublets<-kidney$DF.classifications_0.25_0.26_10813
pdf(file="1.Doublets_TSNE.pdf",width=6,height=5)
DimPlot(kidney,reduction = "tsne",group.by = "Doublets")
dev.off()
pdf(file="1.Doublets_UMAP.pdf",width=6,height=5)
DimPlot(kidney,reduction = "umap",group.by = "Doublets")
dev.off()

kidney_dedoublets <- subset(kidney, subset = Doublets  %in% c("Singlet"))
saveRDS(kidney_dedoublets, "kidney_Capsule_delete_doublets2.rds")
######################################################################################
scRNA_harmony<-readRDS("./kidney_Capsule_delete_doublets2.rds")
table(scRNA_harmony@meta.data$orig.ident)

pdf(file="1.ElbowPlot.pdf",width=6,height=6)
ElbowPlot(scRNA_harmony,ndims = 50)
dev.off()

scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:40, verbose = F)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:40)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony",dims = 1:40)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.4)
scRNA_harmony@meta.data[1:4,]
saveRDS(scRNA_harmony, "kidney_Harmony_Singlet2.rds")
######################################################################################
scRNA_harmony<-readRDS("kidney_Harmony_Singlet2.rds")
scRNA_harmony@meta.data$CellMajor<-"unknow"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==0,"CellMajor"]="T"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==5,"CellMajor"]="T"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==11,"CellMajor"]="T"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==15,"CellMajor"]="T"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==16,"CellMajor"]="T"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==4,"CellMajor"]="NK"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==13,"CellMajor"]="B/Plasma"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==3,"CellMajor"]="Myeloid"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==9,"CellMajor"]="Myeloid"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==10,"CellMajor"]="Myeloid"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==21,"CellMajor"]="Myeloid"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==18,"CellMajor"]="Myeloid"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==7,"CellMajor"]="Fibroblast"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==20,"CellMajor"]="Fibroblast"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==6,"CellMajor"]="Endothelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==1,"CellMajor"]="Epithelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==2,"CellMajor"]="Epithelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==8,"CellMajor"]="Epithelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==17,"CellMajor"]="Epithelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==19,"CellMajor"]="Epithelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==22,"CellMajor"]="Epithelial"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==12,"CellMajor"]="PT"
scRNA_harmony@meta.data[scRNA_harmony@meta.data$seurat_clusters==14,"CellMajor"]="Mast"
table(scRNA_harmony@meta.data$CellMajor,scRNA_harmony@meta.data$seurat_clusters)

saveRDS(scRNA_harmony, "kidney_Harmony_CellMajor.rds")

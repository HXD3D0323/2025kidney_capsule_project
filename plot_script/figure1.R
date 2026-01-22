
library(Seurat)
library(scRNAtoolVis)
library(tidydr)
setwd("/xtdisk/tianchx_group/panyt/18.KIRC_Capsule/Final_Figs")
scRNA_harmony<-readRDS("/xtdisk/tianchx_group/panyt/18.KIRC_Capsule/03.scRNA_reClu/PRAD_Harmony_CellMajor.rds")

mycols <- c("T"="#E95D53","NK"="#F9F871","B/Plasma"="#ECB37B",
"Myeloid"="#D3ACAE","Mast"="#D1D198","Fibroblast"="#4FFBDF",
"Endothelial"="#3EAC71","Epithelial"="#5867AF","PT"="#C5E7A4")

scRNA_harmony_used<-subset(scRNA_harmony,Region %in% c("Capsule"))
dim(scRNA_harmony_used)
Idents(object = scRNA_harmony_used) <- "CellMajor"
		
pdf(file="01.TSNE_Capsule_CellMajor_new.pdf",width=5,height=5)
DimPlot(scRNA_harmony_used, 
        reduction = "tsne",
        pt.size = 0.001,
		label.size = 5,
		cols =  mycols,
        label = T)+theme_dr()+ theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+theme(legend.position = "none")
dev.off()

scRNA_harmony_used<-subset(scRNA_harmony,Region %in% c("Cancer"))
dim(scRNA_harmony_used)
Idents(object = scRNA_harmony_used) <- "CellMajor"
		
pdf(file="01.TSNE_Cancer_CellMajor_new.pdf",width=5,height=5)
DimPlot(scRNA_harmony_used, 
        reduction = "tsne",
        pt.size = 0.001,
		label.size = 5,
		cols =  mycols,
        label = T)+theme_dr()+ theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+theme(legend.position = "none")
dev.off()

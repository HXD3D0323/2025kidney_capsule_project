## figure 1B
library(Seurat)
library(scRNAtoolVis)
library(tidydr)
setwd("/xtdisk/tianchx_group/panyt/18.KIRC_Capsule/Final_Figs")
scRNA_harmony<-readRDS("/xtdisk/tianchx_group/panyt/18.KIRC_Capsule/03.scRNA_reClu/KIRC_Harmony_CellMajor.rds")

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

## figure 1C

from os.path import join
import scanpy as sc
import pandas as pd
from matplotlib import rcParams
sc.set_figure_params(dpi=80, color_map='viridis')
sc.settings.verbosity = 2
sc.logging.print_versions()
file_path = "/xtdisk/tianchx_group/panyt/18.KIRC_Capsule/Final_Figs/Figures/Total_Cell/seurat"
ann = sc.read_10x_mtx(file_path)
medata = pd.read_csv(join(file_path, "metadata.csv"), index_col = 0)
ann.obs = medata
ann.var.index.name = 'index'
cluster_color = ['#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863','#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180']
#Lymphocytes	"CD3D", "CD3E"
#Myeloid	"CD68", "CD14"
#B	"CD79A","CD79B"
#Epithelial	"EPCAM", "KRT19"
#Fibroblast	"DCN","COL1A1"
#Endothelial	"VWF","PECAM1"
#Plasma "IGHG1", "JCHAIN"           
marker_genes = ["IGHA1","MS4A1","CD79A",
"PLVAP","PECAM1","VWF",
"KRT18","KRT8","EPCAM",
"RGS5","TAGLN","ACTA2",
"TPSAB1","CPA3","KIT",
"C1QA","CD14","LYZ",
"GNLY","KLRD1","KLRF1",
"GPX3","GATM","PDZK1IP1",
"CD3D","CD3E","TRAC"] 

###########################################################################
ax = sc.pl.dotplot(ann, var_names = marker_genes, groupby='bulk_labels')
ax = sc.pl.dotplot(ann,marker_genes, groupby='bulk_labels',dot_max=0.5, dot_min=0.3, standard_scale='var')
ax = sc.pl.dotplot(ann,marker_genes, groupby='bulk_labels',dot_max=0.5, dot_min=0.3, standard_scale='var',use_raw=False, save="Total.pdf")
###########################################################################
# do Matrix Plots
gs = sc.pl.matrixplot(ann,marker_genes, groupby='bulk_labels')
gs = sc.pl.matrixplot(ann,marker_genes, groupby='bulk_labels')
gs = sc.pl.matrixplot(ann,marker_genes, groupby='bulk_labels',standard_scale='var',use_raw=False, save="Total.pdf")






















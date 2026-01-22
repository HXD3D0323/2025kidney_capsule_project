##--- figure 1B
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

##--- figure 1C
# (python)
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

##--- figure 1D
# proportion in tumor sample
library(harmony)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(cowplot)
library(data.table)
library(plyr)

data<-read.table("F:/23.肾包膜项目/02.重新聚类/Cell_Number.csv",sep=",",header=T,check.names=F)

data<-subset(data,Region %in% c("Cancer"))

data$Cell <- factor(data$Cell, levels = c("Mast","B/Plasma","PT","Fibroblast","Endothelial","NK","Myeloid","T","Epithelial"))
forestCol=c("T"="#E95D53","NK"="#F9F871","B/Plasma"="#ECB37B",
"Myeloid"="#D3ACAE","Mast"="#D1D198","Endothelial"="#3EAC71","Fibroblast"="#4FFBDF","PT"="#C5E7A4",
"Epithelial"="#5867AF")

data$Proportion<-0

for(i in 1:nrow(data)){
	data[i,"Proportion"]=data[i,"Number"]/sum(data$Number)
}
data$Proportion<-data$Proportion*100
ggplot(data=data,mapping=aes(x=Cell,y=Proportion,fill=Cell))+
theme_bw()+#theme(axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
labs(x="",y="Cell proportion (%)")+scale_y_continuous(limits = c(0, 40),breaks = seq(0, 40,10))+
scale_fill_manual(values=forestCol)+coord_flip()+
#coord_polar()+这是变成旋转的图
 geom_bar(stat="identity")+theme(legend.position="none")+theme(panel.grid.major = element_blank())+theme_classic()+
 theme( axis.text.y = element_text(size = 18,colour = 'black'),#这是更改了y轴的标识的大小
        axis.text.x = element_text(size = 18,colour = 'black'),
		axis.title.x = element_text(size = 15,colour = 'black'),
        axis.title.y = element_text(size = 15,colour = 'black'),
		axis.line = element_line(size = 1.0)
		)+theme (legend.position = 'none')
		
ggsave("02.Region_Cancer_Proportion.pdf",device = "pdf",width = 5.5,height =4.5)

# proportion in capsule sample
library(harmony)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(cowplot)
library(data.table)
library(plyr)

data<-read.table("F:/23.肾包膜项目/02.重新聚类/Cell_Number.csv",sep=",",header=T,check.names=F)

data<-subset(data,Region %in% c("Capsule"))

data$Cell <- factor(data$Cell, levels = c("PT","B/Plasma","Mast","Endothelial","Fibroblast","NK","Myeloid","Epithelial","T"))
forestCol=c("T"="#E95D53","NK"="#F9F871","B/Plasma"="#ECB37B",
"Myeloid"="#D3ACAE","Mast"="#D1D198","Endothelial"="#3EAC71","Fibroblast"="#4FFBDF","PT"="#C5E7A4",
"Epithelial"="#5867AF")

data$Proportion<-0

for(i in 1:nrow(data)){
	data[i,"Proportion"]=data[i,"Number"]/sum(data$Number)
}
data$Proportion<-data$Proportion*100
ggplot(data=data,mapping=aes(x=Cell,y=Proportion,fill=Cell))+
theme_bw()+#theme(axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
labs(x="",y="Cell proportion (%)")+scale_y_continuous(limits = c(0, 50),breaks = seq(0, 50,10))+
scale_fill_manual(values=forestCol)+coord_flip()+
#coord_polar()+这是变成旋转的图
 geom_bar(stat="identity")+theme(legend.position="none")+theme(panel.grid.major = element_blank())+theme_classic()+
 theme( axis.text.y = element_text(size = 18,colour = 'black'),#这是更改了y轴的标识的大小
        axis.text.x = element_text(size = 18,colour = 'black'),
		axis.title.x = element_text(size = 15,colour = 'black'),
        axis.title.y = element_text(size = 15,colour = 'black'),
		axis.line = element_line(size = 1.0))+theme (legend.position = 'none')
		
ggsave("02.Region_Capsule_Proportion.pdf",device = "pdf",width = 5.5,height =4.5)

##--- figure 1E



















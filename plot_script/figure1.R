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
#coord_polar()+
 geom_bar(stat="identity")+theme(legend.position="none")+theme(panel.grid.major = element_blank())+theme_classic()+
 theme( axis.text.y = element_text(size = 18,colour = 'black'),
        axis.text.x = element_text(size = 18,colour = 'black'),
		axis.title.x = element_text(size = 15,colour = 'black'),
        axis.title.y = element_text(size = 15,colour = 'black'),
		axis.line = element_line(size = 1.0))+theme (legend.position = 'none')
		
ggsave("02.Region_Capsule_Proportion.pdf",device = "pdf",width = 5.5,height =4.5)

##--- figure 1E
# here we take the YF-FG84-T-1-0803-0804_add sample as an example in the following code
# region define---
barcode_pos_file = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/0.Raw_data/YF-FG84-T-1-0803-0804_add/BSTViewer_project/subdata/L5_heAuto/barcodes_pos.tsv.gz'
png_path = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/0.Raw_data/YF-FG84-T-1-0803-0804_add/BSTViewer_project/he_roi_small.png'
sp_data_FilePath = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/Spatial_result/YF-FG84-T-1-0803-0804_add/1.Cluster/level5/'
# Capsule_Normal_FilePath = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/out/YF-FG84-T-1-0803-0804_add/subdata/L5_capsulenormal/' #use pyt result
# Capsule_Tumor_FilePath = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/out/YF-FG84-T-1-0803-0804_add/subdata/L5_capsuletumor/' #use pyt result
Capsule_Normal_FilePath = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/out/YF-FG84-T-1-0803-0804_addme/subdata/L5_capsulenormal/' 
Capsule_Tumor_FilePath = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/out/YF-FG84-T-1-0803-0804_addme/subdata/L5_capsuletumor/' 
point_size = 0.56
alpha = 0.6
out_path = '/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/59.空间转录组定义区域3/YF-FG84-T-1-0803-0804_add/'

.libPaths(c("/xtdisk/tianchx_group/gongjy/0.tools/my_R_packages/R4.3.3_mypkg", .libPaths()))
library(Seurat)
library(tidyverse)
library(ggplot2, lib.loc = "/xtdisk/tianchx_group/gongjy/0.tools/my_R_packages/R4.3.3_mypkg") 
library(patchwork)
library(reshape2)
library(ggdark)
library(cluster)
library(png)
library(ggpubr)
library(magrittr)
library(dplyr)
library(Matrix)
library(gridExtra)
library(readr)
library(config)
library(spacexr)
library(future)
library(future.apply)
# conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ r-ggdark
# conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ r-config
# conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/biocon
# da/ r-spacexr

# matrix transfer
Rcpp::sourceCpp(code=' 
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp, NumericVector cp, NumericVector z, int nrows, int ncols){
int k = z.size() ;
IntegerMatrix mat(nrows, ncols);
for (int i = 0; i < k; i++){ mat(rp[i],cp[i]) = z[i];
}
return mat;
}
') 

as_matrix <- function(mat){
	row_pos <- mat@i 
	col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
	tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x, nrows = mat@Dim[1], ncols = mat@Dim[2]) 
	row.names(tmp) <- mat@Dimnames[[1]]
	colnames(tmp) <- mat@Dimnames[[2]] 
	return(tmp)
}

sp_data<-readRDS(paste0(sp_data_FilePath,"object.rds"))
Capsule_Normal_spots<-Read10X(Capsule_Normal_FilePath, cell.column = 1) 
Capsule_Normal_data <- CreateSeuratObject(counts = Capsule_Normal_spots,assay = "Spatial") 
Capsule_Tumor_spots<-Read10X(Capsule_Tumor_FilePath, cell.column = 1) 
Capsule_Tumor_data <- CreateSeuratObject(counts = Capsule_Tumor_spots,assay = "Spatial") 
dim(Capsule_Normal_data)

length(intersect(rownames(Capsule_Tumor_data@meta.data),rownames(Capsule_Normal_data@meta.data)))
Capsule_Region<-intersect(rownames(Capsule_Tumor_data@meta.data),rownames(Capsule_Normal_data@meta.data))
sp_data@meta.data$Region<-"Normal"

sp_data@meta.data[Capsule_Region,"Region"]<-"Capsule"
Tumor_Region<-setdiff(rownames(Capsule_Tumor_data@meta.data),Capsule_Region)
sp_data@meta.data[Tumor_Region,"Region"]<-"Tumor"
table(sp_data@meta.data$Region)

sp_data@meta.data$Barcode<-rownames(sp_data@meta.data)
# read image
cal_zoom_rate <- function(width, height){
std_width = 1000
std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
if (std_width / std_height > width / height){
scale = width / std_width
	}else{
		scale = height / std_height
	}
	return(scale)
}
png <- readPNG(png_path) 

zoom_scale = cal_zoom_rate(dim(png)[2], dim(png)[1]) 
barcode_pos <- read.table(gzfile(barcode_pos_file), header = F)
#barcode_pos <- rename(barcode_pos, Barcode = "V1", pos_w = "V2", pos_h = "V3")
colnames(barcode_pos)<-c("Barcode", "pos_w", "pos_h")
barcode_pos <- mutate(barcode_pos, across(c(pos_w, pos_h), ~ .x * zoom_scale))

Region<-sp_data@meta.data[,c("Barcode","Region")]
barcode_pos_Region <- left_join(Region,barcode_pos,  by = 'Barcode')
barcode_pos_Region[1:4,]	
saveRDS(sp_data,paste0(out_path, 'Region.RDS')) 

p = ggplot(barcode_pos_Region, aes(x = pos_w, y = dim(png)[1] - pos_h))+
		background_image(png)+
		geom_point(shape = 16, size = point_size, alpha = alpha,  aes(color = Region))+
		coord_cartesian(xlim = c(0, dim(png)[2]), y = c(0, dim(png)[1]), expand = FALSE)+
		scale_color_manual(values = c("Capsule"=alpha("#7FABD1", c(1,0.3))[1],
		"Normal"="#89aa7b","Tumor"="#EC6E66"))+
		theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
		guides(color = guide_legend(override.aes = list(size = 2.5, alphe = 0.1)))+theme(legend.position = "none")+
		theme(plot.margin = margin(),axis.ticks.length = unit(0, "pt"))
ggsave(p, file = paste0(out_path, 'Region_Proportion2.png'), width = dim(png)[2]/200, height = dim(png)[1]/200, dpi = 300) 
ggsave(p, file = paste0(out_path, 'Region_Proportion2.pdf'), width = dim(png)[2]/200, height = dim(png)[1]/200, dpi = 300) 

# major cell type deconvolution (SPOTlight)---
# input sample: YF-FG84-T-1-0803-0804_add, YF-T-16-FG73-0710-0711_add, YF-T-17-FG96-20230824-0825_add
.libPaths(c("/xtdisk/tianchx_group/gongjy/0.tools/my_R_packages/R4.3.3_mypkg", .libPaths()))
library(ggplot2, lib.loc = "/xtdisk/tianchx_group/gongjy/0.tools/my_R_packages/R4.3.3_mypkg") 
library(Seurat)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(reshape2)
library(ggdark)
library(cluster)
library(png)
library(magrittr)
library(dplyr)
library(BiocGenerics)
library(SPOTlight, lib.loc = "/xtdisk/tianchx_group/gongjy/0.tools/anaconda3/envs/R4.3.3/lib/R/library") # version ‘0.1.7’

point_size = 0.58
alpha = 0.6

# samples=samples[!samples %in% c("YF-FG84-T-1-0803-0804_add")]
samples=c("YF-FG84-T-1-0803-0804_add", "YF-T-16-FG73-0710-0711_add", "YF-T-17-FG96-20230824-0825_add")
for(j in samples){
    print(paste("start process sample：", j))
    FilePath = paste0('/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/Spatial_result/', j, '/1.Cluster/level5/')
    barcode_pos_file = paste0('/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/0.Raw_data/', j, '/BSTViewer_project/subdata/L5_heAuto/barcodes_pos.tsv.gz')
    if (j %in% c("YF-T-16-FG72-0707-0710", "YF-T-3-FG71-0707-0710"))
        {png_path = paste0('/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/0.Raw_data/', j, '/BSTViewer_project/he_rou_small.png')} else {
            png_path = paste0('/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/0.Raw_data/', j, '/BSTViewer_project/he_roi_small.png')}
    out_path = paste0('/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/Spatial_result/', j, '/2.SPOTlight/level_5/cellproportion3/') ### out

	cal_zoom_rate <- function(width, height){
	std_width = 1000
	std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
	if (std_width / std_height > width / height){
		scale = width / std_width
		} else{
		scale = height / std_height
		}
		return(scale)
	}
    
	png <- readPNG(png_path) 
	zoom_scale = cal_zoom_rate(dim(png)[2], dim(png)[1]) 
	barcode_pos <- read.table(gzfile(barcode_pos_file), header = F)
	#barcode_pos <- rename(barcode_pos, Barcode = "V1", pos_w = "V2", pos_h = "V3")
	colnames(barcode_pos)<-c("Barcode", "pos_w", "pos_h")
	barcode_pos <- mutate(barcode_pos, across(c(pos_w, pos_h), ~ .x * zoom_scale))
	# decon_mtrx$Barcode<-rownames(decon_mtrx)
	# decon_mtrx<-decon_mtrx[,c("Barcode",colnames(decon_mtrx)[1:(ncol(decon_mtrx)-1)])]
	# write.table(decon_mtrx,file = paste0(out_path, 'spotlight_majortype3_0825.csv'),sep=",",row.names=F,quote=F) 
   
    decon_mtrx <- read.csv(paste0(out_path, 'spotlight_majortype0825.csv'),
                        header = TRUE, 
                        sep = ",",
                        stringsAsFactors = FALSE)

	barcode_pos_cluster <- left_join(decon_mtrx,barcode_pos,  by = 'Barcode')

    ## P: major cell type show in plot
    barcode_pos_cluster$`Cell type` <- barcode_pos_cluster$Order1
    barcode_pos_cluster <- barcode_pos_cluster %>%
    mutate(`Cell type` = ifelse(`Cell type` == "Neoplastic", "Malignant\ncells", `Cell type`))

    barcode_pos_cluster$`Cell type` <- factor(barcode_pos_cluster$`Cell type`, 
                                              levels = c(
                                                        "B.Plasma",
                                                        "Endothelial",
                                                        "Epithelial",
                                                        "Fibroblast",
                                                        "Mast",
                                                        "Myeloid",
                                                        "Malignant\ncells",
                                                        "NK",
                                                        "PT",
                                                        "T"
                                                        ))
	barcode_pos_cluster[1:4,]
	col =   c("#7DB954",
            "#96C5D7",
            "#CAA7DD",
            "#2874A6", # "#4682B4"
            "#FA8072",
            "#a15891",
            "#DD3F4E",
            "#64864A",
            "#ff9d5c",
            "#019477"
            )

	p = ggplot(barcode_pos_cluster, aes(x = pos_w, y = dim(png)[1] - pos_h))+
	background_image(png)+
	geom_point(shape = 16, size = point_size, aes(color = `Cell type`))+
	coord_cartesian(xlim = c(0, dim(png)[2]), y = c(0, dim(png)[1]), expand = FALSE)+
	scale_color_manual(values = col)+
	theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
	guides(color = guide_legend(override.aes = list(size = 2.5, alphe = 0.1)))
    ggsave(p, file = paste0(out_path, j, "_", 'SPOTlight_major_0828.png'), width=8, height=7, dpi = 300) 
    ggsave(p, file = paste0(out_path, j, "_",'SPOTlight_major_0828.pdf'), width=8, height=7, dpi = 300) 
    print(paste("P draw plot ending：", j))
}	


##--- figure 1F
setwd("/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/86.计算各个区域的各类细胞的表达差异/")
.libPaths(c("/xtdisk/tianchx_group/gongjy/0.tools/my_R_packages/R4.3.3_mypkg", .libPaths()))

samples<- system("ls -l /xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/59.空间转录组定义区域3 | awk '{print $9}'", intern = TRUE)
samples<-samples[-1]
print(samples)

outaTab<-data.frame()

for(i in samples){
	sample<-readRDS(paste0("/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/59.空间转录组定义区域3/",i,"/Region.RDS"))
	sample@meta.data$Barcode<-rownames(sample@meta.data)
	meta_sample<-sample@meta.data[,c("Barcode","Region")]
	sample_Spotlight<-read.table(paste0("/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/Spatial_result/",i,"/2.SPOTlight/level_5/cellproportion3/spotlight.csv"),sep=",",header=T,check.names=F)
	sample_data<-merge(meta_sample,sample_Spotlight,by="Barcode")
	sample_data$Sample=i
	outaTab<-rbind(outaTab,sample_data)
}

table(outaTab$Sample,outaTab$Region)
outaTab<-outaTab[,c("Barcode","Region","B.Plasma","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","Neoplastic","PT","T","NK","Sample")]
cells<-c("B.Plasma","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","Neoplastic","PT","T","NK")
outTab2<-data.frame()

for(i in levels(factor(outaTab$Sample))){
	for(j in levels(factor(outaTab$Region))){
		for(k in cells){
			outTab2=rbind(outTab2,cbind(Sample=i,Region=j,Cell=k,Ratio=mean(as.numeric(outaTab[outaTab$Sample==i & outaTab$Region==j,k]))))
			}
	}
}

write.table(outTab2,"Ratio_region0824.csv",sep=",",col.names=T,row.names=F,quote=F)

library(ggplot2)
library(tidyverse)
library(ggsignif)
setwd("/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/86.计算各个区域的各类细胞的表达差异/")

data <- read.csv("Ratio_region0824.csv", sep=",", header=TRUE, check.names=FALSE)
data <- data[,-1]
head(data)

data$group<-paste0(data$Cell,"_",data$Region)
colnames(data)<-c("group2","cancer", "values","group")
data<-as.data.frame(data)
levels(factor(data$group))
data$group <- factor(data$group, levels = c("B.Plasma_Normal","B.Plasma_Capsule","B.Plasma_Tumor",
											"T_Normal","T_Capsule","T_Tumor",
											"NK_Normal","NK_Capsule","NK_Tumor",
											"Mast_Normal","Mast_Capsule","Mast_Tumor",
											"Myeloid_Normal","Myeloid_Capsule","Myeloid_Tumor",
											"Endothelial_Normal","Endothelial_Capsule","Endothelial_Tumor",
											"Epithelial_Normal","Epithelial_Capsule","Epithelial_Tumor",
											"Fibroblast_Normal","Fibroblast_Capsule","Fibroblast_Tumor",
											"Neoplastic_Normal","Neoplastic_Capsule","Neoplastic_Tumor",
											"PT_Normal","PT_Capsule","PT_Tumor"))
data<-as.data.frame(data)
p <- ggplot(data)+
  geom_boxplot(aes(group, values, fill = group2, color = group2), 
               # 间距调整：
               width = 0.5,
               outlier.shape = NA
               ) +  labs(y = "Cell ratio in each spots",x="")+
  # 中位数线：
  #geom_line(data = tmp_data, aes(x_value, med_value, group = group), color = "#ffffff")+
  # 顶部灰色方块：
  geom_rect(aes(xmin=0, xmax=31, ymin=0.72, ymax = 0.8), fill="#eaeae0")+
  # 灰色竖线：
  geom_vline(xintercept = c(0.5+seq(3,27,by=3)), color = "#bcbdbf", alpha = 0.8)+
  # x轴标签：
  scale_x_discrete(labels = rep(c("N", "C","T"), 10))+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,0.7,.1),limits = c(0,0.8))+ # 数据范围到 0.8，只显示到 0.7 的刻度
  # 颜色：
#   scale_fill_manual(values = c("#00BFFF","#FFD700",  "#27408B"))+
#   scale_color_manual(values = c("#00BFFF","#FFD700",  "#27408B"))+
#   scale_fill_manual(values = c("#2e8cf0ff", "#46e931ff", "#F28E2BFF"))+
#   scale_color_manual(values = c("#2e8cf0ff", "#46e931ff", "#F28E2BFF"))+
  scale_fill_manual(values = c("#7AC3DF", "#70CDBE", "#EB7E60"))+
  scale_color_manual(values = c("#7AC3DF", "#70CDBE", "#EB7E60"))+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # ---- 坐标轴样式 ----
        axis.text.x = element_text(size = 12, color = "black"),   # 横轴刻度字号
        axis.text.y = element_text(size = 12, color = "black"),   # 纵轴刻度字号
        axis.title.y = element_text(size = 14, color = "black"),  # 纵轴标题
        axis.line = element_line(color = "black"),                # 坐标轴线黑色
        axis.ticks = element_line(color = "black"),               # 坐标轴刻度黑色
        axis.ticks.length = unit(0.18, "cm")                       # 刻度线长度，默认大概是0.15cm
        )+
  # 添加p值：
  geom_signif(aes(group, values),
              comparisons = list(c("B.Plasma_Capsule","B.Plasma_Tumor"),
								c("B.Plasma_Normal","B.Plasma_Tumor"),
								c("T_Normal","T_Capsule"),
								c("T_Capsule","T_Tumor"),
								c("T_Normal","T_Tumor"),
								c("Mast_Capsule","Mast_Tumor"),
								c("Mast_Normal","Mast_Tumor"),
								c("Myeloid_Normal","Myeloid_Capsule"),
								#c("Myeloid_Capsule","Myeloid_Tumor"),
								c("Myeloid_Normal","Myeloid_Tumor"),
								c("Endothelial_Normal","Endothelial_Capsule"),
								c("Endothelial_Capsule","Endothelial_Tumor"),
								c("Endothelial_Normal","Endothelial_Tumor"),
								c("Epithelial_Normal","Epithelial_Capsule"),
								c("Epithelial_Capsule","Epithelial_Tumor"),
								c("Epithelial_Normal","Epithelial_Tumor"),
								c("Fibroblast_Normal","Fibroblast_Capsule"),
								c("Fibroblast_Capsule","Fibroblast_Tumor"),
								c("Fibroblast_Normal","Fibroblast_Tumor"),
								c("Neoplastic_Normal","Neoplastic_Capsule"),
								c("Neoplastic_Capsule","Neoplastic_Tumor"),
								c("Neoplastic_Normal","Neoplastic_Tumor"),
                                c("NK_Capsule","NK_Tumor"), #新增的显著性
								c("PT_Normal","PT_Capsule"),
								c("PT_Capsule","PT_Tumor"),
								c("PT_Normal","PT_Tumor")),
              vjust = 1.7,  # 距离调小，文字会更靠近横线。需设大于1，1为重合
              tip_length = rep(0,8), 
              y_position = c(0.15, 0.20, # B.Plasma
							 0.28, 0.33, 0.38, # T
							 0.06, 0.11, # Mast
							 0.13, 0.18, # Myeloid
							 0.12, 0.17, 0.22, # Endothelial
							 0.56, 0.61, 0.66, # Epithelial
							 0.41, 0.46, 0.51, # Fibroblast
							 0.25, 0.3, 0.35, # Neoplastic
                             0.15, # NK
							 0.12, 0.17, 0.22), # PT
              annotations = c("***", "**", # B.Plasma
							  "*", "***", "***", # T
							  "***", "***", # Mast
							  "*", "***", # Myeloid
							  "**","***","***", # Endothelial
							  "***", "***", "***", # Epithelial
							  "***", "***", "***", # Fibroblast
							  "*", "***", "***", # Neoplastic 
                              "*", # NK
							  "*", "***", "***") # PT
							  )+
  annotate("text", 
            x = seq(2, 31, 3), 
            y = 0.76, 
            size = 3.5, 
            lineheight = 0.8, 
            label = c("B/Plasma",
                        "T",
                        "NK",
                        "Mast",
                        "Myeloid",
                        "Endothelial",
                        "Epithelial",
                        "Fibroblast",
                        "Malignant\ncells", 
                        "PT"))
ggsave("sumcellratio_boxplot0824.png", height = 3.5, width = 9, dpi = 300) 
ggsave('sumcellratio_boxplot0824.pdf', height = 3.5, width = 9, dpi = 300) 


##--- figure 1H
# bulk proteomics pathway analysis
# tumor group GO analysis---
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#install.packages("msigdbr")
#install.packages("shadowtext")
#install.packages("ggstance")
library(clusterProfiler)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(plyr)
library(dplyr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)
setwd("F:/23.肾包膜项目/109.全蛋白质谱的差异分析/11.GO与KEGG功能富集分析的barplot")

data<-read.table("data.txt",header=T,check.names=F,sep="\t")
#data_plot2<-subset(data_plot2,ID %in% pathways)
data[1:4,]
data_plot2<-subset(data,Type %in% c("Tumor"))
data_plot2<-data_plot2[,c("Description","p.adjust")]
colnames(data_plot2)<-c("ID","pvalue")

data_plot2$ID <- factor(data_plot2$ID, levels = data_plot2$ID)
max(data_plot2$pvalue)
min(data_plot2$pvalue)
#data_plot2$label<-"Lump"
data_plot2$label<-c("Yes","Yes","Yes","Yes","Yes")
data_plot2 <- data_plot2[order(data_plot2[, "pvalue"]), ]
data_plot2
data_plot2$ID <- factor(data_plot2$ID, levels = rev(data_plot2$ID))
data_plot2$pvalue<-as.numeric(data_plot2$pvalue)
max(-log2(data_plot2$pvalue))

data_plot2$P_used<-(-log2(data_plot2$pvalue))

ggplot(data_plot2, aes(x = ID, y = P_used, fill = label)) + 
		coord_flip() + 
        #geom_col() +
		geom_col(width = 0.8) +  # 调整条形图宽度
		#scale_fill_manual(values = c(alpha("black", c(1,0.3))[1],alpha("#87CEEB", c(1,0.3))[1]))+
		scale_fill_manual(values = c(alpha("#87CEEB", c(1,0.3))[1]))+
        scale_y_continuous(breaks=seq(0,9,3),limits = c(0,9)) + #两侧留空
        theme_classic() +
		geom_text(data = data_plot2,
            aes(x=ID, y= .1, label= paste0(" ", ID), color = label),#bar跟坐标轴间留出间隙
            size = 6, #字的大小
            hjust = "inward" ) +  #字的对齐方式
		scale_color_manual(values = c("black"))+
		xlab("")+
        ylab("-Log2 adjP value") + 
		geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "grey")+
        theme(axis.title.x = element_text(size = 15)) + 
		theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        theme(axis.text.x = element_text(size = 15))+theme(legend.position = 'none')+
		theme(axis.line = element_line(size = 0.8)) # 设置坐标轴的线宽
		
    #labs(title = "KEGG: Lump malignant cell vs \nno invaded malignant cell")+  # 添加标题
    #theme(plot.title = element_text(size = 20))  # 调整标题字体大小
ggsave("GO_proteomic_Tumor.pdf", width = 4.8, height = 4)

# capsule group GO analysis---
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#install.packages("msigdbr")
#install.packages("shadowtext")
#install.packages("ggstance")
library(clusterProfiler)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(plyr)
library(dplyr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
setwd("F:/23.肾包膜项目/109.全蛋白质谱的差异分析/11.GO与KEGG功能富集分析的barplot")

data<-read.table("data.txt",header=T,check.names=F,sep="\t")
#data_plot2<-subset(data_plot2,ID %in% pathways)
data[1:4,]
data_plot2<-subset(data,Type %in% c("Capsule"))
data_plot2<-data_plot2[,c("Description","p.adjust")]
colnames(data_plot2)<-c("ID","pvalue")

data_plot2$ID <- factor(data_plot2$ID, levels = data_plot2$ID)
max(data_plot2$pvalue)
min(data_plot2$pvalue)
#data_plot2$label<-"Lump"
data_plot2$label<-c("Yes","Yes","Yes","Yes","Yes")
data_plot2 <- data_plot2[order(data_plot2[, "pvalue"]), ]
data_plot2
data_plot2$ID <- factor(data_plot2$ID, levels = rev(data_plot2$ID))
data_plot2$pvalue<-as.numeric(data_plot2$pvalue)
max(-log2(data_plot2$pvalue))

data_plot2$P_used<-(-log2(data_plot2$pvalue))

ggplot(data_plot2, aes(x = ID, y = P_used, fill = label)) + 
		coord_flip() + 
        #geom_col() +
		geom_col(width = 0.8) +  # 调整条形图宽度
		#scale_fill_manual(values = c(alpha("black", c(1,0.3))[1],alpha("#87CEEB", c(1,0.3))[1]))+
		scale_fill_manual(values = c(alpha("#87CEEB", c(1,0.3))[1]))+
        scale_y_continuous(breaks=seq(0,32,8),limits = c(0,32)) + #两侧留空
        theme_classic() +
		geom_text(data = data_plot2,
            aes(x=ID, y= .1, label= paste0(" ", ID), color = label),#bar跟坐标轴间留出间隙
            size = 6, #字的大小
            hjust = "inward" ) +  #字的对齐方式
		scale_color_manual(values = c("black"))+
		xlab("")+
        ylab("-Log2 adjP value") + 
		geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "grey")+
        theme(axis.title.x = element_text(size = 15)) + 
		theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        theme(axis.text.x = element_text(size = 15))+theme(legend.position = 'none')+
		theme(axis.line = element_line(size = 0.8)) # 设置坐标轴的线宽
		
    #labs(title = "KEGG: Lump malignant cell vs \nno invaded malignant cell")+  # 添加标题
    #theme(plot.title = element_text(size = 20))  # 调整标题字体大小
ggsave("GO_proteomic_Capsule.pdf", width = 4.8, height = 4)


##--- figure 1I
setwd("F:/23.肾包膜项目/54.蛋白组学的差异/3.差异蛋白的热图")
ECM<-read.table("F:\\2.研究生课题\\2.重新整合\\26.利用Ecotype的细胞进行gsva\\matrisome_hs_masterlist.csv",sep=",",header=T,check.names=F)
ECM[1:4,1:4]
colnames(ECM)[3]="ID"
ECM<-ECM[,c("ID","Category","Division")]
length_col<-dim(ECM)[2]
rownames(ECM)<-ECM$ID
CellMajor<-read.table("F:/23.肾包膜项目/55.细胞大类的各类基因的差异/Total_CellMajor.csv",sep=",",header=T,check.names=F)
CellMajor<-subset(CellMajor, gene %in% ECM$ID)
levels(factor(CellMajor$cluster))

cell<-subset(CellMajor, cluster %in% c("B/Plasma"))
ECM[,"B_Plasma"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"B_Plasma"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("Endothelial"))
ECM[,"Endothelial"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Endothelial"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("Epithelial"))
ECM[,"Epithelial"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Epithelial"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("Fibroblast"))
ECM[,"Fibroblast"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Fibroblast"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("Mast"))
ECM[,"Mast"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Mast"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("Myeloid"))
ECM[,"Myeloid"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Myeloid"]=cell[i,"avg_log2FC"]
}

ECM[,"Neoplastic"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Neoplastic"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("NK"))
ECM[,"NK"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"NK"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("PT"))
ECM[,"PT"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"PT"]=cell[i,"avg_log2FC"]
}

cell<-subset(CellMajor, cluster %in% c("T"))
ECM[,"T"]=0
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"T"]=cell[i,"avg_log2FC"]
}
Genetype<-ECM[,1:3]
rownames(ECM)<-ECM$ID
ECM<-ECM[,4:ncol(ECM)]

# filter
#ECM_LogFC<-subset(ECM,abs(Num)>0)
ECM_LogFC<-ECM
ECM_LogFC<-cbind(ID=rownames(ECM_LogFC),ECM_LogFC)

##########################################################################
setwd("F:/23.肾包膜项目/54.蛋白组学的差异/3.差异蛋白的热图")
ECM<-read.table("F:\\2.研究生课题\\2.重新整合\\26.利用Ecotype的细胞进行gsva\\matrisome_hs_masterlist.csv",sep=",",header=T,check.names=F)
ECM[1:4,1:4]
colnames(ECM)[3]="ID"
ECM<-ECM[,c("ID","Category","Division")]
length_col<-dim(ECM)[2]
rownames(ECM)<-ECM$ID
CellMajor<-read.table("F:/23.肾包膜项目/55.细胞大类的各类基因的差异/Total_CellMajor.csv",sep=",",header=T,check.names=F)
CellMajor<-subset(CellMajor, gene %in% ECM$ID)
levels(factor(CellMajor$cluster))

cell<-subset(CellMajor, cluster %in% c("B/Plasma"))
ECM[,"B_Plasma"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"B_Plasma"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("Endothelial"))
ECM[,"Endothelial"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Endothelial"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("Epithelial"))
ECM[,"Epithelial"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Epithelial"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("Fibroblast"))
ECM[,"Fibroblast"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Fibroblast"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("Mast"))
ECM[,"Mast"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Mast"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("Myeloid"))
ECM[,"Myeloid"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Myeloid"]=cell[i,"p_val_adj"]
}

ECM[,"Neoplastic"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"Neoplastic"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("NK"))
ECM[,"NK"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"NK"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("PT"))
ECM[,"PT"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"PT"]=cell[i,"p_val_adj"]
}

cell<-subset(CellMajor, cluster %in% c("T"))
ECM[,"T"]=1
for(i in 1:nrow(cell)){
	ECM[cell[i,"gene"],"T"]=cell[i,"p_val_adj"]
}
#Genetype<-ECM[,1:3]
rownames(ECM)<-ECM$ID
ECM<-ECM[,4:ncol(ECM)]

# filter
#ECM_LogFC<-subset(ECM,abs(Num)>0)
ECM_Pvalue<-ECM
ECM_Pvalue<-cbind(ID=rownames(ECM_Pvalue),ECM_Pvalue)
#ECM_Pvalue<-subset(ECM,abs(Num)>0)
#ECM_Pvalue<-cbind(ID=rownames(ECM_Pvalue),ECM_Pvalue)
write.table(ECM_Pvalue,"ECM_Pvalue.csv",quote=F,sep=",",row.names=F,col.names=T)

ECM_LogFC<-subset(ECM_LogFC, ID %in% ECM_Pvalue$ID)
write.table(ECM_LogFC,"ECM_LogFC.csv",quote=F,sep=",",row.names=F,col.names=T)
Genetype<-subset(Genetype, ID %in% ECM_Pvalue$ID)
write.table(Genetype,"GenoType.csv",quote=F,sep=",",row.names=F,col.names=T)

############################################################################################
library(ComplexHeatmap) 
library(circlize) 
#library(ChAMPdata) 
library(data.table)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
setwd("F:/23.肾包膜项目/54.蛋白组学的差异/3.差异蛋白的热图")
scRNA_harmony<-readRDS("F:/23.肾包膜项目/19.NMF鉴定恶性肿瘤细胞亚群/Total_Malignant.rds")
table(scRNA_harmony@meta.data$CellMajor2)
table(scRNA_harmony@meta.data$CellMin2)
Deseq2_used<-read.table("../1.Limma做差异/校正后的差异/ECM_diff_DESeq2_used.xls",sep="\t",header=T,check.names=F)
#scRNA_harmony<-subset(scRNA_harmony,CellMin2 %in% c("iCAF","myCAF","imPVL","dPVL"))
Deseq2_used$Region<-ifelse(Deseq2_used$log2FoldChange>0,"Tumor","Capsule")
#scRNA_harmony<-subset(scRNA_harmony,Region %in% c("Capsule"))

Idents(object = scRNA_harmony) <- "CellMajor2"
heatmap_gene <- c("DCN","COL4A4","CLEC3B","COL14A1","THSD4","MGP","TNXB",
"COL16A1","FRZB","F9","EGLN1","SERPINH1","AGRN","ANXA4","COL21A1","COL5A3")
heatmap.BlWtRd <- c("#1F66AC", "grey90", "#B2192B")
#celltype <- c("Memory","Effector","Exhausted")
expMat <- AverageExpression(scRNA_harmony, assays = "RNA", features = heatmap_gene,verbose = TRUE) %>% .$RNA
expMat

P_Value<-read.table("F:/23.肾包膜项目/55.细胞大类的各类基因的差异/ECM_Pvalue.csv",sep=",",check.names=F,header=T)
ECM_label<-subset(P_Value, ID %in% rownames(expMat))
genetype<-read.table("F:/23.肾包膜项目/55.细胞大类的各类基因的差异/GenoType.csv",sep=",",row.names=1,header=T,check.names=F)
rownames(ECM_label)<-ECM_label$ID
ECM_label<-ECM_label[,-1]

for(i in 1:nrow(ECM_label)){
	for(j in 1:ncol(ECM_label)){
		if(as.numeric(ECM_label[i,j])< 0.05){
			ECM_label[i,j]="*"
		}else{
				ECM_label[i,j]=" "
			}
	}
}

cell_fun <- function(ECM_label,  
                     darkcol = "black", lightcol = "white", digit = 2, fontsize  = 8){
    function(j, i, x, y, width, height, fill){
        if(ECM_label[i,j] == 0){
            grid.text("a", x, y, 
                      gp = gpar(fontsize = 10, col  = darkcol))
        }else{
				if(ECM_label[i,j] == "NaN"){
					grid.text("NA", x, y, 
                    gp = gpar(fontsize = fontsize, col  = darkcol))
				}else{
					if(ECM_label[i,j] == "Inf"){
						grid.text("b", x, y, 
						gp = gpar(fontsize = 10, col  = darkcol))
					}else{
						if(ECM_label[i,j] == "*"){
							grid.text("*", x, y, 
							gp = gpar(fontsize = 12, col  = darkcol))
						}else{
							if(ECM_label[i,j] == " "){
								grid.text(" ", x, y, 
								gp = gpar(fontsize = fontsize, col  = lightcol))
							}
						}
					}
				}
			}
		}
    }

#####################################annRow##########################################
immunomodulator<-genetype[rownames(expMat),]
rownames(Deseq2_used)<-Deseq2_used$ID
immunomodulator<-cbind(immunomodulator,Region=Deseq2_used[rownames(expMat),"Region"])
annRow <- immunomodulator

unique(immunomodulator$Category)
unique(immunomodulator$Division)

annRow<-annRow[order(rownames(annRow)),]

annRow<-annRow[order(annRow[,1]),]
annRow<-annRow[order(annRow[,2]),]
annRow$Category <- factor(annRow$Category, levels = c("ECM Glycoproteins","ECM Regulators",
		"Secreted Factors","Proteoglycans","ECM-affiliated Proteins","Collagens")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$Division <- factor(annRow$Division, levels = c("Core matrisome","Matrisome-associated"))
annRow$Region <- factor(annRow$Region, levels = c("Capsule","Tumor"))
annRowColors <- list("Division" = c("Core matrisome"="#FFA500","Matrisome-associated"="#5D478B"),
			"Category"= c("ECM Glycoproteins"="#ED1450", "ECM Regulators"="#FCCA02", 
			"Secreted Factors" = "#A7CE35", "Proteoglycans" = "#2C92DA", 
			"ECM-affiliated Proteins" = "#228B22", "Collagens" = "#FFF8AD"))
left_anno <- HeatmapAnnotation(df                   = data.frame(Division = annRow$Division,Category = annRow$Category),
                               which                = "row", 
                               gp                   = gpar(col = "grey80"),
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"), 
                               show_annotation_name = F,
                               border               = F)

ECM_label<-ECM_label[rownames(annRow),]
colnames(ECM_label)[1]<-"B/Plasma"
ECM_label<-ECM_label[,colnames(expMat)]

expMat2 <- t(apply(expMat, 1, function(row) {
  min_val <- min(row)
  max_val <- max(row)
  # 校正公式：newValue = (oldValue - minVal) * (newMax - newMin) / (maxVal - minVal) + newMin
  corrected_row <- (row - min_val) * (2 - (-2)) / (max_val - min_val) + (-2)
  return(corrected_row)
}))

col_expr <- colorRamp2(c(min(na.omit(expMat2)),0,max(na.omit(expMat2))), heatmap.BlWtRd) 

ECM.expr <- Heatmap(matrix             = as.matrix(expMat2[rownames(annRow),]),
                   col                = col_expr,
                   border             = NA, 
                   rect_gp = gpar(col = "grey80"), 
                   cluster_rows       = F, 
                   cluster_columns    = F, 
                   show_row_names     = T,
                   row_names_side     = "left", 
				   cell_fun = cell_fun(ECM_label), 
                   row_names_gp       = gpar(fontsize = 10), 
                   show_column_names  = T, 
                   column_names_side  = "bottom", 
                   row_split          = annRow$Region,
				   column_names_rot = 45,
                   #top_annotation     = top_anno, 
                   left_annotation    = left_anno, 
                   name               = "Expression", 
                   width              = ncol(expMat) * unit(6, "mm"),
                   height             = nrow(expMat) * unit(5.5, "mm"))
# do plot
pdf(file = "Protein.pdf", width = 6,height = 7)
draw(ECM.expr,
	 heatmap_legend_side = "right", 
     annotation_legend_side = "right") 
	 #heatmap_legend_side = "bottom", 
     #annotation_legend_side = "bottom") 
invisible(dev.off())

#expMat<-cbind(ID=rownames(expMat),expMat)
#write.table(expMat,"ECM_Log2FC_PanCancer.csv",row.names=F,col.names=T,sep=",",quote=F)



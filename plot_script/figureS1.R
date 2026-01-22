##--- figure S1A
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^7)
setwd("F:\\23.肾包膜项目\\02.重新聚类")
library(harmony)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(cowplot)
library(data.table)
library(plyr)
scRNA_harmony <- readRDS("KIRC_Harmony_CellMajor.rds")

rt1<-scRNA_harmony@meta.data[,c("Sample","Region","CellMajor")]
rt1$Label<-paste0(rt1$Sample,"_",rt1$Region)
rt1<-rt1[,c("Label","CellMajor")]

table(rt1$CellMajor)
table(rt1$Label)
rt1$CellMajor <- factor(rt1$CellMajor, levels = c("T","NK","B/Plasma","Myeloid","Mast","Endothelial","Fibroblast","PT","Epithelial"))
rt1[rt1$Label=="Sample1_Cancer","Label2"]="T1"
rt1[rt1$Label=="Sample2_Cancer","Label2"]="T2"
rt1[rt1$Label=="Sample3_Cancer","Label2"]="T3"
rt1[rt1$Label=="Sample4_Cancer","Label2"]="T4"
rt1[rt1$Label=="Sample1_Capsule","Label2"]="C1"
rt1[rt1$Label=="Sample2_Capsule","Label2"]="C2"
rt1[rt1$Label=="Sample3_Capsule","Label2"]="C3"
rt1[rt1$Label=="Sample4_Capsule","Label2"]="C4"

rt1$Label2 <- factor(rt1$Label2, levels = rev(c("C1","T1","C2","T2","C3","T3","C4","T4")))
rt1<-rt1[,c("Label2","CellMajor")]
df=as.data.frame(table(rt1))
chisq=chisq.test(table(rt1))
pValue=chisq$p.value
if(pValue<0.001){
	pValue="P < 0.001"
	}else{
		pValue=paste0("P = ",sprintf("%.03f",pValue))
	}
	
df=ddply(df, .(Label2), transform, percent = Freq/sum(Freq) * 100)
col = RColorBrewer::brewer.pal(n = 11, name = 'Paired')
#names(col) = c('Frame_Shift_Del','Missense_Cell', 'Nonsense_Cell', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Cell','Translation_Start_Site','Multi_Hit')
	
df=ddply(df, .(Label2), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label2=paste0(sprintf("%.0f", df$percent), "%")
	
# do plot
p=ggplot(df, aes(x = factor(Label2), y = percent, fill = CellMajor)) +
    geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
	scale_fill_manual(values=c("T"="#E95D53","NK"="#F9F871","B/Plasma"="#ECB37B",
"Myeloid"="#D3ACAE","Mast"="#D1D198","Endothelial"="#3EAC71","Fibroblast"="#4FFBDF","PT"="#C5E7A4",
"Epithelial"="#5867AF"))+
	ggtitle(paste0("")) +coord_flip()+
	#xlab("")+ ylab("Percent")+  #guides(fill=guide_legend(title=""))+
	xlab("")+ ylab("")+
	#geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
	theme_bw()+#+theme(axis.text.x = element_text(angle = 45, hjust = 1))
	theme(panel.border = element_blank(),
		panel.grid.major = element_blank())+theme(legend.position = "none")
pdf(file=paste0("CellMajor_Ratio",".pdf"), width=6, height=4)
print(p)
dev.off()


##--- figure S1B
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#install.packages("epitools")
library(epitools)
library(ComplexHeatmap)
library(circlize)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
setwd("F:/23.肾包膜项目/整理所有的图")
scRNA_harmony <- readRDS("F:/23.肾包膜项目/11.总结细胞亚类/Total_Sub_combined.rds")
table(scRNA_harmony@meta.data$CellMajor)
cells<-c("B/Plasma","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","NK","PT","T")
scRNA_harmony<-subset(scRNA_harmony,CellMajor %in% cells)
observe.data<-as.matrix(as.data.frame.matrix(table(scRNA_harmony@meta.data$CellMajor, scRNA_harmony@meta.data$Region)))
expected.data <- expected(observe.data)
plot.data <- observe.data/expected.data
#plot.data[plot.data>3] = 3

col.order <- c("Cancer", "Capsule") 
plot.data <- plot.data[, col.order] 
#colnames(plot.data)<-c("T","C")

# filter: odd ratio > 1
capsule <- plot.data[plot.data[, "Capsule"] > 1, ]
capsule <- capsule[order(capsule[, "Capsule"]), ]
cancer <- plot.data[plot.data[, "Cancer"] > 1, ]
cancer <- cancer[order(cancer[, "Cancer"]), ]

plot.data <- rbind(capsule, cancer)
colnames(plot.data)<-c("Tumor","Capsule")
#plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] 
#rownames(plot.data) <- substr(rownames(plot.data), 5, 6)
ECM_label<-plot.data

cols <- setNames(object = c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
                 nm = c("1", "1.5", "3", ">3"))
cell_fun <- function(ECM_label,  
                     darkcol = "black", lightcol = "white", fontsize  = 8){
    function(j, i, x, y, width, height, fill){
        if(ECM_label[i,j] < as.numeric(names(cols)[1])){
            grid.text("\u00B1", x, y, 
                      gp = gpar(fontsize = 15, col  = darkcol))
        }else{
				if(ECM_label[i,j] < as.numeric(names(cols)[2])){
					grid.text("+", x, y, 
                    gp = gpar(fontsize = 15, col  = darkcol))
				}else{
					if(ECM_label[i,j] < as.numeric(names(cols)[3])){
						grid.text("++", x, y, 
						gp = gpar(fontsize = 15, col  = darkcol))
					}else{
						if(ECM_label[i,j] < as.numeric(names(cols)[4])){
							grid.text("+++", x, y, 
							gp = gpar(fontsize = 15, col  = darkcol))
						}else{
								grid.text("++++", x, y, 
								gp = gpar(fontsize = fontsize, col  = darkcol))
						}
					}
				}
			}
		}
    }

#heatmap.BlWtRd <- c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D")
library(ComplexHeatmap)
library(circlize)
heatmap.BlWtRd <- c("#317cb7", "white", "#b72230")
col_expr <- colorRamp2(c(min(na.omit(plot.data)),1,max(na.omit(plot.data))), heatmap.BlWtRd)

# do plot
hm = Heatmap(plot.data, cell_fun = cell_fun(ECM_label),
			 col = col_expr,rect_gp = gpar(col = "grey80"),
             cluster_rows = F, cluster_columns = F,column_names_rot = 0, #肿瘤类型呈45度
             width = unit(4, "cm"), height = unit(7.5, "cm"), name  = "R[o/e]", # 热图颜色图例的名称
             show_heatmap_legend = T)
#lgd = Legend(labels = names(cols), title = expression(R[o/e]), legend_gp = gpar(fill = cols))
draw(hm)

pdf("04.Total_RatioHeatmap.pdf", width = 4, height = 4)
#draw(hm, annotation_legend_list = lgd)
draw(hm)
invisible(dev.off())


##--- figure S1C


















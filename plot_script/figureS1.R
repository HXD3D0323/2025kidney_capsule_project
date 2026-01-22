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
	
#计算高低风险组中突变基因百分率
df=ddply(df, .(Label2), transform, percent = Freq/sum(Freq) * 100)

col = RColorBrewer::brewer.pal(n = 11, name = 'Paired')
#names(col) = c('Frame_Shift_Del','Missense_Cell', 'Nonsense_Cell', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Cell','Translation_Start_Site','Multi_Hit')
	
#百分比位置
df=ddply(df, .(Label2), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label2=paste0(sprintf("%.0f", df$percent), "%")
	
#绘制图形
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

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
             cluster_rows = F, cluster_columns = F,column_names_rot = 0,
             width = unit(4, "cm"), height = unit(7.5, "cm"), name  = "R[o/e]", 
             show_heatmap_legend = T)
#lgd = Legend(labels = names(cols), title = expression(R[o/e]), legend_gp = gpar(fill = cols))
draw(hm)

pdf("04.Total_RatioHeatmap.pdf", width = 4, height = 4)
#draw(hm, annotation_legend_list = lgd)
draw(hm)
invisible(dev.off())


##--- figure S1C



##--- figure S1D



##--- figure S1F
library(harmony)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(cowplot)
library(data.table)
library(plyr)

setwd("/xtdisk/tianchx_group/gongjy/geyuze/19.KIRC_Capsule_ST/86.计算各个区域的各类细胞的表达差异")
data <- read.csv("Ratio_region0826.csv")
dim(data)
data <- na.omit(data)
dim(data)
outTab<-data.frame()
for(i in levels(factor(data$Region))){
	for(j in levels(factor(data$Cell))){
		outTab<-rbind(outTab,cbind(Region=i,Cell=j,Mean=mean(as.numeric(data[data$Region==i & data$Cell==j,"Ratio"]))))
	}
}
outTab$Mean<-as.numeric(outTab$Mean)*100
write.table(outTab,"Ratio_region_mean.csv",sep=",",col.names=T,row.names=F,quote=F)

##----- draw capsule sample ratio plot
data <- read.table("Ratio_region_mean.csv", sep=",", header=T, check.names=F)
data <- subset(data, Region %in% c("Capsule"))

# descent order by ratio
data <- data %>% arrange(Mean) #desc(Mean)
data$Cell <- factor(data$Cell, levels = data$Cell)

# set color
forestCol = c(
  "T"="#E95D53","NK"="#F9F871","B.Plasma"="#ECB37B",
  "Myeloid"="#D3ACAE","Mast"="#D1D198","Endothelial"="#3EAC71",
  "Fibroblast"="#4FFBDF","PT"="#C5E7A4",
  "Epithelial"="#5867AF","Neoplastic"="#465886"
)

# do plot
p <- ggplot(data, aes(x=Cell, y=Mean, fill=Cell)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=forestCol) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x="", y="Mean cell proportion in Capsule (%)") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 18, colour = 'black'),
    axis.text.x = element_text(size = 18, colour = 'black'),
    axis.title.x = element_text(size = 15, colour = 'black'),
    axis.title.y = element_text(size = 15, colour = 'black'),
    axis.line = element_line(size = 1.0),
    legend.position = "none"
  )

ggsave("02.Region_Capsule_Proportion.pdf", plot=p, device="pdf", width=5.5, height=4.5)

##----- draw tymor sample ratio plot
data <- read.table("Ratio_region_mean.csv", sep=",", header=T, check.names=F)
data <- subset(data, Region %in% c("Tumor"))

# descent order by ratio
data <- data %>% arrange(Mean) #desc(Mean)
data$Cell <- factor(data$Cell, levels = data$Cell)

# set color
forestCol = c(
  "T"="#E95D53","NK"="#F9F871","B.Plasma"="#ECB37B",
  "Myeloid"="#D3ACAE","Mast"="#D1D198","Endothelial"="#3EAC71",
  "Fibroblast"="#4FFBDF","PT"="#C5E7A4",
  "Epithelial"="#5867AF","Neoplastic"="#465886"
)

# do plot
p <- ggplot(data, aes(x=Cell, y=Mean, fill=Cell)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=forestCol) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x="", y="Mean cell proportion in Tumor (%)") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 18, colour = 'black'),
    axis.text.x = element_text(size = 18, colour = 'black'),
    axis.title.x = element_text(size = 15, colour = 'black'),
    axis.title.y = element_text(size = 15, colour = 'black'),
    axis.line = element_line(size = 1.0),
    legend.position = "none"
  )

ggsave("02.Region_Tumor_Proportion.pdf", plot=p, device="pdf", width=5.5, height=4.5)

##----- draw normal sample ratio plot
data <- read.table("Ratio_region_mean.csv", sep=",", header=T, check.names=F)
data <- subset(data, Region %in% c("Normal"))

# descent order by ratio
data <- data %>% arrange(Mean) #desc(Mean)
data$Cell <- factor(data$Cell, levels = data$Cell)

# set color
forestCol = c(
  "T"="#E95D53","NK"="#F9F871","B.Plasma"="#ECB37B",
  "Myeloid"="#D3ACAE","Mast"="#D1D198","Endothelial"="#3EAC71",
  "Fibroblast"="#4FFBDF","PT"="#C5E7A4",
  "Epithelial"="#5867AF","Neoplastic"="#465886"
)

# do plot
p <- ggplot(data, aes(x=Cell, y=Mean, fill=Cell)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=forestCol) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 10)) +
  labs(x="", y="Mean cell proportion in Normal (%)") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 18, colour = 'black'),
    axis.text.x = element_text(size = 18, colour = 'black'),
    axis.title.x = element_text(size = 15, colour = 'black'),
    axis.title.y = element_text(size = 15, colour = 'black'),
    axis.line = element_line(size = 1.0),
    legend.position = "none"
  )

ggsave("02.Region_Normal_Proportion.pdf", plot=p, device="pdf", width=5.5, height=4.5)


##--- figure S1H
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("plyr")
library(ggplot2)
library(dplyr)
library(plyr)
setwd("F:/23.肾包膜项目/54.蛋白组学的差异/5.蛋白组做主成分分析")

expr_df <- read.csv(file='easy_input_expr.csv',row.names = 1, 
                    header = TRUE, sep=",", stringsAsFactors = FALSE)
ECM <-	read.table(file='F:/23.肾包膜项目/54.蛋白组学的差异/1.Limma做差异/校正后的差异/ECM_diff_DESeq2_used.xls',sep="\t",row.names = 1, 
                    header = TRUE)				
meta_df <- read.csv(file='easy_input_meta.csv', row.names = 1,
                    header = TRUE, sep=",",stringsAsFactors = FALSE)
expr_df[1:3,1:4]
expr_df<-expr_df[,rownames(ECM)]
head(meta_df, n=3)
pca.results <- prcomp(expr_df, center = TRUE, scale. = FALSE)

mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
	
pca.rotation <- pca.results$rotation
pca.rotation
pca.pv <- summary(pca.results)$importance[2,]
pca.pv
low_dim_df <- as.data.frame(pca.results$x[,c(1,2)])
low_dim_df$group <- meta_df$group
low_dim_df[1:3,]
add_ellipase <- function(p, x="PC1", y="PC2", group="group",
                         ellipase_pro = 0.95,
                         linetype="dashed",
                         colour = "black",
                         lwd = 2,...){
  obs <- p$data[,c(x, y, group)]
  colnames(obs) <- c("x", "y", "group")
  ellipse_pro <- ellipase_pro
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(obs, 'group', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$x, x$y))
    mu <- c(mean(x$x), mean(x$y))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
    })
  names(ell)[2:3] <- c('x', 'y')
  
  ell <- ddply(ell, .(group) , function(x) x[chull(x$x, x$y), ])
  p <- p + geom_polygon(data = ell, aes(x=x,y=y,group = group), 
                   colour = colour,
                   alpha = 1,fill = NA,
                   linetype=linetype,
                   lwd =lwd)
  return(p)
}
pc1.pv <- paste0(round(pca.pv['PC1'],digits = 3) * 100, "%")
pc2.pv <- paste0(round(pca.pv['PC2'],digits = 3) * 100, "%")

p <- ggplot(low_dim_df) + 
  geom_point(aes(x=PC1, y=PC2, color=group), size=2, #点的大小
             shape=20,
             alpha=0.5) +
  scale_color_manual(values = mycol[1:length(unique(meta_df$group))]) +
  #scale_colour_hue(l=45) + 
  theme_bw() +
  theme(panel.grid =element_blank()) + 
  annotate("text",x=(-1.5e+07),y=-1.2,label = "Capsule",color = mycol[1]) +
  annotate("text",x=1.35e+07,y=5.6,label = "Tumor",color = mycol[2]) +
  #annotate("text",x=5,y=-2.5,label = "OV",color = mycol[3]) +
  guides(color=guide_legend(title = NULL)) +
  theme(legend.background = element_blank(), 
        legend.position = c(0,1),legend.justification = c(0,1),
        legend.text = element_text(size=12)) + #字体大小
  xlab(paste0("PC1 ( ", pc1.pv," variance )")) + 
  ylab(paste0("PC2 ( ", pc2.pv," variance )")) 
p
p1 <- add_ellipase(p,ellipase_pro = 0.95,colour = "dimgrey",linetype=2,lwd=1)
p1
ggsave('PCA_DIY1.pdf',width = 4.7,height = 4)


##--- figure S1I



##--- figure S1J
library(ggplot2)
library(ggrepel)
library(dplyr)
setwd("F:/23.肾包膜项目/54.蛋白组学的差异/8.做二维")
rt=read.delim("Protein.txt",sep="\t",header=T,check.names=F)
rt<-rt[,c("PG.Genes","[1] C5521.raw.PG.Quantity","[2] C7412.raw.PG.Quantity","[3] C7441.raw.PG.Quantity",
"[4] C7469.raw.PG.Quantity","[5] C7968.raw.PG.Quantity","[6] T5521.raw.PG.Quantity","[7] T7412.raw.PG.Quantity",
"[8] T7441.raw.PG.Quantity","[9] T7469.raw.PG.Quantity","[10] T7968.raw.PG.Quantity")]
colnames(rt)<-c("ID",paste0("C",rep(1:5)),paste0("T",rep(1:5)))
rt[is.na(rt)] <- 0
			   
ECM<-read.table("F:\\2.研究生课题\\2.重新整合\\26.利用Ecotype的细胞进行gsva\\matrisome_hs_masterlist.csv",sep=",",header=T,check.names=F)
ECM[1:4,1:4]
colnames(ECM)[3]="ID"
ECM<-ECM[,c("ID","Category","Division")]

for(i in 1:nrow(rt)){
	rt[i,"ID"] <- strsplit(rt[i,"ID"], ";")[[1]][1]
}
rt<-subset(rt,ID %in% ECM$ID)
dim(rt)
rt <- rt[complete.cases(rt), ]
dim(rt)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(rt)
norm_level<-mean(colSums(rt))
rt <- apply(rt, 2, function(x) x / sum(x))
rt[1:4,1:4]
#rt<-rt*norm_level
#rt[1:4,1:4]
#rt<-round(rt)
rt<-cbind(rt,Capsule=0,Tumor=0)
for(i in 1:nrow(rt)){
	rt[i,"Capsule"]<-mean(as.numeric(rt[i,1:5]))
	rt[i,"Tumor"]<-mean(as.numeric(rt[i,6:10]))
}

dat2<-rt[,c("Capsule","Tumor")]
dat2<-dat2*100
dat2<-cbind(ID=rownames(dat2),dat2)
			
ECM<-read.table("F:\\2.研究生课题\\2.重新整合\\26.利用Ecotype的细胞进行gsva\\matrisome_hs_masterlist.csv",sep=",",header=T,check.names=F)
ECM[1:4,1:4]
colnames(ECM)[3]="ID"
ECM<-ECM[,c("ID","Category","Division")]
			
dat2<-merge(ECM,dat2,by="ID")
head(dat2)
colnames(dat2)[1]<-"Gene"
write.table(dat2,"data_used.csv",sep=",",col.names=T,row.names=F,quote=F)
dat2[1:4,]
dat2$Category<-factor(dat2$Category,levels=c("Collagens","ECM Glycoproteins","Proteoglycans","ECM Regulators","ECM-affiliated Proteins","Secreted Factors"))
dim(dat2)
			
per2<-subset(dat2,Capsule > 0 & Tumor > 0)
dim(per2)
max(-log2(as.numeric(per2$Capsule)))
max(-log2(as.numeric(per2$Tumor)))
per2$Tumor<-as.numeric(per2$Tumor)
per2$Capsule<-as.numeric(per2$Capsule)
#per2<-as.data.frame(per2)
per2<-as.data.frame(per2)
rownames(per2)<-per2$Gene
table(per2$Category)
write.table(per2,"data_used2.csv",sep=",",col.names=T,row.names=F,quote=F)
			
# per2<-read.table("data_used2.csv",sep=",",header=T,check.names=F)
rownames(per2)<-per2$Gene
per2$Tumor_Capsule<-per2$Tumor/per2$Capsule
per2$Capsule_Tumor<-per2$Capsule/per2$Tumor
# do plot			
p<-ggplot(per2, aes(x=-log2(Tumor),y=-log2(Capsule)))+theme_bw()+
	geom_point(aes(fill=Category,shape=Category,color=Category))+ 
	labs(x="-log2 Protein abundance in Tumor (mean%)",y="-log2 Protein abundance in Capsule (mean%)")+
	theme(axis.title = element_text(size = 16),
		legend.title = element_blank(),
		legend.position = c(0.22,.83),
		legend.background = element_rect(fill = "transparent"),
		panel.grid=element_blank(),
		legend.text = element_text(size=12,face="bold"))+
	scale_shape_manual(values = c(21:25,21))+
	scale_color_manual(values = c("Collagens"='#cc0000',"ECM-affiliated Proteins"="#B8860B","ECM Glycoproteins"="#1c4587",
	"ECM Regulators"="#2C92DA","Proteoglycans"="#228B22","Secreted Factors"="#6529c0"))+
	scale_fill_manual(values = c("Collagens"='#cc0000',"ECM-affiliated Proteins"="#B8860B","ECM Glycoproteins"="#1c4587",
	"ECM Regulators"="#2C92DA","Proteoglycans"="#228B22","Secreted Factors"="#6529c0"))+
	scale_y_continuous(breaks=seq(0,25,5))+
	scale_x_continuous(breaks=seq(0,25,5))+
	geom_blank(aes(y = 15))+
	geom_blank(aes(x = 15))+
	#theme(legend.text = element_text(size=20,face="bold"))+#改变legend的大小
	#slope控制斜率，intercept控制截距
	geom_abline(slope = 1,intercept = 0,lwd=0.8,linetype="dashed",color="grey")+
	geom_abline(slope = 1,intercept = log2(6),lwd=0.5,linetype="dashed",color="grey")+
	geom_abline(slope = 1,intercept = -log2(6),lwd=0.5,linetype="dashed",color="grey")+
	geom_text_repel(
	data = per2[-log2(per2$Tumor_Capsule) >log2(6) | -log2(per2$Capsule_Tumor)>log2(6),],
	aes(label = per2[-log2(per2$Tumor_Capsule) >log2(6) | -log2(per2$Capsule_Tumor)>log2(6),1]),
	size = 5,colour="black",box.padding=unit(0.5,"lines"),
	point.padding=unit(0.8,"lines"),
	show.legend=F,max.overlaps=19,
	force=1)+
	guides(shape = guide_legend(override.aes = list(size = 4)))  # Add this line to adjust legend symbols
	
ggsave("Tumor vs Capsule.pdf",p,width = 6, height = 6)














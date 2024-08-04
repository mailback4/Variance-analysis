#####加载R包
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
cran_packages <- c('tidyr',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   'factoextra',
                   'FactoMineR',
                   'devtools',
                   'cowplot',
                   'patchwork',
                   'basetheme',
                   'paletteer',
                   'AnnoProbe',
                   'ggthemes',
                   'VennDiagram',
                   'tinyarray') 
Biocductor_packages <- c('GEOquery',
                         'hgu133plus2.db',
                         'ggnewscale',
                         "limma",
                         "impute",
                         "GSEABase",
                         "GSVA",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "preprocessCore",
                         "enrichplot",
                         "ggplotify")

for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}

rm(list = ls())


#####获取GEO数据
library(GEOquery)
gse_number = "GSE26459"
eSet26459 <- getGEO(gse_number, destdir = '.', AnnotGPL = F,getGPL = F)
class(eSet26459)
length(eSet26459)
eSet26459 = eSet26459[[1]]

#(1)提取表达矩阵exp
exp26459 <- exprs(eSet26459)
dim(exp26459)
exp26459[1:4,1:4]
#检查矩阵是否正常，如果是空的就会报错，空的和有负值的、有异常值的矩阵需要处理原始数据。
#自行判断是否需要log
ex <- exp26459
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exp26459 <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
boxplot(exp26459)
# 将原始数据标准化（针对画出的图水平不一致）
exp26459 =normalizeBetweenArrays(exp26459)
boxplot(exp26459)

#(2)提取临床信息
pd26459 <- pData(eSet26459)
#(3)让exp列名与pd的行名顺序完全一致
p = identical(rownames(pd26459),colnames(exp26459));p
if(!p) exp26459 = exp26459[,match(rownames(pd26459),colnames(exp26459))]
#(4)提取芯片平台编号
gpl_number <- eSet26459@annotation;gpl_number
save(gse_number,pd26459,exp26459,gpl_number,file = "step1output26459.Rdata")

# Group(实验分组)和ids(探针注释)
load(file = "step1output26459.Rdata")
library(stringr)
# 生成Group向量的三种常规方法，三选一，选谁就把第几个逻辑值写成T，另外两个为F。如果三种办法都不适用，可以继续往后写else if
if(F){
  # 1.Group----
  # 第一种方法，有现成的可以用来分组的列
  Group = pd$`disease state:ch1` 
}else if(F){
  # 第二种方法，自己生成
  Group = c(rep("RA",times=13),
            rep("control",times=9))
  Group = rep(c("RA","control"),times = c(13,9))
}else if(T){
  # 第三种方法，使用字符串出理的函数获取分组
  Group=ifelse(str_detect(pd26459$source_name_ch1,"Sensitive"),
               "Sensitive",
               "Resistant")
  
}

# 需要把Group转换成因子，并设置参考水平，指定levels，对照组在前，处理组在后
Group = factor(Group,levels = c("Sensitive","Resistant"))
Group

#2.探针注释的获取-----------------
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
if(F){
  #注：表格读取参数、文件列名不统一，活学活用，有的表格里没有symbol列，也有的GPL平台没有提供注释表格
  b = read.table("GPL570-55999.txt",header = T,
                 quote = "\"",sep = "\t",check.names = F)
  colnames(b)
  ids2 = b[,c("ID","Gene Symbol")]
  colnames(ids2) = c("probe_id","symbol")
  ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]
}

load(file = "step1output26459.Rdata")

# 1.PCA 图----
dat26459=as.data.frame(t(exp26459))
library(FactoMineR)
library(factoextra) 
Group
dat.pca26459 <- PCA(dat26459, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca26459,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)

pca_plot
save(pca_plot,file = "pca_plot26459.Rdata")

# 2.top 1000 sd 热图---- 
cg26459=names(tail(sort(apply(exp26459,1,sd)),10000))
n=exp26459[cg26459,]

# 直接画热图，对比不鲜明
library(pheatmap)
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(n) 
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col
)

# 按行标准化
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100)
) 
dev.off()

load(file = "step1output26459.Rdata")
#需要表达矩阵和Group，不需要改
library(limma)
Group
design=model.matrix(~0+factor(Group))
colnames(design)=levels(factor(Group))
rownames(design)=colnames(exp26459)
contrast.matrix<-makeContrasts(paste0(unique(Group),collapse = '-'),levels = design)
fit=lmFit(exp26459,design)
fit2=contrasts.fit(fit,contrast.matrix)
fit2=eBayes(fit2)
deg=topTable(fit2,coef=1,number = Inf)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
#2.加上探针注释
ids <- toTable(hgu133plus2SYMBOL)
ids = ids[!duplicated(ids$symbol),]
deg <- inner_join(deg,ids,by="probe_id")
nrow(deg)



########糖酵解
file.path1<-"C:\\Users\\lenovo\\Desktop\\Variance analysis\\糖酵解.txt"
exp_genelist1<-read.table(file.path1,header = TRUE,sep='\t');exp_genelist1[,2]
deg1=deg[match(exp_genelist1[,2],deg[,8]),]

#3.加change列,标记上下调基因
logFC_t=1
P.Value_t = 0.05
k1 = (deg1$P.Value < P.Value_t)&(-deg1$logFC < -logFC_t)
k2 = (deg1$P.Value < P.Value_t)&(-deg1$logFC > logFC_t)
deg1 <- mutate(deg1,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg1$change)
#4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg1$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
deg1 <- inner_join(deg1,s2e,by=c("symbol"="SYMBOL"))
save(Group,deg1,logFC_t,P.Value_t,gse_number,file = "step4output26459deg1.Rdata")

load(file = "step1output26459.Rdata")
load(file = "step4output26459deg1.Rdata")
#1.火山图----
library(dplyr)
library(ggplot2)
dat1  = deg1[!duplicated(deg1$symbol),]

p <- ggplot(data = dat1, 
            aes(x = -logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c('blue', "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

for_label <- dat1%>% 
  filter(symbol %in% c("HK1","BPGM"))

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot

#2.差异基因热图----

load(file = 'step1output26459.Rdata')
# 表达矩阵行名替换
exp26459 = exp26459[dat1$probe_id,]
rownames(exp26459) = dat1$symbol
if(F){
  #全部差异基因
  cg = dat1$symbol[dat1$change !="stable"]
  length(cg)
}else{
  #取前10上调和前10下调
  library(dplyr)
  dat2 = dat1 %>%
    filter(change!="stable") %>%
    arrange(logFC)
  cg = c(head(dat2$symbol,1),
         tail(dat2$symbol,1))
}
n=exp26459[cg,]
dim(n)

#差异基因热图
library(pheatmap)
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(n) 
heatmap_plot <- pheatmap(n,show_colnames =F,
                         scale = "row",
                         #cluster_cols = F, 
                         annotation_col=annotation_col,
                         color = colorRampPalette(c('blue','white','red'))(100),
                         breaks = seq(-3,3,length.out = 100)
) 
heatmap_plot
#拼图
library(patchwork)
library(ggplotify)
volcano_plot +as.ggplot(heatmap_plot)

#(1)输入数据
gene_up = deg1$ENTREZID[deg1$change == 'up'] 
gene_down = deg1$ENTREZID[deg1$change == 'down'] 
gene_diff1 = c(gene_up,gene_down)



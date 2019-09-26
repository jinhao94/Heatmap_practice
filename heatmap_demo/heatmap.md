---
title: "Complexheat_map"
author: "jinh"
date: "2019年9月26日"
output: html_document
---

# 长注释热图做法
### 加载包， 导入数据
### 数据格式打印如下，加减乘除不能做行名或者列名，所以把数据放到第一行


```r
library(printr)
library(ComplexHeatmap)
```

```
## Loading required package: grid
```

```
## ========================================
## ComplexHeatmap version 2.1.0
## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
## Github page: https://github.com/jokergoo/ComplexHeatmap
## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
## 
## If you use it in published research, please cite:
## Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
##   genomic data. Bioinformatics 2016.
## ========================================
```

```r
library(ggplotify)
library(circlize)
```

```
## ========================================
## circlize version 0.4.6
## CRAN page: https://cran.r-project.org/package=circlize
## Github page: https://github.com/jokergoo/circlize
## Documentation: http://jokergoo.github.io/circlize_book/book/
## 
## If you use it in published research, please cite:
## Gu, Z. circlize implements and enhances circular visualization 
##   in R. Bioinformatics 2014.
## ========================================
```

```r
#Data format
dd<-read.table("input_file.txt",header = T,sep = "\t")
knitr::kable(head(dd,10),digits =4)
```



|ID                                               |ID2 |     A37|     A38|     A39|     A40|     A41|     A42|     A43|     A44|      B9|     B10|     B11|     B12|     B13|     B14|     B15|
|:------------------------------------------------|:---|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|
|Acetylcarnitine                                  |s1  |  1.2538|  2.1964|  1.4251|  2.3590|  1.7572|  1.6408|  1.7745|  1.6270|  3.2280|  2.4749|  3.1506|  3.0083|  2.6531|  2.2494|  2.2943|
|Butyryl carnitine (isomer of 920); PlaSMA ID-919 |s2  |  1.0449|  0.3529|  0.1377|  1.1599|  1.1984|  0.2263|  0.6495|  1.1082|  2.2031|  1.8724|  2.1125|  2.0442|  1.7482|  1.2391|  2.0947|
|Dodecanoylcarnitine                              |s3  | -3.0112| -2.2904| -3.8954| -2.7518| -4.0316| -4.1208| -3.4053| -3.7464| -1.7686| -1.8831| -1.4768| -1.6826| -2.2131| -2.0643| -1.8410|
|Palmitoyl-L-carnitine                            |s4  | -0.9348| -0.7688| -1.9140| -0.8696| -1.1714| -1.3621| -1.0338| -1.2786|  0.1190| -0.1074|  0.2525|  0.3024| -0.0379|  0.1411| -0.0652|
|LPC 18:3; PlaSMA ID-2744                         |s5  | -0.0519|  0.5608| -0.1659|  0.5176| -0.6468| -0.4148|  0.4534| -0.1001|  1.5960|  1.3621|  1.3110|  1.7454|  2.0605|  0.6712|  0.7129|
|Hexosyl LPE 18:2; PlaSMA ID-3108                 |s6  | -2.4077| -2.5625| -2.3920| -1.8129| -2.8995| -3.1411| -2.5776| -3.0955| -0.9546| -1.2184| -0.4552|  0.2288| -0.5047| -1.4606| -1.7965|
|3-Hydroxypropanoic acid                          |s7  |  2.5304|  1.3820|  1.6771|  1.3348|  1.9632|  1.3621|  1.5910|  1.1688|  2.6637|  2.5142|  2.6416|  3.0459|  2.6345|  2.5187|  2.6605|
|Threonate                                        |s8  | -3.8552| -2.4511| -1.3975| -0.4506| -4.0397| -0.7038| -1.6417| -3.0933| -5.0684| -3.8880| -6.1825| -5.2094| -4.4661| -5.9432| -3.7890|
|O-Acetyl-L-serine                                |s9  | -3.1036| -2.8572| -2.4330| -2.0327| -2.8791| -1.8111| -2.8626| -2.7310| -4.3562| -3.4767| -4.5315| -5.0137| -3.6088| -3.9949| -3.3100|
|3-Hydroxycapric acid                             |s10 | -5.2478| -5.1319| -5.5401| -4.5807| -5.9431| -5.8706| -5.3124| -5.7603| -4.5654| -4.8522| -4.2781| -4.4977| -4.6012| -4.0542| -4.8035|
#### step1

```r
Heatmap(dd[,3:ncol(dd)])
```

```
## Warning: The input is a data frame, convert it to the matrix.
```

![plot of chunk unnamed-chunk-2](/figure/./heatmap/unnamed-chunk-2-1.png)
####step2 add annotations

```r
dd_2<-dd[,3:ncol(dd)]
row.names(dd_2)<-dd$ID2
row_labels=structure(as.character(dd$ID),names=as.character(dd$ID2))
Heatmap(dd_2,row_labels=row_labels[rownames(dd_2)], show_column_names = T, cluster_rows = F,show_row_names=T)
```

```
## Warning: The input is a data frame, convert it to the matrix.
```

![plot of chunk unnamed-chunk-3](/figure/./heatmap/unnamed-chunk-3-1.png)
####step3 modified annotations

```r
Heatmap(dd_2,row_labels=row_labels[rownames(dd_2)], show_column_names = T, cluster_rows = F,show_row_names=T,row_names_gp = gpar(fontsize = 5))
```

```
## Warning: The input is a data frame, convert it to the matrix.
```

![plot of chunk unnamed-chunk-4](/figure/./heatmap/unnamed-chunk-4-1.png)
#### set a beautiful colour

```r
Heatmap(dd_2,row_labels=row_labels[rownames(dd_2)], show_column_names = T, cluster_rows = F,show_row_names=T,row_names_gp = gpar(fontsize = 5),col = colorRamp2(c(5,0,-5,10,-15),c( "white","#fdbe85","#fd8d3c","#e6550d","#a63603")))
```

```
## Warning: The input is a data frame, convert it to the matrix.
```

![plot of chunk unnamed-chunk-5](/figure/./heatmap/unnamed-chunk-5-1.png)

***
## 参数总结
1. 要使row_labels生效，需先使用structure函数对名字进行处理
2. 使用colorRamp2之前要加载circlize
3. 这个例子也展示了长注释的写法

#复杂热图绘制
### KEGG SCFA 模块

```r
library(ggpubr)
```

```
## Loading required package: ggplot2
```

```
## Loading required package: magrittr
```

```r
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(ggplotify)
```


```r
df<-read.table("452_kegg.out",sep = "\t",header = T,row.names=1)
ko<-read.table("ko_others",sep="\t") #ko_others
df<-data.frame(t(df))
bin_info<-read.table("all_452_mag_info.txt",sep="\t",header = T)

#sub ko
sub_ko_df<-df[,colnames(df) %in% ko[,1]]
sub_ko_df_ant<-data.frame(Bin=row.names(sub_ko_df),rawrank=c(1:nrow(sub_ko_df)))
sub_ko_df_plot<-sub_ko_df
sub_ko_df_ant_info<-merge(sub_ko_df_ant,bin_info,by="Bin")

#order mapping 
sub_ko_df_ant_info_ordered<-sub_ko_df_ant_info[order(sub_ko_df_ant_info$rawrank),]
#order raw data
sub_ko_df_ant_info_ordered<-sub_ko_df_ant_info[order(sub_ko_df_ant_info$rawrank),]



c1<-c("#76c5ad","#e89f67") #Cultured uncultured
c2<-c("#6a92a7","#c0e9f2") #High abundance Low abundance
c3<-rev(c("#4B78A5","#A8514B","#8BA156","#6A6599FF","#4697AB","#CE8844","#5A8CBF","#A20056B2","#008B45B2","#BB0021B2","#54B0C5","#D5928F","#90ACD2","#F39B7FB2","#BBCE95","#AC9CC0","#84D7E1FF","#FBB98B","#C3CFE3","#E4C4C2","#D6E1C5"))
#This function design to get col set
get_col<-function(cha_grp,cc){

	cha_grp<-as.character(cha_grp)
	cha_grp<-unique(cha_grp)
	a=character()
	for (i in 1:length(cha_grp)){
	a[cha_grp[i]]=cc[i]
	}
	return(a)
}

phylum_col<-get_col(sub_ko_df_ant_info_ordered$Phylum, c3)
# Function_type_col<-get_col(sub_ko_df_ant_info_ordered$Function_type)
Cultured_uncultured<-get_col(sub_ko_df_ant_info_ordered$Culture_or_Unculture,c1)
Abundance_col<-get_col(sub_ko_df_ant_info_ordered$Abundance_v3, c2)

ha_row = rowAnnotation(df = data.frame(Phylum =sub_ko_df_ant_info_ordered$Phylum, Cultured_uncultured = sub_ko_df_ant_info_ordered$Culture_or_Unculture, Abundance=sub_ko_df_ant_info_ordered$Abundance_v3), col = list(Phylum = phylum_col, Cultured_uncultured=Cultured_uncultured, Abundance=Abundance_col), width = unit(0.5, "cm"))

colant<-data.frame(ko=colnames(sub_ko_df_plot))
colant$rank<-1:nrow(colant)
colant<-merge(colant,ko,by.x="ko",by.y="V1")
colant<-colant[order(colant$V2),]
colant$frank<-1:nrow(colant)
# colant<-colant[order(colant$rank),]

sub_ko_df_plot<-sub_ko_df_plot[,as.numeric(colant$rank)]


 # top_annotation = HeatmapAnnotation(Pathway = colant$V2)
#右边的注释可以直接加入到图中，变为一个图， rect_gp 控制线的粗细
h1=Heatmap(sub_ko_df_plot, name = "Gene number",col = colorRamp2(c(0,1,10,100),c( "white","orange","#d73027","#762a83")),cluster_rows = T,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray",lwd = 0.01,alpha=0.2), row_names_gp = gpar(fontsize = 1) ,column_names_gp = gpar(fontsize = 2),column_names_side = "top",top_annotation = HeatmapAnnotation(Pathway = colant$V2),column_split =colant$V2 ,show_row_names=T,row_split=sub_ko_df_ant_info_ordered$Sample,right_annotation=ha_row)
```

```
## Warning: The input is a data frame, convert it to the matrix.
```

```r
h1
```

![plot of chunk unnamed-chunk-7](/figure/./heatmap/unnamed-chunk-7-1.png)

## 总结
1.上图由于太密集，所以显示不全，可使用如下参数导出合适的图


```r
h1_out<-as.ggplot(h1)
ggsave(h1_out, file="KEGG_scfa_plot.pdf", height = 16,width = 20)
```

2.rect_gp 调节每个网格之间线的粗细，如果太密集可以设置alpha参数
3.top_annotation是column注释，后续可添加其它花样
4.right_annotation的用法


```r
## 绘制CAZY图
library(ggplotify)
library(ggpubr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
df<-read.table("452_cazy.out",sep = "\t",header = T,row.names=1)
df<-data.frame(t(df))
ko<-read.table("cazy_list.txt",sep="\t") #ko_others
bin_info<-read.table("all_452_mag_info.txt",sep="\t",header = T)
sub_ko_df<-df[,colnames(df) %in% ko[,1]]
sub_ko_df_ant<-data.frame(Bin=row.names(sub_ko_df),rawrank=c(1:nrow(sub_ko_df)))
sub_ko_df_plot<-sub_ko_df
sub_ko_df_ant_info<-merge(sub_ko_df_ant,bin_info,by="Bin")

#order mapping 
sub_ko_df_ant_info_ordered<-sub_ko_df_ant_info[order(sub_ko_df_ant_info$rawrank),]
#order raw data
sub_ko_df_ant_info_ordered<-sub_ko_df_ant_info[order(sub_ko_df_ant_info$rawrank),]

c1<-c("#76c5ad","#e89f67") #Cultured uncultured
c2<-c("#6a92a7","#c0e9f2") #High abundance Low abundance
c3<-rev(c("#4B78A5","#A8514B","#8BA156","#6A6599FF","#4697AB","#CE8844","#5A8CBF","#A20056B2","#008B45B2","#BB0021B2","#54B0C5","#D5928F","#90ACD2","#F39B7FB2","#BBCE95","#AC9CC0","#84D7E1FF","#FBB98B","#C3CFE3","#E4C4C2","#D6E1C5"))
# This function design to get col set
# get_col<-function(cha_grp,cc){
# 
#   cha_grp<-as.character(cha_grp)
#   cha_grp<-unique(cha_grp)
#   a=character()
#   for (i in 1:length(cha_grp)){
#     a[cha_grp[i]]=cc[i]
#   }
#   return(a)
# }

sub_ko_df_ant_info_ordered$Phylum<-factor(sub_ko_df_ant_info_ordered$Phylum, levels=c("p_Firmicutes","p_Proteobacteria","p_Bacteroidetes","p_Actinobacteria","p_Fusobacteria","p_Verrucomicrobia","p_Euryarchaeota"))

#用于第二轮颜色对齐
a=character()
cc=c("#d6ffc5","#E4C4C2","#c3cfe3","#fbb98b","#84d7e1","#ac9cc0","#bbce95")
cha_grp<-c("p_Firmicutes","p_Proteobacteria","p_Bacteroidetes","p_Actinobacteria","p_Fusobacteria","p_Verrucomicrobia","p_Euryarchaeota")
for (i in 1:length(cha_grp)){
    a[cha_grp[i]]=cc[i]
  }
phylum_col<-a


Cultured_uncultured<-c("Y"="#e89f67","N"="#76c5ad")
Abundance_col<-c("Low"="#c0e9f2","High"="#6a92a7")

colant<-data.frame(ko=colnames(sub_ko_df_plot))
colant$rank<-1:nrow(colant)
colant<-merge(colant,ko,by.x="ko",by.y="V1")
colant<-colant[order(colant$rank),]

#列名排序，很重要
colant$V2<-factor(colant$V2, levels = c("Animal Carbohydrates","Plant Cell Wall Carbohydrates","Mucin","pectin","Sucrose/Fructans","Starch","xylan","Cellulose","Other"))

sub_ko_df_ant_info_ordered$Phylum<-factor(sub_ko_df_ant_info_ordered$Phylum, levels=c("p_Firmicutes","p_Proteobacteria","p_Bacteroidetes","p_Actinobacteria","p_Fusobacteria","p_Verrucomicrobia","p_Euryarchaeota"))


ha_row2 = rowAnnotation(df = data.frame(Phylum =sub_ko_df_ant_info_ordered$Phylum, Cultured_uncultured = sub_ko_df_ant_info_ordered$Culture_or_Unculture, Abundance=sub_ko_df_ant_info_ordered$Abundance_v3), col = list(Phylum = phylum_col, Cultured_uncultured=Cultured_uncultured, Abundance=Abundance_col), width = unit(0.5, "cm"))

h2=Heatmap(sub_ko_df_plot, name = "Gene number",col = colorRamp2(c(0,1,10,100),c( "white","orange","#d73027","#762a83")),cluster_rows = T,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray",lwd = 0.01,alpha=0.2), row_names_gp = gpar(fontsize = 1) ,column_names_gp = gpar(fontsize = 2),column_names_side = "top",top_annotation = HeatmapAnnotation(Pathway = colant$V2),column_split =colant$V2 ,show_row_names=T,row_split=sub_ko_df_ant_info_ordered$Sample,right_annotation=ha_row2)
```

```
## Warning: The input is a data frame, convert it to the matrix.
```

```r
h2
```

![plot of chunk unnamed-chunk-9](/figure/./heatmap/unnamed-chunk-9-1.png)

```r
#save data
h2_out<-as.ggplot(h2)
ggsave(h2_out, file="CAZYme_plot.pdf", height = 16,width = 20)
```

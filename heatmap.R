dt = read.table("../00.data/cazy.profile.1277.float",sep="\t", header=T, row.names=1,check.names = F)
map = read.table("../../../00.data/strain.taxo.txt", sep="\t", header=T, check.names=F)
cazy_map = read.table("./cazy.map", sep="\t", row.names=1, header=T, quote = "")

#map = map[which(map$phylum == "Bacteroidetes" ),]
#map = map[which(map$merged_genus == "Prevotella" ),]
map = map[which(map$phylum != "Firmicutes" & map$complete >= 90 & map$enriched != "na"),]
map = map[which(map$complete >= 90),]
dt = dt[,map$user_genome]

dt[dt>0]=1
dt = dt[rowSums(dt)!=0 , ]

#write.table(dt, "clipboard-1290", sep="\t")
kk = heatmap(as.matrix(dt))
roworder = rownames(dt)[kk$rowInd]
colorder = colnames(dt)[kk$colInd]



library(ComplexHeatmap)
library(dplyr)
library(circlize)

dt = dt[roworder, colorder]
#dt = dt[rev(order(rowSums(dt))), rev(order(colSums(dt)))]
bb = cazy_map[rownames(dt),]

rownames(map) = map$user_genome
map = map[colnames(dt),]

rsd = data.frame(levelA=bb[,2]) #row_split_data
rsd$levelA = factor(rsd$levelA, unique(rsd$levelA))

csd = data.frame(enriched=map$enriched)
csd$enriched = ifelse(csd$enriched=="ESRD-enriched", "na", csd$enriched)
csd$enriched = factor(csd$enriched, unique(csd$enriched))

mycols = colorRamp2(breaks=c(0, 50, 51, 300),
                    colors=c("#ffffff", "#e41a1c", "#377eb8", "#4daf4a"))


dt = as.matrix(dt)
cazy_map = cazy_map[rownames(dt),]
rownames(dt) = rownames(cazy_map)
map = map[colnames(dt),]
colnames(dt) = map$merged_genus
Heatmap(dt,
        cluster_rows=FALSE, cluster_columns=FALSE,
        row_split=rsd, 
        column_split=csd,
        row_title_side='right',column_title_side = 'bottom',
        row_title_rot = 0, column_title_rot=90,
        border=T,
        show_row_names=T, show_column_names = T,
        row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
        #col=mycols
        )



# dt = 矩阵文件 （每一行为物种， 每一列为样本）
# map = 样本分组文件
# kk = heatmap（dt）
# 获得聚类的行列顺序

# row_order = rownames(dt)[kk$rowInd]
# col_order = colnames(dt)[kk$colInd]
# 对dt 重新排序  dt = dt[row_order, col_order]
# 然后就可以直接画图了


# ------分组
# 新建一个dataframe， 储存分组，顺序按照

# map = map[colnames(dt),] #将map按照dt顺序排列
# csd = data.frame(enriched=map$family)
# csd$enriched = factor(csd$enriched, unique(csd$enriched))
# 
# mycols = colorRamp2(breaks=c(0, 50, 51, 300),
#                     colors=c("#ffffff", "#e41a1c", "#377eb8", "#4daf4a"))
# 
# Heatmap(as.matrix(dt),
#         cluster_rows=FALSE, cluster_columns=FALSE, # 是否启用聚类
#         row_split=rsd, column_split=csd, # 行/列分割
#         row_title_side='right',column_title_side = 'bottom', # 文字显示方位
#         row_title_rot = 0, column_title_rot=90, # 文字旋转角度
#         border=T,
#         show_row_names=F, show_column_names = F, 
#         row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
#         #col=mycols
# )
# 



#====================
dt = read.table("../00.data/kegg.profile.01",sep="\t", header=T, row.names=1,check.names = F)
map = read.table("../../../00.data/strain.taxo.txt", sep="\t", header=T)
module_map = read.table("../../02.kegg/00.data/pathway.levelB.desc", sep="\t", row.names=1, header=T, quote = "")
map = map[which(map$phylum == "Firmicutes" & map$enriched != "na"),]
se = read.table("./select.ko.list",sep="\t", header=F)
se$V2 = factor(se$V2, levels=unique(se$V2))
#map = map[which(map$complete >= 90),]
dt = dt[se$V1,map$user_genome]
dt = dt[rowSums(dt)!=0, ]
dtt = t(dt)
names = rownames(dtt)




map = map[which(map$phylum == "Firmicutes" & map$enriched != "na"),]
map = map[which(map$complete >= 90),]
dt = dt[,map$user_genome]
dt = dt[rowSums(dt)!=0, ]

kk = heatmap(as.matrix(dt))
roworder = rownames(dt)[kk$rowInd]
colorder = colnames(dt)[kk$colInd]



library(ComplexHeatmap)
library(dplyr)
library(circlize)

dt = dt[roworder, colorder]
dt[dt>0] = 1

#dt = dt[rev(order(rowSums(dt))), rev(order(colSums(dt)))]
se$V2 = factor(se$V2, levels=unique(se$V2))
rownames(se) = se$V1
bb = se[rownames(dt),]
bb$V2 = factor(bb$V2, levels=unique(se$V2))
rownames(map) = map$user_genome
map = map[colnames(dt),]

rsd = data.frame(levelA=bb[,2]) #row_split_data
rsd$levelA = factor(rsd$levelA, unique(rsd$levelA))

csd = data.frame(enriched=map$enriched)
csd$enriched = factor(csd$enriched, unique(csd$enriched))

mycols = colorRamp2(breaks=c(0, 50, 51, 300),
                    colors=c("#ffffff", "#e41a1c", "#377eb8", "#4daf4a"))

#map = map[colnames(dt),]
#colnames(dt) = map$merged_genus
enrich = HeatmapAnnotation(bar = bb[,3], col=list(bar=c("HC"="green", "ESRD"="red")))
Heatmap(as.matrix(dt),
        cluster_rows=FALSE, cluster_columns=FALSE,
        row_split=rsd, 
        column_split=csd,
        row_title_side='right',column_title_side = 'bottom',
        row_title_rot = 0, column_title_rot=90,
        border=T,
        show_row_names=T, show_column_names = T,
        row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
        right_annotation = enrich
        #col=mycols
)


#========================
#       remove Firm
#========================

dt = read.table("../00.data/cazy.profile.1277.float",sep="\t", header=T, row.names=1,check.names = F)
map = read.table("../../../00.data/strain.taxo.txt", sep="\t", header=T, check.names=F)
cazy_map = read.table("./cazy.map", sep="\t", row.names=1, header=T, quote = "")


map = map[which(map$phylum != "Firmicutes" & map$complete >= 90 & map$enriched != "na"),]
dt = dt[,map$user_genome]
dt = dt[rowSums(dt)!=0, ]
dtt = t(dt)
names = rownames(dtt)


dt[dt>0]=1
dt = dt[rowSums(dt)!=0 & rowSums(dt)!= 66, ]
#write.table(dt, "clipboard-1290", sep="\t")

kk = heatmap(as.matrix(dt))
roworder = rownames(dt)[kk$rowInd]
colorder = colnames(dt)[kk$colInd]



library(ComplexHeatmap)
library(dplyr)
library(circlize)

dt = dt[roworder, colorder]
#dt = dt[rev(order(rowSums(dt))), rev(order(colSums(dt)))]
bb = cazy_map[rownames(dt),]

rownames(map) = map$user_genome
map = map[colnames(dt),]



csd = data.frame(enriched=map$enriched)
csd$enriched = factor(csd$enriched, unique(csd$enriched))

mycols = colorRamp2(breaks=c(0, 1),
                    colors=c("#ffffff", "#375623"))


dt = as.matrix(dt)
cazy_map = cazy_map[rownames(dt),]
rownames(dt) = rownames(cazy_map)
map = map[colnames(dt),]
colnames(dt) = map$merged_genus
Heatmap(dt,
        cluster_rows=FALSE, cluster_columns=FALSE,
        column_split=csd,
        row_title_side='right',column_title_side = 'bottom',
        row_title_rot = 0, column_title_rot=90,
        border=T,
        show_row_names=T, show_column_names = T,
        row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
        col=mycols
)



#-----------------植物糖动物糖
dt = read.table("../00.data/cazy.profile.1277.float.glucide",sep="\t", header=T, row.names=1,check.names = F)
map = read.table("../../../00.data/strain.taxo.txt", sep="\t", header=T, check.names=F)


map = map[which(map$phylum != "Firmicutes" & map$complete >= 90 & map$enriched != "na"),]
dt = dt[,map$user_genome]
dt = dt[rowSums(dt)!=0, ]
dtt = t(dt)
names = rownames(dtt)


dt[dt>0]=1
dt = dt[rowSums(dt)!=0 & rowSums(dt)!= 66, ]
#write.table(dt, "clipboard-1290", sep="\t")

kk = heatmap(as.matrix(dt))
roworder = rownames(dt)[kk$rowInd]
colorder = colnames(dt)[kk$colInd]



library(ComplexHeatmap)
library(dplyr)
library(circlize)

dt = dt[roworder, colorder]
#dt = dt[rev(order(rowSums(dt))), rev(order(colSums(dt)))]
bb = cazy_map[rownames(dt),]

rownames(map) = map$user_genome
map = map[colnames(dt),]



csd = data.frame(enriched=map$enriched)
csd$enriched = factor(csd$enriched, unique(csd$enriched))

mycols = colorRamp2(breaks=c(0, 1),
                    colors=c("#ffffff", "#375623"))


dt = as.matrix(dt)
cazy_map = cazy_map[rownames(dt),]
rownames(dt) = rownames(cazy_map)
map = map[colnames(dt),]
colnames(dt) = map$merged_genus
Heatmap(dt,
        cluster_rows=FALSE, cluster_columns=FALSE,
        column_split=csd,
        row_title_side='right',column_title_side = 'bottom',
        row_title_rot = 0, column_title_rot=90,
        border=T,
        show_row_names=T, show_column_names = T,
        row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
        col=mycols
)


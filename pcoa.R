dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)


library(vegan)
library(ggplot2)
library(ggpubr)



otu.dist = vegdist(as.matrix(t(dt)), method='bray')

pcoa = cmdscale(otu.dist, k=4, eig=T)
eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
point = pcoa$points
data = merge(point, map, by.x='row.names', by.y='sample')

xlab = paste("PCoA 1 ( ", eigs[1], "% )", sep="")
ylab = paste("PCoA 2 ( ", eigs[2], "% )", sep="")

p = ggscatter(data, x="V2", y="V3", color = 'group3',
          ellipse = T,shape='group2',
          ellipse.level = 0.8,
          star.plot = T,
          xlab=xlab,ylab=ylab)+
  scale_color_brewer(palette='Set1')+
  theme_bw()
p

#ggsave("abundance.pdf", p, height=5,width=6)




# ================================================
#                   除去CKD， 根据地区比较ESRD和HC
# ================================================

library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsignif)

dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)

map = map[which(map$group3 != "CKD"), ]
dt = dt[,map$sample]


otu.dist = vegdist(as.matrix(t(dt)), method='bray')

pcoa = cmdscale(otu.dist, k=8, eig=T)
eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
point = pcoa$points
data = merge(point, map, by.x='row.names', by.y='sample')

# 这几个距离是否有差异
t.test(data$V1 ~ data$group2)$p.value
t.test(data$V2 ~ data$group2)$p.value
t.test(data$V3 ~ data$group2)$p.value
t.test(data$V4 ~ data$group2)$p.value
t.test(data$V5 ~ data$group2)$p.value
t.test(data$V6 ~ data$group2)$p.value


xn=2
yn=3

xlab = paste("PCoA ", xn, " ( ", eigs[xn], "% )", sep="")
ylab = paste("PCoA ", yn, " ( ", eigs[yn], "% )", sep="")
title = "PCoA 	Bray-Curtis"
x = paste("V", xn, sep="")
y = paste("V", yn, sep="")
p = ggscatter(data, x=x, y=y, color = 'group1',
              ellipse = T,
              ellipse.level = 0.8,
              star.plot = T,
              ellipse.alpha = 0.1,
              title = title,
              mean.point = T,
              ellipse.border.remove = F,
              xlab=xlab,ylab=ylab)+
  scale_color_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))+
  theme_bw()
p

#s.class(point[,c(2,3)], as.factor(map$group1), col=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"), cstar=0)
rownames(data) = data$Row.names
compar = list(c("BJ_ESRD", "BJ_HC"), 
              c("SH_ESRD", "SH_HC"),
              c("BJ_ESRD", "SH_ESRD"),
              c("BJ_HC", "SH_HC")
              )

compar = list(c("BJ", "SH"))
#--------
title="V1"
p1 = ggplot(data=data, aes(y=V1, x=group2, group=group2, fill=group2))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = compar, 
              step_increase = 0.1,
              test=t.test,
              #map_signif_level=T,
              )+
  ggtitle(title)+
  scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))
  



title="V2"
p2 = ggplot(data=data, aes(y=V2, x=group1, group=group1, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = compar, 
              step_increase = 0.1,
              test=t.test,
              #map_signif_level=T,
  )+
  ggtitle(title)+
  scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))



title="V3"
p3 = ggplot(data=data, aes(y=V3, x=group1, group=group1, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = compar, 
              step_increase = 0.1,
              test=t.test,
              #map_signif_level=T,
  )+
  ggtitle(title)+
  scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))



title="V4"
p4 = ggplot(data=data, aes(y=V4, x=group1, group=group1, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = compar, 
              step_increase = 0.1,
              test=t.test,
              #map_signif_level=T,
  )+
  ggtitle(title)+
  scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))



title="V5"
p5 = ggplot(data=data, aes(y=V5, x=group1, group=group1, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = compar, 
              step_increase = 0.1,
              test=t.test,
              #map_signif_level=T,
  )+
  ggtitle(title)+
  scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))



title="V6"
p6 = ggplot(data=data, aes(y=V6, x=group1, group=group1, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = compar, 
              step_increase = 0.1,
              test=t.test,
              #map_signif_level=T,
  )+
  ggtitle(title)+
  scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))


p = p1+p2+p3+p4+p5+p6

p

ggsave("distance.boxplot.t.test.pdf", p,  width=14, height=9)





# =======================================================
#    除去CKD， 根据地区比较ESRD和HC, 查看他们的距离的差异
# =======================================================

library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsignif)

dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)

map = map[which(map$group3 != "CKD"), ]
dt = dt[,map$sample]


otu.dist = vegdist(as.matrix(t(dt)), method='bray')

library(reshape2)
bb = as.matrix(otu.dist)
dl = merge(bb)
dl = dl[duplicated(t(apply(dl, 1, sort))),]

map = map[,c(1,2)]
dm1 = merge(dl, map, by.x='Var1', by.y='sample')
dm2 = merge(dm1, map, by.x='Var2', by.y='sample')
dm2 = dm2[which(dm2$Var2 != dm2$Var1),]
dm2$group = paste(dm2$group1.x, dm2$group1.y, sep=".vs.")

ggplot(dm2, aes(x=group, y=value, fill=group))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))


# =======================================================
#    原发病
# =======================================================

rm(list=ls())

dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)

zy_pcoa = function(current_group = NA, select_region =NA){
  index_g = which(colnames(map) == current_group)
  map = map[which(map$group3!="CKD" & map[,current_group] != "na" & map$group2 %in% select_region),] # 筛选
  dt = dt[,map$sample]
  
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  
  
  
  #----adonis
  ado = adonis(t(dt)~map[,current_group], method='bray')
  r2 = ado$aov.tab$R2[1]
  p = ado$aov.tab$`Pr(>F)`[1]
  
  #----pcoa
  otu.dist = vegdist(as.matrix(t(dt)), method='bray')
  pcoa = cmdscale(otu.dist, k=4, eig=T)
  eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
  point = pcoa$points
  data = merge(point, map, by.x='row.names', by.y='sample')
  xlab = paste("PCoA 1 ( ", eigs[1], "% )", sep="")
  ylab = paste("PCoA 2 ( ", eigs[2], "% )", sep="")
  
  
  #-----正常的pcoa图
  p1 = ggscatter(data, x="V1", y="V2", color = current_group,
                ellipse = T,
                #shape='group2',
                ellipse.level = 0.3,
                xlab=xlab,ylab=ylab,
                title =paste(select_region, current_group, "\nR2=",r2, '\npvalue=', p, collapse = "_and_")
                )+
    scale_color_npg()+
    theme_bw()
  p1
  
  
  #-----只画mean point 和 std
  se  <- function(x) sqrt(var(x)/length(x))
  std <- function(x) sd(x)/sqrt(length(x))
  dat = aggregate(data[,2:5], list(group=data[,current_group]), mean)
  dat2 = aggregate(data[,2:5], list(group=data[,current_group]), std)
  
  dat$Low.1 = dat$V1-dat2$V1 
  dat$Up.1 = dat$V1+dat2$V1 
  dat$Low.2 = dat$V2-dat2$V2 
  dat$Up.2 = dat$V2+dat2$V2 
  
  p2 = ggplot(data=dat,aes(x=V1,y=V2,color=group))+
    geom_hline(yintercept = 0,linetype=2,color='grey')+
    geom_vline(xintercept = 0,linetype=2,color='grey')+
    geom_point(size=10)+
    geom_errorbar(aes(ymin=Low.2,ymax=Up.2),width = 0.002)+
    geom_errorbar(aes(xmin = Low.1, xmax = Up.1),width = 0.005)+
    geom_text(aes(label=group),color='black',size=5)+
    theme_bw()+
    scale_color_npg()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    )+
    labs(title=paste(select_region, current_group, "\nR2=",r2, '\npvalue=', p, collapse = "_and_"))+
    xlab(xlab)+ylab(ylab)
  
  p2
  #list(p1,p2)
  #ggsave("new_pcoa.bj_protopathy.pdf", p2, width=5,height=5)
}

region = "BJ"
region = "SH"
region = c("BJ","SH")

# - 原发病
p = zy_pcoa("protopathy", region)

# - 并发症
p = zy_pcoa("cvd", region)
p = zy_pcoa("Diabetes", region)
p = zy_pcoa("Dialysis period(level)", region)
p = zy_pcoa("high_blood_pressure", region)

ggsave(paste("new_pcoa.", region, "_protopathy.pdf",sep=""), p, width=5,height=5)




#-------------------------------------------------
#                 esrd 根据毒素浓度分型
#-------------------------------------------------


library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsignif)
rm(list=ls())
dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)
#toxin = read.table("../../00.data/scale.toxin.txt",sep="\t", header=T, check.names=F, row.names=1)
toxin = read.table("../../00.data/esrd_scale.toxin.txt",sep="\t", header=T, check.names=F, row.names=1)
map = map[which(map$group3 != "CKD"), ]
dt = dt[,map$sample]


otu.dist = vegdist(as.matrix(t(dt)), method='bray')

pcoa = cmdscale(otu.dist, k=8, eig=T)
eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
point = pcoa$points
data = merge(point, map, by.x='row.names', by.y='sample')

xn=1
yn=2

xlab = paste("PCoA ", xn, " ( ", eigs[xn], "% )", sep="")
ylab = paste("PCoA ", yn, " ( ", eigs[yn], "% )", sep="")
title = "PCoA 	Bray-Curtis"
x = paste("V", xn, sep="")
y = paste("V", yn, sep="")

dm = merge(data,toxin,by.x='Row.names', by.y='row.names', all.x=T)




dm = merge(data,toxin,by.x='Row.names', by.y='row.names', all.x=T)

tox = 'PAG'
xx = summary(dm[,tox])
q3 = xx["3rd Qu."]
q1 = xx["1st Qu."]
q2 = xx["Median"]

dm[,tox] = ifelse(dm[,tox] >  q3, "l4", dm[,tox])
dm[,tox] = ifelse(dm[,tox] >  q2 & dm[,tox] <=  q3, "l3", dm[,tox])
dm[,tox] = ifelse(dm[,tox] >  q1 & dm[,tox] <=  q2, "l2", dm[,tox])
dm[,tox] = ifelse(dm[,tox] <= q1, "l1",dm[,tox])

dm[,tox] = as.vector(dm[,tox])

ggplot(dm, aes(x=.data[[x]], y=.data[[y]], z=.data[[tox]]))+
  geom_point(aes(color=group3, alpha=.data[[tox]]),size=3)+
  scale_color_manual(values=c("#e41a1c", "#4daf4a"))+
  scale_alpha_continuous(range = c(0.1,1), n.breaks = 5)+
  stat_density2d()


ggscatter(dm, x=x, y=y,color = 'PAG',
          ellipse = T,
          ellipse.level = 0.8,
          star.plot = F,
          ellipse.alpha = 0.1,
          title = title,
          mean.point = T,
          alpha="PAG",
          ellipse.border.remove = F,
          xlab=xlab,ylab=ylab)
#scale_color_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))+
#scale_alpha_continuous(range=c(0.1,1))+
theme_bw()






# =======================================================
#   boxplot 单看上海的所有ckd control esrd
# =======================================================

library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggsignif)
library(reshape2)

dt = read.table("../../00.data/new_mag.clust.profile.rename", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)

map = map[which(map$group2 == "SH"), ]
dt = dt[,map$new_sample]


otu.dist = vegdist(as.matrix(t(dt)), method='bray')

library(reshape2)
bb = as.matrix(otu.dist)
dl = melt(bb)
dl = dl[duplicated(t(apply(dl, 1, sort))),]

map = map[,c(5,2)]
dm1 = merge(dl, map, by.x='Var1', by.y='new_sample')
dm2 = merge(dm1, map, by.x='Var2', by.y='new_sample')
dm2 = dm2[which(dm2$Var2 != dm2$Var1),]
dm2$group = paste(dm2$group1.x, dm2$group1.y, sep=".vs.")

ggplot(dm2, aes(x=group, y=value, fill=group))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))



# =======================================================
#   pcoa 单看上海的所有ckd control esrd  
# =======================================================

dt = read.table("../..//00.data/new_mag.clust.profile.rename", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)


library(vegan)
library(ggplot2)
library(ggpubr)

map = map[which(map$group2 == "SH"),]

otu.dist = vegdist(as.matrix(t(dt)), method='bray')

pcoa = cmdscale(otu.dist, k=4, eig=T)
eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
point = pcoa$points
data = merge(point, map, by.x='row.names', by.y='new_sample')

xlab = paste("PCoA 1 ( ", eigs[1], "% )", sep="")
ylab = paste("PCoA 2 ( ", eigs[2], "% )", sep="")

p = ggscatter(data, x="V2", y="V3", color = 'group3',
              ellipse = T,shape='group2',
              ellipse.level = 0.8,
              star.plot = T,
              xlab=xlab,ylab=ylab)+
  scale_color_brewer(palette='Set1')+
  theme_bw()
p



# =======================================================
#   pcoa新方式 单看上海的所有ckd control esrd  
# =======================================================

library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)

rm(list=ls())

dt = read.table("../../00.data/new_mag.clust.profile.rename", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)
map = map[which(map$group2 == "SH"),]

zy_pcoa = function(current_group = NA, select_region =NA){
  index_g = which(colnames(map) == current_group)
  dt = dt[,map$new_sample]
  
  #----adonis
  ado = adonis(t(dt)~map[,current_group], method='bray')
  r2 = ado$aov.tab$R2[1]
  p = ado$aov.tab$`Pr(>F)`[1]
  
  #----pcoa
  otu.dist = as.matrix(vegdist(as.matrix(t(dt)), method='bray'))
  pcoa = cmdscale(otu.dist, k=4, eig=T)
  eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
  point = pcoa$points
  data = merge(point, map, by.x='row.names', by.y='new_sample')
  xlab = paste("PCoA 1 ( ", eigs[1], "% )", sep="")
  ylab = paste("PCoA 2 ( ", eigs[2], "% )", sep="")
  
  
  #--------barplot distance
  dl = melt(otu.dist)
  dl = dl[duplicated(t(apply(dl, 1,sort))),]
  temp_map = map[,c('new_sample', current_group)]
  dm1 = merge(dl, temp_map, by.x='Var1',by.y='new_sample')
  dm2 = merge(dm1, temp_map, by.x='Var2', by.y='new_sample')
  dm2 = dm2[which(dm2$Var2 != dm2$Var1),]
  dm2$group = paste(dm2$group4.x, dm2$group4.y, sep=".vs.")
  dm2$group = gsub("^CKD-5nD.vs.CKD$", "CKD.vs.CKD-5nD", dm2$group)
  comp = combn(unique(dm2$group),2,list)
  dm2$group = factor(dm2$group, level=c("Control.vs.Control","CKD.vs.Control","CKD-5nD.vs.Control","Control.vs.ESRD","CKD.vs.CKD","CKD.vs.CKD-5nD","CKD.vs.ESRD","CKD-5nD.vs.CKD-5nD","CKD-5nD.vs.ESRD","ESRD.vs.ESRD"))
  
  
  p0 = ggplot(data=dm2, aes(x=group, y=value,))+
    geom_boxplot(fill='#66ccff')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  
  #-----正常的pcoa图
  p1 = ggscatter(data, x="V1", y="V2", color = current_group,
                 ellipse = T,
                 #shape='group2',
                 ellipse.level = 0.3,
                 xlab=xlab,ylab=ylab,
                 title =paste(select_region, current_group, "\nR2=",r2, '\npvalue=', p, collapse = "_and_")
  )+
    scale_color_npg()+
    theme_bw()
  p1
  
  
  #-----只画mean point 和 std
  se  <- function(x) sqrt(var(x)/length(x))
  std <- function(x) sd(x)/sqrt(length(x))
  dat = aggregate(data[,2:5], list(group=data[,current_group]), mean)
  dat2 = aggregate(data[,2:5], list(group=data[,current_group]), std)
  
  dat$Low.1 = dat$V1-dat2$V1 
  dat$Up.1 = dat$V1+dat2$V1 
  dat$Low.2 = dat$V2-dat2$V2 
  dat$Up.2 = dat$V2+dat2$V2 
  
  p2 = ggplot(data=dat,aes(x=V1,y=V2,color=group))+
    geom_hline(yintercept = 0,linetype=2,color='grey')+
    geom_vline(xintercept = 0,linetype=2,color='grey')+
    geom_point(size=10)+
    geom_errorbar(aes(ymin=Low.2,ymax=Up.2),width = 0.002)+
    geom_errorbar(aes(xmin = Low.1, xmax = Up.1),width = 0.005)+
    geom_text(aes(label=group),color='black',size=5)+
    theme_bw()+
    scale_color_npg()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    )+
    labs(title=paste(select_region, current_group, "\nR2=",r2, '\npvalue=', p, collapse = "_and_"))+
    xlab(xlab)+ylab(ylab)
  
  #-----dbrda
  dbrda = capscale(t(dt)~map[,current_group], distance='bray')
  point = scores(dbrda)$sites
  data = merge(point, map, by.x='row.names', by.y='new_sample')
  #xlab = paste("PCoA 1 ( ", eigs[1], "% )", sep="")
  #ylab = paste("PCoA 2 ( ", eigs[2], "% )", sep="")
  
  p3 = ggscatter(data, x="CAP1", y="CAP2", color = current_group,
                 ellipse = T,
                 star.plot = T,
                 #shape='group2',
                 ellipse.level = 0.3,
                 #xlab=xlab,ylab=ylab,
                 title =paste(select_region, current_group, "\nR2=",r2, '\npvalue=', p, collapse = "_and_")
  )+
    scale_color_npg()+
    theme_bw()
  
  p0
  #list(p1,p2)
  #ggsave("new_pcoa.bj_protopathy.pdf", p2, width=5,height=5)
}


########上海没有cvd ，所以不要用了

p = zy_pcoa('group4', 'SH')

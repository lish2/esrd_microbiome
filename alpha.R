library(vegan)
library(ggplot2)

library(patchwork)
library(ggsignif)

dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names=F, row.names=1)
map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)

map = map[which(map$group3 != "CKD"), ]
dt = dt[,map$sample]

div_zy = function(index=NA){
  
  sha = diversity(t(dt), index=index)
  compar = list(c("BJ_ESRD", "BJ_HC"), 
                c("SH_ESRD", "SH_HC"),
                c("BJ_ESRD", "SH_ESRD"),
                c("BJ_HC", "SH_HC"),
                c("SH_CKD", "SH_HC")
  )
  
  data = merge(sha, map, by.x='row.names', by.y='sample')
  p = ggplot(data=data, aes(x=group3,y=x,fill=group3))+
    geom_boxplot()+
    theme_bw() +
    scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a", "#ff7f00"))+
    geom_signif(comparisons = compar, step_increase = 0.1,test=t.test,map_signif_level=T,)+
    ggtitle(index)
  p
}

p1 = div_zy('shannon')
p2 = div_zy('simpson')
p3 = div_zy('invsimpson')
p = p1+p2+p3
p

ggsave("alpha.diversity.g5.pdf", width=14, height=5)



#=================

div_zy = function(xx=NA, title=NA){
  
  dd = colSums(dt > xx)
  

  compar = list(c("BJ_ESRD", "BJ_HC"), 
                c("SH_ESRD", "SH_HC"),
                c("BJ_ESRD", "SH_ESRD"),
                c("BJ_HC", "SH_HC")
  )
  
  data = merge(dd, map, by.x='row.names', by.y='sample')
  p = ggplot(data=data, aes(x=group1,y=x,fill=group1))+
    geom_boxplot()+
    theme_bw() +
    scale_fill_manual(values=c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a"))+
    geom_signif(comparisons = compar, 
                step_increase = 0.1,
                test=t.test,
                #map_signif_level=T,
    )+
  ggtitle(title)
  p
}

div_zy(xx=20)


p1 = div_zy(20, "1e-6")
p2 = div_zy(200, "1e-5")
p3 = div_zy(2000, "1e-4")
p4 = div_zy(20000, "1e-3")


p  = p1 + p2 + p3 + p4
p


ggsave("xxx.pdf", p, width=10, height=8)



compar = list(c("HC","CGN"),c("HC","DKD"),c("HKD","HC"),c("HC","other"))
map = map[which(map$group3 != "CKD" & map$protopathy != "na"), ]






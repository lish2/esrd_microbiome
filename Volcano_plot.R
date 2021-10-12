#===================================
#--------菌株差异
dt= read.table("./kk.txt",sep="\t", header=T, check.names=F)
library(ggplot2)
dt1 = dt
dt1$lfdr = -log10(dt1$lfdr)

dt1$color = ifelse(dt1$fold_change<1.2 & dt1$fold_change > -1.2, "ns", dt1$color)
dt1$fold_change = ifelse(dt1$enriched=="ESRD-enriched", 0-log2(dt1$fold_change), log2(dt1$fold_change))

ggplot(dt1, aes(x=fold_change,
               y=lfdr,
               colour=color,
               alpha=color,
               shape=color))+
  geom_point()+
  theme_bw()+
  geom_vline(xintercept=log(1.2), linetype=5)+
  geom_vline(xintercept=-log(1.2), linetype=5)+
  geom_hline(yintercept = -log10(0.2), linetype=5)+
  geom_hline(yintercept = -log10(0.05), linetype=5)+
  facet_grid(region~., scales='free')+
  scale_color_manual(values=c("#d53e4f", "#3288bd", "gray"))+
  scale_alpha_manual(values=c(1,1,0.5))+
  scale_shape_manual(values=c(16,16,1))+
  xlab("log2(fold_change)")+
  ylab("-log10(lfdr)")

#----------------------
dt= read.table("./module.txt",sep="\t", header=T, check.names=F)
library(ggplot2)
dt1 = dt
dt1$lfdr = -log10(dt1$lfdr)

dt1$color = ifelse(dt1$fold<1.2 & dt1$fold > -1.2, "ns", dt1$color)
dt1$fold = ifelse(dt1$enrich=="ESRD", 0-log2(dt1$fold), log2(dt1$fold))

ggplot(dt1, aes(x=fold,
                y=lfdr,
                colour=color,
                alpha=color,
                shape=color))+
  geom_point()+
  theme_bw()+
  geom_vline(xintercept=log2(1.2), linetype=5)+
  geom_vline(xintercept=-log2(1.2), linetype=5)+
  geom_hline(yintercept = -log10(0.2), linetype=5)+
  geom_hline(yintercept = -log10(0.05), linetype=5)+
  facet_grid(region~., scales='free')+
  scale_color_manual(values=c("#d53e4f", "#3288bd", "gray"))+
  scale_alpha_manual(values=c(1,1,0.5))+
  scale_shape_manual(values=c(16,16,1))+
  xlab("log2(fold_change)")+
  ylab("-log10(lfdr)")


#----------------------ko
dt= read.table("./ko.txt",sep="\t", header=T, check.names=F)
library(ggplot2)
dt1 = dt
dt1$lfdr = -log10(dt1$lfdr)

dt1$color = ifelse(dt1$fold<1.2 & dt1$fold > -1.2, "ns", dt1$color)
dt1$fold = ifelse(dt1$enrich=="ESRD", 0-log2(dt1$fold), log2(dt1$fold))

ggplot(dt1, aes(x=fold,
                y=lfdr,
                colour=color,
                alpha=color,
                shape=color))+
  geom_point()+
  theme_bw()+
  geom_vline(xintercept=log2(1.2), linetype=5)+
  geom_vline(xintercept=-log2(1.2), linetype=5)+
  geom_hline(yintercept = -log10(0.2), linetype=5)+
  geom_hline(yintercept = -log10(0.05), linetype=5)+
  facet_grid(region~., scales='free')+
  scale_color_manual(values=c("#d53e4f", "#3288bd", "gray"))+
  scale_alpha_manual(values=c(1,1,0.5))+
  scale_shape_manual(values=c(16,16,1))+
  xlab("log2(fold_change)")+
  ylab("-log10(lfdr)")





#---------------------------------------------------
#       单看上海的菌株，在四个分组中，是否有差异
#---------------------------------------------------
rm(list=ls())
dt= read.table("./x.txt",sep="\t", header=T, check.names=F)
library(ggplot2)
library(patchwork)
library(grid)

dt1 = dt

# color_dict = c("#3288bd","green", "#c51b7d", "#d53e4f", 'gray', "blue", "#e08214") # 当前差异菌也涂色
color_dict = c("gray","gray", "gray", "gray", 'gray', "#3288bd", "#d53e4f")
names(color_dict) = c("Control","CKD","CKD-5nD","ESRD", 'ns',"HC-enriched", "ESRD-enriched")

position_dict = structure( c("CKD-5nD","CKD","CKD-5nD","ESRD","ESRD","ESRD"),
                           names=c("CKD .vs CKD-5nD","Control .vs CKD","Control .vs CKD-5nD","Control .vs ESRD","ESRD .vs CKD","ESRD .vs CKD-5nD"))

y_label = structure(c("CKD","Control","Control","Control","CKD","CKD-5nD"),
          names= c("CKD .vs CKD-5nD","Control .vs CKD","Control .vs CKD-5nD","Control .vs ESRD","ESRD .vs CKD","ESRD .vs CKD-5nD"))

alpha_dict = c(0.2, 0.2, 0.2, 0.2, 1, 1, 0.2)
names(alpha_dict) = c("Control","CKD","CKD-5nD","ESRD", "HC-enriched", "ESRD-enriched", 'ns')
shape_dict = c(1,1,1,1,16,16,1)
names(shape_dict) = c("Control","CKD","CKD-5nD","ESRD", "HC-enriched", "ESRD-enriched", 'ns')



grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=2)))
defint_region <- function(row,col){
  viewport(layout.pos.row=row, layout.pos.col=col)
}

x=2 # 控制列
y=3 #控制行
z = 1

for (i in seq(2,17,3)){
  
  temp_dt = dt1[,c(1,i,i+1,i+2,20)]
  temp_dt = na.omit(temp_dt)
  title = colnames(temp_dt)[3]
  title = gsub("lfdr \\(|\\)","",title)
  colnames(temp_dt) = c("name",'fold','lfdr','enriched','old_enriched')

  temp_dt$fold_c = ifelse(temp_dt$enriched==position_dict[title], 0-log2(temp_dt$fold), log2(temp_dt$fold))
  
  # 这个用来计算一致性
  #temp_dt$enriced_count = ifelse(abs(temp_dt$fold) > 1.2 & temp_dt$lfdr < 0.2,  temp_dt$enriched, "ns")  # 这个是看差异倍数和显著性
  temp_dt$enriced_count = ifelse(abs(temp_dt$fold) > 1.2,  temp_dt$enriched, "ns")   # 这个是只看差异倍数是否一致 
  
  l = paste("^",position_dict[title],"$",sep="")
  r = paste("^",y_label[title],"$",sep="")
  
  temp_dt$enriced_count = gsub(l,"ESRD-enriched",temp_dt$enriced_count)
  temp_dt$enriced_count = gsub(r,"HC-enriched",temp_dt$enriced_count)
  

  mm = temp_dt[which(temp_dt$old_enriched != "na"),]
  same = sum(mm$old_enriched == mm$enriced_count)
  percent = round(same/348*100,digits=2)
  title = paste(title,"\nconsistency: ", percent,"%",sep="")
  # 这个enriched用来画图
  temp_dt$enriched = ifelse(temp_dt$old_enriched != "na", temp_dt$old_enriched, ifelse(abs(temp_dt$fold) > 1.2 & temp_dt$lfdr < 0.2,  temp_dt$enriched, "ns"))
  
  temp_dt$lfdr = -log10(temp_dt$lfdr)
  temp_dt = na.omit(temp_dt)
  
  p = ggplot(temp_dt, aes(x=fold_c,
                  y=lfdr,
                  colour=enriched,
                  alpha=enriched,
                  shape=enriched
                  ))+
    geom_point()+
    theme_bw()+
    ggtitle(label=title)+
    geom_vline(xintercept=log2(1.2), linetype=5)+
    geom_vline(xintercept=-log2(1.2), linetype=5)+
    geom_hline(yintercept = -log10(0.2), linetype=5)+
    geom_hline(yintercept = -log10(0.05), linetype=5)+
    scale_alpha_manual(values=alpha_dict)+
    scale_shape_manual(values=shape_dict)+
    scale_color_manual(values=color_dict)+
    xlab("log2(fold_change)")+
    ylab("-log10(lfdr)")
  print(p, vp=defint_region(y%%3+1, x%%2+1))
  x=x+1
  if(z%%2 == 0){y=y+1}
  z=z+1
  
}



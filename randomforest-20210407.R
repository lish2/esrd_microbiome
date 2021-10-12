rm(list=ls())

library(randomForest)
library('pROC')

#---------------
ntree_plot <- function(out=NA,p=NA){
  pdf(out)
  plot(p)
  dev.off()
}
importance_plot <- function(out=NA,p=NA){
  pdf(out)
  varImpPlot(p)
  dev.off()
}
roc_plot <- function(out=NA,p=NA){
  pdf(out)
  plot.roc(p,add=F, 
           reuse.auc=TRUE,col="blue", partial.auc=c(1, 0.8), print.auc = F,
           print.auc.cex=2, print.auc.col='Black',
           max.auc.polygeon=T,
           axes=TRUE,
           grid = c(0.2,0.2),grid.col = 'grey')
  ciauc = round(ci.auc(p), digits = 3)
  textauc = paste("AUC: ", ciauc[2], " (", ciauc[1], " ~ ", ciauc[3],")", sep="" )
  text(x=0.4 ,y=0.4, labels=textauc, cex=2)
  dev.off()
}

ntree = 2000
dt = read.table("../../00.data/new_mag.clust.profile.rename",sep="\t", check.names=F, row.names = 1,header=T)
sample_map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)
sample_map = sample_map[which(sample_map$group2 != "BJ"), c("new_sample", "group3")]
sample_map$rf = ifelse(sample_map$group3 == "Control","Control","kd")
sample_map$rf = as.factor(sample_map$rf)

dtt = t(dt)
dtt = dtt[sample_map$new_sample,]


# ckd -> test;    hc/esrd -> model
group_train = sample_map[which(sample_map$group3 != "CKD"),]
dtt_train = dtt[group_train$new_sample,]

group_test = sample_map[which(sample_map$group3 != "ESRD"),]
dtt_test = dtt[group_test$new_sample,]

set.seed(123)
fit = randomForest(dtt_train, group_train$rf, ntree=ntree,proximity=TRUE)
pred = as.data.frame(predict(fit, dtt_test, type='prob'))
roc_p = roc(as.character(group_test$rf), pred[,2])
roc_plot("roc.only_SH_esrd2ckd.pdf",roc_p)
#roc_plot("roc.only_SH_ckd2esrd.pdf",roc_p)










#----------------------------------------------------------
#               肾病分期，只看上海的
#----------------------------------------------------------
rm(list=ls())

library(randomForest)
library('pROC')

#---------------
ntree_plot <- function(out=NA,p=NA){
  pdf(out)
  plot(p)
  dev.off()
}
importance_plot <- function(out=NA,p=NA){
  pdf(out)
  varImpPlot(p)
  dev.off()
}
roc_plot <- function(out=NA,p=NA, color='blue', title=NA,y=0.4, x=x, add=FALSE){
  #pdf(out)
  plot.roc(p,add=add, 
           reuse.auc=TRUE,col=color, partial.auc=c(1, 0.8), print.auc = F,
           print.auc.cex=2, print.auc.col='Black',
           max.auc.polygeon=T,
           axes=TRUE, #legacy.axes=T,
           grid = c(0.2,0.2),grid.col = 'grey')
  ciauc = round(ci.auc(p), digits = 3)
  textauc = paste(title,"AUC: ", ciauc[2], " (", ciauc[1], " ~ ", ciauc[3],")", sep="" )
  text(x=x ,y=y, labels=textauc, cex=1, col=color)
  #dev.off()
}

ntree = 2000
dt = read.table("../../00.data/new_mag.clust.profile.rename",sep="\t", check.names=F, row.names = 1,header=T)
sample_map = read.table("../../00.data/mapping.file", sep="\t", header=T, check.names=F)
sample_map = sample_map[which(sample_map$group2 != "BJ"), ]

sample_map$rf = ifelse(sample_map$group3 == "Control","Control","kd")
sample_map$rf = as.factor(sample_map$rf)

dtt = t(dt)
dtt = dtt[sample_map$new_sample,]
group_control = sample_map[which(sample_map$group4 == "Control"),]
dtt_control = dtt[group_control$new_sample,]


plot_zy = function(a="CKD", b='ESRD', color='blue',y=0.4,x=0.4, add=add){
  
  # 分割健康人数据集 0.7 -> train
  set.seed(123)
  control = sample(seq(1, nrow(group_control)), nrow(group_control)*0.7)
  group_train_con = group_control[control, ]
  group_test_con = group_control[-control, ]
  dtt_train_con = dtt_control[group_train_con$new_sample, ]
  dtt_test_con = dtt_control[group_test_con$new_sample, ]
  
  # 删选病人的
  group_train = sample_map[which(sample_map$group4 == a),]
  group_test = sample_map[which(sample_map$group4 == b),]
  
  group_train = rbind(group_train, group_train_con)
  group_test = rbind(group_test, group_test_con)
  dtt_train = dtt[group_train$new_sample, ]
  dtt_test = dtt[group_test$new_sample, ]

  
  fit = randomForest(dtt_train, group_train$rf, ntree=ntree,proximity=TRUE)
  pred = as.data.frame(predict(fit, dtt_test, type='prob'))
  roc_p = roc(as.character(group_test$rf), pred[,2])
  roc_plot(paste("roc.only_SH_",a,"2",b,".pdf",sep=""),roc_p, color=color, title=paste(a,"->",b,spe=""),y=y, x=x, add=add)
  
}




plot_zy("CKD-5nD","ESRD",color='red',y=0.6,add=FALSE)
plot_zy("ESRD","CKD-5nD", color='orange',y=0.5, add=T)

plot_zy("CKD","ESRD", color='blue', y=0.4, add=T)
plot_zy("ESRD","CKD", color='green', y=0.3, add=T)

plot_zy("CKD-5nD","CKD", color='#a65628', y=0.2, add=T)
plot_zy("CKD","CKD-5nD", color='#984ea3', y=0.1, add=T)

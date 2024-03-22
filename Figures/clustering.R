rm(list = ls())
#args <- commandArgs(trailingOnly = TRUE)
#strain = args[1]
#print(strain)
library(ape)
library(cluster)


#strain="/finaloutput/Malawi_bov_rem/"
strain="scripts/finaloutput/Malawi_bov_rem_rerun/"
diverg=read.csv(paste(strain,"diverg_id.csv",sep=''))
print(diverg)
avector <- as.character(diverg[, "X0"])
avector_f <- sapply(avector, function(x) substr(x, 1, 1))
print(names(data)[2])
data=read.csv(paste(strain,"kldiverg.csv",sep=''))
theta=read.csv(paste(strain,"Prior.csv",sep=''))
names(data)<-avector
matrixdata = as.matrix(data)
print(matrixdata)
matrixdata=as.dist(matrixdata)
fit <- hclust(matrixdata, method="ward.D2")
png(paste(strain,"cluster_grouping.png",sep=""), width = 3000, height = 2800, res = 300)
#plot(fit,cex=2,cex.lab=2,cex.main=2,cex.axis=2,lwd=2)
plot(as.phylo(fit), type = "fan",cex=1.5,cex.lab=2,cex.main=2,cex.axis=2)
dev.off()
final_prop<- matrix(0, nrow = nrow(theta), ncol = 2)
final_strain<- matrix(0, nrow = nrow(theta), ncol = 2)
final_format_strain<- matrix(0, nrow = nrow(theta), ncol = 1)
for (i in 1:nrow(theta)){
for(j in 1:length(fit$height)){
cut <-cutree(fit,h=fit$height[j])
cluster_prop<-numeric(max(cut))
strain_p<-numeric(max(cut))
for(k in 1:max(cut)){
cutvector<-as.matrix(cut)
cluster_indices<-which(cutvector==k)
cluster_prop[k]<-sum(theta[i,cluster_indices])
uni<-unique(avector_f[cluster_indices])
strain_p[k]<-paste0(uni, collapse = "/")
}
non_zero_prop<-which(cluster_prop>0)
if(length(non_zero_prop)<=2){
ind<-1
for(k in non_zero_prop)
{
final_prop[i,ind]<-cluster_prop[k]
final_strain[i,ind]<-paste("L",strain_p[k],sep="")
if(ind==2)
{
  if(final_prop[i,ind-1]<final_prop[i,ind])
  {
    temp<-final_strain[i,ind-1]
    final_strain[i,ind-1]<-final_strain[i,ind]
    final_strain[i,ind]<-temp
  }
  if(final_strain[i,ind]!="0")
    final_format_strain[i,1]<-paste(final_strain[i,ind-1],final_strain[i,ind],sep="/")
  else
  {
    final_format_strain[i,1]<-final_strain[i,ind]
  }
    
}
else
{
  final_format_strain[i,1]<-final_strain[i,ind] 
}
ind<-ind+1
}
break;
}
}
}
write.csv(final_prop,file=paste(strain,"final_proportion.csv",sep=""), row.names = FALSE)
write.csv(final_format_strain,file=paste(strain,"final_lineage.csv",sep=""), row.names = FALSE)
library(ggplot2)
#library(pheatmap)
library(viridisLite)

library(ComplexHeatmap)

theta=read.csv(paste(strain,"Prior.csv",sep=''))
row_index <- c(1644,1850,1881,883,885,906)
selected_rows <- theta[row_index,] 
highlighted_column<-41
cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(10),
          colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(50))


#,labels_col=avector,
png(paste(strain,"heatmap.png",sep=""), width = 1200, height = 500, res = 300)

ht1<-ComplexHeatmap::pheatmap(  selected_rows,color=cols,fontsize = 2, cluster_rows = FALSE,fontsize_row = 4,  
           fontsize_col = 4,labels_col=paste0("T", 1:42),labels_row=c("ERR181980","ERR212052","ERR212091","ERR161044","ERR161048","ERR161071"),border_color = "black" ,cluster_cols = FALSE,legend = TRUE,fontsize_number = 4,
           angle_col = c("90"),  row_names_side = c("left"),heatmap_legend_param = list(
             legend_direction = "horizontal", 
             legend_width = unit(5, "cm"),title = "Strain proportion",fontsize=1
           ))
draw(ht1, heatmap_legend_side = "bottom")
dev.off()


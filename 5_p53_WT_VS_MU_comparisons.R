
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

wt_table = readRDS("wt_table.rds")
mu_table = readRDS("mu_table.rds")

library(VennDiagram)

v <- venn.diagram(list(P53_WT=rownames(wt_table)[which(wt_table$mean.diff<(-.25) & wt_table$diffmeth.p.adj.fdr<0.05)], 
                       P53_MU=rownames(mu_table)[which(mu_table$mean.diff<(-.25) & mu_table$diffmeth.p.adj.fdr<0.05)]),
                  fill = c("orange", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 1, cex=1.5,
                  filename=NULL)

pdf("venn_25diff.pdf")
grid.draw(v)
dev.off()

v <-  venn.diagram(list(P53_WT=rownames(wt_table)[which(wt_table$mean.diff<(-.2) & wt_table$diffmeth.p.adj.fdr<0.05)], 
                       P53_MU=rownames(mu_table)[which(mu_table$mean.diff<(-.2) & mu_table$diffmeth.p.adj.fdr<0.05)]),
                  fill = c("orange", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 1, cex=1.5,
                  filename=NULL)

pdf("venn_20diff.pdf")
grid.draw(v)
dev.off()

pdf("meth_changes_Scatter_fdr5.pdf")
ix = (wt_table$mean.diff<(-.2) & wt_table$diffmeth.p.adj.fdr<0.05) | (mu_table$mean.diff<(-.2) & mu_table$diffmeth.p.adj.fdr<0.05)
plot(wt_table$mean.diff[ix],mu_table$mean.diff[ix],xlim=c(-.5,.2),ylim=c(-.5,.2),col=alpha("grey",.6),pch=20)
abline(lm(mu_table$mean.diff[ix]~wt_table$mean.diff[ix]),col="red")
abline(0,1,col="black",lty=2)
dev.off()

pdf("meth_changes_diff_of_diff_20diff_fdr5.pdf")
ix = (wt_table$mean.diff<(-.2) & wt_table$diffmeth.p.adj.fdr<0.05) | (mu_table$mean.diff<(-.2) & mu_table$diffmeth.p.adj.fdr<0.05)
plot(density(wt_table$mean.diff[ix]-mu_table$mean.diff[ix]),lwd=2,main="(P53 WT: Treated - Untreated ) - (P53 MU: Treated - Untreated )")
legend("topright",legend=c("TP53 WT"), border=T, bty="n" )
legend("topleft",legend=c("TP53 MU"), border=T, bty="n" )
abline(v=0,lty=2)
dev.off()

pdf("meth_changes_Scatter.pdf")
plot(wt_table$mean.diff,mu_table$mean.diff,xlim=c(-.5,.2),ylim=c(-.5,.2),col=alpha("grey",.6),pch=20,
     xlab="P53 WT: Treated - Untreated",ylab="P53 MU: Treated - Untreated")
abline(lm(mu_table$mean.diff~wt_table$mean.diff),col="red")
abline(0,1,col="black",lty=2)
dev.off()

########################################
track=as.character(rnb.set.filtered@pheno$TP53)
track[track=="WT"]=1
track[track=="MU"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])

track=as.character(rnb.set.filtered@pheno$Treatment)
track[track=="Untreated"]=1
track[track=="Treated"]=2
track=as.numeric(track)
colores=c("black","red")
clab2=as.character(colores[track])

ix = (wt_table$mean.diff<(-.25) & wt_table$diffmeth.p.adj.fdr<0.05) | (mu_table$mean.diff<(-.25) & mu_table$diffmeth.p.adj.fdr<0.05)
meth.filtered.sig = meth.filtered[ix,]
meth.filtered.sig = meth.filtered.sig[complete.cases(meth.filtered.sig),]
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))

###########
meth.norm = meth.filtered[ix, ]
meth.norm.centered = meth.norm
for(ix in 1:dim(meth.norm)[1]){ 
           meth.norm.centered[ix,1] = meth.norm[ix,1]-mean(meth.norm[ix,1:2])
           meth.norm.centered[ix,2] = meth.norm[ix,2]-mean(meth.norm[ix,1:2])
           meth.norm.centered[ix,3] = meth.norm[ix,3]-mean(meth.norm[ix,3:4])
           meth.norm.centered[ix,4] = meth.norm[ix,4]-mean(meth.norm[ix,3:4])
           meth.norm.centered[ix,5] = meth.norm[ix,5]-mean(meth.norm[ix,5:6])
           meth.norm.centered[ix,6] = meth.norm[ix,6]-mean(meth.norm[ix,5:6])
           meth.norm.centered[ix,7] = meth.norm[ix,7]-mean(meth.norm[ix,7:8])
           meth.norm.centered[ix,8] = meth.norm[ix,8]-mean(meth.norm[ix,7:8])
           meth.norm.centered[ix,9] = meth.norm[ix,9]-mean(meth.norm[ix,9:10])
           meth.norm.centered[ix,10] = meth.norm[ix,10]-mean(meth.norm[ix,9:10])
           meth.norm.centered[ix,11] = meth.norm[ix,11]-mean(meth.norm[ix,11:12])
           meth.norm.centered[ix,12] = meth.norm[ix,12]-mean(meth.norm[ix,11:12])
           meth.norm.centered[ix,13] = meth.norm[ix,13]-mean(meth.norm[ix,13:14])
           meth.norm.centered[ix,14] = meth.norm[ix,14]-mean(meth.norm[ix,13:14])
           meth.norm.centered[ix,15] = meth.norm[ix,15]-mean(meth.norm[ix,15:16])
           meth.norm.centered[ix,16] = meth.norm[ix,16]-mean(meth.norm[ix,15:16])
           meth.norm.centered[ix,17] = meth.norm[ix,17]-mean(meth.norm[ix,17:18])
           meth.norm.centered[ix,18] = meth.norm[ix,18]-mean(meth.norm[ix,17:18])
           meth.norm.centered[ix,19] = meth.norm[ix,19]-mean(meth.norm[ix,19:20])
           meth.norm.centered[ix,20] = meth.norm[ix,20]-mean(meth.norm[ix,19:20])
           meth.norm.centered[ix,21] = meth.norm[ix,21]-mean(meth.norm[ix,21:22])
           meth.norm.centered[ix,22] = meth.norm[ix,22]-mean(meth.norm[ix,21:22])
           meth.norm.centered[ix,23] = meth.norm[ix,23]-mean(meth.norm[ix,23:24])
           meth.norm.centered[ix,24] = meth.norm[ix,24]-mean(meth.norm[ix,23:24])
           meth.norm.centered[ix,25] = meth.norm[ix,25]-mean(meth.norm[ix,25:26])
           meth.norm.centered[ix,26] = meth.norm[ix,26]-mean(meth.norm[ix,25:26])
           meth.norm.centered[ix,27] = meth.norm[ix,27]-mean(meth.norm[ix,27:28])
           meth.norm.centered[ix,28] = meth.norm[ix,28]-mean(meth.norm[ix,27:28])
           meth.norm.centered[ix,29] = meth.norm[ix,29]-mean(meth.norm[ix,29:30])
           meth.norm.centered[ix,30] = meth.norm[ix,30]-mean(meth.norm[ix,29:30])
           meth.norm.centered[ix,31] = meth.norm[ix,31]-mean(meth.norm[ix,31:32])
           meth.norm.centered[ix,32] = meth.norm[ix,32]-mean(meth.norm[ix,31:32])
           meth.norm.centered[ix,33] = meth.norm[ix,33]-mean(meth.norm[ix,33:34])
           meth.norm.centered[ix,34] = meth.norm[ix,34]-mean(meth.norm[ix,33:34])
           meth.norm.centered[ix,35] = meth.norm[ix,35]-mean(meth.norm[ix,35:36])
           meth.norm.centered[ix,36] = meth.norm[ix,36]-mean(meth.norm[ix,35:36])
           meth.norm.centered[ix,37] = meth.norm[ix,37]-mean(meth.norm[ix,37:38])
           meth.norm.centered[ix,38] = meth.norm[ix,38]-mean(meth.norm[ix,37:38])
           meth.norm.centered[ix,39] = meth.norm[ix,39]-mean(meth.norm[ix,39:40])
           meth.norm.centered[ix,40] = meth.norm[ix,40]-mean(meth.norm[ix,39:40])
           meth.norm.centered[ix,41] = meth.norm[ix,41]-mean(meth.norm[ix,41:42])
           meth.norm.centered[ix,42] = meth.norm[ix,42]-mean(meth.norm[ix,41:42])
           meth.norm.centered[ix,43] = meth.norm[ix,43]-mean(meth.norm[ix,43:44])
           meth.norm.centered[ix,44] = meth.norm[ix,44]-mean(meth.norm[ix,43:44])
           meth.norm.centered[ix,45] = meth.norm[ix,45]-mean(meth.norm[ix,45:46])
           meth.norm.centered[ix,46] = meth.norm[ix,46]-mean(meth.norm[ix,45:46])
           meth.norm.centered[ix,47] = meth.norm[ix,47]-mean(meth.norm[ix,47:48])
           meth.norm.centered[ix,48] = meth.norm[ix,48]-mean(meth.norm[ix,47:48])
           meth.norm.centered[ix,49] = meth.norm[ix,49]-mean(meth.norm[ix,49:50])
           meth.norm.centered[ix,50] = meth.norm[ix,50]-mean(meth.norm[ix,49:50])
           meth.norm.centered[ix,51] = meth.norm[ix,51]-mean(meth.norm[ix,51:52])
           meth.norm.centered[ix,52] = meth.norm[ix,52]-mean(meth.norm[ix,51:52])
           meth.norm.centered[ix,53] = meth.norm[ix,53]-mean(meth.norm[ix,53:54])
           meth.norm.centered[ix,54] = meth.norm[ix,54]-mean(meth.norm[ix,53:54])
           meth.norm.centered[ix,55] = meth.norm[ix,55]-mean(meth.norm[ix,55:56])
           meth.norm.centered[ix,56] = meth.norm[ix,56]-mean(meth.norm[ix,55:56])
           meth.norm.centered[ix,57] = meth.norm[ix,57]-mean(meth.norm[ix,57:58])
           meth.norm.centered[ix,58] = meth.norm[ix,58]-mean(meth.norm[ix,57:58])
           meth.norm.centered[ix,59] = meth.norm[ix,59]-mean(meth.norm[ix,59:60])
           meth.norm.centered[ix,60] = meth.norm[ix,60]-mean(meth.norm[ix,59:60])
           meth.norm.centered[ix,61] = meth.norm[ix,61]-mean(meth.norm[ix,61:62])
           meth.norm.centered[ix,62] = meth.norm[ix,62]-mean(meth.norm[ix,61:62])
           meth.norm.centered[ix,63] = meth.norm[ix,63]-mean(meth.norm[ix,63:64])
           meth.norm.centered[ix,64] = meth.norm[ix,64]-mean(meth.norm[ix,63:64])
           meth.norm.centered[ix,65] = meth.norm[ix,65]-mean(meth.norm[ix,65:66])
           meth.norm.centered[ix,66] = meth.norm[ix,66]-mean(meth.norm[ix,65:66])
           meth.norm.centered[ix,67] = meth.norm[ix,67]-mean(meth.norm[ix,67:68])
           meth.norm.centered[ix,68] = meth.norm[ix,68]-mean(meth.norm[ix,67:68])
           meth.norm.centered[ix,69] = meth.norm[ix,69]-mean(meth.norm[ix,69:70])
           meth.norm.centered[ix,70] = meth.norm[ix,70]-mean(meth.norm[ix,69:70])
           meth.norm.centered[ix,71] = meth.norm[ix,71]-mean(meth.norm[ix,71:72])
           meth.norm.centered[ix,72] = meth.norm[ix,72]-mean(meth.norm[ix,71:72])
           meth.norm.centered[ix,73] = meth.norm[ix,73]-mean(meth.norm[ix,73:74])
           meth.norm.centered[ix,74] = meth.norm[ix,74]-mean(meth.norm[ix,73:74])
           meth.norm.centered[ix,75] = meth.norm[ix,75]-mean(meth.norm[ix,75:76])
           meth.norm.centered[ix,76] = meth.norm[ix,76]-mean(meth.norm[ix,75:76])
           meth.norm.centered[ix,77] = meth.norm[ix,77]-mean(meth.norm[ix,77:78])
           meth.norm.centered[ix,78] = meth.norm[ix,78]-mean(meth.norm[ix,77:78])
           meth.norm.centered[ix,79] = meth.norm[ix,79]-mean(meth.norm[ix,79:80])
           meth.norm.centered[ix,80] = meth.norm[ix,80]-mean(meth.norm[ix,79:80])
}

meth.norm.centered = meth.norm.centered[complete.cases(meth.norm.centered),]

###########ColSideColors=clab
track=as.character(rnb.set.filtered.pheno$TP53)
track[track=="WT"]=1
track[track=="MU"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])

track=as.character(rnb.set.filtered.pheno$Response)
track[track=="NO"]=1
track[track=="YES"]=2
track=as.numeric(track)
colores=c("black","red")
clab2=as.character(colores[track])

png("heatmap_FDR5e-3_DIF25_MT_OR_WT_centeredOnPatient.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.3(as.matrix(meth.norm.centered),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),
srtCol=90,labRow = FALSE,xlab="", ylab="11808 CpGs",key.title="Methylation lvl",ColSideColors=cbind(TP53=clab,Response=clab2))
dev.off()


###########
ix = (wt_table$mean.diff<(-.25) & wt_table$diffmeth.p.adj.fdr<0.05) | (mu_table$mean.diff<(-.25) & mu_table$diffmeth.p.adj.fdr<0.05)

meth.filtered_diff = meth.filtered[ix, which(rnb.set.filtered.pheno$Treatment=="Treated")] - 
                     meth.filtered[ix, which(rnb.set.filtered.pheno$Treatment=="Untreated")]

meth.filtered_diff = meth.filtered_diff[complete.cases(meth.filtered_diff),]

track=as.character(rnb.set.filtered.pheno$TP53)[which(rnb.set.filtered.pheno$Treatment=="Treated")]
track[track=="WT"]=1
track[track=="MU"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])

track=as.character(rnb.set.filtered.pheno$Response)[which(rnb.set.filtered.pheno$Treatment=="Treated")]
track[track=="NO"]=1
track[track=="YES"]=2
track=as.numeric(track)
colores=c("black","red")
clab2=as.character(colores[track])

###########ColSide
png("heatmap_FDR5e-3_DIF25_MT_OR_WT_difference.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.3(as.matrix(meth.filtered_diff),col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),
srtCol=90,labRow = FALSE,xlab="", ylab="11772 CpGs",key.title="Methylation lvl",ColSideColors=cbind(TP53=clab,Response=clab2))
dev.off()

pdf("Labels.pdf")
plot(NULL)
legend("topright",legend=c("TP53 WT","TP53 MT","NoResponse","Response"),fill=c("#ffb3ba","#baffc9","black","red"), border=T, bty="n" )
dev.off()

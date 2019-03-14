### STD PIPE ###
library(RnBeads)
options(bitmapType="cairo")
save.rnb.set(rnb.set.filtered,path="/home/rtm/NEJM/RnBeads/rnb.set.filtered.RData")

rnb.set.untreated=remove.samples(rnb.set.filtered,samples(rnb.set.filtered)[which(rnb.set.filtered@pheno$Treatment=="Treated")])
rnb.set.treated=remove.samples(rnb.set.filtered,samples(rnb.set.filtered)[which(rnb.set.filtered@pheno$Treatment=="Untreated")])

dmc_untreated_byTP53 <- rnb.execute.computeDiffMeth(rnb.set.untreated,pheno.cols=c("TP53"))
dmc_treated_byTP53 <- rnb.execute.computeDiffMeth(rnb.set.treated,pheno.cols=c("TP53"))

comparison <- get.comparisons(dmc_untreated_byTP53)[1]
dmc_untreated_byTP53_table <-get.table(dmc_untreated_byTP53, comparison, "sites", return.data.frame=TRUE)
saveRDS(dmc_untreated_byTP53_table,"dmc_untreated_byTP53_table.rds")

comparison <- get.comparisons(dmc_treated_byTP53)[1]
dmc_treated_byTP53_table <-get.table(dmc_treated_byTP53, comparison, "sites", return.data.frame=TRUE)
saveRDS(dmc_treated_byTP53_table,"dmc_treated_byTP53_table.rds")
######################################################################################################################################
meth.filtered<-meth(rnb.set.filtered)
saveRDS(meth.filtered,"meth.filtered.rds")
saveRDS(rnb.set.filtered@pheno,"rnb.set.filtered.pheno.rds")
######################################################################################################################################
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(10))
######################################################################################################################################
meth.filtered.sig=meth.filtered[which(dmc_untreated_byTP53_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_untreated_byTP53_table$mean.diff)>.25),
                               which(rnb.set.filtered@pheno$Treatment=="Untreated")]

meth.filtered.sig = meth.filtered.sig[complete.cases(meth.filtered.sig),]

track=as.character(rnb.set.untreated@pheno$TP53)
track[track=="WT"]=1
track[track=="MU"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])

png("heatmap_FDR5e-3_DIF25_dmc_untreated_byTP53.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(meth.filtered.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),
srtCol=90,labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TP53 WT","TP53 MT"),fill=c("#ffb3ba","#baffc9"), border=T, bty="n" )
dev.off()
######################################################################################################################################
length(which(dmc_treated_byTP53_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_treated_byTP53_table$mean.diff)>.25))

meth.filtered.sig=meth.filtered[which(dmc_untreated_byTP53_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_untreated_byTP53_table$mean.diff)>.25),
                               which(rnb.set.filtered@pheno$Treatment=="Treated")]

meth.filtered.sig = meth.filtered.sig[complete.cases(meth.filtered.sig),]

track=as.character(rnb.set.treated@pheno$TP53)
track[track=="WT"]=1
track[track=="MU"]=2
track=as.numeric(track)
colores=c("#ffb3ba","#baffc9")
clab=as.character(colores[track])

png("heatmap_FDR5e-3_DIF25_dmc_in_Untreated_projected_treated_byTP53.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(meth.filtered.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),
srtCol=90,labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TP53 WT","TP53 MT"),fill=c("#ffb3ba","#baffc9"), border=T, bty="n" )
dev.off()
######################################################################################################################################
dmc_untreated_byTP53_table = readRDS("dmc_untreated_byTP53_table.rds")
meth.filtered = readRDS("meth.filtered.rds")
rnb.set.filtered.pheno = readRDS("rnb.set.filtered.pheno.rds")

meth.filtered.sig.untreated=meth.filtered[which(dmc_untreated_byTP53_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_untreated_byTP53_table$mean.diff)>.25),
                               which(rnb.set.filtered.pheno$Treatment=="Untreated")]

meth.filtered.sig.treated=meth.filtered[which(dmc_untreated_byTP53_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_untreated_byTP53_table$mean.diff)>.25),
                               which(rnb.set.filtered.pheno$Treatment=="Treated")]


pdf("primed_change_scatter.pdf")
ix = rnb.set.filtered.pheno$TP53=="WT"
plot(meth.filtered.sig.untreated[ix], meth.filtered.sig.treated[ix], xlab= "Untreated", ylab= "Decitabine Treated")
abline(0,1,col="blue")
ix = rnb.set.filtered.pheno$TP53=="MU"
points(meth.filtered.sig.untreated[ix], meth.filtered.sig.treated[ix],col="red")
legend("topleft",legend=c("TP53 WT","TP53 MT"),fill=c("black","red"), border=T, bty="n" )
abline(0,1,col="blue")
dev.off()

pdf("primed_change_density.pdf")
ix = rnb.set.filtered.pheno$TP53=="WT"
plot(density(meth.filtered.sig.untreated[ix]- meth.filtered.sig.treated[ix],na.rm=T))
ix = rnb.set.filtered.pheno$TP53=="MU"
lines(density(meth.filtered.sig.untreated[ix]- meth.filtered.sig.treated[ix],na.rm=T),col="red")
legend("topright",legend=c("TP53 WT","TP53 MT"),fill=c("black","red"), border=T, bty="n" )
dev.off()


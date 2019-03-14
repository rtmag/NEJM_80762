#######################################################################
# Paired Analysis
rnb.options("columns.pairing"=c("Treatment"="Patient"))
rnb.options("differential.variability"=FALSE)
rnb.options("covariate.adjustment.columns"=c("Treatment"))

dmc_p53 <- rnb.execute.computeDiffMeth(rnb.set.filtered,pheno.cols=c("TP53"))
comparison <- get.comparisons(dmc_p53)[1]
dmc_p53_table <-get.table(dmc_p53, comparison, "sites", return.data.frame=TRUE)

meth.filtered.sig = meth.filtered[which(dmc_p53_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_p53_table$mean.diff)>.25),]
meth.filtered.sig = meth.filtered.sig[complete.cases(meth.filtered.sig),]


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

source("https://raw.githubusercontent.com/rtmag/tumor-meth-pipe/master/heatmap3.R")
png("heatmap_FDR5e-3_DIF25_dmc_TP53_covarTREATMENT.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.3(as.matrix(meth.filtered.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),
srtCol=90,labRow = FALSE,xlab="", ylab="146 CpGs",key.title="Methylation lvl",ColSideColors=cbind(TP53=clab,Treatment=clab2))

legend("topright",legend=c("TP53 WT","TP53 MT","Untreated","Treated"),fill=c("#ffb3ba","#baffc9","black","red"), border=T, bty="n" )
dev.off()

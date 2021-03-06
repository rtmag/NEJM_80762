###############################################################################################################################
### STD PIPE ###
library(RnBeads)
options(bitmapType="cairo")

## Options and Parameters ##

#idat files
idat.dir <- file.path("/home/rtm/NEJM/GSE80762_RAW/idats")

# Sample annotation
sample.annotation <- file.path("/home/rtm/NEJM/NEJM_sample_sheet.csv")
rnb.options(import.table.separator=",")

# Report directory
system("rm -fr /home/rtm/NEJM/RnBeads/RnBeads_normalization")
report.dir <- file.path("/home/rtm/NEJM/RnBeads/RnBeads_normalization")

# Vanilla parameters
rnb.options(identifiers.column="Sample_ID")

# Multiprocess
num.cores <- 20
parallel.setup(num.cores)

data.source <- c(idat.dir, sample.annotation)
result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)
rnb.set.norm <- rnb.execute.normalization(result$rnb.set, method="swan",bgcorr.method="methylumi.noob")

save.rnb.set(rnb.set.norm,path="/home/rtm/NEJM/RnBeads/rnb.set.norm.RData")
###############################################################################################################################

meta = read.csv("../450k_patient_annotation.csv")
Patient = gsub("Patient_","",as.character(rnb.set.norm@pheno$Sample_ID))
Patient = gsub("_Treated","",Patient)
Patient = gsub("_Untreated","",Patient)
rnb.set.norm@pheno = cbind(rnb.set.norm@pheno,Patient=Patient)

rnb.set.norm@pheno = merge(rnb.set.norm@pheno,meta,by="Patient")

Treatment = gsub(".+\\_","",as.character(rnb.set.norm@pheno$Sample_ID),perl=TRUE)
rnb.set.norm@pheno = cbind(rnb.set.norm@pheno,Treatment=Treatment)
rnb.set.norm@pheno = rnb.set.norm@pheno[,c(2:10,1)]
###############################################################################################################################
# Filters
rnb.set.filtered <- rnb.execute.context.removal(rnb.set.norm)$dataset
rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.filtered)$dataset
rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
#######################################################################
# Paired Analysis
rnb.options("columns.pairing"=c("Treatment"="Patient"))
rnb.options("differential.variability"=FALSE)

dmc_treatment <- rnb.execute.computeDiffMeth(rnb.set.filtered,pheno.cols=c("Treatment"))

comparison <- get.comparisons(dmc_treatment)[1]
dmc_table <-get.table(dmc_treatment, comparison, "sites", return.data.frame=TRUE)
#######################################################################
#######################################################################

meth.norm<-meth(rnb.set.norm)

options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

track=as.character(rnb.set.norm@pheno$TP53)
track[track=="WT"]=1
track[track=="MU"]=2


track=as.numeric(track)
colores=c("#ffb3ba","#baffc9","#bae1ff")
clab=as.character(colores[track])

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))


meth.norm.sig=meth.norm[which(dmc_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_table$mean.diff)>.25),]
meth.norm.sig = meth.norm.sig[complete.cases(meth.norm.sig),]
##
# CEntering
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

##

meth.norm.sig=meth.norm.centered[which(dmc_table$diffmeth.p.adj.fdr<0.05 & abs(dmc_table$mean.diff)>.25),]
meth.norm.sig = meth.norm.sig[complete.cases(meth.norm.sig),]

png("heatmap_FDR5e-3_DIF25_Centered.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
x = heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
legend("topright",legend=c("TP53 WT","TP53 MT"),fill=c("#ffb3ba","#baffc9"), border=T, bty="n" )
dev.off()
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
# P53-WT TREATED VS UNTREATED
# 1) filtering out P53 MU samples
rnb.set.filtered_WT=remove.samples(rnb.set.filtered,samples(rnb.set.filtered)[which(rnb.set.filtered@pheno$TP53=="MU")])

# 2) filtering out P53 WT samples
rnb.set.filtered_MU=remove.samples(rnb.set.filtered,samples(rnb.set.filtered)[which(rnb.set.filtered@pheno$TP53=="WT")])

# 3) Differential methylation
rnb.options("columns.pairing"=c("Treatment"="Patient"))
rnb.options("differential.variability"=FALSE)

dmc_WT <- rnb.execute.computeDiffMeth(rnb.set.filtered_WT,pheno.cols=c("Treatment"))
dmc_MU <- rnb.execute.computeDiffMeth(rnb.set.filtered_MU,pheno.cols=c("Treatment"))

comparison <- get.comparisons(dmc_WT)[1]
wt_table <-get.table(dmc_WT, comparison, "sites", return.data.frame=TRUE)
saveRDS(wt_table,"wt_table.rds")

comparison <- get.comparisons(dmc_MU)[1]
mu_table <-get.table(dmc_MU, comparison, "sites", return.data.frame=TRUE)
saveRDS(mu_table,"mu_table.rds")
# 4) stats
wt_table = readRDS("wt_table.rds")
mu_table = readRDS("mu_table.rds")

options(scipen=999)
library(DESeq2)
library(gplots)
library(factoextra)
library(RColorBrewer)
library(graphics)

png("Volcano_wt.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)

plot(wt_table$mean.diff,-log10(wt_table$diffmeth.p.adj.fdr),xlab="Methylation Difference",
              ylab=expression('-Log'[10]*' q-values'),col=alpha("grey",.01),pch=20 )

  abline(v=-.2,lty = 2,col="grey")
  abline(v=.2,lty = 2,col="grey")
  abline(h=-log10(0.05),lty = 2,col="grey")
  points(wt_table$mean.diff[abs(wt_table$mean.diff)>.2 & wt_table$diffmeth.p.adj.fdr<0.05],
       -log10(wt_table$diffmeth.p.adj.fdr)[abs(wt_table$mean.diff)>.2 & wt_table$diffmeth.p.adj.fdr<0.05],
      col=alpha("red",.01),pch=20)
  legend("topright", paste("Gained Methylation",":",length(which(wt_table$mean.diff>.2 & wt_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
  legend("topleft", paste("Lost Methylation",":",length(which(wt_table$mean.diff<(-.2) & wt_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
dev.off()
############################
png("Volcano_mu.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
plot(mu_table$mean.diff,-log10(mu_table$diffmeth.p.adj.fdr),xlab="Methylation Difference",
              ylab=expression('-Log'[10]*' q-values'),col=alpha("grey",.01),pch=20 )

  abline(v=-.2,lty = 2,col="grey")
  abline(v=.2,lty = 2,col="grey")
  abline(h=-log10(0.05),lty = 2,col="grey")
  points(mu_table$mean.diff[abs(mu_table$mean.diff)>.2 & mu_table$diffmeth.p.adj.fdr<0.05],
       -log10(mu_table$diffmeth.p.adj.fdr)[abs(mu_table$mean.diff)>.2 & mu_table$diffmeth.p.adj.fdr<0.05],
      col=alpha("red",.01),pch=20)
  legend("topright", paste("Gained Methylation",":",length(which(mu_table$mean.diff>.2 & mu_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
  legend("topleft", paste("Lost Methylation",":",length(which(mu_table$mean.diff<(-.2) & mu_table$diffmeth.p.adj.fdr<0.05))), bty="n") 
dev.off()
#########################
pdf("density.pdf")
par(mfrow=c(1,1))
ix = mu_table$diffmeth.p.adj.fdr<0.05
plot(density(mu_table$mean.diff[ix]),col="red",main="Methylation Difference (FDR 5%)")

ix = wt_table$diffmeth.p.adj.fdr<0.05
lines(density(wt_table$mean.diff[ix]),col="blue")

legend("topright",c("p53 MT","p53 WT"),fill=c('red','blue'))
dev.off()

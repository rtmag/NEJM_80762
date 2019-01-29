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

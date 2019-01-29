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

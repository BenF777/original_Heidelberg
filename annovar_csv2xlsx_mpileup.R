options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]

require(openxlsx)

csv <- read.csv(file, stringsAsFactors=FALSE)

if(nrow(csv) != 0) {
otherInfo <- paste("DP=", unlist(lapply(strsplit(csv$Otherinfo, "DP="), '[[', 2)), sep="")

info <- unlist(lapply(strsplit(otherInfo, "\t"), '[[', 1))
infoSplit <-strsplit(info, ";")
dp <- unlist(lapply(strsplit(unlist(lapply(strsplit(info, "DP="), '[[', 2)), ";"), '[[', 1))
dp4 <-unlist(lapply(strsplit(unlist(lapply(strsplit(info, "DP4="), '[[', 2)), ";"), '[[', 1)) 
mq <- unlist(lapply(strsplit(unlist(lapply(strsplit(info, "MQ="), '[[', 2)), ";"), '[[', 1)) 
fq <- unlist(lapply(strsplit(unlist(lapply(strsplit(info, "FQ="), '[[', 2)), ";"), '[[', 1))
af1 <- unlist(lapply(strsplit(unlist(lapply(strsplit(info, "AF1="), '[[', 2)), ";"), '[[', 1))
ac1 <- unlist(lapply(strsplit(unlist(lapply(strsplit(info, "AC1="), '[[', 2)), ";"), '[[', 1))
rest <-unlist(lapply(strsplit(unlist(lapply(strsplit(info, "AC1="), '[[', 2)), ";"), '[[', 2))
format <- unlist(lapply(strsplit(otherInfo, "\t"), '[[', 2))
sample <- unlist(lapply(strsplit(otherInfo, "\t"), '[[', 3))
sample <- gsub("nan,nan,nan", ".", sample)

out <- csv
out$DP <- dp
out$DP4 <- dp4
out$MQ <- mq
out$FQ <- fq
out$AF1 <- af1
out$AC1 <- ac1
out$Other <- rest
out$format <-format
out$sample <- sample
out <- out[,!(names(out) %in% c("Otherinfo"))]
out[is.na(out)] <- "NA"

out <- out[,c(1:11,23:33,12:22)]

out.seq <- out
} else
{
    cols <-  c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", "cytoBand", "snp138", "LJB_PhyloP", "LJB_PhyloP_Pred", "LJB_SIFT", "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred", "LJB_GERP..", "X1000g2012apr_all", "cosmic68", "Otherinfo")
    out.seq <- as.data.frame(t(cols))[-1,]
    colnames(out.seq) <- cols
    
}

outFile <- paste(substr(file, 1, nchar(file)-3), "xlsx",sep="")
outFileCsv <- paste(substr(file, 1, nchar(file)-3), "_sort.csv",sep="")
write.csv2(out.seq, outFileCsv, row.names=FALSE)
write.xlsx(out.seq, file=outFile, rowNames=FALSE)


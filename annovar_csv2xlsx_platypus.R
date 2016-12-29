options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]

require(openxlsx)

library(GenomicFeatures)
library(ChIPpeakAnno)
#library(BSgenome.Hsapiens.ensemble.grch37)

csv <- read.csv(file, stringsAsFactors=FALSE)


if(nrow(csv) != 0) {
    
otherInfo <- csv$Otherinfo

info <- strsplit(otherInfo, ";")
dp.col <- grepl("TC=", info[[1]])
nf.col <- grepl("NF=", info[[1]])
nr.col <- grepl("NR=", info[[1]])
otherSample.col <- grepl("GT:", info[[1]])

dp <- unlist(lapply(strsplit(unlist(lapply(info, function(x) x[dp.col])), "TC="), function(x) x[2]))
nf <- unlist(lapply(strsplit(unlist(lapply(info, function(x) x[nf.col])), "NF="), function(x) x[2]))
nr <- unlist(lapply(strsplit(unlist(lapply(info, function(x) x[nr.col])), "NR="), function(x) x[2]))
dp4 <-paste("-,-,",nf,",",nr, sep="")

sampleInfo <-unlist(lapply(info, function(x) x[otherSample.col]))
format <- unlist(lapply(strsplit(sampleInfo, "\t"), '[[', 2))
sample <- unlist(lapply(strsplit(sampleInfo, "\t"), '[[', 3))
sample <- gsub("nan,nan,nan", ".", sample)

out <- csv
out$DP <- dp
out$DP4 <- dp4
out$MQ <- ""
out$FQ <- ""
out$AF1 <- ""
out$AC1 <- ""
out$Other <- otherInfo
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


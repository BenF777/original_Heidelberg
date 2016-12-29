options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

minDP <- as.integer(args[1])
minAF <- as.numeric(args[2])
minStrand <- as.numeric(args[3])
file <- args[4]

library(VariantAnnotation)

pre <- FilterRules(list(isLowCoverageExomeSnp = function(x) {
    !grepl("./.", x, fixed=TRUE) & grepl("PASS", x, fixed=TRUE)
}
))


 fileExtention <- paste(paste(paste(paste(paste(paste("_filtered_DP",minDP,sep=""), sep=""), minAF, sep=""), "", sep=""),"", sep=""), ".vcf", sep="")
 destination <- paste( substr(file, 1, nchar(file)-7), fileExtention,sep="")
 filt <- FilterRules(list(
     isDP = function(x) {
        info(x)$TC >= minDP      
     },
    isGQ = function(x) {
    (!is.na(geno(x)$GQ) & geno(x)$GQ >= 99)
    }
   ))

 filtered <- filterVcf(file, "Grch37d5", destination, filters=filt,  prefilters=pre)

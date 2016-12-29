options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
minDP <- as.integer(args[1])
file <- args[2]

library(VariantAnnotation)
                                        #library(cgdv17)
#setwd("/home/neuro/sequenzierung/NGS_DATA/tmp_L0002/out/var_gatk/")
#setwd("/home/neuro/sequenzierung/NGS_DATA/141007_NS500122_0005_AH128NBGXX/out/var_gatk")

#flist.Inter <- list.files(pattern = "_SNP_inter.vcf.gz$")
#flist.Compl <- list.files(pattern = "_SNP_complete.vcf.gz$")

pre <- FilterRules(list(isLowCoverageExomeSnp = function(x) {
    !grepl("./.", x, fixed=TRUE)
}
))

 destination <- paste( substr(file, 1, nchar(file)-7), "_filtered.vcf",sep="")#tempfile()
 print(file)
 filt <- FilterRules(list(
     isSNP = function(x) {
     geno(x)$GT != "0/0"
     },
     isDP = function(x) {
        info(x)$DP >= minDP        
     },
    isGQ80 = function(x) {
    (!is.na(geno(x)$GQ) & geno(x)$GQ >= 99)
    }
   ))

 filtered <- filterVcf(file, "Grch37d5", destination, filters=filt,  prefilters=pre)

    checkStrand = function(x){
        pv4 <- info(x)$PV4
        pValue <- unlist(lapply(pv4, '[[', 1))
        (geno(x)$GT == "0/1" & !is.na(pValue)) | (geno(x)$GT != "0/1")
    }

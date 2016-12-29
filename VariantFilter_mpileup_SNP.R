options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

minDP <- as.integer(args[1])
minAF <- as.numeric(args[2])
minStrand <- as.numeric(args[3])
file <- args[4]

library(VariantAnnotation)


pre <- FilterRules(list(isLowCoverageExomeSnp = function(x) {
    !grepl("./.", x, fixed=TRUE)
}
))


 fileExtention <- paste(paste(paste(paste(paste(paste("_filtered_DP",minDP,sep=""), "_AF", sep=""), minAF, sep=""), "", sep=""),"", sep=""), ".vcf", sep="")
 destination <- paste( substr(file, 1, nchar(file)-7), fileExtention,sep="")#tempfile()

filt <- FilterRules(list(
     isSNP = function(x) {
     geno(x)$GT != "0/0"
     },
     isDP = function(x) {
        info(x)$DP >= minDP      
     },
    isGQ80 = function(x) {
    (!is.na(geno(x)$GQ) & geno(x)$GQ >= 80)
    },
     checkAllelFrequency = function(x){
        dp4 <- info(x)$DP4
        ref.for <- unlist(lapply(dp4, '[[', 1))
        ref.rev <- unlist(lapply(dp4, '[[', 2))
        alt.for <- unlist(lapply(dp4, '[[', 3))
        alt.rev <- unlist(lapply(dp4, '[[', 4))
        ratioRef.for <- ref.for/(ref.for+ref.rev)
        ratioRef.rev <- ref.for/(ref.for+ref.rev)
        ratioAlt.for <- alt.for/(alt.for+alt.rev)
        ratioAlt.rev <- alt.rev/(alt.for+alt.rev)
      # fisher.pValues <- unlist(lapply(dp4, function(x){fisher.test(matrix(x,2))$p.value}))
        af <- (alt.for + alt.rev) /(ref.for+ref.rev+alt.for + alt.rev)
        ratioRef <- if((ref.for+ref.rev)>=10) {
            ratioAlt.for >= minStrand & ratioAlt.rev >= minStrand}
            else {
                TRUE
            }
       (geno(x)$GT == "0/1" & !is.na(af) & af >= minAF & ratioAlt.for >= minStrand & ratioAlt.rev >= minStrand) | (geno(x)$GT != "0/1")
#        (geno(x)$GT == "0/1" & !is.na(af) & af >= minAF) | (geno(x)$GT != "0/1")
     }
   ))

file.remove(destination)
 filtered <- filterVcf(file, "Grch37d5", destination, filters=filt)

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


fileExtention <- paste(paste(paste(paste(paste(paste("_filtered_DP",minDP,sep=""), "_AF", sep=""), minAF, sep=""), "", sep=""),"", sep=""), ".vcf", sep="")
destination <- paste( substr(file, 1, nchar(file)-7), fileExtention,sep="")
filt <- FilterRules(list(
    isSNP = function(x) {
        geno(x)$GT != "0/0"
    },
    isDP = function(x) {
        info(x)$TC >= minDP      
    },
    isGQ = function(x) {
        (!is.na(geno(x)$GQ) & geno(x)$GQ >= 80)
    },
    checkAllelFrequency = function(x){
        dp <- info(x)$TC
        alt.for <- unlist(info(x)$NF)
        alt.rev <- unlist(info(x)$NR)
      
        ratioAlt.for <- alt.for/(alt.for+alt.rev)
        ratioAlt.rev <- alt.rev/(alt.for+alt.rev)
      
        af <- (alt.for + alt.rev) / dp

       (geno(x)$GT == "0/1" & !is.na(af) & af >= minAF & ratioAlt.for >= minStrand & ratioAlt.rev >= minStrand) | (geno(x)$GT != "0/1")
    }
    ))

filtered <- filterVcf(file, "Grch37d5", destination, filters=filt,  prefilters=pre)

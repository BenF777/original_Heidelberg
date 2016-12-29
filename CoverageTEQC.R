library(TEQC)
library(Rsamtools)
library(GenomicAlignments)
data(genesymbol, package = "biovizBase")

args <- commandArgs(trailingOnly = TRUE)
print(args)

file <- args[3]

out <- paste(args[1], "/", sep="")

bedTargetFile <- args[2]
offset <- 50
targets <- get.targets(targetsfile=bedTargetFile)

outname <- tail(unlist(strsplit(file, "/")), n=1)
outname <- substr(outname, 1, nchar(outname)-4)

sink(file = paste(paste(out, outname, sep=""), ".log", sep=""), append = FALSE, type = c("output", "message"), split = FALSE)
outname

bam <-  as(readGAlignments(file), "GRanges")
bamPos <- bam[which(bam@strand == "+")]
bamNeg <- bam[which(bam@strand == "-")]

reads <- as(bam, "RangedData")
readsPos <- as(bamPos, "RangedData")
readsNeg <- as(bamNeg, "RangedData")

fr <- fraction.reads.target(reads, targets, Offset=offset)
print(paste(paste("Fraction=", round(fr*100, digits=2), sep=""), "%", sep=""))

Coverage <- coverage.target(reads, targets, perTarget=T, perBase=T)
CoveragePos <- coverage.target(readsPos, targets, perTarget=T, perBase=T)
CoverageNeg <- coverage.target(readsNeg, targets, perTarget=T, perBase=T)

print("avgTargetCoverage")
Coverage$avgTargetCoverage
print("avgTargetCoverageSD")
Coverage$targetCoverageSD
print("avgTargetCoverageQuantiles")
Coverage$targetCoverageQuantiles

print("avgTargetCoverage_Pos")
CoveragePos$avgTargetCoverage
print("avgTargetCoverageSD_Pos")
CoveragePos$targetCoverageSD
print("avgTargetCoverageQuantiles_Pos")
CoveragePos$targetCoverageQuantiles

print("avgTargetCoverage_Neg")
CoverageNeg$avgTargetCoverage
print("avgTargetCoverageSD_Neg")
CoverageNeg$targetCoverageSD
print("avgTargetCoverageQuantiles_Neg")
CoverageNeg$targetCoverageQuantiles

coverK <- covered.k(Coverage$coverageTarget, k=c(10, 50, 80, 120))
coverK.Pos <- covered.k(CoveragePos$coverageTarget, k=c(10, 50, 80, 120))
coverK.Neg <- covered.k(CoverageNeg$coverageTarget, k=c(10, 50, 80, 120))
coverK.DF <- data.frame(coverK)
coverK.DF <- cbind(coverK.DF, coverK.Pos)
coverK.DF <- cbind(coverK.DF, coverK.Neg)

write.table(coverK, file=paste(paste(out, outname, sep=""), "_fraction.csv", sep=""), sep="\t", row.names=T, quote=F, col.names=T)

png(file.path(out, paste(outname, "_hist.png", sep="")), width = 1000, height = 1000, units = "px")
coverage.hist(Coverage$coverageTarget, covthreshold=80)
dev.off()

png(file.path(out, paste(outname, "_uniformity.png", sep="")),  width = 1000, height = 1000, units = "px")
coverage.uniformity(Coverage)
dev.off()

coverageTargets.DF <- as.data.frame(Coverage$targetCoverages)
covPos.DF <- as.data.frame(CoveragePos$targetCoverages)
covNeg.DF <- as.data.frame(CoverageNeg$targetCoverages)

coverageTargets.DF$avgCoverage.Pos <- covPos.DF$avgCoverage
coverageTargets.DF$coverageSD.Pos <- covPos.DF$coverageSD

coverageTargets.DF$avgCoverage.Neg <- covNeg.DF$avgCoverage
coverageTargets.DF$coverageSD.Neg <- covNeg.DF$coverageSD

write.table(coverageTargets.DF, file=paste(paste(out, outname, sep=""), "_targetCoverage.csv", sep=""), sep="\t", row.names=F, quote=F)

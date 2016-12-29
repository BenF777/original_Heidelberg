library(seqCNA)
library(openxlsx)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
print(args)

result.path <- paste(args[1],"/",sep="")

out.path <- paste(args[2],"/",sep="")

window <- 100
type <- args[3]

sample <- args[4]
sample.file <- args[5]
lib <- args[6]

out.filter.path <- paste(out.path, sample,sep="")

dir.create(out.filter.path,showWarning=FALSE)

#out.filter.path <- paste(out.filter.path,"/",sep="")

#TODO own reference 
references <- data.frame(type=c("NPHD2015A","MNG2015A","EXOM_V5"), id= c("69001-DNA-BLUT","69001-DNA-BLUT","21823-DNA-BLUT"), lib=c("Agilent_L0016","Agilent_L0016","Agilent_L0029"))
references$window <- window
references$path <- paste(result.path, references$id,"/",references$lib,"/cnv/", references$id,"_",references$type,"_",references$lib,"_w",references$window,".txt" ,sep="")

reference <- references[grepl(type,references$type),]
                          
out <- paste(out.path, sample,"_",type,"_",lib,"_w",window,sep="" )

print(paste("Process ",sample, "-", type, " ", lib, sep=""))

out.file <- paste(out,".txt",sep="")
pdf.file <- paste(out,"_ref_",reference$id,".pdf",sep="")
png.file <- paste(out,"_ref_",reference$id,".png",sep="")

if(!file.exists(out.file)) {
  runSeqsumm(summ.win=window, file=sample.file, folder=NULL, output.file=out.file, samtools.path="samtools")
}

rco <- readSeqsumm(build="",tumour.data=read.table(out.file, header=TRUE),normal.data=read.table(reference$path, header=TRUE))

rco <- applyFilters(rco, trim.filter=1, mapq.filter=4, plots=FALSE, folder=out.filter.path)
rco <- runSeqnorm(rco, plots=TRUE, folder=out.filter.path)
rco <- runGLAD(rco, nproc=8)

testout <- rco@output
testout$chrom <- as.character(testout$chrom)
testout$chrom <- gsub("chr","",testout$chrom)
testout$win <- as.numeric(as.character(testout$win.start))
testout$x <- 1:length(testout$win)

testout <- testout[!testout$chrom %in% c("X","Y"),]

linesChr <- unlist(lapply(split(testout,testout$chrom),function(x) nrow(x)))
ordChr <- order(as.numeric(names(linesChr)))
linesChr <- linesChr[ordChr]
linesChrom <- linesChr
for (i in 1:length(linesChrom)) {
    if (i>1) {
        linesChrom[i] <- linesChrom[i-1]+linesChr[i]
    }else {
        linesChrom[i] <- linesChr[i]
    }
}

breaks <- linesChrom
for (i in 1:length(linesChr)) {
    if (i>1) {
        breaks[i] <- linesChrom[i-1]+round(linesChr[i]/2)
            }else {
                breaks[i] <- round(linesChr[i]/2)
            }
}

testout$normalized[testout$normalized > 4] <- NA
title <- paste(sample,"   Type: ", type, "   Lib: ", lib, "   Window: ", window, "   Ref: ", reference$id ,sep="")
plotout <-
    ggplot(testout,aes(x, normalized)) + geom_point(aes(color=normalized))+ scale_colour_gradient2(low = "red", high="green", mid="lightgrey",midpoint=0) +  geom_point(aes(x,segmented),color="#FF0000") + theme_bw() + scale_y_continuous(limit=c(-2.25,2.25)) +  scale_x_continuous(expand=c(0,0), breaks=breaks, labels=c(1:22)) + geom_vline(xintercept=linesChrom) + xlab("Chromosom") + ylab("CNV") + theme(axis.ticks.x=element_blank(), panel.grid=element_blank()) + ggtitle(title) + geom_hline(yintercept=0, color="darkgrey")

ggsave(pdf.file,plotout,width=20,height=10)
ggsave(png.file,plotout,width=20,height=10)  

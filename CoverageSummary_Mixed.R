args <- commandArgs(trailingOnly = TRUE)
print(args)

require(openxlsx)

targetLocation <- args[1]
targetFiles <- read.csv(paste(targetLocation,"/","targetFiles.csv", sep=""), stringsAsFactors=FALSE)
targetFiles$TargetFile <- paste(targetLocation,"/",targetFiles$TargetFile,sep="")

folder <-  paste(args[2], "/", sep="")
                                  
samples <- read.csv(args[3], header=FALSE)
colnames(samples) <- c("ID","Barcode", "Panel")
                                  
filesComplete <- list.files(folder, pattern ="*_targetCoverage.csv")

samplesByPanel <- split(samples, samples$Panel)

for (panelName in names(samplesByPanel)) {

    files <- filesComplete[grepl(paste(samplesByPanel[[panelName]]$ID,collapse="|"), filesComplete)]
    targetFile <- targetFiles[targetFiles$Panel == panelName,]$TargetFile


    data <- read.table(paste(folder,files[1], sep=""), stringsAsFactors=FALSE, skip=1)
    colnames(data) <- c("chr","start","end","width","x","y", "x.Pos", "y.Pos", "x.Neg", "y.Neg")
    data <- data[,!(names(data) %in% c("x","y", "x.Pos","y.Pos","x.Neg","y.Neg"))]
    header <- data
    for (file in files) { 
        table <- read.table(paste(folder,file, sep=""), stringsAsFactors=FALSE, skip=1)
        avgname <- paste("AVG", file, sep="_")
        sdname <- paste("SD", file, sep="_")
        avgname.Pos <- paste("AVG.Pos", file, sep="_")
        sdname.Pos <- paste("SD.Pos", file, sep="_")
        avgname.Neg <- paste("AVG.Neg", file, sep="_")
        sdname.Neg <- paste("SD.Neg", file, sep="_")
        
        colnames(table) <- c("chr","start","end","width",avgname,sdname,avgname.Pos,sdname.Pos,avgname.Neg,sdname.Neg)
        data[avgname] <- table[avgname]
        data[sdname] <- table[sdname]
        data[avgname.Pos] <- table[avgname.Pos]
        data[sdname.Pos] <- table[sdname.Pos]
        data[avgname.Neg] <- table[avgname.Neg]
        data[sdname.Neg] <- table[sdname.Neg]
    }

    str(data)
    f <- function(x) {
        chr <- x[1]
        start <- x[2]
        end <- x[3]
        width <- x[4]
        avg <- as.numeric(x[grep("AVG_",names(x))])
        avgSum <- summary(avg)
        avgSum <- avgSum[c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")]
        avgPosSum <- summary(as.numeric(x[grep("AVG.Pos_",names(x))]))
        avgPosSum <- avgPosSum[c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")]
        avgNegSum <- summary(as.numeric(x[grep("AVG.Neg_",names(x))]))
        avgNegSum <- avgNegSum[c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")]
        res <- c(avgSum, avgPosSum, avgNegSum)
        names(res) <- c("Comp_Min.", "Comp_1st_Qu.", "Comp_Median", "Comp_Mean", "Comp_3rd_Qu.", "Comp_Max", "Pos_Min.", "Pos_1st_Qu.", "Pos_Median", "Pos_Mean", "Pos_3rd_Qu.", "Pos_Max",  "Neg_Min.", "Neg_1st_Qu.", "Neg_Median", "Neg_Mean", "Neg_3rd_Qu.", "Neg_Max")
        data.frame(as.list(res))
    }

    data.summary <- apply(data, 1, f)
    data.summary <- do.call("rbind", data.summary)

    header <- cbind(data[,1:4], data.summary)

    print("summary")

    print(targetFile)
    targets <- read.table(targetFile)
    print("read")
                                       
    colnames(targets) <- c("chr", "start", "end", "ID")#, "nr")
    targets$start <- targets$start+1
    out <- data.frame()
    out <- merge(x = header, y = targets, by = c("chr", "start", "end"), all.x=TRUE)
                                   

    write.csv(out, file=paste(folder,"/",panelName,"_", "SummaryTargetPositions.csv", sep=""), row.names=FALSE)
    

    filesLog <- list.files(folder, pattern ="*.log")
    filesLog <- filesLog[grepl(paste(samplesByPanel[[panelName]]$ID,collapse="|"), filesLog)]
    filesLog <- paste(folder, filesLog, sep="")
    filesLog <- data.frame(filesLog)

    filesearch <- function (x) {

        f <- readLines(x)

        sample <- gsub("\"", "", f[1])
        sample <- substring(sample, 5)
        sample <- gsub("_MERGED","", sample)

        fr <- gsub("\"", "", f[2])
        fr <- substring(fr, 5)
        fr <- gsub("Fraction=","", fr)

        avgTarget <- gsub("\"", "", f[4])
        avgTarget <- substring(avgTarget, 5)
        avgTarget <- gsub("[[1]] ","", avgTarget)

        avgSD <- gsub("\"", "", f[6])
        avgSD <- substring(avgSD, 5)
        avgSD <- gsub("[[1]] ","", avgSD)

        avgQuantiles <- unlist(strsplit(f[9]," "))
        avgQuantiles <- avgQuantiles[avgQuantiles!=""]

        avgQuantiles.Pos <- unlist(strsplit(f[16]," "))
        avgQuantiles.Pos <- avgQuantiles.Pos[avgQuantiles.Pos!=""]

        avgQuantiles.Neg <- unlist(strsplit(f[23]," "))
        avgQuantiles.Neg <- avgQuantiles.Neg[avgQuantiles.Neg!=""]

        avgTarget.Pos <- gsub("\"", "", f[11])
        avgTarget.Pos <- substring(avgTarget.Pos, 5)
        avgTarget.Pos <- gsub("[[1]] ","", avgTarget.Pos)

        avgSD.Pos <- gsub("\"", "", f[13])
        avgSD.Pos <- substring(avgSD.Pos, 5)
        avgSD.Pos <- gsub("[[1]] ","", avgSD.Pos)

        avgTarget.Neg <- gsub("\"", "", f[18])
        avgTarget.Neg <- substring(avgTarget.Neg, 5)
        avgTarget.Neg <- gsub("[[1]] ","", avgTarget.Neg)

        avgSD.Neg <- gsub("\"", "", f[20])
        avgSD.Neg <- substring(avgSD.Neg, 5)
        avgSD.Neg <- gsub("[[1]] ","", avgSD.Neg)

        c(sample, fr, avgTarget, avgSD, avgTarget.Pos, avgSD.Pos, avgTarget.Neg, avgSD.Neg, avgQuantiles, avgQuantiles.Pos,avgQuantiles.Neg )
    }

    overall.out  <- apply(filesLog, 1, filesearch)
    overall.out <- t(overall.out)
    colnames(overall.out) <- c("Sample", "OnTarget", "Comp_Mean_coverage", "Comp_Mean SD", "Pos_Mean_coverage", "Pos_Mean SD","Neg_Mean_coverage", "Neg_Mean SD",  "Com_Min", "Com_1st_Qu.", "Com_Median", "Com_3rd_Qu.", "Com_Max", "Pos_Min", "Pos_1st_Qu.", "Pos_Median", "Pos_3rd_Qu.", "Pos_Max", "Neg_Min", "Neg_1st_Qu.", "Neg_Median", "Neg_3rd_Qu.", "Neg_Max")

    write.csv(overall.out, paste(folder,"/",panelName,"_", "SummarySamples.csv", sep=""), row.names=FALSE)
    
    write.xlsx(out, file=paste(folder,"/",panelName,"_", "SummaryTargetPositions.xlsx", sep=""), row.names=FALSE)
    write.xlsx(overall.out, paste(folder,"/",panelName,"_", "SummarySamples.xlsx", sep=""), row.names=FALSE)
}

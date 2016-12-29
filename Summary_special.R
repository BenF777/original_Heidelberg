library(openxlsx)

options(echo=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

call.path <- args[1]
out.path <- args[2]
panel <- args[3]
out.file <- args[4]

#TODO run for samples (split, lapply)
files <- list.files(call.path,pattern="*vcf")
files <- paste(call.path,files,sep="")

files <- files[grepl(out.file,files)]

results <- data.frame(sample="",lib="",chr="",pos="", ref="",alt="",dp="",dp.rf="",dp.rr="",dp.af="", dp.ar="", special.pos="",stringsAsFactors=FALSE)
for (i in files) {
    file <- try(read.table(i,stringsAsFactors=FALSE))
    if (!inherits(file, 'try-error')){
            dp.info <- strsplit(file$V8,";")
        dp <- unlist(lapply(dp.info,function(x) {as.numeric(substring(x[grepl("DP=",x)],4))}))
            dp4 <- lapply(dp.info,function(x) { strsplit(substring(x[grepl("DP4", x)],5),",")})
            dp4.rf <- unlist(lapply(dp4, function(x) x[[1]][[1]]))
            dp4.rr <-  unlist(lapply(dp4, function(x) x[[1]][[2]]))
             dp4.af <- unlist(lapply(dp4, function(x) x[[1]][[3]]))
            dp4.ar <-  unlist(lapply(dp4, function(x) x[[1]][[4]]))

            path.parts <- unlist(strsplit(i,"/"))
            path.parts <- path.parts[length(path.parts)]
            sample <- strsplit(path.parts, "_")[[1]][1]
             specialPos <- gsub(paste(panel,"_",sep=""),"",i)
            specialPos <- gsub(".vcf","",specialPos)
            specialPos <- gsub(paste(sample,"_",sep=""),"",specialPos)
            dat <- data.frame(sample=sample,lib=panel,file$V1, file$V2, file$V4, file$V5,dp,dp4.rf,dp4.rr,dp4.af,dp4.ar,specialPos)
            colnames(dat) <- colnames(results)
        results <- rbind(results,dat)
    }else print(paste("Error: " ,i))
}
results <- results[-1,]
results$dp <- as.numeric(results$dp)
results$dp.rf <- as.numeric(results$dp.rf)
results$dp.rr <- as.numeric(results$dp.rr)
results$dp.af <- as.numeric(results$dp.af)
results$dp.ar <- as.numeric(results$dp.ar)

write.xlsx(results,paste(out.path,"/",results$sample[1],"_",panel,"_special.xlsx",sep=""))


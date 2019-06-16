########## maxPeak

options(stringsAsFactors=FALSE)

# Define folder containing .wig files and genome annotations
# start in the working directory this script is located in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder = "./"

# Set location and naming of .wig files to be processed, containing raw data
inFiles = c("Fkh1_log"=paste0(folder,"Fkh1_log.wig"),
            "Fkh1_stat"=paste0(folder,"Fkh1_stat.wig"),
            "Fkh2_log"=paste0(folder,"Fkh2_log.wig"),
            "Fkh2_stat"=paste0(folder,"Fkh2_stat.wig"))


# Import annotations and chromosomes
genome = list()
genome[["ORF"]] = read.delim(paste0(folder,"sacCer3.bed"),sep="\t",header=TRUE)
colnames(genome[["ORF"]]) = c("Chrom","Strand","ORF_st","ORF_end","Gene","ComName")
genome[["chromSize"]] = read.delim(paste0(folder,"sacCer3.chrom.sizes"),sep="\t",header=FALSE)


# Import .wig data
wigDat <- list()
for (datIn in names(inFiles)) {
  readsIn <- read.delim(as.character(inFiles[datIn]),sep=" ",header=TRUE)
  chroms <- unique(readsIn[grep("chrom=",readsIn[,2]),2])
  readsVarstep <- rownames(readsIn[readsIn[,"track"]=="variableStep",])
  for (y in 1:length(chroms)) {
    if (y < length(chroms)) {
      chromSpan <- as.character(seq(as.numeric(rownames(readsIn[readsIn[,2]==chroms[y],][1,]))+1,
                                    as.numeric(rownames(readsIn[readsIn[,2]==chroms[y+1],][1,]))-1))
    } else if (y == length(chroms)) {
      chromSpan <- as.character(seq(as.numeric(rownames(readsIn[readsIn[,2]==chroms[y],][1,]))+1,nrow(readsIn)))
    }
    chromSpan <- chromSpan[!(chromSpan %in% readsVarstep)]
    wigDat[[datIn]][[unlist(strsplit(chroms[y],"="))[2]]] = readsIn[chromSpan,]
  }
}

# Create promoter positions genome-wide, defined as 1000bp upstream of ORF start
distPromoter <- 1000
for (gen in unique(genome[["ORF"]][,"Gene"])) {
  genTmp <- genome[["ORF"]][genome[["ORF"]][,"Gene"]==gen,]
  for (chr in unique(genTmp[,"Chrom"])) {
    chrSize <- genome[["chromSize"]][genome[["chromSize"]][,"V1"]==chr,"V2"]
    genTmpChr <- genTmp[genTmp[,"Chrom"]==chr,]
    genome[[chr]][[gen]][["prom"]] <- c()
    for (x in 1:nrow(genTmpChr)) {
      if (genTmpChr[x,"Strand"]=="+") {
        genome[[chr]][[gen]][["prom"]] <- c(genome[[chr]][[gen]][["prom"]],
                                                   seq(max(1,genTmpChr[x,"ORF_st"]-distPromoter),genTmpChr[x,"ORF_st"]))
      } else if (genTmpChr[x,"Strand"]=="-") {
        genome[[chr]][[gen]][["prom"]] <- c(genome[[chr]][[gen]][["prom"]],
                                                   seq(genTmpChr[x,"ORF_end"],min(genTmpChr[x,"ORF_end"]+distPromoter,chrSize)))
      }
    }
  }
}


# Extract peak maxima for each promoter
maxPeak = list()
for (datIn in names(inFiles)) {
  for (gen in unique(genome[["ORF"]][,"Gene"])) {
    genTmp <- genome[["ORF"]][genome[["ORF"]][,"Gene"]==gen,]
    for (chr in unique(genTmp[,"Chrom"])) {
      # get read counts upstream of gene ORF: dataframe with rows for each base pair
      countsTmp = wigDat[[datIn]][[chr]][wigDat[[datIn]][[chr]][,"track"] %in% 
                                           as.character(genome[[chr]][[gen]][["prom"]]),,drop=F]
      if (nrow(countsTmp)>0) {
        # find maximum read count bp, round and calculate the median to obtain a single location if multiple bp have the same number of reads
        peakPos = round(median(as.numeric(countsTmp[as.numeric(countsTmp[,2])==max(as.numeric(countsTmp[,2])),1])))
        maxPeak[[datIn]] = rbind(maxPeak[[datIn]],cbind("Chromosome"=chr,
                                                        "peak-5"=peakPos-5,
                                                        "peak+5"=peakPos+5,
                                                        "Gene"=gen,
                                                        "S/N"=max(as.numeric(countsTmp[,2]))))
      } else {
        # apparently no reads were detected so set peakPos to NA and SNR = 0
        maxPeak[[datIn]] = rbind(maxPeak[[datIn]],cbind("Chromosome"=chr,
                                                        "peak-5"=NA,
                                                        "peak+5"=NA,
                                                        "Gene"=gen,
                                                        "S/N"=0))
      }
    }
  }
}

# Calculate quantiles and export maxPeak value relative to selected quantile (signal/noise) value to table
quantiles = list()
relQuant = "65%"
for (datIn in names(inFiles)) {
  quantiles[[datIn]] = quantile(as.numeric(maxPeak[[datIn]][maxPeak[[datIn]][,"S/N"]>0,"S/N"]),probs=seq(0.1,1,0.05))
  maxOut = maxPeak[[datIn]]
  maxOut[,"S/N"] = as.numeric(maxOut[,"S/N"]) / as.numeric(quantiles[[datIn]][relQuant])
  write.table(maxOut,
              file=paste0(folder,datIn,"_maxPeak.bed"),sep="\t",row.names=F, col.names=F, quote=F)
}
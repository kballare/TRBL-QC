args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript endDNA.R [path to flagstats files] [path to seqprep files] [path to lgdist folders]\n")
  cat("Example: Rscript endDNA.R ./BWA_analyses/ ./SeqPrep_output/ ./MapDamage_output/\n")
  quit()
}
#paths....
filepathfs=args[1]#c("../../alisa_horses/beringia_horses_1aug16/BWA_analyses/")
filesfs <- list.files(path=filepathfs, pattern="*.sorted.flagstats.*", full.names=T, recursive=FALSE)
#complexity info
filesRMDUP <- list.files(path=filepathfs, pattern="*.sorted.rmdup.flagstats.*", full.names=T, recursive=FALSE)
filepathsp=args[2]#c("../../alisa_horses/beringia_horses_1aug16/SeqPrep_output/")
filessp <- list.files(path=filepathsp, pattern="*SeqPrep_output.txt", full.names=T, recursive=FALSE)
folderpathlgdist=args[3]#c("../../alisa_horses/beringia_horses_1aug16/MapDamage_output/")#args[3] #lgdist folders here
lgdistpath <- list.dirs(path = folderpathlgdist, full.names = TRUE, recursive = TRUE)
#filepathMIA=args[4]#c("../../alisa_horses/beringia_horses_1aug16/MIA_analyses/")
#filesMIA <- list.files(path=filepathMIA, pattern="*.mia_stats.*", full.names=T, recursive=FALSE)
#AV148_S89_all_seqprep.complexity_filtered.duplicates_removed.EqCab_mt.20M.maln.F.mia_stats.txt
cat ("checking number of found files\n")
if (length(filessp) == length(filesfs) && length(filesfs) == length(lgdistpath)-1)
{
  cat (paste0("OK seqpreps:", length(filessp), ", flagstats:", length(filesfs),", lgdistfolders:", length(lgdistpath)-1, "\n"))
  #create output matrix
  output <- matrix(ncol=20, nrow=length(filesfs))
  #go through each sample
  i = 0
  for (f in filesfs)
  {
    i = i + 1
    #obtain sample name
    filesname <- unlist(strsplit(f, "/"))
    samplename <- unlist(strsplit(unlist(filesname[length(filesname)]), "_"))
    samplen = paste(samplename[1], samplename[2],samplename[3],samplename[4],sep="_")
    cat(paste0("sample name:", samplen, "\n"))
    #just to make sure on correct pairing of flagstats and seqprep files:
    cat(paste0("flagstats filename:", filesfs[grep(samplen, filesfs)], "\n"))
    cat(paste0("SeqPrep filename:", filessp[grep(samplen, filessp)], "\n"))
    #cat(paste0("MIA filename:", filesMIA[grep(samplen, filesMIA)], "\n"))
    cat(paste0("RMDUP filename:", filesRMDUP[grep(samplen, filesRMDUP)], "\n"))
    #look up files with sample name and read them
    seqprep <- read.table(filessp[grep(samplen, filessp)], fill=TRUE, header = FALSE, sep = ":")
    flagstat <- read.table(filesfs[grep(samplen, filesfs)], fill = TRUE, header = FALSE, sep = "+")
    #mia <- read.csv(filesMIA[grep(samplen, filesMIA)], fill=TRUE, header = FALSE, sep = ":")
    rmdup <- read.table(filesRMDUP[grep(samplen, filesRMDUP)], fill = TRUE, header = FALSE, sep = "+")
    #convert flagstats first column to numeric since the last row has chars...
    flagstat$V1 <- as.numeric(as.character(flagstat$V1))
    #same for mia
    #mia$V2 <- suppressWarnings(as.numeric(as.character(mia$V2)))
    #same for RMDUP
    rmdup$V1 <- as.numeric(as.character(rmdup$V1))

    #same for mia
    #mia$V2 <- as.numeric(as.character(mia$V2))
    #look up for lgdistribution.txt
    #lgdistpath <- c("/home/alex/Documents/university/research/alisa_horses/beringia_horses_1aug16/MapDamage_output/mapDamage_AV34_S4/lgdistribution.txt") 
    #lgdistpath <- list.dirs(path = "/home/alex/Documents/university/research/alisa_horses/beringia_horses_1aug16/MapDamage_output/", full.names = TRUE, recursive = TRUE)
    lgdistfilename <- paste0(lgdistpath[grep(samplen, lgdistpath)],"/lgdistribution.txt")
    
    #alisa's lgdistribution script
    dat <- read.table(lgdistfilename, header=T)
    mean = sum(dat$Length * dat$Occurences) / sum(dat$Occurences)
    SD = sqrt(sum((dat$Length - mean)**2 * dat$Occurences) / (sum(dat$Occurences)-1))
    
    
    #alisa's script part
    RawReads <- seqprep$V2[2]
    Discarded <- seqprep$V2[5]
    NumUse <- (RawReads-Discarded)
    NumMer <- seqprep$V2[3]
    PropMer <-  NumMer/NumUse
    PropUnmer <- ((NumUse-NumMer)/NumUse)
    AllMapped <- flagstat$V1[5] #corrected to 3
    PEmapped <- flagstat$V1[6] #corrected to 4
    NumMapMer <-(AllMapped-PEmapped)
    NumMapUnmer <- PEmapped/2
    PropMapMer <- NumMapMer/NumMer
    PropMapUnmer <- (NumMapUnmer/ (NumUse-NumMer))
    MerCont <- PropMer*PropMapMer
    UnmerCont <- PropUnmer*PropMapUnmer
    Endog <- MerCont+UnmerCont
    #MiaCov <- mia$V2[5]
    #MiaNumMap <- mia$V2[3]
    DedupMap <-  rmdup$V1[5]
    ComplexProp <- DedupMap/AllMapped
    
    #print(c(samplen, RawReads, Discarded, NumUse, NumMer, PropMer, PropUnmer, AllMapped, PEmapped, NumMapMer, NumMapUnmer, PropMapMer, PropMapUnmer, MerCont, UnmerCont, Endog))
    cat(paste0("Endog:", Endog,"\n"))
    #put results into the matrix
    output[i,] <- c(samplen, RawReads, Endog, ComplexProp, Discarded, NumUse, NumMer, PropMer, PropUnmer, AllMapped, PEmapped, DedupMap, NumMapMer, NumMapUnmer, PropMapMer, PropMapUnmer, MerCont, UnmerCont, mean, SD)
  }
} else {
  cat ("Error, different number of files/folders found\n")
  cat (paste0("seqpreps:", length(filessp), ", flagstats:", length(filesfs),", lgdistfolders:", length(lgdistpath)-1, "\n"))
  quit()
}
#transform matrix to dataframe
output <- data.frame(output)
#label columns
colnames(output) <- c("Sample", "RawReads", "Endog", "ComplexProp", "Discarded", "NumUse", "NumMer", "PropMer", "PropUnmer", "AllMapped", "PEmapped", "DedupMap", "NumMapMer", "NumMapUnmer", "PropMapMer", "PropMapUnmer", "MerCont", "UnmerCont", "lgdist mean", "lgdist SD")
#write ouput
wd <- basename(getwd())
write.csv(output, file=paste(wd, "summary.csv", sep="."))
#write.table(utms, file=paste(x, ".mean", sep=""))
cat("done\n")

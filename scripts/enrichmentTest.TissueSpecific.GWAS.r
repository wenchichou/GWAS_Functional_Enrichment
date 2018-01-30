args <- commandArgs(trailingOnly=TRUE)

functionScript <-args[1]
resultDIR <- args[2]
outputTableFile <- args[3]
inputSNP <- readLines(args[4])
iteration <- as.numeric(as.character(args[5]))
cellTableFile <- args[6]
source(functionScript)

#resultDIR = "/home/unix/wcchou/gsapWenChi/tmp"
#iteration = 10000
#inputSNPList <-read.line("/home/unix/wcchou/gsapWenChi/tmp/SNPs_leanMass2ndPaper.txt")
#inputSNP <-readLines("/home/unix/wcchou/gsapWenChi/tmp/SNPs_leanMass2ndPaper.txt")
#source("/home/unix/wcchou/gsapWenChi/tmp/enrichmentTest.TissueSpecific.GWAS.function.r")
#allCellTable = read.table("./allCell.group.ID.description.txt", header=F, sep="\t")
allCellTable = read.table(cellTableFile, header=F, sep="\t")
colnames(allCellTable)=c("SN","cellGroup","cellID","cellDescription")

allChromHMMCells <- list()
for(i in 1:length(inputSNP)){
	#allChromHMMCells[i] <- read.table("./rs733381/ChromHMM.Cell", header=F, sep="\t")
	allChromHMMCells[i] <- read.table(paste(resultDIR,"/",inputSNP[i],"/ChromHMM.Cell",sep=""), header=F, sep="\t")
}
names(allChromHMMCells) <- inputSNP

pollAllChromCells <- data.frame(matrix(unlist(allChromHMMCells), ncol=1, byrow=T),stringsAsFactors=FALSE)

countCutoff.cell <- as.numeric(summary(as.numeric(table(pollAllChromCells))))[2] + 1


# 112817 separate skeletal MUS and smooth MUS
# replace "SM.MUS" to "sm.mus"
levels(allCellTable$cellID)=c(levels(allCellTable$cellID), sub("SM.MUS","sm.mus",grep("SM.MUS",allCellTable$cellID,value=T)))
allCellTable$cellID[grep("SM.MUS",allCellTable$cellID)]=sub("SM.MUS","sm.mus",grep("SM.MUS",allCellTable$cellID,value=T))


targetedCellGroup <- c("MUS","sm.mus","BRN","BLD","GI","FET","SKIN","FAT","LNG","BRST")

resultTable<-data.frame()
for(cell in 1:length(targetedCellGroup)){
	pollAllChromCells.Table <- prepareTable(pollAllChromCells, allCellTable, countCutoff.cell, targetedCellGroup[cell])
	all.permutedTestStatistics = replicate(iteration, onePermutation(pollAllChromCells.Table, allCellTable, countCutoff.cell, targetedCellGroup[cell]))
	#summary(all.permutedTestStatistics)
	observedPvalue=NULL
	finalPvalue=NULL
	countCutoff = 1
	for( SNP in 1:length(allChromHMMCells)){
		observedTestStatistics = prepareTableNgeneratePermutationStat(allChromHMMCells[SNP], allCellTable, countCutoff, targetedCellGroup[cell])
		observedPvalue=c(observedPvalue, observedTestStatistics)
		p = mean(c(all.permutedTestStatistics,observedTestStatistics) <= observedTestStatistics)
		finalPvalue = c(finalPvalue, p)
		cat(names(allChromHMMCells)[SNP],"in",targetedCellGroup[cell],"cells has enrichment p-value:",p,"\n")
		resultTable[cell,SNP]=p
	}
}
row.names(resultTable) <- targetedCellGroup
colnames(resultTable) <- names(allChromHMMCells)
write.table(resultTable, outputTableFile, sep="\t", quote=F, col.names=NA)

# Checking the permutation results
#plot(density(all.permutedTestStatistics))
#cabline(v=observedTestStatistics, col=2)
#hist(all.permutedTestStatistics)
#abline(v=observedTestStatistics, col=2)
#dev.off()

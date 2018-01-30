onePermutation = function(func_ObservedTable, allCellTable, cutoffHit, targetedCellGroup){
	targetedCellDescription=allCellTable$cellDescription[grep(targetedCellGroup, allCellTable$cellID)]
	permutedTable = func_ObservedTable[sample.int(length(func_ObservedTable), replace=FALSE)]
    names(permutedTable) = names(func_ObservedTable)
	"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    targetedCellTypeHit = sum(!is.na(match(targetedCellDescription, names(which(permutedTable > cutoffHit))))*1) 
	nontargetedCellTypeHit = length(names(which(permutedTable > cutoffHit)) %w/o% targetedCellDescription) 
    targetedCellTypeNoHit = sum(!is.na(match(targetedCellDescription, names(which(permutedTable <= cutoffHit))))*1) 
	nontargetedCellTypeNoHit = length(names(which(permutedTable <= cutoffHit)) %w/o% targetedCellDescription) 
	#cat(targetedCellTypeHit, targetedCellTypeNoHit, nontargetedCellTypeHit, nontargetedCellTypeNoHit,"\n")
    permutedTestStatistics=1-phyper(targetedCellTypeHit-1, targetedCellTypeHit+targetedCellTypeNoHit, nontargetedCellTypeHit+nontargetedCellTypeNoHit, targetedCellTypeHit+nontargetedCellTypeHit)
#	cat(permutedTestStatistics,"\n")
    return(permutedTestStatistics)
}



prepareTableNgeneratePermutationStat <-function(chromcells, allCellTable, cutoffHit, targetedCellGroup){
    observedTable = table(chromcells)

	#if there are cell lines not listed in input cell line table, I give 0 to those cell lines.
    cellNameNotFound = allCellTable$cellDescription[!(allCellTable$cellDescription %in% names(observedTable))]
	if(length(cellNameNotFound) >0){
    	additionalTable = rep(0,length(cellNameNotFound))
	    names(additionalTable) = cellNameNotFound
    	observedTable = c(observedTable, additionalTable)
	}
	targetedCellDescription=allCellTable$cellDescription[grep(targetedCellGroup, allCellTable$cellID)]
    #countTargetedCellType = sum(observedTable[grep("uscle",names(observedTable))])
	countTargetedCellType = sum(observedTable[match(targetedCellDescription,names(observedTable))])
    countNonTargetedCellType = sum(observedTable) - countTargetedCellType

	
	"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    targetedCellTypeHit = sum(!is.na(match(targetedCellDescription, names(which(observedTable > cutoffHit))))*1) 
	nontargetedCellTypeHit = length(names(which(observedTable > cutoffHit)) %w/o% targetedCellDescription) 
    targetedCellTypeNoHit = sum(!is.na(match(targetedCellDescription, names(which(observedTable <= cutoffHit))))*1) 
	nontargetedCellTypeNoHit = length(names(which(observedTable <= cutoffHit)) %w/o% targetedCellDescription) 
#    cat("#",length(observedTable),"Total cell lines or tissues in the observed table including cell lines with 0 count.\n")
#    cat("#",countTargetedCellType,"\tTotal Targeted cell lines or tissues with > 0 counts\n")
#    cat("#",countNonTargetedCellType,"\tTotal nonTargeted cell lines or tissues with > 0 counts\n")
#    cat("cutoff for a hit is:",cutoffHit,"\n")
#   cat("targetedCellTypeHit",targetedCellTypeHit,"\t")
#   cat("targetedCellTypeNoHit",targetedCellTypeNoHit,"\t")
#   cat("nontargetedCellTypeHit",nontargetedCellTypeHit,"\t")
#   cat("nontargetedCellTypeNoHit",nontargetedCellTypeNoHit,"\t")

    observedTestStatistics = 1-phyper(targetedCellTypeHit-1, targetedCellTypeHit+targetedCellTypeNoHit, nontargetedCellTypeHit+nontargetedCellTypeNoHit, targetedCellTypeHit+nontargetedCellTypeHit)
#   cat("observedTestStatistics:",observedTestStatistics,"\n")
    return(observedTestStatistics)
}

#chromcells=pollAllChromCells
#allCellTable=allCellTable
#cutoffHit=20
#targetedCellGroup="MUS"
prepareTable <-function(chromcells, allCellTable, cutoffHit, targetedCellGroup){
    observedTable = table(chromcells)

	#if there are cell lines not listed in input cell line table, I give 0 to those cell lines.
    cellNameNotFound = allCellTable$cellDescription[!(allCellTable$cellDescription %in% names(observedTable))]
	if(length(cellNameNotFound) >0){
    	additionalTable = rep(0,length(cellNameNotFound))
	    names(additionalTable) = cellNameNotFound
    	observedTable = c(observedTable, additionalTable)
	}
#   cat("#",length(observedTable),"Total cell lines or tissues in the observed table including cell lines with 0 count.\n")
	targetedCellDescription=allCellTable$cellDescription[grep(targetedCellGroup, allCellTable$cellID)]
    #countTargetedCellType = sum(observedTable[grep("uscle",names(observedTable))])
	countTargetedCellType = sum(observedTable[match(targetedCellDescription,names(observedTable))])
    countNonTargetedCellType = sum(observedTable) - countTargetedCellType

#   cat("#",countTargetedCellType,"\tTotal Targeted cell lines or tissues counts\n")
#   cat("#",countNonTargetedCellType,"\tTotal nonTargeted cell lines or tissues counts\n")
	
	"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    targetedCellTypeHit = sum(!is.na(match(targetedCellDescription, names(which(observedTable > cutoffHit))))*1) 
	nontargetedCellTypeHit = length(names(which(observedTable > cutoffHit)) %w/o% targetedCellDescription) 
    targetedCellTypeNoHit = sum(!is.na(match(targetedCellDescription, names(which(observedTable <= cutoffHit))))*1) 
	nontargetedCellTypeNoHit = length(names(which(observedTable <= cutoffHit)) %w/o% targetedCellDescription) 

#   cat("cutoff for a hit is:",cutoffHit,"\n")
#   cat("targetedCellTypeHit",targetedCellTypeHit,"\n")
#   cat("targetedCellTypeNoHit",targetedCellTypeNoHit,"\n")
#   cat("nontargetedCellTypeHit",nontargetedCellTypeHit,"\n")
#   cat("nontargetedCellTypeNoHit",nontargetedCellTypeNoHit,"\n")
	return(observedTable)
}


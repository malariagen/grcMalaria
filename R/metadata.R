##############################################################################
# Sample metadata management
################################################################################
meta.setDatasetMeta <- function (ctx, datasetName, meta, store=TRUE) {
    dataset <- ctx[[datasetName]]
    dataset$meta <- meta
    if (store) {
        metaFile <- meta.getMetaDataFile(ctx, datasetName)
        writeRdaSampleData(meta, metaFile)
    }
}

meta.getMetaDataFile <- function (ctx, datasetName) {
    metaFile <- getContextCacheFile(ctx, datasetName, "meta", "sampleMeta")
    metaFile
}

###############################################################################
# Metadata Column Statistics - General
################################################################################
#
# Returns a tabulated array of alleles sample counts in a column, sorted by decreasing sample count
# The names of the elements are the allele labels. Missing and <NA> are discarded.
#
meta.getColumnValueCounts <- function (sampleMeta, colName, excludeMultiValues=TRUE) {
    vals <- sampleMeta[,colName]			#; print(vals)
    vals <- vals[which(!(vals %in% c("-","<NA>")))]	#; print(vals)
    if (excludeMultiValues) {
        multi <- which(grepl(",", vals))		#; print(multi)
        if (length(multi) > 0) {
            vals <- vals[-multi]			#; print(vals)
        }
    }
    counts <- table(vals)				#; print(counts)
    counts <- sort(counts, decreasing=TRUE)
    counts
}
#
# Returns an allele sample count summary string for the columns specified (one line per column).
# The value sample counts are semicolon-separated.
#
meta.getValueCounts <- function (sampleMeta, colNames, params=NULL, excludeMultiValues=TRUE) {
    if (is.null(colNames) || (length(colNames) == 0)) {
        return (c());
    }
    result <- c()
    for (mIdx in 1:length(colNames)) {
        colName <- colNames[mIdx]			#; print(colName)
        counts <- meta.getColumnValueCounts (sampleMeta, colName, excludeMultiValues)
        nVals <- names(counts)
        valPairs <- paste(nVals, counts, sep=":")
        valPairStr <- paste(valPairs, collapse="; ")
        result <- c(result, valPairStr)			#; print(valPairStr)
    }
    result <- as.character(result)
    names(result) <- colNames
    result
}
#
###############################################################################
# Metadata Statistics - Resistance Predictions
################################################################################
#
#
#
meta.getResistancePrevalence <- function (ctx, sampleMeta, drugNames, params=NULL) {
    if (is.null(drugNames) || (length(drugNames) == 0)) {
        return (c());
    }
    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- param.getParam ("map.aggregateCountMin", params)
    }  							#; print(aggregateCountMin)
    result <- c()
    for (mIdx in 1:length(drugNames)) {
        drugName <- drugNames[mIdx]			#; print(drugName)
        phenos <- sampleMeta[,drugName]			#; print(phenos)
        rCount <- length(which(phenos=="Resistant"))	#; print(rCount)
        sCount <- length(which(phenos=="Sensitive"))	#; print(sCount)
        totalCount <- rCount + sCount
        if (totalCount >= aggregateCountMin) {
            preval <- rCount / totalCount
        } else {
            preval <- NA
        }
        result <- c(result, preval)			#; print(preval)
    }
    result <- as.numeric(result)			#; print(result)
    names(result) <- drugNames				#; print(result)
    result
}
#
###############################################################################
# Metadata Statistics - Amino alleles
################################################################################
#
#
#
meta.getPositionOfMutation <- function (mutation) {
    positionName <- substring(mutation, 1, nchar(mutation)-1)
    positionName
}
#
meta.getAlleleOfMutation <- function (mutation) {
    allele <- substring(mutation, nchar(mutation), nchar(mutation))
    allele
}
#
meta.getPositionData <- function (config, positionName) {
    posData <- list()
    mParts <- unlist(strsplit(positionName, "_"))
    posData$gene <- mParts[1]
    posData$ref  <- substring(mParts[2], 1, 1)
    posData$pos  <- as.integer(substring(mParts[2], 2))
    #
    # FUTURE
    #columnName <- paste("M", positionName, sep="_")
    #
    # RIGHT NOW, BODGE: "crt_C72" must give column name "PfCRT:72"
    #
    posData$columnName <- paste0(config$species, toupper(posData$gene), ":", posData$pos)
    posData
}
#
###############################################################################
# Metadata Statistics - Prevalence
################################################################################
#
meta.getMutationPrevalence <- function (ctx, sampleMeta, mutationNames, params) {
    positionNames <- c()
    alleles <- c()
    for (mIdx in 1:length(mutationNames)) {
        mut <- mutationNames[mIdx]
        pos <- meta.getPositionOfMutation (mut)
        all <- meta.getAlleleOfMutation (mut)
        positionNames <- c(positionNames, pos)
        alleles <- c(alleles, all)
    }											#; print(positionNames); print(alleles)
    prevData <- meta.getAllelePrevalence (ctx, sampleMeta, positionNames, alleles)	#; print(prevData)
    prevData
}
#
meta.getAllelePrevalence <- function (ctx, sampleMeta, positionNames, alleles, params=NULL) {
    if (is.null(positionNames) || (length(positionNames) == 0)) {
        return (c());
    }
    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- param.getParam ("map.aggregateCountMin", params)
    }  							#; print(aggregateCountMin)
    result <- c()
    for (mIdx in 1:length(positionNames)) {
        positionName <- positionNames[mIdx]		#; print(positionName)
        allele <- alleles[mIdx]				#; print(allele)

        posData <- meta.getPositionData (ctx$config, positionName)
        columnName <- posData$columnName	 	#; print(columnName)
        
        # Remove hets and missing
        genos <- sampleMeta[,columnName]		#; print(genos)
        genos <- genos[which(!(genos %in% c("*","-")))]	#; print(genos)
        genos <- genos[which(nchar(genos)==1)]	    	#; print(genos)

        totalCount <- length(genos)
        preval <- NA
        if (totalCount >= aggregateCountMin) {
            alleleCount <- length(which(genos==allele))	#; print(alleleCount)
            preval <- alleleCount / totalCount
        }
        result <- c(result, preval)			#; print(preval)
    }
    result <- as.numeric(result)
    names(result) <- paste(positionNames, alleles, sep="")
    result
}
#
###############################################################################
# Metadata processing
################################################################################
#
# Select sample metadata as specified in the config
#
meta.select <- function(sampleMeta, selectDefs) {
    for (sIdx in 1:length(selectDefs)) {
        selDef <- selectDefs[[sIdx]]
        selField <- selDef$field
        selValues <- selDef$values
        #print (paste("Select field: ",selField))
        #print (paste("Select values: ",selValues))
        metaValues <- sampleMeta[,selField]
        sampleMeta <- sampleMeta[which (metaValues %in% selValues),]
        #print (paste("Selected samples: ",nrow(sampleMeta)))
    }
    sampleMeta
}
#
# Select sample metadata within a given time period (start and end date inclusive)
#
meta.defaultDateFormat <- "%d/%m/%Y"
meta.naDate <- "01/01/1000"		# Date value to use instead of NA or "-"
#
meta.filterByDate <- function(sampleMeta, startDate, endDate, format=meta.defaultDateFormat) {
    metaDates <- meta.getSampleDates (sampleMeta, format)		#; print(head(metaDates))
    selectIdx <- !is.na(metaDates)
    if (!is.null(startDate)) {
        selectIdx <- selectIdx & (metaDates>=startDate)
    }								#; print(selectIdx[1:10])
    if (!is.null(endDate)) {
        selectIdx <- selectIdx & (metaDates<=endDate)
    }								#; print(selectIdx[1:10])
    filterMeta <- NULL						#; print(length(which(selectIdx==TRUE)))
    if (length(which(selectIdx==TRUE)) > 0) {
         filterMeta <- sampleMeta[selectIdx,]
    }								#; print(nrow(filterMeta))
    filterMeta
}
#
#
#
meta.getSampleDates <- function (sampleMeta, format=meta.defaultDateFormat) {
    .naDate <- 
    sampleCount <- nrow(sampleMeta)
    if ("CollectionDate" %in% colnames(sampleMeta)) {
        metaDatesIn <- as.character(sampleMeta$CollectionDate)
    } else {
        #message("Your GRC data does not contain the CollectionDate Field. Using the Year field instead- see the documentation.")
        metaDatesIn <- rep("-", sampleCount)
    }									#; print(metaDatesIn[1:500])
    #
    # If date is missing, use the year which shoul be there
    #
    missingDateIndexes <- which (metaDatesIn == "-")			#; print(missingDateIndexes)
    if (length(missingDateIndexes) > 0) {
        years <- suppressWarnings(as.integer(sampleMeta$Year))		#; print(years[1:500])
        yearDates <- paste0("31/12/", years)				#; print(yearDates[1:500])
        for (mdIdx in missingDateIndexes) {
            if (is.na(years[mdIdx])) {
                #
                # Silly workaround. Placing a "-" can cause as.Date() errors if it is in the first element of the vector!!!
                # So we put a valid but improbably date instead. Not ideal, but avoids errors
                #
                #metaDatesIn[mdIdx] <- "-" 
                metaDatesIn[mdIdx] <- meta.naDate 
            } else {
                metaDatesIn[mdIdx] <- yearDates[mdIdx]
            }
        }								#; print(metaDatesIn[1:500])
    }
    metaDates <- as.Date(metaDatesIn, tryFormats=format, optional=TRUE)		#; print(metaDates[1:500])
    metaDates[which(metaDates==meta.naDate)] <- NA
    metaDates
}

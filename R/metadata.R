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
meta.getColumnValueCounts <- function (sampleMeta, colName, excludeMultiValues=TRUE, excludeHets=TRUE) {
    vals <- sampleMeta[,colName]			#; print(vals)
    vals <- vals[which(!(vals == "<NA>"))]		#; print(vals)
    vals <- vals[which(!grepl("-", vals, fixed=TRUE))]	#; print(vals)
    if (excludeMultiValues) {
        multi <- which(grepl(",", vals))		#; print(multi)
        if (length(multi) > 0) {
            vals <- vals[-multi]			#; print(vals)
        }
    }
    if (excludeHets) {
        het <- which(grepl("\\*", vals))			#; print(het)
        if (length(het) > 0) {
            vals <- vals[-het]				#; print(vals)
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
        drugCol <- setup.getFeatureColumn (ctx$config, drugName)	#; print(drugCol)
        phenos <- sampleMeta[,drugCol]			#; print(phenos)
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
# Metadata Statistics - V2
################################################################################
#
# Mutation Prevalence
#
meta.getMutationPrevalence <- function (ctx, sampleMeta, mutationNames, params) {	#; print(mutationNames)

    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- param.getParam ("map.aggregateCountMin", params)
    }  							#; print(aggregateCountMin)
    prevalenceData <- c()
    for (mIdx in 1:length(mutationNames)) {
        mutName <- mutationNames[mIdx]				#; print(mutName)
        mutData <- meta.getPositionData (ctx$config, mutName)
        if (is.null(mutData)) {
            stop(paste("Invalid mutation specified:",mutName))
        }
        columnName <- mutData$columnName			#; print(columnName)
        allele  <- mutData$alt
        #
        # Get genotypes and Remove hets and missing
        genos <- sampleMeta[,columnName]			#; print(genos)
        genos <- genos[which(!(genos %in% c("*","-","X","N")))]	#; print(genos)
        genos <- genos[which(nchar(genos)==1)]	    		#; print(genos)
        #
        # Calculate the fraction of the desirted allele
        #
        totalCount <- length(genos)
        prevalence <- NA
        if (totalCount >= aggregateCountMin) {
            alleleCount <- length(which(genos==allele))		#; print(alleleCount)
            prevalence <- alleleCount / totalCount
        }
        prevalenceData <- c(prevalenceData, prevalence)		#; print(prevalence)
    }
    prevalenceData <- as.numeric(prevalenceData)
    names(prevalenceData) <- paste(mutationNames)
    prevalenceData
}
#
# Statistics - Amino alleles
#
meta.getPositionData <- function (config, posName) {		#; print(posName)
    feat <- config$features					#; print(rownames(feat))
    if (!(posName %in% rownames(feat))) {
        return (NULL)
    }
    posData <- list()
    posData <- list(
        featureName=feat[posName,"FeatureName"],
        columnName=feat[posName,"ColumnName"],
        ref=feat[posName,"Reference"],
        alt=feat[posName,"Alternative"]
    )
    posData
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

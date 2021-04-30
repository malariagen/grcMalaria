#
meta.colName.sample  <- "SampleId"
meta.colName.barcode <- "GenBarcode"
#
#
###############################################################################
# Declarations: sample data and metadata
################################################################################
folder.data.meta <- getDataFolder ("meta")

sampleMetaFname <- "sampleMeta.tab"

meta.initializeMetadata <- function (ctx, store=TRUE) {
    if (exists("sampleMeta.xls.file") && !is.null(sampleMeta.xls.file)) {
        sampleMeta <- meta.readExcelSampleMeta ()
    } else if (exists("sampleMeta.tab.file") && !is.null(sampleMeta.tab.file)) {
        sampleMeta <- meta.readTabSampleMeta ()
    } else {
        stop("No metadata file specified in configuration")
    }
    ctx <- meta.setContextSampleMeta(ctx, sampleMeta, store=store)
    ctx
}

meta.setContextSampleMeta <- function (ctx, newMeta, store=TRUE) {
    ctx$meta <- newMeta
    if (store) {
        writeSampleData(newMeta, 
                        analysis.getDataFile(ctx, folder.data.meta, sampleMetaFname))
    }
    ctx
}

###############################################################################
# Data input routines
################################################################################
#
# Read sample metadata file
#
meta.readExcelSampleMeta <- function () {
    sampleMetaFile <- paste(sampleMeta.folder, sampleMeta.xls.file, sep="/")
    sampleMeta <- data.frame(read_excel(sampleMetaFile, sheet=sampleMeta.xls.sheet, col_names=FALSE, col_types="text"))
    colnames(sampleMeta) <- sampleMeta[1,]			#; print(str(sampleMeta))
    sampleMeta <- sampleMeta[-1,]
    sampleNames <- sampleMeta[,meta.colName.sample]	#; print (sampleNames)
    rownames(sampleMeta) <- sampleNames				#; print(str(sampleMeta))
    return (sampleMeta)
}

meta.readTabSampleMeta <- function () {
    sampleMetaFile <- paste(sampleMeta.folder, sampleMeta.tab.file, sep="/")
    sampleMeta <- read.table(sampleMetaFile, as.is=TRUE, stringsAsFactors=TRUE, quote="", header=FALSE, sep="\t")
    colnames(sampleMeta) <- sampleMeta[1,]			#; print(str(sampleMeta))
    sampleMeta <- sampleMeta[-1,]
    sampleNames <- sampleMeta[,meta.colName.sample]	#; print (sampleNames)
    rownames(sampleMeta) <- sampleNames				#; print(str(sampleMeta))
    return (sampleMeta)
}

###############################################################################
# Metadata Statistics
################################################################################
meta.getResistancePrevalence <- function (sampleMeta, drugNames, params=NULL) {
    if (is.null(drugNames) || (length(drugNames) == 0)) {
        return (c());
    }
    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- analysis.getParam ("map.aggregateCountMin", params, default.map.aggregateCountMin)
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
    names(result) <- drugNames
    result
}
#
meta.getAllelePrevalence <- function (sampleMeta, positionNames, alleles, params=NULL) {
    if (is.null(positionNames) || (length(positionNames) == 0)) {
        return (c());
    }
    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- analysis.getParam ("map.aggregateCountMin", params, default.map.aggregateCountMin)
    }  							#; print(aggregateCountMin)
    result <- c()
    for (mIdx in 1:length(positionNames)) {
        positionName <- positionNames[mIdx]		; print(positionName)
        allele <- alleles[mIdx]				; print(allele)

        # FUTURE
        #columnName <- paste("M", positionName, sep="_")
        # RIGHT NOW, BODGE: "crt_C72" must give column name "PfCRT:72"
        mParts <- unlist(strsplit(positionName, "_"))
        gene <- toupper(mParts[1])
        ref <- substring(mParts[2], 1, 1)
        pos <- as.integer(substring(mParts[2], 2))
        columnPrefix <- config.species
        columnName <- paste(columnPrefix, gene, ":", pos, sep="")	; print(columnName)

        # Remove hets and missing
        genos <- sampleMeta[,columnName]		#; print(genos)
        genos <- genos[which(!(genos %in% c("*","-")))]	#; print(genos)
        genos <- genos[which(nchar(genos)==1)]		#; print(genos)
        
        totalCount <- length(genos)
        preval <- NA
        if (totalCount >= aggregateCountMin) {
            alleleCount <- length(which(genos==allele))	#; print(alleleCount)
            preval <- alleleCount / totalCount
        }
        result <- c(result, preval)			#; print(preval)
    }
    names(result) <- paste(positionNames, alleles, sep="")
    result
}
#
meta.getValueCounts <- function (sampleMeta, colNames, params=NULL) {
    if (is.null(colNames) || (length(colNames) == 0)) {
        return (c());
    }
    aggregateCountMin <- 0
    if (!is.null(params)) {
        aggregateCountMin <- analysis.getParam ("map.aggregateCountMin", params, default.map.aggregateCountMin)
    }  							#; print(aggregateCountMin)
    result <- list()
    for (mIdx in 1:length(colNames)) {
        colName <- colNames[mIdx]			#; print(colName)
        vals <- sampleMeta[,colName]			#; print(vals)
        vals <- table(vals[which(vals!="-")])		#; print(vals)
        vals <- sort(vals, decreasing=TRUE)
        nVals <- names(vals)
        valPairs <- paste(nVals, vals, sep=":")
        valPairStr <- paste(valPairs, collapse = "\n")
        result <- c(result, valPairStr)			#; print(valPairStr)
    }
    names(result) <- colNames
    result
}
#
###############################################################################
# Metadata processing
################################################################################
#
# Add fields generated from the merger of multiple fields, as specified in the config
#
meta.addMergedFields <- function(sampleMeta, mergeFieldsDefs) {
    for (mIdx in 1:length(mergeFieldsDefs)) {
        fieldList <- mergeFieldsDefs[[mIdx]]
    
        outFieldName <- fieldList[[1]]
        outValues <- sampleMeta[,outFieldName]
    
        for (fIdx in 2:length(fieldList)) {
            srcFieldName <- fieldList[[fIdx]]
            outFieldName <- paste(outFieldName, srcFieldName, sep="__")
            srcValues <- sampleMeta[,srcFieldName]
            outValues <- paste(outValues, srcValues, sep="__")
        }
        #print (paste("Merged field: ",outFieldName))
        #print (outValues[1:10])
        headers <- colnames(sampleMeta)
        sampleMeta <- cbind(sampleMeta, outValues)
        colnames(sampleMeta) <- c(headers, outFieldName)
    }
    sampleMeta
}
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

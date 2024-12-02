# #####################################################################################
# GRC-specific setup configurations
# #####################################################################################
setup.getSpeciesConfig <- function (species, version, grcFeatures) {
    chrLengths <- species.chrLengths[[species]]
    chromosomes <- names(chrLengths)
    #
    config <- list (
        species = species,
        version = version,
        #
        chromosomes = chromosomes,
        chrCount = length(chromosomes),
        chrLengths = chrLengths
    )
    #
    #
    #
    config$drugs <- grcFeatures[which(grcFeatures$Class == "Drug Resistance Prediction"),]
    #
    # Drug Resistant Positions
    #
    posFeat <- grcFeatures[which(grcFeatures$Class == "Amino Variant"),]	#; print(posFeat)
    config$drugResistancePositions <- posFeat
    #
    # We create the Drug Resistant Mutations features based on the alternative alleles of the Drug Resistant Positions
    #
    fNames <- c(); fCols <- c(); fRefs <- c(); fAlts <- c()
    for (i in 1:nrow(posFeat)) {						#;print(rownames(posFeat)[i])
        alt <- posFeat[i,"Alternative"]
        if ((alt == "-") | (alt == "")) next
        alt <- trimws(unlist(strsplit(alt, ",", fixed=TRUE)))
        for (j in 1:length(alt)) {
            fNames <- c(fNames, paste0(posFeat[i,"FeatureName"],alt[j]))
            fCols <- c(fCols, posFeat[i,"ColumnName"])
            fRefs <- c(fRefs, posFeat[i,"Reference"])
            fAlts <- c(fAlts, alt[j])
        }
    }
    mutFeat <- data.frame(FeatureName=fNames, ColumnName=fCols, Class=rep("Amino Mutation",length(fNames)), 
                          DataType=rep("AminoAcidMutation",length(fNames)), Reference=fRefs, Alternative=fAlts)
    rownames(mutFeat) <- fNames
    config$drugResistanceMutations <- mutFeat					#; print(mutFeat)
    #
    #
    #
    config$countColumns <- grcFeatures[which(grcFeatures$DataType %in% c("MutationList","AminoAcidSequence")),]
    #
    #
    #
    config$amplificationColumns <- grcFeatures[which(grcFeatures$DataType == "Amplification"),]
    #
    #
    #
    config$barcodeMeta <- grcFeatures[which(grcFeatures$Class == "Genetic Barcode SNP"),]
    #
    #
    #
    config$cluster.stats.drugs        <- config$drugs
    config$cluster.stats.mutations    <- config$drugResistanceMutations
    config$cluster.stats.alleleCounts <- config$countColumns
    #
    #
    # Attach the newly created features to the GRC features before sticking them in the config
    #
    config$features <- rbind (grcFeatures, mutFeat)
    config
}

# #####################################################################################
# Species-specific stuff
# #####################################################################################
species.chrLengths <- list(
    Pf = list(
            Pf3D7_01_v3=640851,
            Pf3D7_02_v3=947102,
            Pf3D7_03_v3=1067971,
            Pf3D7_04_v3=1200490,
            Pf3D7_05_v3=1343557,
            Pf3D7_06_v3=1418242,
            Pf3D7_07_v3=1445207,
            Pf3D7_08_v3=1472805,
            Pf3D7_09_v3=1541735,
            Pf3D7_10_v3=1687656,
            Pf3D7_11_v3=2038340,
            Pf3D7_12_v3=2271494,
            Pf3D7_13_v3=2925236,
            Pf3D7_14_v3=3291936
    ),
    Pv = list(
            PvP01_00_v1=4802351,
            PvP01_01_v1=1021664,
            PvP01_02_v1=956327,
            PvP01_03_v1=896704,
            PvP01_04_v1=1012024,
            PvP01_05_v1=1524814,
            PvP01_06_v1=1042791,
            PvP01_07_v1=1652210,
            PvP01_08_v1=1761288,
            PvP01_09_v1=2237066,
            PvP01_10_v1=1548844,
            PvP01_11_v1=2131221,
            PvP01_12_v1=3182763,
            PvP01_13_v1=2093556,
            PvP01_14_v1=3153402
    )
)

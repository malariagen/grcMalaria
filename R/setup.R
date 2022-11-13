# #####################################################################################
# Build configuration
# #####################################################################################
setup.getConfig <- function (grc, dir, minSnpTypability, minSampleTypability) {
    data <- grc$data
    speciesConfig <- grc$speciesConfig

    print(paste("Root folder:", dir))
    config <- list (
        version=grc$version,
        folder.root=dir,
        folder.data=getSubFolder (dir, "data"),
        folder.out=getSubFolder (dir, "out"),
        minSnpTypability=minSnpTypability,
        minSampleTypability=minSampleTypability,
        defaultPalette=setup.getDefaultPalette(),
        defaultTextPalette=graphics.makeTextPalette(setup.getDefaultPalette())
    )

    # Merge the config with the species config
    config <- c(config, speciesConfig)
    config
}

# #####################################################################################
# Species-specific setup configurations
# #####################################################################################
setup.getSpeciesConfig <- function (species, version) {
    if (species == "Pf") {
        spConfig <- setup.getPfConfig(version)
    } else if (species == "Pf") {
        spConfig <- setup.getPvConfig(version)
    } else {
        stop (paste("Invalid species specified: ", species))
    }

    versionLabel <- paste("v",version,sep="")
    barcodeMeta <- barcode.metadata[[species]][[versionLabel]]

    config <- list (
        species=species,
        version=version,
        barcodeMeta=barcodeMeta
    )
    config <- c (config, spConfig)
    config
}

# #####################################################################################
# Plasmodium falciparum
#
setup.getPfConfig <- function (version) {
    config <- setup.buildSpeciesConfig (
        species="Pf", version=version,
        chrLengths = list(
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
        drugs = c(
            "Artemisinin",
            "Piperaquine",
            "DHA-PPQ",
            "Mefloquine",
            "AS-MQ",
            "Chloroquine",
            "Pyrimethamine",
            "Sulfadoxine",
            "S-P",
            "S-P-IPTp"
        ),
        drugResistancePositions = c(
            "crt_C72", "crt_M74",  "crt_N75", "crt_K76",
            "crt_T93",  "crt_H97", "crt_I218", "crt_A220", "crt_Q271",
            "crt_N326", "crt_T333", "crt_I356", "crt_R371",
            "dhfr_N51", "dhfr_C59", "dhfr_S108", "dhfr_I164",
            "dhps_S436", "dhps_A437", "dhps_K540", "dhps_A581", "dhps_A613",
            "mdr1_N86", "mdr1_Y184", "mdr1_S1034", "mdr1_F1226",
            "mdr1_D1246", "arps10_V127", "arps10_D128", "fd_D193", "mdr2_T484", "exo_E415"
        ),
        drugResistanceMutations = c(
            "crt_C72S", "crt_M74I",  "crt_N75E",  "crt_N75D",  "crt_K76T",
            "crt_T93S",  "crt_H97Y", "crt_H97L", "crt_I218F", "crt_A220S", "crt_Q271E",
            "crt_N326S", "crt_N326D", "crt_T333S", "crt_I356T", "crt_I356L", "crt_R371I",
            "dhfr_N51I", "dhfr_C59R", "dhfr_C59H", "dhfr_S108N", "dhfr_S108T", "dhfr_I164L",
            "dhps_S436A", "dhps_S436F", "dhps_A437G", "dhps_K540E", "dhps_K540N", "dhps_A581G", "dhps_A613T", "dhps_A613S",
            "mdr1_N86Y", "mdr1_Y184F", "mdr1_S1034I", "mdr1_F1226Y",
            "mdr1_D1246Y", "arps10_V127M", "arps10_D128H", "arps10_D128Y", "fd_D193Y", "mdr2_T484I", "exo_E415G"
        ),
        countColumns = c(
            "Pfkelch13"
        ),
        amplificationColumns = c(
            "pm23-Amp", "mdr1-Amp"
        )
    )
    config
}

# #####################################################################################
# Plasmodium vivax
#
setup.getPvConfig <- function (version) {
    config <- setup.buildSpeciesConfig (
        species="Pv", version=version,
        chrLengths = list(
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
        ),
        drugs = c(
        ),
        drugResistancePositions = c(
            "dhfr_F57", "dhfr_R58", "dhfr_T61", "dhfr_N117",
            "dhps_E380", "dhps_S382", "dhps_G383", "dhps_Y385", "dhps_A553",
            "mdr1_Y976"
        ),
        drugResistanceMutations = c(
            "dhfr_F57L", "dhfr_F57I", "dhfr_R58S", "dhfr_R58K", "dhfr_T61M", "dhfr_N117S", "dhfr_N117T",
            "dhps_S382C", "dhps_S382A", "dhps_G383A", "dhps_A553G",
            "mdr1_Y976F"
        ),
        countColumns = c(
        ),
        amplificationColumns = c(
        )
    )
    config
}

# #####################################################################################
# Species Config Buildup
# #####################################################################################
setup.buildSpeciesConfig <- function (species, version, chrLengths, drugs,
                     drugResistancePositions, drugResistanceMutations, countColumns, amplificationColumns) {
    chromosomes = names(chrLengths)
    config <- list (
        species = species,
        version = version,
        #
        # Chromosomes in order for this species
        #
        chromosomes = chromosomes,
        chrCount = length(chromosomes),
        chrLengths = chrLengths,
        #
        # Things we can calculate the prevalence of for this species
        #
        drugs = drugs,
        drugResistancePositions = drugResistancePositions,
        drugResistanceMutations = drugResistanceMutations,
        countColumns = countColumns,
        amplificationColumns = amplificationColumns,
        cluster.stats.drugs = drugs,
        cluster.stats.mutations = drugResistanceMutations,
        cluster.stats.alleleCounts = countColumns
    )
    config
}

# #####################################################################################
# Palette for pie charts.
# #####################################################################################
#
setup.getDefaultPalette <- function () {
    # This one was submitted by Kevin Wright on https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
    # Other similar palettes are shown in the same article
    c25Palette <- c(
        "dodgerblue2",
        "#E31A1C",	# red
        "green4",
        "#6A3D9A",	# purple
        "#FF7F00",	# orange
        "black",
        "gold1",
        "skyblue2",
        "#FB9A99",	# lt pink
        "palegreen2",
        "#CAB2D6",	# lt purple
        "#FDBF6F",	# lt orange
        "gray70",
        "khaki2",
        "maroon",
        "orchid1",
        "deeppink1",
        "blue1",
        "steelblue4",
        "darkturquoise",
        "green1",
        "yellow4",
        "yellow3",
        "darkorange4",
        "brown"
    )
    c25Palette
}

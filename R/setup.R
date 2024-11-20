#
# Override warning message hiding (i.e. those cause by packages that are invoked)
# This is set to TRUE in the version we distribute, but we can set it to FALSE to make sure we see the warnings.
#
OVERRIDE_EXTERNAL_WARNINGS <- TRUE
options(warn=1)
#
# #####################################################################################
# Build configuration
# #####################################################################################
setup.getConfig <- function (grc, dir, outDir, minSnpTypability, minSampleTypability, maxImputedProportion) {
    data <- grc$data
    speciesConfig <- grc$speciesConfig

    print(paste("Root folder:", dir))
    config <- list (
        version=grc$version,
        folder.root=dir,
        folder.data=getSubFolder (dir, "data"),
        folder.out=getSubFolder (dir, outDir),
        minSnpTypability=minSnpTypability,
        minSampleTypability=minSampleTypability,
        maxImputedProportion=maxImputedProportion,
        defaultPalette=setup.getDefaultPalette(),
        defaultTextPalette=graphics.makeTextPalette(setup.getDefaultPalette())
    )

    # Merge the config with the species config
    config <- c(config, speciesConfig)
    config
}

# #####################################################################################
# 
# #####################################################################################
setup.getFeatureNames <- function (features) {
    names <- rownames(features)
    names
}

setup.getFeatureColumn <- function (config, featureName) {
    cNames <- config$features[featureName,"ColumnName"]
    cNames
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

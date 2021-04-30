config.species <- "Pv"

# #####################################################################################
# Chromosomes in order for this species
# #####################################################################################
chrLengths <- list(
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
chromosomes <- names(chrLengths)
chrCount <- length(chromosomes)

# #####################################################################################
# Things we can calculate the prevalence of for this species
# #####################################################################################

allDrugs <- c()

allDrugResistancePositions <- c(
"dhfr_F57", "dhfr_R58", "dhfr_T61", "dhfr_N117", "dhps_E380", "dhps_S382", "dhps_G383", "dhps_Y385", "dhps_A553", "mdr1_Y976")

allDrugResistanceMutations <- c(
    "dhfr_F57L", "dhfr_F57I", "dhfr_R58S", "dhfr_R58K", "dhfr_T61M", "dhfr_N117S", "dhfr_N117T", 
    "dhps_S382C", "dhps_S382A", "dhps_G383A", "dhps_A553G", "mdr1_Y976F")

cluster.prevalenceColumns <- c()
cluster.countColumns <- c()

# #####################################################################################
# The rest of the setup should be the same for all species.
# #####################################################################################
source(paste(folder.code, "setup.common.R", sep="/"), local=FALSE)

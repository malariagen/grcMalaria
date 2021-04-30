config.species <- "Pf"

# #####################################################################################
# Chromosomes in order for this species
# #####################################################################################
chrLengths <- list(
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
)
chromosomes <- names(chrLengths)
chrCount <- length(chromosomes)

# #####################################################################################
# Things we can calculate the prevalence of for this species
# #####################################################################################
allDrugs <- c("Artemisinin", "Piperaquine", "DHA-PPQ", "Chloroquine", "Pyrimethamine", "Sulfadoxine", "S-P", "S-P-IPTp")

allDrugResistancePositions <- c(
    "crt_C72", "crt_M74",  "crt_N75", "crt_K76",  "crt_T93",  "crt_H97", "crt_I218", "crt_A220", "crt_Q271", 
                "crt_N326", "crt_T333", "crt_I356", "crt_R371", 
    "dhfr_N51", "dhfr_C59", "dhfr_S108", "dhfr_I164", 
    "dhps_S436", "dhps_A437", "dhps_K540", "dhps_A581", "dhps_A613",
    "mdr1_N86", "mdr1_Y184", "mdr1_S1034", "mdr1_F1226", 
    "mdr1_D1246", "arps10_V127", "arps10_D128", "fd_D193", "mdr2_T484", "exo_E415")

allDrugResistanceMutations <- c(
    "crt_C72S", "crt_M74I",  "crt_N75E",  "crt_N75D",  "crt_K76T",  "crt_T93S",  "crt_H97Y", "crt_H97L", "crt_I218F", "crt_A220S", "crt_Q271E", 
                "crt_N326S", "crt_N326D", "crt_T333S", "crt_I356T", "crt_I356L", "crt_R371I", 
    "dhfr_N51I", "dhfr_C59R", "dhfr_C59H", "dhfr_S108N", "dhfr_S108T", "dhfr_I164L", 
    "dhps_S436A", "dhps_S436F", "dhps_A437G", "dhps_K540E", "dhps_K540N", "dhps_A581G", "dhps_A613T", "dhps_A613S",
    "mdr1_N86Y", "mdr1_Y184F", "mdr1_S1034I", "mdr1_F1226Y", 
    "mdr1_D1246Y", "arps10_V127M", "arps10_D128H", "arps10_D128Y", "fd_D193Y", "mdr2_T484I", "exo_E415G")

cluster.prevalenceColumns <- allDrugs
cluster.countColumns <- c("K13")

# #####################################################################################
# The rest of the setup should be the same for all species.
# #####################################################################################
source(paste(folder.code, "setup.common.R", sep="/"), local=FALSE)

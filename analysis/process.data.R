## Process EPT Data

# Bulk on Hyak -----------------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
#source("analysis/fx.R")

# ( fn <- list.files("data/followup", full.names = TRUE) )
# # fn <- list.files(pattern = "data/followup/n[3-4][0-9][0-9][0-9].rda")
# for (i in fn) {
#     load(i)
#     sim <- truncate_sim(sim, at = 2600)
#     save(sim, file = i, compress = TRUE)
#     cat("*")
# }

# truncate and limit sim 3000/4000 files on Hyak
# cd /gscratch/csde/kweiss2/sti/data
# module load r_3.2.4
# R
fn <- list.files(pattern = "n[8-9][0-9][0-9][0-9].rda")
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 5200)
  vars.needed <- c(# HIV
    "incid", "ir100", "hivtests.nprep", "i.prev",

    # GC
    "incid.gc", "ir100.gc", "ir100.rgc", "ir100.ugc",
    "prev.gc", "prev.uct", "prev.rct",
    "GCasympttests", "GCsympttests",
    "txGC","txGC_asympt",

    # CT
    "incid.ct", "ir100.ct", "ir100.rct", "ir100.uct",
    "prev.ct", "prev.uct", "prev.rct",
    "CTasympttests", "GCsympttests",
    "txCT","txCT_asympt",

    # Combined
    "incid.gcct", "ir100.gcct",
    "prev.STI", "prev.sti.tttraj1", "prev.sti.tttraj2",
    "stiasympttests", "stisympttests",
    "stiasympttests.tttraj1", "stiasympttests.tttraj2",
    "stisympttests.tttraj1", "stisympttests.tttraj2",
    "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
    "txSTI_asympt",
    "tt.traj.sti1", "tt.traj.sti2",

    # Other
    "num",
    'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
    'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',

    # EPT
    "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
    "eptTx", "propindexeptElig",
    "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
    "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
    "eptgcctinfectundiaghiv",
    "gc.timesInf", "ct.timesInf", "sti.timesInf",
    "eptgcinfecthiv", "eptctinfecthiv",
    "eptgcctinfecthiv",
    "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
    "eptgcctinfecthiv_inst",
    "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
    "eptgcctinfectundiaghiv_inst",
    "eptindexprovided_gc", "eptindexprovided_ct",
    "eptpartprovided_gc", "eptpartprovided_ct",
    "eptpartprovided_main", "eptpartprovided_pers",
    "eptpartprovided_inst", "eptpartuptake_main",
    "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
    "eptpartuptake_pers", "eptpartuptake_inst",
    "eptpartuptake_gc", "eptpartuptake_ct")

  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  out.fn <- paste0("followup/", i)
  save(sim, file = out.fn, compress = "gzip")
  file.remove(i)
  cat(i, "\n")
}


### 1 by 1 processing on Hyak - cd /gscratch/csde/kweiss2/sti/data -------------
rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
load("sim.n3037.rda")
sim <- truncate_sim(sim, at = 5200)
vars.needed <- c(# HIV
  "incid", "ir100", "hivtests.nprep", "i.prev",

  # GC
  "incid.gc", "ir100.gc", "ir100.rgc", "ir100.ugc",
  "prev.gc", "prev.uct", "prev.rct",
  "GCasympttests", "GCsympttests",
  "txGC","txGC_asympt",

  # CT
  "incid.ct", "ir100.ct", "ir100.rct", "ir100.uct",
  "prev.ct", "prev.uct", "prev.rct",
  "CTasympttests", "GCsympttests",
  "txCT","txCT_asympt",

  # Combined
  "incid.gcct", "ir100.gcct",
  "prev.STI", "prev.sti.tttraj1", "prev.sti.tttraj2",
  "stiasympttests", "stisympttests",
  "stiasympttests.tttraj1", "stiasympttests.tttraj2",
  "stisympttests.tttraj1", "stisympttests.tttraj2",
  "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
  "txSTI_asympt",
  "tt.traj.sti1", "tt.traj.sti2",

  # Other
  "num",
  'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
  'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',

  # EPT
  "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
  "eptTx", "propindexeptElig",
  "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
  "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
  "eptgcctinfectundiaghiv",
  "gc.timesInf", "ct.timesInf", "sti.timesInf",
  "eptgcinfecthiv", "eptctinfecthiv",
  "eptgcctinfecthiv",
  "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
  "eptgcctinfecthiv_inst",
  "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
  "eptgcctinfectundiaghiv_inst",
  "eptindexprovided_gc", "eptindexprovided_ct",
  "eptpartprovided_gc", "eptpartprovided_ct",
  "eptpartprovided_main", "eptpartprovided_pers",
  "eptpartprovided_inst", "eptpartuptake_main",
  "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
  "eptpartuptake_pers", "eptpartuptake_inst",
  "eptpartuptake_gc", "eptpartuptake_ct")
i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "followup/sim.n3037.rda", compress = "gzip")


## Locally merge files --------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

sim <- merge_simfiles(simno = 3118, indir = "data/", ftype = "min")
sim <- truncate_sim(sim, at = 5200)
vars.needed <-  c(# HIV
  "incid", "ir100", "hivtests.nprep", "i.prev",

  # GC
  "incid.gc", "ir100.gc", "ir100.rgc", "ir100.ugc",
  "prev.gc", "prev.uct", "prev.rct",
  "GCasympttests", "GCsympttests",
  "txGC","txGC_asympt",

  # CT
  "incid.ct", "ir100.ct", "ir100.rct", "ir100.uct",
  "prev.ct", "prev.uct", "prev.rct",
  "CTasympttests", "GCsympttests",
  "txCT","txCT_asympt",

  # Combined
  "incid.gcct", "ir100.gcct",
  "prev.STI", "prev.sti.tttraj1", "prev.sti.tttraj2",
  "stiasympttests", "stisympttests",
  "stiasympttests.tttraj1", "stiasympttests.tttraj2",
  "stisympttests.tttraj1", "stisympttests.tttraj2",
  "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
  "txSTI_asympt",
  "tt.traj.sti1", "tt.traj.sti2",

  # Other
  "num",
  'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
  'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',

  # EPT
  "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
  "eptTx", "propindexeptElig",
  "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
  "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
  "eptgcctinfectundiaghiv",
  "gc.timesInf", "ct.timesInf", "sti.timesInf",
  "eptgcinfecthiv", "eptctinfecthiv",
  "eptgcctinfecthiv",
  "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
  "eptgcctinfecthiv_inst",
  "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
  "eptgcctinfectundiaghiv_inst",
  "eptindexprovided_gc", "eptindexprovided_ct",
  "eptpartprovided_gc", "eptpartprovided_ct",
  "eptpartprovided_main", "eptpartprovided_pers",
  "eptpartprovided_inst", "eptpartuptake_main",
  "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
  "eptpartuptake_pers", "eptpartuptake_inst",
  "eptpartuptake_gc", "eptpartuptake_ct")

i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "data/followup/sim.3118.rda", compress = "gzip")


# 1 by 1 burnin on Hyak---------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
sim <- merge_simfiles(simno = 1000, indir = "data/", ftype = "max")
save(sim, file = "data/sim.n1000.rda", compress = "gzip")


#### Merge 1 by 1 on Hyak ------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

sim <- merge_simfiles(simno = 6041, indir = "data/", ftype = "min")
sim <- truncate_sim(sim, at = 5200)
vars.needed <- c(# HIV
  "incid", "ir100", "hivtests.nprep", "i.prev",

  # GC
  "incid.gc", "ir100.gc", "ir100.rgc", "ir100.ugc",
  "prev.gc", "prev.uct", "prev.rct",
  "GCasympttests", "GCsympttests",
  "txGC","txGC_asympt",

  # CT
  "incid.ct", "ir100.ct", "ir100.rct", "ir100.uct",
  "prev.ct", "prev.uct", "prev.rct",
  "CTasympttests", "GCsympttests",
  "txCT","txCT_asympt",

  # Combined
  "incid.gcct", "ir100.gcct",
  "prev.STI", "prev.sti.tttraj1", "prev.sti.tttraj2",
  "stiasympttests", "stisympttests",
  "stiasympttests.tttraj1", "stiasympttests.tttraj2",
  "stisympttests.tttraj1", "stisympttests.tttraj2",
  "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
  "txSTI_asympt",
  "tt.traj.sti1", "tt.traj.sti2",

  # Other
  "num",
  'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
  'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',

  # EPT
  "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
  "eptTx", "propindexeptElig",
  "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
  "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
  "eptgcctinfectundiaghiv",
  "gc.timesInf", "ct.timesInf", "sti.timesInf",
  "eptgcinfecthiv", "eptctinfecthiv",
  "eptgcctinfecthiv",
  "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
  "eptgcctinfecthiv_inst",
  "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
  "eptgcctinfectundiaghiv_inst",
  "eptindexprovided_gc", "eptindexprovided_ct",
  "eptpartprovided_gc", "eptpartprovided_ct",
  "eptpartprovided_main", "eptpartprovided_pers",
  "eptpartprovided_inst", "eptpartuptake_main",
  "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
  "eptpartuptake_pers", "eptpartuptake_inst",
  "eptpartuptake_gc", "eptpartuptake_ct")
i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "data/sim.n3000.rda", compress = "gzip")

### More Hyak merge------------------------------------------------------------------

# sims <- c(3001:3009, 3018, 3027, 3036, 3045, 3054, 3063, 3072, 3081, 3090, 3099, 3108, 3117, 3126, 3135, 3144, 3153, 3162, 3171, 3180, 3189:3198, 3221:3513)
#sims <- c(4001:4009, 4018, 4027, 4036, 4045, 4054, 4063, 4072, 4081, 4090, 4099, 4108, 4117, 4126, 4135, 4144, 4153, 4162, 4171, 4180, 6001:6009)
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

sims <- c(8000:8070)
for (i in sims) {

  sim <- merge_simfiles(simno = i, indir = "data/", ftype = "min")
  sim <- truncate_sim(sim, at = 2600)
  vars.needed <- c(# HIV
    "incid", "ir100", "hivtests.nprep", "i.prev",

    # GC
    "incid.gc", "ir100.gc", "ir100.rgc", "ir100.ugc",
    "prev.gc", "prev.uct", "prev.rct",
    "GCasympttests", "GCsympttests",
    "txGC","txGC_asympt",

    # CT
    "incid.ct", "ir100.ct", "ir100.rct", "ir100.uct",
    "prev.ct", "prev.uct", "prev.rct",
    "CTasympttests", "GCsympttests",
    "txCT","txCT_asympt",

    # Combined
    "incid.gcct", "ir100.gcct",
    "prev.STI", "prev.sti.tttraj1", "prev.sti.tttraj2",
    "stiasympttests", "stisympttests",
    "stiasympttests.tttraj1", "stiasympttests.tttraj2",
    "stisympttests.tttraj1", "stisympttests.tttraj2",
    "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
    "txSTI_asympt",
    "tt.traj.sti1", "tt.traj.sti2",

    # Other
    "num",
    'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
    'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',

    # EPT
    "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
    "eptTx", "propindexeptElig",
    "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
    "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
    "eptgcctinfectundiaghiv",
    "gc.timesInf", "ct.timesInf", "sti.timesInf",
    "eptgcinfecthiv", "eptctinfecthiv",
    "eptgcctinfecthiv",
    "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
    "eptgcctinfecthiv_inst",
    "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
    "eptgcctinfectundiaghiv_inst",
    "eptindexprovided_gc", "eptindexprovided_ct",
    "eptpartprovided_gc", "eptpartprovided_ct",
    "eptpartprovided_main", "eptpartprovided_pers",
    "eptpartprovided_inst", "eptpartuptake_main",
    "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
    "eptpartuptake_pers", "eptpartuptake_inst",
    "eptpartuptake_gc", "eptpartuptake_ct")
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  filename <- paste0("data/sim.n", i, ".rda")
  save(sim, file = filename, compress = "gzip")

}

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
library("EpiModel")


## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
eptcov <- as.numeric(Sys.getenv("EPTCOV"))
eptint <- as.numeric(Sys.getenv("EPTINT"))
prov.main.ong <- as.numeric(Sys.getenv("PROVMAINONG"))
prov.pers.ong <- as.numeric(Sys.getenv("PROVPERSONG"))
prov.main.end <- as.numeric(Sys.getenv("PROVMAINEND"))
prov.pers.end <- as.numeric(Sys.getenv("PROVPERSEND"))
prov.inst <- as.numeric(Sys.getenv("PROVINST"))
uptake.main <- as.numeric(Sys.getenv("UPTAKEMAIN"))
uptake.pers <- as.numeric(Sys.getenv("UPTAKEPERS"))
uptake.inst <- as.numeric(Sys.getenv("UPTAKEINST"))
gctxsuccess <- as.numeric(Sys.getenv("GCTXSUCCESS"))
cttxsuccess <- as.numeric(Sys.getenv("CTTXSUCCESS"))

cat("Array number is ", jobno)
cat("\n fsimno is ", fsimno)

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   # AI Scale
                   ai.scale = 1.061338,
                   ai.scale.pospos = 1.061338,

                   # Probability of rectal testing
                   tst.rect.sti.rr = 1,

                   # Correlation
                   sti.correlation.time = 0,

                   # STI acquisition
                   rgc.tprob = 0.5425, #0.5364416,
                   ugc.tprob = 0.4405, #0.434692,
                   rct.tprob = 0.2480, #0.2493814,
                   uct.tprob = 0.1940, #0.1944415,

                   # HIV acquisition
                   hiv.rgc.rr = 2.175918,
                   hiv.ugc.rr = 1.564797,
                   hiv.rct.rr = 2.175918,
                   hiv.uct.rr = 1.564797,

                   # STI symptom probability
                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.80,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.48,

                   # EPT
                   ept.provision.partner.main.ong = prov.main.ong,
                   ept.provision.partner.pers.ong = prov.pers.ong,
                   ept.provision.partner.main.end = prov.main.end,
                   ept.provision.partner.pers.end = prov.pers.end,
                   ept.provision.partner.inst = prov.inst,
                   ept.uptake.partner.main = uptake.main,
                   ept.uptake.partner.pers = uptake.pers,
                   ept.uptake.partner.inst = uptake.inst,
                   ept.gc.success = gctxsuccess,
                   ept.ct.success = cttxsuccess,
                   ept.coverage = eptcov,
                   ept.risk.int = eptint,

                   # Partner cutoff
                   partnercut = 1,

                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0,
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0,
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,

                   # Other
                   prep.coverage = 0,


                   prep.start = 7000,
                   stitest.start = 1,
                   ept.start = 2601,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 3120,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       prev.FUN = prevalence_msm_ept,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/sti.tnt.burnin.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 truncate.at = 2600,
                 vars =
                   c(# HIV
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
                     "eptpartuptake_gc", "eptpartuptake_ct"))


## EPT analyses
## Combined NG/CT
## Coverage, window, provision, uptake
## Incidence, PIA, Total Doses taken by partners, per-capita times infected
## Table 2 ---------------------------------------------------------

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
#source("analysis/fx.R")

# Base - No EPT
# Reference scenario here
#load("data/followup/EPT/sim.n8000.rda")
#load("data/followup/sim.n8000.rda")
load("data/sim.n8000.rda")
sim.base <- sim

incid.base.gcct <- unname(colSums(sim.base$epi$incid.ct)) +
  unname(colSums(sim.base$epi$incid.gc))

sims <- c(8000:8015, 8018:8022, 8016:8017)

qnt.low <- 0.25
qnt.high <- 0.75

eptcov <- rep(NA, length(sims))
eptint <- rep(NA, length(sims))
mainuptake <- rep(NA, length(sims))
persuptake <- rep(NA, length(sims))
instuptake <- rep(NA, length(sims))
mainongprov <- rep(NA, length(sims))
mainendprov <- rep(NA, length(sims))
persongprov <- rep(NA, length(sims))
persendprov <- rep(NA, length(sims))
instprov <- rep(NA, length(sims))
gctxsuccess <- rep(NA, length(sims))
cttxsuccess <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
gcct.nia <- rep(NA, length(sims))
gcct.nnt <- rep(NA, length(sims))

sti.timesInf <- rep(NA, length(sims))
hiv.undiag <- rep(NA, length(sims))
eptuninfectedprovided <- rep(NA, length(sims))

df <- data.frame(eptcov, eptint, mainuptake, persuptake, instuptake,
                 mainongprov, mainendprov, persongprov, persendprov, instprov,
                 gctxsuccess, cttxsuccess,

                 gcct.incid, gcct.pia, gcct.nia, gcct.nnt, sti.timesInf,
                 hiv.undiag, eptuninfectedprovided
)

for (i in seq_along(sims)) {

  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  #fn <- list.files("data/followup/EPT/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$eptcov[i] <- sim$param$ept.coverage
  df$eptint[i] <- sim$param$ept.risk.int
  df$mainuptake[i] <- sim$param$ept.uptake.partner.main
  df$persuptake[i] <- sim$param$ept.uptake.partner.pers
  df$instuptake[i] <- sim$param$ept.uptake.partner.inst
  df$mainongprov[i] <- sim$param$ept.provision.partner.main.ong
  df$mainendprov[i] <- sim$param$ept.provision.partner.main.end
  df$persongprov[i] <- sim$param$ept.provision.partner.pers.ong
  df$persendprov[i] <- sim$param$ept.provision.partner.pers.end
  df$instprov[i] <- sim$param$ept.provision.partner.inst
  df$gctxsuccess[i] <- sim$param$ept.gc.success
  df$cttxsuccess[i] <- sim$param$ept.ct.success

  # Incidence Rate over last year
  vec.ir.gcct <- unname(colMeans(tail(sim$epi$ir100.ct, 52))) + unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             ", ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")

  # PIA (Cumulative)
  incid.gcct <- unname(colSums(sim$epi$incid.ct)) +
    unname(colSums(sim$epi$incid.gc))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  df$gcct.nia[i] <- paste0(round(quantile(vec.nia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(vec.nia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           ", ", round(quantile(vec.nia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           ", ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  # Doses
  # Number of doses provided to index partners
  eptdoses.gcct <- unname(colSums(sim$epi$eptindexprovided_gc, na.rm = TRUE)) +
    unname(colSums(sim$epi$eptindexprovided_ct, na.rm = TRUE))

  # NNT
  vec.gcct.nnt <- (eptdoses.gcct) / (incid.base.gcct - incid.gcct)

  if (is.nan(mean(vec.gcct.nnt))) {
    vec.gcct.nnt <- rep(0, 256)
  }

  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           ", ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  # Per-capita number of times infected with an STI
  sti.timesInf <- unname(colMeans(tail(sim$epi$sti.timesInf, 52)))
  df$sti.timesInf[i] <- paste0(round(quantile(sti.timesInf, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(sti.timesInf, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               ", ", round(quantile(sti.timesInf, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")

  # Prevalence of undiagnosed HIV infection among partners
  hiv.undiag <- unname(colMeans(sim$epi$eptgcctinfectundiaghiv / sim$epi$eptpartprovided, na.rm = TRUE))
  df$hiv.undiag[i] <- paste0(round(quantile(100*hiv.undiag, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(100*hiv.undiag, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               ", ", round(quantile(100*hiv.undiag, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")

  # Percent of eligible parners provided EPT who did not have any STI
  vec.eptuninfectedprovided <- unname(colMeans(sim$epi$eptuninfectedprovided / sim$epi$eptpartprovided, na.rm = TRUE))
  df$eptuninfectedprovided[i] <- paste0(round(quantile(100*vec.eptuninfectedprovided, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(100*vec.eptuninfectedprovided, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        ", ", round(quantile(100*vec.eptuninfectedprovided, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")

  cat("*")

}

View(df)

write.csv(df, "analysis/EPT Table 2.csv")


library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))
suppressMessages(library("EasyABC"))

f <- function(x) {

  set.seed(x[1])

  suppressMessages(library("EpiModelHIV"))

  data(st)
  param <- param_msm(nwstats = st,

                     rgc.dur.asympt = x[2],
                     ugc.dur.asympt = x[2],
                     rct.dur.asympt = x[3],
                     uct.dur.asympt = x[3],

                     hiv.rgc.rr = x[4],
                     hiv.ugc.rr = x[5],
                     hiv.rct.rr = x[4],
                     hiv.uct.rr = x[5])

  init <- init_msm(nwstats = st)

  control <- control_msm(simno = 1,
                         nsteps = 2600,
                         nsims = 1, ncores = 1,
                         verbose = FALSE)

  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 52)

  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.prev <- mean(df$i.prev)

  out <- c(gc.incid, ct.incid, hiv.prev)

  return(out)
}

priors <- list(c("unif", 32, 38),
               c("unif", 42, 48),
               c("unif", 2.5, 3),
               c("unif", 1.5, 2))

targets <- c(4.2, 6.6, 0.26)

( nsim <- as.numeric(Sys.getenv("NSIM")) )
( pacc <- as.numeric(Sys.getenv("PACC")) )

a <- ABC_sequential(method = "Lenormand",
                    model = f,
                    prior = priors,
                    nb_simul = nsim,
                    summary_stat_target = targets,
                    p_acc_min = pacc,
                    progress_bar = TRUE,
                    n_cluster = 16,
                    use_seed = TRUE,
                    verbose = FALSE)

fn <- paste0("data/smc3.", pacc*100, "pct.", nsim, "sim.rda")
save(a, file = fn)

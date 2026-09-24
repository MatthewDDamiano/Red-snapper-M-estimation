# Function to write dat file for BAM

BAM.write.dat <- function(fname, styr, styr.Rdevs,
                          nyr, nages,
                          spawn.parms = spawn.parms,
                          lcomp.config = lcomp.config,
                          dat.survey, dat.L1, dat.L2, dat.D,
                          parms, dev.tseries, prop.f, mat.age,
                          wgt.parms, fec.parms, len.cm, M.lor.scale) {
  sink(fname)

  ########## GENERAL HEADER INPUT############################################################
  cat(as.character(paste("#######BAM SIM (BAM MSE & General Simulation) INPUT#######", collapse = "\t")), sep = "\n")
  # Starting and ending year of the model (year data starts)
  cat(as.character(paste("# Model start year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(styr, collapse = "\t")), sep = "\n")
  cat(as.character(paste("# Number of years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(nyr, collapse = "\t")), sep = "\n", append = TRUE)

  # Starting and ending year to estimate recruitment deviation from S-R curve
  cat(as.character(paste("# S-R start year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(styr.Rdevs, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# S-R end year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(nyr, collapse = "\t")), sep = "\n", append = TRUE)

  # Possible 3 phases of constraints on recruitment deviations
  cat(as.character(paste("# Possible 3 phases of recdev constraint", collapse = "\t")), sep = "\n")
  cat(as.character(paste(5, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste((nyr - 3), collapse = "\t")), sep = "\n", append = TRUE)

  # ending year for selectivity blocks (not used in this sim, but architecture left in place)
  cat(as.character(paste("# End year for selectivity time blocks (if used)", collapse = "\t")), sep = "\n")
  cat(as.character(paste((nyr - 10), collapse = "\t")), sep = "\n", append = TRUE)

  # SR switch and spr proxy (as a ratio) # currently sets Fproxy to F40
  cat(as.character(paste("# SPR switch", collapse = "\t")), sep = "\n")
  cat(as.character(paste(spawn.parms$SR.switch, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# SPR proxy", collapse = "\t")), sep = "\n")
  cat(as.character(paste(spawn.parms$SPR.proxy, collapse = "\t")), sep = "\n", append = TRUE)

  # nages and vector of ages in the population
  cat(as.character(paste("# Number of ages in the population", collapse = "\t")), sep = "\n")
  cat(as.character(paste(nages, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Vector of ages in population", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1:nages, collapse = "\t")), sep = "\n", append = TRUE)

  # nages and vector of ages in the age comps
  cat(as.character(paste("# Number of ages in age comps data", collapse = "\t")), sep = "\n")
  cat(as.character(paste(nages, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Vector ages in age comps data", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1:nages, collapse = "\t")), sep = "\n", append = TRUE)

  # nlength bins, width of bins, and vector of length bins in length comps
  cat(as.character(paste("# Number of length bins for length comps data", collapse = "\t")), sep = "\n")
  cat(as.character(paste(lcomp.config$nlenbins, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Bin width", collapse = "\t")), sep = "\n")
  cat(as.character(paste(lcomp.config$binw, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Vector of length bins in length comps data", collapse = "\t")), sep = "\n")
  cat(as.character(paste(lcomp.config$lenbins, collapse = "\t")), sep = "\n", append = TRUE)

  # Max F used in spr and msy calcs
  cat(as.character(paste("# Max F used in SPR and MSY calcs", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1, collapse = "\t")), sep = "\n", append = TRUE)
  # Total number of iterations for spr calcs
  cat(as.character(paste("# Number of iterations for SPR calcs", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1001, collapse = "\t")), sep = "\n", append = TRUE)
  # Number years at end of time series over which to average sector F's, for weighted selectivities
  cat(as.character(paste("# Number of years at end of time series to avg fleet Fs (for weighted selex)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(3, collapse = "\t")), sep = "\n", append = TRUE)
  # Bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
  cat(as.character(paste("# S-R bias correction", collapse = "\t")), sep = "\n")
  cat(as.character(paste(spawn.parms$bam.bc.switch, collapse = "\t")), sep = "\n", append = TRUE)

  ########## SURVEY DATA############################################################
  cat(as.character(paste("#######SURVEY DATA#######", collapse = "\t")), sep = "\n")
  # number of years
  cat(as.character(paste("# Survey/index number of years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$nyr, collapse = "\t")), sep = "\n", append = TRUE)
  # survey years
  cat(as.character(paste("# Survey/index years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$yrs, collapse = "\t")), sep = "\n", append = TRUE)
  # observed index
  cat(as.character(paste("# Survey/index observed values", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$vals.obs, collapse = "\t")), sep = "\n", append = TRUE)
  # cv's of observed index
  cat(as.character(paste("# Survey/index cv of observed values", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$cv, collapse = "\t")), sep = "\n", append = TRUE)

  # number of years of age comps
  cat(as.character(paste("# Survey/index number of age comp years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$nyr.ages, collapse = "\t")), sep = "\n", append = TRUE)
  # years of age comps
  cat(as.character(paste("# Survey/index age comp years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$yrs.age, collapse = "\t")), sep = "\n", append = TRUE)
  # sample size of age comps
  cat(as.character(paste("# Survey/index age comp sample size", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$nsamp, collapse = "\t")), sep = "\n", append = TRUE)
  # nfish of age comps
  cat(as.character(paste("# Survey/index number of fish in age comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.survey$nfish, collapse = "\t")), sep = "\n", append = TRUE)
  # age comps
  cat(as.character(paste("# Survey/index observed age comps", collapse = "\t")), sep = "\n")
  for (i in 1:dat.survey$nyr.ages) {
    cat(as.character(paste(dat.survey$acomp[i, ], collapse = "\t")), sep = "\n", append = TRUE)
  }

  ########## LANDINGS DATA############################################################
  ### FLEET 1
  cat(as.character(paste("#######FLEET 1 DATA#######", collapse = "\t")), sep = "\n")
  # start year of landings
  cat(as.character(paste("# Fleet 1 landings start year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$styr, collapse = "\t")), sep = "\n", append = TRUE)
  # end year of landings
  cat(as.character(paste("# Fleet 1 landings end year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$endyr, collapse = "\t")), sep = "\n", append = TRUE)
  # observed L
  cat(as.character(paste("# Fleet 1 observed landings", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$vals.obs, collapse = "\t")), sep = "\n", append = TRUE)
  # cv's of observed L
  cat(as.character(paste("# Fleet 1 cv of observed landings", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$cv, collapse = "\t")), sep = "\n", append = TRUE)

  # number of years of age comps
  cat(as.character(paste("# Fleet 1 number of years of age comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$nyr.ages, collapse = "\t")), sep = "\n", append = TRUE)
  # years of age comps
  cat(as.character(paste("# Fleet 1 age comps years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$yrs.age, collapse = "\t")), sep = "\n", append = TRUE)
  # sample size of age comps
  cat(as.character(paste("# Fleet 1 age comps sample size", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$nsamp, collapse = "\t")), sep = "\n", append = TRUE)
  # nfish of age comps
  cat(as.character(paste("# Fleet 1 number of fish in age comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L1$nfish, collapse = "\t")), sep = "\n", append = TRUE)
  # age comps
  cat(as.character(paste("# Fleet 1 observed age comps", collapse = "\t")), sep = "\n")
  for (i in 1:dat.L1$nyr.ages) {
    cat(as.character(paste(dat.L1$acomp[i, ], collapse = "\t")), sep = "\n", append = TRUE)
  }

  ### FLEET 2
  cat(as.character(paste("#######FLEET 2 DATA#######", collapse = "\t")), sep = "\n")
  # start year of landings
  cat(as.character(paste("# Fleet 2 start year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$styr, collapse = "\t")), sep = "\n", append = TRUE)
  # end year of landings
  cat(as.character(paste("# Fleet 2 end year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$endyr, collapse = "\t")), sep = "\n", append = TRUE)
  # observed L
  cat(as.character(paste("# Fleet 2 observed landings", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$vals.obs, collapse = "\t")), sep = "\n", append = TRUE)
  # cv's of observed L
  cat(as.character(paste("# Fleet 2 cv of observed landings", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$cv, collapse = "\t")), sep = "\n", append = TRUE)

  # number of years of age comps
  cat(as.character(paste("# Fleet 2 number of years of age comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$nyr.ages, collapse = "\t")), sep = "\n", append = TRUE)
  # years of age comps
  cat(as.character(paste("# Fleet 2 age comps years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$yrs.age, collapse = "\t")), sep = "\n", append = TRUE)
  # sample size of age comps
  cat(as.character(paste("# Fleet 2 age comps sample size", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$nsamp, collapse = "\t")), sep = "\n", append = TRUE)
  # nfish of age comps
  cat(as.character(paste("# Fleet 2 number of fish in age comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.L2$nfish, collapse = "\t")), sep = "\n", append = TRUE)
  # age comps
  cat(as.character(paste("# Fleet 2 observed age comps", collapse = "\t")), sep = "\n")
  for (i in 1:dat.L2$nyr.ages) {
    cat(as.character(paste(dat.L2$acomp[i, ], collapse = "\t")), sep = "\n", append = TRUE)
  }
  ### FLEET 3
  cat(as.character(paste("#######FLEET 3 DATA#######", collapse = "\t")), sep = "\n")
  # start year of discards
  cat(as.character(paste("# Fleet 3 start year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$styr, collapse = "\t")), sep = "\n", append = TRUE)
  # end year of discards
  cat(as.character(paste("# Fleet 3 end year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$endyr, collapse = "\t")), sep = "\n", append = TRUE)
  # observed discards
  cat(as.character(paste("# Fleet 3 observed discards", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$vals.obs, collapse = "\t")), sep = "\n", append = TRUE)
  # cv's of observed D
  cat(as.character(paste("# Fleet 3 cv of observed discards", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$cv, collapse = "\t")), sep = "\n", append = TRUE)

  # number of years of length comps
  cat(as.character(paste("# Fleet 3 number of years of length comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$D.lc.nyr, collapse = "\t")), sep = "\n", append = TRUE)
  # years of length comps
  cat(as.character(paste("# Fleet 3 length comps years", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$yrs.len, collapse = "\t")), sep = "\n", append = TRUE)
  # sample size of length comps
  cat(as.character(paste("# Fleet 3 length comps sample size", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$nsamp, collapse = "\t")), sep = "\n", append = TRUE)
  # nfish of length comps
  cat(as.character(paste("# Fleet 3 number of fish in length comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dat.D$nfish, collapse = "\t")), sep = "\n", append = TRUE)
  # length comps
  cat(as.character(paste("# Fleet 3 observed length comps", collapse = "\t")), sep = "\n")
  for (i in 1:dat.D$D.lc.nyr) {
    cat(as.character(paste(dat.D$lcomp[i, ], collapse = "\t")), sep = "\n", append = TRUE)
  }
  
  ########## CKMR Data ############################################################
  
  # These were experimental at the time and not used in the current model
  cat(as.character(paste("#######CKMR DATA#######", collapse = "\t")), sep = "\n")
  # first ckmr sample (1)
  cat(as.character(paste(ckmr.dat$ckmr_fsamp, collapse = "\t")), sep = "\n", append = TRUE)
  # last ckmr sample
  cat(as.character(paste(ckmr.dat$ckmr_lsamp, collapse = "\t")), sep = "\n", append = TRUE)
  # cohort vector 1
  cat(as.character(paste(ckmr.dat$c1, collapse = "\t")), sep = "\n", append = TRUE)
  # cohort vector 2
  cat(as.character(paste(ckmr.dat$c2, collapse = "\t")), sep = "\n", append = TRUE)
  # sampling year vector
  cat(as.character(paste(ckmr.dat$s, collapse = "\t")), sep = "\n", append = TRUE)
  # hsp vector
  cat(as.character(paste(ckmr.dat$hsp, collapse = "\t")), sep = "\n", append = TRUE)
  # pop vector
  cat(as.character(paste(ckmr.dat$pop, collapse = "\t")), sep = "\n", append = TRUE)
  # # 'not' vector
  cat(as.character(paste(ckmr.dat$not, collapse = "\t")), sep = "\n", append = TRUE)

  ########## PARAMETER SECTION############################################################
  cat(as.character(paste("#######PARAMETER SECTION#######", collapse = "\t")), sep = "\n")
  nparm <- nrow(parms)
  for (i in 1:nparm) {
    cat(as.character(paste(parms[i, 1:7], collapse = "\t")), sep = "\n", append = TRUE)
  }

  ########## DEV CONSTRAINTS and DEVS############################################################
  cat(as.character(paste("#######DEV CONSTRAINTS and DEVS#######", collapse = "\t")), sep = "\n")
  cat(as.character(paste("# Fleet 1 F devs (lower bound, upper bound, phase)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(c(-10, 5, dev.tseries$F1.phase), collapse = "\t")), sep = "\n", append = TRUE)
  # Fleet2 F devs (lower bound, upper bound, phase)
  cat(as.character(paste("# Fleet 2 F devs (lower bound, upper bound, phase)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(c(-10, 5, dev.tseries$F2.phase), collapse = "\t")), sep = "\n", append = TRUE)
  # Fleet3 F devs (lower bound, upper bound, phase)
  cat(as.character(paste("# Fleet 3 F devs (lower bound, upper bound, phase)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(c(-10, 5, dev.tseries$F3.phase), collapse = "\t")), sep = "\n", append = TRUE)
  # Recruitment devs (lower bound, upper bound, phase)
  cat(as.character(paste("# Rec devs (lower bound, upper bound, phase)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(c(-5, 5, dev.tseries$R.phase), collapse = "\t")), sep = "\n", append = TRUE)
  # Nage devs (lower bound, upper bound, phase)
  cat(as.character(paste("# N at age devs (lower bound, upper bound, phase)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(c(-5, 5, -3), collapse = "\t")), sep = "\n", append = TRUE)

  # Fleet1 F devs (initial guesses)
  cat(as.character(paste("# Fleet 1 F devs (Initial guesses)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dev.tseries$logF1.dev, collapse = "\t")), sep = "\n", append = TRUE)
  # Fleet2 F devs (initial guesses)
  cat(as.character(paste("# Fleet 2 F devs (Initial guesses)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dev.tseries$logF2.dev, collapse = "\t")), sep = "\n", append = TRUE)
  # Fleet3 F devs (initial guesses)
  cat(as.character(paste("# Fleet 3 F devs (Initial guesses)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dev.tseries$logF3.dev, collapse = "\t")), sep = "\n", append = TRUE)
  # Recruitment devs (initial guesses)
  cat(as.character(paste("# Rec devs (Initial guesses)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(dev.tseries$logR.dev, collapse = "\t")), sep = "\n", append = TRUE)
  # Nage devs (initial guesses)
  cat(as.character(paste("# N ate age devs (Initial guesses)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(rep(0, (nages - 1)), collapse = "\t")), sep = "\n", append = TRUE)

  ########## LIKELIHOOD WEIGHTS############################################################
  cat(as.character(paste("#######LIKELIHOOD WEIGHTS#######", collapse = "\t")), sep = "\n")
  # Landings (all fleets), discards, and survey index
  cat(as.character(paste("# Fleet 1, Fleet 2 Landings LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fleet 3 Discards LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE) # added for discard likelihood weight
  cat(as.character(paste("# Survey/index values LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  # Age comps for fleet1,fleet2, length comps for fleet 3, age comps for survey
  cat(as.character(paste("# Fleet 1 age comps LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fleet 2 age comps LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fleet 3 length comps LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Survey/Index age comps LL weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  # Additional weights: 1.Nage_init, 2.rec devs, 3.rec devs early, 4.rec devs late, 5.Ftune
  cat(as.character(paste("# Additional weights", collapse = "\t")), sep = "\n")
  cat(as.character(paste("# Nage_init weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(0.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# rec devs weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# rec devs early weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(0.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# rec devs late weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(0.0, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Ftune weight", collapse = "\t")), sep = "\n")
  cat(as.character(paste(0.0, collapse = "\t")), sep = "\n", append = TRUE)

  ########## MISCELLANEOUS STUFF SECTION############################################################
  cat(as.character(paste("#######MISCELLANEOUS STUFF SECTION#######", collapse = "\t")), sep = "\n")
  # length-weight (TL-whole wgt) coefficients a and b, W=aL^b (W assumed in kg)
  cat(as.character(paste("# Weight-length a par", collapse = "\t")), sep = "\n")
  cat(as.character(paste(wgt.parms$a.lw, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Weight-length b par", collapse = "\t")), sep = "\n")
  cat(as.character(paste(wgt.parms$b.lw, collapse = "\t")), sep = "\n", append = TRUE)
  # prop female, maturity, fecundity, and M at age
  cat(as.character(paste("# Proportion female", collapse = "\t")), sep = "\n")
  cat(as.character(paste(prop.f, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Maturity at age", collapse = "\t")), sep = "\n")
  cat(as.character(paste(mat.age, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fecundity parameters: batch size", collapse = "\t")), sep = "\n")
  cat(as.character(paste(fec.parms$batch, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fecundity parameters: a", collapse = "\t")), sep = "\n")
  cat(as.character(paste(fec.parms$fec.a, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fecundity parameters: b", collapse = "\t")), sep = "\n")
  cat(as.character(paste(fec.parms$fec.b, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fecundity parameters: c", collapse = "\t")), sep = "\n")
  cat(as.character(paste(fec.parms$fec.c, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fecundity parameters: thresh", collapse = "\t")), sep = "\n")
  cat(as.character(paste(fec.parms$fec.thresh, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Fecundity parameters: fecmin", collapse = "\t")), sep = "\n")
  cat(as.character(paste(fec.parms$fec.min, collapse = "\t")), sep = "\n", append = TRUE)
  # Natural mortality parameters
  cat(as.character(paste("# Natural mortality lengths (for Lorenzen fcn)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(len.cm, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Natural mortality scalar (for Lorenzen fcn)", collapse = "\t")), sep = "\n")
  cat(as.character(paste(M.lor.scale, collapse = "\t")), sep = "\n", append = TRUE)
  # cat(as.character(paste("# Natural mortality at age vector", collapse = "\t")), sep = "\n")
  # cat(as.character(paste(M.age, collapse = "\t")), sep = "\n", append = TRUE)
  # spawn time fraction
  cat(as.character(paste("# Spawn fraction", collapse = "\t")), sep = "\n")
  cat(as.character(paste(spawn.parms$sp.frac, collapse = "\t")), sep = "\n", append = TRUE)
  # tuning year, tuning F, min SS for age comps, min SS for length comps
  cat(as.character(paste("# Tuning year", collapse = "\t")), sep = "\n")
  cat(as.character(paste(10, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Tuning F", collapse = "\t")), sep = "\n")
  cat(as.character(paste(0.2, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Minimum sample size for age comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Minimum sample size for length comps", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1, collapse = "\t")), sep = "\n", append = TRUE)
  # ageing error (identity matrix for none)
  # cat(as.character(paste("# Indentity matrix (for ageing error)", collapse = "\t")), sep = "\n")
  I <- diag(nages)
  cat(as.character(paste("# Ageing error matrix (if Identity then none)", collapse = "\t")), sep = "\n")
  for (i in 1:nages) {
    cat(as.character(paste(I[i, ], collapse = "\t")), sep = "\n", append = TRUE)
  }
  # projections: endyr, styr regs, Fproj type, Fproj mult
  cat(as.character(paste("# Projection end year", collapse = "\t")), sep = "\n")
  cat(as.character(paste((nyr + 10), collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# Projection start year", collapse = "\t")), sep = "\n")
  cat(as.character(paste((nyr + 1), collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# F Projection type", collapse = "\t")), sep = "\n")
  cat(as.character(paste(2, collapse = "\t")), sep = "\n", append = TRUE)
  cat(as.character(paste("# F Projection mult", collapse = "\t")), sep = "\n")
  cat(as.character(paste(1, collapse = "\t")), sep = "\n", append = TRUE)

  ########## END OF DATA FILE############################################################
  cat(as.character(paste("#######END OF DATA FILE#######", collapse = "\t")), sep = "\n")
  cat(as.character(paste(999, collapse = "\t")), sep = "\n", append = TRUE, fill = TRUE)
  sink()
  closeAllConnections()
}

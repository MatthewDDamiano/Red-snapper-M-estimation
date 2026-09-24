##################################################################
# Run multiple red snapper simulations and summarize relative error and LL profiles
# 11/12/2024
# Last update: 9/24/2026 MDD
##################################################################
graphics.off()
rm(list = ls(all = TRUE))
library(here)

dir <- here()
setwd(dir)
n.iter = 2 # number of model iterations
iter = 1 # starting iteration
n.vals = 7 # set to the number of values across which the model is profiled,
           # which depends on which parameter is being profiled
years = 70 # when measuring time series quantities over multiple iterations

# storage for relative error
mscale <- array(NA,dim=c(n.iter))
mb <- array(NA,dim=c(n.iter))
r0 <- array(NA,dim=c(n.iter))
# steep <- array(NA,dim=c(n.iter))
max.grad <- array(NA,dim=c(n.iter))
true.f40 <- array(NA,dim=c(n.iter))
pred.f40 <- array(NA,dim=c(n.iter))
true.ssbf40 <- array(NA, dim=c(n.iter))
pred.ssbf40 <- array(NA, dim=c(n.iter))
true.N <- array(NA, dim=c(n.iter,years))
pred.N <- array(NA, dim=c(n.iter,years+1))

# storage for LL profiles
lk.dat.mat <- array(NA, dim=c(n.iter,n.vals))
lk.tot.mat <- array(NA, dim=c(n.iter,n.vals))
lk.sr.mat <- array(NA, dim=c(n.iter,n.vals))
lk.cpue.mat <- array(NA, dim=c(n.iter,n.vals))
lk.sagec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.fl1agec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.fl2agec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.allagec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.lenc.mat <- array(NA, dim=c(n.iter,n.vals))
# non-convergence counter
nonconv <- 0

# READ THIS BEFORE RUNNING
# set the working directory to the folder for whatever model configuration you're testing 
# and remember to TURN OFF any old code clearing at the top of 1RunSim.R file within that folder

while(iter <= n.iter){
  setwd("Results/Model 1") 
  unlink("bam-sim.std") # remove old std file
  # run simulation model with stochastic error
  # ptype = NULL
  source("1RunSim.R") 
  # check for convergence and 
  fileNames=c("bam-sim.std")
  # if (file.exists(fileNames)==TRUE){ # check model convergence by looking for .std file
  #   unlink("bam-sim.std") # remove old std file
  # }
  # 
  while(file.exists(fileNames)==FALSE){ 
    nonconv = nonconv + 1 # record nonconvergence if it occurs 
    source("1RunSim.R") # re-run 
  }
  
  # record results from run
  # steep[iter] <- bam$parms$BH.steep
  mscale[iter] <- bam$parms$M.scale
  mb[iter] <- bam$parms$M.b
  r0[iter] <- bam$parms$R0 #/exp(bam$parms$R.sigma.par^2/2)
  # steep[iter] <- bam$parms$BH.steep
  max.grad[iter] <- bam$like[28]
  true.f40[iter] <- sim$spr$F40
  pred.f40[iter] <- bam$parms$Fproxy
  true.ssbf40[iter] <- sim$spr$SSB.Fproxy
  pred.ssbf40[iter] <- bam$parms$SSB.Fproxy
  true.N[iter,] <- sim$abundance
  pred.N[iter,] <- bam$t.series$N

  # rename files for profiling
  datname <- 'BAM-Sim.dat'
  rdatname <- 'BAM-Sim.rdat'
  newrdatname <- paste0("Run_", iter, "-set.rdat")
  newdatname <- paste0("Run_",iter,"-set.dat")
  file.rename(rdatname, newrdatname)
  file.rename(datname, newdatname)

  # run profile code
  setwd(dir)
  source("01run_profile_iter.r") # remember that set.seed is used here for profiling
  source("01compare_runs_iter_R0.r")

  # record LLs
  lk.dat.mat[iter,] <- lk.dat
  lk.tot.mat[iter,] <- lk.tot
  lk.sr.mat[iter,] <- lk.SR
  lk.cpue.mat[iter,] <- lk.cpue
  lk.sagec.mat[iter,] <- lk.age.survey1
  lk.fl1agec.mat[iter,] <- lk.age.fleet1
  lk.fl2agec.mat[iter,] <- lk.age.fleet2
  lk.allagec.mat[iter,] <- lk.age
  lk.lenc.mat[iter,] <- lk.len

  # update iteration for the loop
  set.seed(NULL) # unset the seed from LL profiling
  iter = iter+1
}

# save outputs - add steepness back in when using BH
save(n.vals,vals,
     lk.dat.mat,lk.tot.mat,lk.sr.mat,lk.cpue.mat,lk.sagec.mat,lk.fl1agec.mat,lk.fl2agec.mat,
     lk.allagec.mat,lk.lenc.mat,
     mscale, mb, r0, max.grad, true.f40, pred.f40, true.ssbf40, pred.ssbf40,
     true.N, pred.N, file = "P1R0NoMEst_100iter_Results.Rdata")

# load or re-load as needed
rm(list=ls())
setwd("Results/")
# these were the names of results I saved to the Results folder 
# load("P1_100iter_Results.Rdata")
# load("P2a_100iter_Results.Rdata")
# load("P2b_100iter_Results.Rdata")
# load("P2c_100iter_Results.Rdata")
# load("P2d_100_Results.Rdata")
# load("P3_100iter_NoMest_Results.Rdata")

# load("S1a_100iter_Results.Rdata")
# load("S1b_100iter_Results.Rdata")
# load("S2_100iter_Results.Rdata")
# load("S3_100iter_Results.Rdata")
# load("S4_100iter_Results.Rdata")
# load("S5_100iter_Results.Rdata")
# load("S6_100iter_Results.Rdata")
# load("S7a_100iter_Results.Rdata")
# load("S7b_100iter_Results.Rdata")
# load("S8_100iter_Results.Rdata")
# load("S9_mscale_100iter_Results.Rdata")
# load("P1Nonstat_100iter_Results.Rdata")
# load("P2cNonstat_100iter_Results.Rdata")

# load("P2cAltObs_100iter_Results.Rdata")
# load("P2cOrigObs_100iter_Results.Rdata")
# load('P2cNewErr_100iter_Results.Rdata')

# load("P1R0_100iter_Results.Rdata")
# load("P1R0NoMEst_100iter_Results.Rdata")
# med.true.N.p2d <- apply(true.N, 2, median)
# med.pred.N.p2d <- apply(pred.N, 2, median)
# 
# setwd("C:/Users/mddamiano/Desktop/NOAA projects/Red Snapper/RS M sim_Damiano mod_LorenzenM/Profiling/Results/Model S3")
# years <- 70
# load("P3_100iter_NoMest_Results.Rdata")
# med.true.N.p3 <- apply(true.N, 2, median)
# med.pred.N.p3 <- apply(pred.N, 2, median)
# 
# # windows(height = 8, width = 10, record = T)
# tiff(filename="Fig-3.tiff", width=6, height=8, units="in", compression="none", res=1000)
# mat <- matrix(1:2, ncol = 1, nrow = 2)
# layout(mat = mat, widths = rep.int(1, ncol(mat)), heights = rep.int(1, nrow(mat)))
# par(las = 1, mar = c(5, 5, 0.5, 0.5), cex = 0.75)
# plot(1:years, med.true.N.p2d / 1000,
#      type = "o", lwd = 2, ylab = "P2d Median Abundance (1000 fish)",
#      xlab = "", panel.first = grid(lty = 1), col = "blue",
#      ylim = c(0, max(c(med.true.N.p2d, med.pred.N.p2d[1:years])) / 1000)
# )
# points(1:years, med.pred.N.p3[1:years] / 1000, pch = 3, col = "black", cex = 1.2)
# legend("topright", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))
# plot(1:years, med.true.N.p3 / 1000,
#      type = "o", lwd = 2, ylab = "P3 Median Abundance (1000 fish)",
#      xlab = "Timestep", panel.first = grid(lty = 1), col = "blue",
#      ylim = c(0, max(c(med.true.N.p3, med.pred.N.p3[1:years])) / 1000)
# )
# points(1:years, med.pred.N.p3[1:years] / 1000, pch = 3, col = "black", cex = 1.2)
# legend("topright", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))
# dev.off()
# plot(mscale, r0, xlab = "M_scale", ylab = "R0", lwd = 1.5, cex.main = 1.5, cex.lab =1.5, cex.axis = 1.5)

n.iter = 100
n.vals = 13
profile.vals <- seq(-3, 3, by = 0.5)
vals <- 2 * exp(profile.vals) / (1 + exp(profile.vals))
# profile.vals <- seq(10.0, 16.0, by = 1.0)
# vals <- profile.vals
# profile.vals <- seq(-2.0,-0.1, by = 0.2) # and these for M.b
# profile.vals <- seq(0.5,1.0, by = 0.05) # for steepness

par.name <- "Mscale"
M.b <- -1.0
M.mult <- 1.0
R0 <- 383000
h <- 0.75
# 
# percent coverage of true value
cover <- c()
for(i in 1:100){
  if(mscale[i] > M.mult*0.9 & mscale[i] < M.mult*1.1){
    cover[i] <- 1
  }
  else
    cover[i] <- 0
}
sum(cover)
range(mscale)
# View(lk.tot.mat)
# empty vectors for quantiles
lk.dat.5 <- c()
lk.dat.95 <- c()
lk.tot.5 <- c()
lk.tot.95 <- c()
lk.sr.5 <- c()
lk.sr.95 <- c()
lk.cpue.5 <- c()
lk.cpue.95 <- c()
lk.sagec.5 <- c()
lk.sagec.95 <- c()
lk.fl1agec.5 <- c()
lk.fl1agec.95 <- c()
lk.fl2agec.5 <- c()
lk.fl2agec.95 <- c()
lk.allagec.5 <- c()
lk.allagec.95 <- c()
lk.lenc.5 <- c()
lk.lenc.95 <- c()
# fill quantile vectors using likelihood results from profiling
for(i in 1:n.vals){
  lk.dat.5[i] <- quantile(lk.dat.mat[,i], probs = c(0.05))
  lk.dat.95[i] <- quantile(lk.dat.mat[,i], probs = c(0.95))
  lk.tot.5[i] <- quantile(lk.tot.mat[,i], probs = c(0.05))
  lk.tot.95[i] <- quantile(lk.tot.mat[,i], probs = c(0.95))
  lk.sr.5[i] <- quantile(lk.sr.mat[,i], probs = c(0.05))
  lk.sr.95[i] <- quantile(lk.sr.mat[,i], probs = c(0.95))
  lk.cpue.5[i] <- quantile(lk.cpue.mat[,i], probs = c(0.05))
  lk.cpue.95[i] <- quantile(lk.cpue.mat[,i], probs = c(0.95))
  lk.sagec.5[i] <- quantile(lk.sagec.mat[,i], probs = c(0.05))
  lk.sagec.95[i] <- quantile(lk.sagec.mat[,i], probs = c(0.95))
  lk.fl1agec.5[i] <- quantile(lk.fl1agec.mat[,i], probs = c(0.05))
  lk.fl1agec.95[i] <- quantile(lk.fl1agec.mat[,i], probs = c(0.95))
  lk.fl2agec.5[i] <- quantile(lk.fl2agec.mat[,i], probs = c(0.05))
  lk.fl2agec.95[i] <- quantile(lk.fl2agec.mat[,i], probs = c(0.95))
  lk.allagec.5[i] <- quantile(lk.allagec.mat[,i], probs = c(0.05))
  lk.allagec.95[i] <- quantile(lk.allagec.mat[,i], probs = c(0.95))
  lk.lenc.5[i] <- quantile(lk.lenc.mat[,i], probs = c(0.05))
  lk.lenc.95[i] <- quantile(lk.lenc.mat[,i], probs = c(0.95))
}

# subtract minimum from each iteration so y-axes are comparable; skip if running a single profile
lk.dat.mat <- apply(lk.dat.mat,2,mean)
lk.tot.mat <- apply(lk.tot.mat,2,mean)
lk.sr.mat <- apply(lk.sr.mat,2,mean)
lk.cpue.mat <- apply(lk.cpue.mat,2,mean)
lk.sagec.mat <- apply(lk.sagec.mat,2,mean)
lk.fl1agec.mat <- apply(lk.fl1agec.mat,2,mean)
lk.fl2agec.mat <- apply(lk.fl2agec.mat,2,mean)
lk.allagec.mat <- apply(lk.allagec.mat,2,mean)
lk.lenc.mat <- apply(lk.lenc.mat,2,mean)
# empty vectors for the means
lk.dat.mean <-c()
lk.tot.mean <-c()
lk.sr.mean <-c()
lk.cpue.mean <-c()
lk.sagec.mean <-c()
lk.fl1agec.mean <-c()
lk.fl2agec.mean <-c()
lk.allagec.mean <-c()
lk.lenc.mean <-c()
# standardize the profiles by substracting the minimum
for(i in 1:length(vals)){
  lk.dat.mean[i] <- lk.dat.mat[i]-min(lk.dat.mat) 
  lk.tot.mean[i] <- lk.tot.mat[i]-min(lk.tot.mat)
  lk.sr.mean[i] <- lk.sr.mat[i]-min(lk.sr.mat)
  lk.cpue.mean[i] <- lk.cpue.mat[i]-min(lk.cpue.mat)
  lk.sagec.mean[i] <- lk.sagec.mat[i]-min(lk.sagec.mat)
  lk.fl1agec.mean[i] <- lk.fl1agec.mat[i]-min(lk.fl1agec.mat)
  lk.fl2agec.mean[i] <- lk.fl2agec.mat[i]-min(lk.fl2agec.mat)
  lk.allagec.mean[i] <- lk.allagec.mat[i]-min(lk.allagec.mat)
  lk.lenc.mean[i] <- lk.lenc.mat[i]-min(lk.lenc.mat)
}
# Get rid of nonconverging runs as needed
# c(44,45,74) # remove for S7a
# c(4,49,51, 63) # remove for 8
# c(25,45) # remove for S9 mscale

# lk.tot.mat <- lk.tot.mat[-c(25,45),]
# lk.dat.mat <- lk.dat.mat[-c(25,45),]
# lk.sr.mat <- lk.sr.mat[-c(25,45),]
# lk.cpue.mat <- lk.cpue.mat[-c(25,45),]
# lk.sagec.mat <- lk.sagec.mat[-c(25,45),]
# lk.fl1agec.mat <- lk.fl1agec.mat[-c(25,45),]
# lk.fl2agec.mat <- lk.fl1agec.mat[-c(25,45),]
# lk.allagec.mat <- lk.allagec.mat[-c(25,45),]
# lk.lenc.mat <- lk.lenc.mat[-c(25,45),]

# Redo after removing bad runs
# take mean of centered LL profiles; technically mean delta-NLL
# lk.dat.mean <- apply(lk.dat.mat,2,mean)
# lk.tot.mean <- apply(lk.tot.mat,2,mean)
# lk.sr.mean <- apply(lk.sr.mat,2,mean)
# lk.cpue.mean <- apply(lk.cpue.mat,2,mean)
# lk.sagec.mean <- apply(lk.sagec.mat,2,mean)
# lk.fl1agec.mean <- apply(lk.fl1agec.mat,2,mean)
# lk.fl2agec.mean <- apply(lk.fl2agec.mat,2,mean)
# lk.allagec.mean <- apply(lk.allagec.mat,2,mean)
# lk.lenc.mean <- apply(lk.lenc.mat,2,mean)

# # get 5th and 95th percentiles for plotting uncertainty

# plot single LL profile with all sources
# plot(lk.tot.mean, type = "l", col = "black", ylim = c(0,40))
# lines(lk.dat.mean, type = "l", col = "blue")
# lines(lk.sr.mean, type = "l", col = "red")
# lines(lk.allagec.mean, col = "green")
# lines(lk.fl1agec.mean, col = "pink")
# lines(lk.fl2agec.mean, col = "orange")
# lines(lk.sagec.mean, col = "yellow")
# lines(lk.lenc.mean, col = "brown")

# lk.tot.mean.mat <- matrix(NA, nrow = 13, ncol = 13)
# for(i in 1:13){
#   lk.tot.mean.mat[i,] <- lk.tot.mean
# }
# # bivariate prpfile plot
# persp(vals, profile.vals2, lk.tot.mean.mat, xlab = "Mscale", ylab = "b", zlab = "total NLL", 
#       axes = TRUE, ticktype = "detailed")

# plot mean LL profiles
# windows(10,10)
# par(mfrow=c(2,5))
# plot(vals, lk.tot.mean, type = "o", xlab = par.name, ylab = "mean nLL (penalized)")
# plot(vals, lk.dat.mean, type = "o", xlab = par.name, ylab = "mean nLL data")
# plot(vals, lk.sr.mean, type = "o", xlab = par.name, ylab = "mean nLL SR")
# plot(vals, lk.cpue.mean, type = "o", xlab = par.name, ylab = "mean nLL index")
# plot(vals, lk.allagec.mean, type = "o", xlab = par.name, ylab = "mean nLL ages")
# plot(vals, lk.sagec.mean, type = "o", xlab = par.name, ylab = "mean nLL ages survey1")
# plot(vals, lk.fl1agec.mean, type = "o", xlab = par.name, ylab = "mean nLL ages fleet1")
# plot(vals, lk.fl2agec.mean, type = "o", xlab = par.name, ylab = "mean nLL ages fleet2")
# plot(vals, lk.lenc.mean, type="o", xlab=par.name, ylab="mean nLL lengths")

windows(10,10)
par(mfrow=c(3,3))
plot(vals, lk.tot.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL (total)", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.tot.5), max(lk.tot.95)))
lines(vals, lk.tot.5, type = "c")
lines(vals, lk.tot.95, type = "c")
plot(vals, lk.dat.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL data", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.dat.5),max(lk.dat.95)))
lines(vals, lk.dat.5, type = "c")
lines(vals, lk.dat.95, type = "c")
plot(vals, lk.sr.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL SRR", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.sr.5),max(lk.sr.95)))
lines(vals, lk.sr.5, type = "c")
lines(vals, lk.sr.95, type = "c")
plot(vals, lk.cpue.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL survey", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.cpue.5), max(lk.cpue.95)))
lines(vals, lk.cpue.5, type = "c")
lines(vals, lk.cpue.95, type = "c")
plot(vals, lk.allagec.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL all age comps", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.allagec.5), max(lk.allagec.95)))
lines(vals, lk.allagec.5, type = "c")
lines(vals, lk.allagec.95, type = "c")
plot(vals, lk.sagec.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL survey age comps", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.sagec.5), max(lk.sagec.95)))
lines(vals, lk.sagec.5, type = "c")
lines(vals, lk.sagec.95, type = "c")
plot(vals, lk.fl1agec.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL rec age comps", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.fl1agec.5), max(lk.fl2agec.95)))
lines(vals, lk.fl1agec.5, type = "c")
lines(vals, lk.fl1agec.95, type = "c")
plot(vals, lk.fl2agec.mean, type = "o", xlab = par.name, ylab = "Delta mean NLL comm age comps", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.fl2agec.5), max(lk.fl2agec.95)))
lines(vals, lk.fl2agec.5, type = "c")
lines(vals, lk.fl2agec.95, type = "c")
plot(vals, lk.lenc.mean, type="o", xlab=par.name, ylab="Mean delta NLL disc len comps", lwd = 2.0, cex.main = 2.0, cex.lab =2.0, cex.axis = 2.0, ylim = c(min(lk.lenc.5), max(lk.lenc.95)))
lines(vals, lk.lenc.5, type = "c")
lines(vals, lk.lenc.95, type = "c")
# dev.off()

# RMSE stuff
# r0 <- r0*exp(bam$rec_sigma)
# root mean square error function
# rmse.single <- function(nobs, true, pred){
#   sqdiff <- c()
#   for(i in 1:nobs){
#     sqdiff[i] <- (pred[i]-true)^2
#   }
#   rmse <- sqrt(sum(sqdiff)/nobs)
#   return(rmse)
# }
# 
# rmse.annual <- function(nobs, true, pred){
#   sqdiff <- c()
#   for(i in 1:nobs){
#     sqdiff[i] <- (pred[i]-true[i])^2
#   }
#   rmse <- sqrt(sum(sqdiff)/nobs)
#   return(rmse)
# }
# # M.mult = 1.5
# rmse.single(70, R0, r0)
# rmse.single(70, M.mult, mscale)
# rmse.annual(70, true.f40, pred.f40)
# rmse.annual(70, true.ssbf40, pred.ssbf40)
# 
# Median relative error
# RE.R0 <- (median(r0) - R0) / R0
# RE.Mscale <- (median(mscale) - M.mult) / M.mult
# RE.Mb <- (median(mb) - M.b) / M.b 
# RE.F40 <- (median(pred.f40) - median(true.f40)) / median(true.f40)
# RE.SSBF40 <- (median(pred.ssbf40) - median(true.ssbf40)) / median(true.ssbf40)
# # RE.h <- (median(steep)-h)/h
# 
# RE <- c(RE.R0, RE.Mscale, RE.Mb, RE.F40, RE.SSBF40)
# RE
# # RE*100
# 
# 
# 
# RE.names <-c ("R0", "Mscale", "Mb", "F40", "SSBF40")
# windows(height = 5, width = 6, record = T)
# par(mfrow=c(1,1))
# plot(1:length(RE), RE,
#      type = "p", pch = 15, col = "blue", ylim = c(-max(abs(RE)), max(abs(RE))), xaxt = "n",
#      ylab = "Relative Error", xlab = ""
# )
# axis(1, at = 1:length(RE), labels = RE.names)
# abline(h = 0)
# # 
# years <- 70
# windows(10,10)
# med.true.N <- apply(true.N, 2, median)
# med.pred.N <- apply(pred.N, 2, median)
# plot(1:years, med.true.N / 1000,
#      type = "o", lwd = 2, ylab = "Median Abundance (1000 fish)",
#      xlab = "", panel.first = grid(lty = 1), col = "blue",
#      ylim = c(0, max(c(med.true.N, med.pred.N[1:years])) / 1000)
# )
# points(1:years, med.pred.N[1:years] / 1000, pch = 3, col = "black", cex = 1.2)
# legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))
# text(65, 0.98 * max(med.pred.N)/1000, "(B)", cex = 1.1)
# # # 

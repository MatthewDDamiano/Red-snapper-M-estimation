##################################################################
# Run multiple red snapper simulations and summarize relative error and LL profiles
# 11/12/2024
# Last update: 9/24/2026 MDD
##################################################################
graphics.off()
rm(list = ls(all = TRUE))
library(here())

dir <- here()
setwd(dir)
n.iter = 3 # number of iterations
iter = 1 # starting iteration
n.vals = 13 # set to the number of values across which the model is profiled
n.vals.b = 7 # 11 if it's steepness
years = 70 # when measuring time series quantities over multiple iterations

# storage for relative error
mscale <- array(NA,dim=c(n.iter))
# mb <- array(NA,dim=c(n.iter))
r0 <- array(NA,dim=c(n.iter))
steep <- array(NA,dim=c(n.iter))
max.grad <- array(NA,dim=c(n.iter))
true.f40 <- array(NA,dim=c(n.iter))
pred.f40 <- array(NA,dim=c(n.iter))
true.ssbf40 <- array(NA, dim=c(n.iter))
pred.ssbf40 <- array(NA, dim=c(n.iter))
true.N <- array(NA, dim=c(n.iter,years))
pred.N <- array(NA, dim=c(n.iter,years+1))

# storage for mscale LL profiles
lk.dat.mat <- array(NA, dim=c(n.iter,n.vals))
lk.tot.mat <- array(NA, dim=c(n.iter,n.vals))
lk.sr.mat <- array(NA, dim=c(n.iter,n.vals))
lk.cpue.mat <- array(NA, dim=c(n.iter,n.vals))
lk.sagec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.fl1agec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.fl2agec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.allagec.mat <- array(NA, dim=c(n.iter,n.vals))
lk.lenc.mat <- array(NA, dim=c(n.iter,n.vals))
# for other par LL profiles
lk.lenc.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.dat.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.tot.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.sr.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.cpue.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.sagec.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.fl1agec.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.fl2agec.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.allagec.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter))
lk.lenc.arr.b <- array(NA, dim=c(n.vals, n.vals.b, n.iter)) # save total NLL for bivariate profile (to start, can add more later)
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
  # mb[iter] <- bam$parms$M.b
  r0[iter] <- bam$parms$R0
  max.grad[iter] <- bam$like[28]
  true.f40[iter] <- sim$spr$F40
  pred.f40[iter] <- bam$parms$Fproxy
  true.ssbf40[iter] <- sim$spr$SSB.Fproxy
  pred.ssbf40[iter] <- bam$parms$SSB.Fproxy
  true.N[iter,] <- sim$abundance
  pred.N[iter,] <- bam$t.series$N

  # rename files for profiling # need distinctions between inner and outer loop for bivariate profile
  datname <- 'BAM-Sim.dat'
  rdatname <- 'BAM-Sim.rdat'
  newrdatname <- paste0("Run_", iter, "-set.rdat")
  newdatname <- paste0("Run_",iter,"-set.dat")
  file.rename(rdatname, newrdatname)
  file.rename(datname, newdatname)
  #

  # run profile code
  setwd(dir)
  source("01run_bivar_profile_iter.r") # remember that set.seed is used here for profiling
  # # record total NLL for mb
  lk.dat.arr.b[,,iter] <- ll.dat.b
  lk.tot.arr.b[,,iter] <- ll.tot.b
  lk.sr.arr.b[,,iter] <- ll.sr.b
  lk.cpue.arr.b[,,iter] <- ll.cpue.b
  lk.sagec.arr.b[,,iter] <- ll.sagec.b
  lk.fl1agec.arr.b[,,iter] <- ll.fl1agec.b
  lk.fl2agec.arr.b[,,iter] <- ll.fl2agec.b
  lk.allagec.arr.b[,,iter] <- ll.allagec.b
  lk.lenc.arr.b[,,iter] <- ll.fl3lenc.b
  
  source("01compare_runs_iter.r")
  # record LLs for mscale
  lk.dat.mat[iter,] <- lk.dat
  lk.tot.mat[iter,] <- lk.tot
  lk.sr.mat[iter,] <- lk.SR
  lk.cpue.mat[iter,] <- lk.cpue
  lk.sagec.mat[iter,] <- lk.age.survey1
  lk.fl1agec.mat[iter,] <- lk.age.fleet1
  lk.fl2agec.mat[iter,] <- lk.age.fleet2
  lk.allagec.mat[iter,] <- lk.age
  lk.lenc.mat[iter,] <- lk.len
  #

  # update iteration for the loop
  set.seed(NULL) # unset the seed from LL profiling
  iter = iter+1
}

# # save outputs - add steepness back in when using BH
save(n.vals,vals,
     lk.dat.mat,lk.tot.mat,lk.sr.mat,lk.cpue.mat,lk.sagec.mat,lk.fl1agec.mat,lk.fl2agec.mat,
     lk.allagec.mat,lk.lenc.mat,
     mscale, r0, max.grad, true.f40, pred.f40, true.ssbf40, pred.ssbf40,
     true.N, pred.N, lk.dat.arr.b,lk.tot.arr.b,lk.sr.arr.b,lk.cpue.arr.b,lk.sagec.arr.b,
     lk.fl1agec.arr.b,lk.fl2agec.arr.b,lk.allagec.arr.b,lk.lenc.arr.b,
     file = "M_R0_bivar_3iter_Results.Rdata")

# load or re-load as needed
rm(list=ls())
setwd("/Results")

load("S3_bivar_100iter_Results.Rdata")
# load("S9_bivar_100iter_Results.Rdata")
# load("M_R0_bivar_100iter_Results.Rdata")
# plot(steep,mscale)
# plot(mscale,mb)

# plot(mscale, r0, xlab = "M_scale", ylab = "R0", lwd = 1.5, cex.main = 1.5, cex.lab =1.5, cex.axis = 1.5)

n.iter = 100
n.vals = 13
# profile.vals <- seq(-3, 3, by = 0.5)
# vals <- 2 * exp(profile.vals) / (1 + exp(profile.vals))
# profile.vals <- seq(-2.0,-0.1, by = 0.2) # and these for M.b
# # profile.vals <- seq(0.5,1.0, by = 0.05) # for steepness
# vals <- profile.vals2
par.name <- "mscale"
M.b <- -1.0
M.mult <- 1.0
R0 <- 383000
h <- 0.75
# 
# View(lk.tot.mat)
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

lk.dat.mean <-c()
lk.tot.mean <-c()
lk.sr.mean <-c()
lk.cpue.mean <-c()
lk.sagec.mean <-c()
lk.fl1agec.mean <-c()
lk.fl2agec.mean <-c()
lk.allagec.mean <-c()
lk.lenc.mean <-c()

for(i in 1:length(vals)){
  lk.dat.mean[i] <- lk.dat.mat[i]-min(lk.dat.mat) # I think this will work...
  lk.tot.mean[i] <- lk.tot.mat[i]-min(lk.tot.mat)
  lk.sr.mean[i] <- lk.sr.mat[i]-min(lk.sr.mat)
  lk.cpue.mean[i] <- lk.cpue.mat[i]-min(lk.cpue.mat)
  lk.sagec.mean[i] <- lk.sagec.mat[i]-min(lk.sagec.mat)
  lk.fl1agec.mean[i] <- lk.fl1agec.mat[i]-min(lk.fl1agec.mat)
  lk.fl2agec.mean[i] <- lk.fl2agec.mat[i]-min(lk.fl2agec.mat)
  lk.allagec.mean[i] <- lk.allagec.mat[i]-min(lk.allagec.mat)
  lk.lenc.mean[i] <- lk.lenc.mat[i]-min(lk.lenc.mat)
}


# c(44,45,74) # remove for S7a
# c(4,49,51, 63) # remove for 8
# c(25,45) # remove for S9 mscale

# in case wonky runs need to be removed
# lk.tot.mat <- lk.tot.mat[-c(49),]
# lk.dat.mat <- lk.dat.mat[-c(49),]
# lk.sr.mat <- lk.sr.mat[-c(49),]
# lk.cpue.mat <- lk.cpue.mat[-c(49),]
# lk.sagec.mat <- lk.sagec.mat[-c(49),]
# lk.fl1agec.mat <- lk.fl1agec.mat[-c(49),]
# lk.fl2agec.mat <- lk.fl1agec.mat[-c(49),]
# lk.allagec.mat <- lk.allagec.mat[-c(49),]
# lk.lenc.mat <- lk.lenc.mat[-c(49),]

# # get 5th and 95th percentiles for plotting uncertainty
# lk.dat.5 <- c()
# lk.dat.95 <- c()
# lk.tot.5 <- c()
# lk.tot.95 <- c()
# lk.sr.5 <- c()
# lk.sr.95 <- c()
# lk.cpue.5 <- c()
# lk.cpue.95 <- c()
# lk.sagec.5 <- c()
# lk.sagec.95 <- c()
# lk.fl1agec.5 <- c()
# lk.fl1agec.95 <- c()
# lk.fl2agec.5 <- c()
# lk.fl2agec.95 <- c()
# lk.allagec.5 <- c()
# lk.allagec.95 <- c()
# lk.lenc.5 <- c()
# lk.lenc.95 <- c()
# 
# for(i in 1:n.vals){
#   lk.dat.5[i] <- quantile(lk.dat.mat[,i], probs = c(0.05))
#   lk.dat.95[i] <- quantile(lk.dat.mat[,i], probs = c(0.95))
#   lk.tot.5[i] <- quantile(lk.tot.mat[,i], probs = c(0.05))
#   lk.tot.95[i] <- quantile(lk.tot.mat[,i], probs = c(0.95))
#   lk.sr.5[i] <- quantile(lk.sr.mat[,i], probs = c(0.05))
#   lk.sr.95[i] <- quantile(lk.sr.mat[,i], probs = c(0.95))
#   lk.cpue.5[i] <- quantile(lk.cpue.mat[,i], probs = c(0.05))
#   lk.cpue.95[i] <- quantile(lk.cpue.mat[,i], probs = c(0.95))
#   lk.sagec.5[i] <- quantile(lk.sagec.mat[,i], probs = c(0.05))
#   lk.sagec.95[i] <- quantile(lk.sagec.mat[,i], probs = c(0.95))
#   lk.fl1agec.5[i] <- quantile(lk.fl1agec.mat[,i], probs = c(0.05))
#   lk.fl1agec.95[i] <- quantile(lk.fl1agec.mat[,i], probs = c(0.95))
#   lk.fl2agec.5[i] <- quantile(lk.fl2agec.mat[,i], probs = c(0.05))
#   lk.fl2agec.95[i] <- quantile(lk.fl2agec.mat[,i], probs = c(0.95))
#   lk.allagec.5[i] <- quantile(lk.allagec.mat[,i], probs = c(0.05))
#   lk.allagec.95[i] <- quantile(lk.allagec.mat[,i], probs = c(0.95))
#   lk.lenc.5[i] <- quantile(lk.lenc.mat[,i], probs = c(0.05))
#   lk.lenc.95[i] <- quantile(lk.lenc.mat[,i], probs = c(0.95))
# }
# 
# # take mean of centered LL profiles
# lk.dat.mean <- apply(lk.dat.mat,2,mean)
# lk.tot.mean <- apply(lk.tot.mat,2,mean)
# lk.sr.mean <- apply(lk.sr.mat,2,mean)
# lk.cpue.mean <- apply(lk.cpue.mat,2,mean)
# lk.sagec.mean <- apply(lk.sagec.mat,2,mean)
# lk.fl1agec.mean <- apply(lk.fl1agec.mat,2,mean)
# lk.fl2agec.mean <- apply(lk.fl2agec.mat,2,mean)
# lk.allagec.mean <- apply(lk.allagec.mat,2,mean)
# lk.lenc.mean <- apply(lk.lenc.mat,2,mean)

# lk.tot.mean.mat <- matrix(NA, nrow = 13, ncol = 13)
# for(i in 1:13){
#   lk.tot.mean.mat[i,] <- lk.tot.mean
# }
wonky <- lk.tot.arr.b[which(lk.tot.arr.b>100000)] # record number of wonky runs

# get total NLL averaged over 100 iterations to produce 10x13 matrix for persp plot
# for steepness/mscale bivar LL profiles, need to replace out of bounds runs with avg
# to get shape of profile

profile.vals <- seq(-3, 3, by = 0.5) # These values are Mscale in logit space
profile.vals1 <- 2 * exp(profile.vals) / (1 + exp(profile.vals))
profile.vals2 <- seq(-2.0,-0.1, by = 0.2) # and these for M.b
# profile.vals2 <- seq(10.0, 16.0, by = 1.0) # or R0
# profile.vals2 <- seq(0.5,1.0, by = 0.05) # use these for profiling steepness


for(k in 1:length(profile.vals1)){
for(j in 1:length(profile.vals2)){
for(i in 1:n.iter){
  if(lk.tot.arr.b[k,j,i]>100000){
    lk.tot.arr.b[k,j,i] <- 75000
      }
    }
  }
}
delta.lk.tot.arr <- array(NA, dim=c(length(profile.vals1),length(profile.vals2),n.iter))
for(i in 1:n.iter){
      delta.lk.tot.arr[,,i] <- lk.tot.arr.b[,,i]-min(lk.tot.arr.b[,,i])
}

lk.tot.arr.mean <- apply(delta.lk.tot.arr,c(1,2),mean)

# bivariate prpfile plot
# tiff(filename="M_R0_bivar.tiff", width=12, height=10, units="in", compression="none", res=1000)
# par(mgp=c(3,2,0),las = 1, mar = c(5, 5, 0.5, 0.5), cex = 1.3)
# # windows(10,10)
# persp(profile.vals1, profile.vals2, lk.tot.arr.mean, xlab = "Mscale", ylab = "R0", zlab = "Mean delta NLL", 
#       axes = TRUE, ticktype = "detailed", theta = 140, phi = 10)
# dev.off()
# contour(profile.vals1, profile.vals2, lk.tot.arr.mean)
#alt approach using plotly
# rename profile.vals 1and 2 and lk.tot.arr.mean

library(htmlwidgets)
library(webshot)
M_scale <- profile.vals1
M_b <- profile.vals2
total_delta_mean_NLL <- lk.tot.arr.mean
library(plotly)
plot_ly(x = ~ M_b, y = ~ M_scale, z = ~ total_delta_mean_NLL) %>% 
  add_surface()
p <- plot_ly(x = ~ M_b, y = ~ M_scale, z = ~ total_delta_mean_NLL) %>% 
  add_surface()
saveWidget(p, "temp.html")
# webshot("temp.html", "Mscale_R0_bivar.png")

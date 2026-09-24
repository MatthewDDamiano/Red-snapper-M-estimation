# some code for evaluating multiple runs

# Last updated 4-2025

# graphics.off()
# rm(list = ls(all = TRUE))

# setwd("C:/Users/matt.damiano/Desktop/Red Snapper project/Iterations/Results/full model full data/runs19")

ptype <- "pdf"

folder <- "runs"

# par.name <- "M b"
# profile.vals <- seq(-2.0,-0.1, by = 0.2)
# vals <- profile.vals

par.name <- "R0"
profile.vals <- seq(10.0,16.0,by=1.0)
vals <- profile.vals


# dget(paste(newname, "-set.rdat", sep = ""))


spp1 <- dget(paste(OPdir, "/", folder, iter, "/BAM-Sim-", as.character(iboot1),"-","1_b.rdat", sep = ""))
spp2 <- dget(paste(OPdir, "/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","2_b.rdat", sep = ""))
spp3 <- dget(paste(OPdir, "/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","3_b.rdat", sep = ""))
spp4 <- dget(paste(OPdir, "/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","4_b.rdat", sep = ""))
spp5 <- dget(paste(OPdir, "/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","5_b.rdat", sep = ""))
spp6 <- dget(paste(OPdir, "/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","6_b.rdat", sep = ""))
spp7 <- dget(paste(OPdir, "/", folder, iter, "/BAM-Sim-", as.character(iboot1),"-","7_b.rdat", sep = ""))
# spp8 <- dget(paste(OPdir,"/", folder, iter, "/BAM-Sim-", as.character(iboot1),"-","8_b.rdat", sep = ""))
# spp9 <- dget(paste(OPdir,"/", folder, iter, "/BAM-Sim-", as.character(iboot1),"-","9_b.rdat", sep = ""))
# spp10 <- dget(paste(OPdir,"/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","10_b.rdat", sep = ""))
# spp11 <- dget(paste(OPdir,"/",folder, iter, "/BAM-Sim-", as.character(iboot1),"-","11_b.rdat", sep = ""))
# spp12 <- dget(paste(folder, iter, "/BAM-Sim-12.rdat", sep = ""))
# spp13 <- dget(paste(folder, iter, "/BAM-Sim-13.rdat", sep = ""))
# spp14=dget(paste(folder,"/BAM-Sim-14.rdat",sep=""))
# spp15=dget(paste(folder,"/BAM-Sim-15.rdat",sep=""))
# spp16=dget(paste(folder,"/BAM-Sim-16.rdat",sep=""))
# spp17=dget(paste(folder,"/BAM-Sim-17.rdat",sep=""))
# spp18=dget(paste(folder,"/BAM-Sim-18.rdat",sep=""))
# spp19=dget(paste(folder,"/BAM-Sim-19.rdat",sep=""))
# spp20=dget(paste(folder,"/BAM-Sim-20.rdat",sep=""))
# spp21=dget(paste(folder,"/BAM-Sim-21.rdat",sep=""))
# spp22=dget(paste(folder,"/BAM-Sim-22.rdat",sep=""))
# spp23=dget(paste(folder,"/BAM-Sim-23.rdat",sep=""))
# spp24=dget(paste(folder,"/BAM-Sim-24.rdat",sep=""))
# spp25=dget(paste(folder,"/BAM-Sim-25.rdat",sep=""))
# spp26=dget(paste(folder,"/BAM-Sim-26.rdat",sep=""))
# spp27=dget(paste(folder,"/BAM-Sim-27.rdat",sep=""))
# spp28=dget(paste(folder,"/BAM-Sim-28.rdat",sep=""))
# spp29=dget(paste(folder,"/BAM-Sim-29.rdat",sep=""))
# spp30=dget(paste(folder,"/BAM-Sim-30.rdat",sep=""))
# spp31=dget(paste(folder,"/BAM-Sim-31.rdat",sep=""))
# spp32=dget(paste(folder,"/BAM-Sim-32.rdat",sep=""))
# spp33=dget(paste(folder,"/BAM-Sim-33.rdat",sep=""))
# spp34=dget(paste(folder,"/BAM-Sim-34.rdat",sep=""))
# spp35=dget(paste(folder,"/BAM-Sim-35.rdat",sep=""))
# spp36=dget(paste(folder,"/BAM-Sim-36.rdat",sep=""))
# spp37=dget(paste(folder,"/BAM-Sim-37.rdat",sep=""))
# spp38=dget(paste(folder,"/BAM-Sim-38.rdat",sep=""))
# spp39=dget(paste(folder,"/BAM-Sim-39.rdat",sep=""))
# spp40=dget(paste(folder,"/BAM-Sim-40.rdat",sep=""))
# spp41=dget(paste(folder,"/BAM-Sim-41.rdat",sep=""))


runnames <- as.character(vals)
nruns <- length(runnames) # set number of runs for formatting output

.SavedPlots <- NULL


description <- runnames
###### Obj fcn values
out2 <- matrix(data = NA, nrow = nruns, ncol = 11) # ncol is number of ll components 
label <- rep("", times = nrow(out2))


# Runs
for (i in 1:nruns) {
  fname <- paste("spp", i, sep = "")
  x <- get(fname)
  attach(x)
  out2[i, ] <- c(
    like["lk.unwgt.data"], like["lk.total"], like["lk.SRfit"],
    like["lk.U.survey1"],
    like['lk.L.fleet1'], like['lk.L.fleet2'], like['lk.D.fleet3'],
    like['lk.lenc.fleet3'],
    like["lk.agec.survey1"], like["lk.agec.fleet1"], like["lk.agec.fleet2"]
  )

  detach(x)

  label[i] <- paste("Run", i, sep = "")
}

out2 <- round(out2, 3)
tab <- cbind(label, description, out2)
col.names <- c(
  "label", par.name, "nLL(data)", "nLL(penalized)", "nLL(SR)",
  "U.survey1",
  "L.fleet1","L.fleet2","D.fleet3",
  "len.fleet3",
  "age.survey1", "age.fleet1", "age.fleet2"
)
colnames(tab) <- col.names

# write.csv(tab, file = paste("01Table-Likeprof-lkhd-", par.name, ".csv", sep = ""), quote = F)

# windows(width = 9, height = 8, record = T)

lk.dat.b <- out2[, 1]
lk.tot.b <- out2[, 2]
lk.SR.b <- out2[, 3]
lk.cpue.b <- out2[, 4] # rowSums(out2[,4:6])
# lk.cpue.fleet1 <- out2[, 5]
lk.len.b <- (out2[, 8])
lk.age.b <- rowSums(out2[, 9:11])
lk.age.survey1.b <- out2[, 9]
lk.age.fleet1.b <- out2[, 10]
lk.age.fleet2.b <- out2[, 11]
# 

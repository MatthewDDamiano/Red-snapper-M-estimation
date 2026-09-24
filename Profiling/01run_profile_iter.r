##################################################################
# Run lkhd profiles on red snapper simulation model
# 6/18/11
# Last update: 4/2/2025 MDD
##################################################################
# rm(list = ls(all = TRUE))
set.seed(4886)

library(foreach)
library(doParallel)
library(parallel)

totCores <- detectCores()
numCores <- 8
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Set root level working directory
OPdir <- getwd() # the operating directory
setwd(OPdir)

# admb file to run, without tpl extension
filename <- "BAM-Sim"
newname <- paste0("Run_",iter) # updated here and throughout to use iteration-specific files

# switches applied during admb execution, list with space in between
admb.switch <- "-est -nox"

# folder name for bootstrap results
boot.folder <- paste0("runs",iter) #updated to save profile runs by iteration of the model

# create ADMB input and output files, to be filled in later
# newname.dat <- paste(newname, ".dat", sep = "") # updated each iteration
# newname.rdat <- paste(newname, ".rdat", sep = "") # updated each iteration
# newname.fix.rdat <- paste(newname, ".rdat", sep = "") # same each iteration

# shell(paste("del ",filename,".exe", sep=""))

dat <- scan(file = paste(newname, "-set.dat", sep = ""), what = "character", sep = "&", blank.lines.skip = T, na.strings = "$", comment.char = "#")
write(dat, file = "countlines.dat") # countlines.dat used only for development, does not include blank lines or comments
spp <- dget(paste(newname, "-set.rdat", sep = ""))

##################################################################################
### Run    ####################################################

dir.create(boot.folder, showWarnings = FALSE)
# shell(paste("admb ",filename, sep="")) #This creates the executable using ADMB


# profile.vals1 <- seq(-3, 3, by = 0.5) # These values are Mscale in logit space
profile.vals <- seq(10.0, 16.0, by = 1.0)# for R0
# profile.vals2 <- seq(-2.0,0.5, by = 0.2) # and these for M.b
# profile.vals <- seq(0.5,1.0, by = 0.05) # use these for profiling steepness
nboot <- length(profile.vals)

dat.tmp <- dat

bamexe <- paste(filename, ".exe", sep = "")
bamsource <- paste(OPdir, "/", bamexe, sep = "")
bootout <- paste(OPdir, "/", boot.folder, sep = "")

# for bivariate profile, need to run at each combination of the two vectors of parameter values 
# so you need a nested loop: one that says, for each value of Mscale, profile over b or h, then move to the next value of Mscale
foreach(iboot = 1:nboot) %dopar% { # do in parallel
  # for (iboot in 1:nboot) {
  process.dir <- paste(OPdir, "/", iboot, sep = "") # created for each iteration
  dir.create(process.dir)
  setwd(process.dir)
  ###### Write new admb input file for each bootstrap replicate
  newname.dat <- paste(newname, "-", as.character(iboot), ".dat", sep = "")
  newname.rdat <- paste(filename, "-", as.character(iboot), ".rdat", sep = "") 

  ###############################################################
  ###### Create data set ###############################
  # vec.tmp <- c(profile.vals1[iboot], -10, 10, -4, 0, -0.25, 1) # for Mscale
  # vec.tmp <- c(profile.vals2[iboot], -2.0, -0.1, -4, M.b, -0.25, 1) # for Mb
  # vec.tmp <- c(profile.vals[iboot], 0.21, 1.0, -3, 0.999, -0.25, 1) # for steepness
  vec.tmp <- c(profile.vals[iboot], 10.0, 16.0, -1, 12.86, -0.25, 1) # for R0
  
  # dat.tmp[164] <- as.character(paste(vec.tmp, collapse = "\t")) # use for full model restricted data styr1975 (S4)
  # dat.tmp[345] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P1 or P3; NOTE, this line number found manually in countlines.dat, 337 for mscale, 338 for mb
  # dat.tmp[346] <- as.character(paste(vec.tmp, collapse = "\t")) # for m.b profiling only
  # dat.tmp[347] <- as.character(paste(vec.tmp, collapse = "\t")) # for profiling steepness
  dat.tmp[348] <- as.character(paste(vec.tmp, collapse = "\t")) # for profiling R0
  # dat.tmp[285] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2a
  # dat.tmp[210] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2b
  # dat.tmp[170] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2c
  # dat.tmp[186] <- as.character(paste(vec.tmp, collapse = "\t")) # use for S10 (P2c with 4 more years)
  # dat.tmp[255] <- as.character(paste(vec.tmp, collapse = "\t")) # use for B1 
  ###############################################################
  ####### Run admb. -ind switch changes the name of the data input file each bootstrap iteration
  write(file = newname.dat, dat.tmp)
  bamboot <- paste(process.dir, "/", basename(filename), "-", as.character(iboot), ".exe", sep = "")
  file.copy(bamsource, bamboot, overwrite = TRUE)
  bamrun <- paste(basename(filename), "-", as.character(iboot), ".exe", sep = "")
  run.command <- paste(bamrun, admb.switch, "-ind", newname.dat, sep = " ")
  shell(run.command)

  ########## Copy data files to boot.folder
  file.copy(from = newname.dat, to = bootout, overwrite = TRUE) 
  file.copy(from = newname.rdat, to = bootout, overwrite = TRUE)  

  ####### Remove individual processing folders
  setwd(OPdir)
  unlink(process.dir, recursive = T)
} # end nboot
stopCluster(cl)

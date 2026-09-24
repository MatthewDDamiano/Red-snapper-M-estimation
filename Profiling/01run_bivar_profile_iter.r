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
# OPdir <- getwd() # the operating directory
OPdir <- "C:/Users/mddamiano/Desktop/NOAA projects/Red Snapper/RS M sim_Damiano mod_LorenzenM/Profiling/Results/Model S3 bivar"
setwd(OPdir)

# iter = 1
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
spp <- dget(paste(newname, "-set.rdat", sep = "")) # the -set files are the permanent ones that are copied and renamed as .dat and .rdats

##################################################################################
### Run    ####################################################

dir.create(boot.folder, showWarnings = FALSE)
# shell(paste("admb ",filename, sep="")) #This creates the executable using ADMB


profile.vals1 <- seq(-3, 3, by = 0.5) # These values are Mscale in logit space
# profile.vals2 <- seq(-2.0,-0.1, by = 0.2) # and these for M.b
profile.vals2 <- seq(10.0,16.0, by = 1.0) # use these for profiling steepness

# create matrices to store LL profiles for inner par
ll.tot.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.dat.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.sr.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.cpue.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.allagec.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.sagec.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.fl1agec.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.fl2agec.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
ll.fl3lenc.b <- array(NA,dim=c(length(profile.vals1),length(profile.vals2)))
  
nboot <- length(profile.vals1)
nboot2 <- length(profile.vals2)

dat.tmp <- dat

bamexe <- paste(filename, ".exe", sep = "")
bamsource <- paste(OPdir, "/", bamexe, sep = "")
bootout <- paste(OPdir, "/", boot.folder, sep = "")
# iboot1 = 1
# for bivariate profile, need to run at each combination of the two vectors of parameter values 
# so you need a nested loop: one that says, for each value of Mscale, profile over b or h, then move to the next value of Mscale

for(iboot1 in 1:nboot){ # do in parallel - now just a for loop so internal loop can run parallel
  # for (iboot in 1:nboot) {
  process.dir <- paste(OPdir, "/", iboot1, sep = "") # created for each iteration
  dir.create(process.dir)
  setwd(process.dir)
  ###### Write new admb input file for each bootstrap replicate
  newname.dat <- paste(newname, "-", as.character(iboot1), ".dat", sep = "")
  newname.rdat <- paste(filename, "-", as.character(iboot1), ".rdat", sep = "") 

  ###############################################################
  ###### Create data set ###############################
  vec.tmp1 <- c(profile.vals1[iboot1], -10, 10, -4, 0, -0.25, 1) # for Mscale
  # vec.tmp2 <- c(profile.vals2[iboot], -2.0, -0.1, -4, M.b, -0.25, 1) # for Mb
  # vec.tmp <- c(profile.vals[iboot], 0.21, 1.0, -3, 0.999, -0.25, 1) # for steepness
  
  # dat.tmp[164] <- as.character(paste(vec.tmp, collapse = "\t")) # use for full model restricted data styr1975 (S4)
  dat.tmp[345] <- as.character(paste(vec.tmp1, collapse = "\t")) # use for P1 or P3; NOTE, this line number found manually in countlines.dat, 337 for mscale, 338 for mb
  # dat.tmp[346] <- as.character(paste(vec.tmp2, collapse = "\t")) # for m.b profiling only
  # dat.tmp[347] <- as.character(paste(vec.tmp, collapse = "\t")) # for profiling steepness
  # dat.tmp[285] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2a
  # dat.tmp[210] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2b
  # dat.tmp[170] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2c
  # dat.tmp[255] <- as.character(paste(vec.tmp, collapse = "\t")) # use for B1 
  ###############################################################
  ####### Run admb. -ind switch changes the name of the data input file each bootstrap iteration
  write(file = newname.dat, dat.tmp)
  bamboot <- paste(process.dir, "/", basename(filename), "-", as.character(iboot1), ".exe", sep = "")
  file.copy(bamsource, bamboot, overwrite = TRUE)
  bamrun <- paste(basename(filename), "-", as.character(iboot1), ".exe", sep = "")
  run.command <- paste(bamrun, admb.switch, "-ind", newname.dat, sep = " ")
  shell(run.command)

  ########## Copy data files to boot.folder
  file.copy(from = newname.dat, to = bootout, overwrite = TRUE) # these go to the "runs[iter]" folder, which we keep
  file.copy(from = newname.rdat, to = bootout, overwrite = TRUE)  
  
  #inner loop
  
  #set directory to the temporary boot folder run profile of inner par
  # Idir <- getwd() # the operating directory
  # setwd(Idir)
  # set profile vals for inner par
  # profile.vals2 <- seq(-2.0,-0.1, by = 0.2) # and these for M.b
  # set length for loop
  # you can keep the old dat <- dat.tmp since you're using the same dat file

  foreach(iboot2 = 1:nboot2) %dopar% { # do in parallel
    # for (iboot in 1:nboot) {
    process.dir2 <- paste(OPdir, "/", iboot1, "/", iboot2, sep = "") # created for each iteration within the inner par's temp boot folder
    dir.create(process.dir2)
    setwd(process.dir2)
    ###### Write new admb input file for each bootstrap replicate
    newname.b.dat <- paste(newname, "-", as.character(iboot1),"-",as.character(iboot2), "_b", ".dat", sep = "")
    newname.b.rdat <- paste(filename, "-", as.character(iboot1),"-",as.character(iboot2), "_b", ".rdat", sep = "") 
    
    ###############################################################
    ###### Create data set ###############################
    # vec.tmp1 <- c(profile.vals1[iboot], -10, 10, -4, 0, -0.25, 1) # for Mscale
    # M.b = -1.0
    # vec.tmp2 <- c(profile.vals2[iboot2], -2.0, -0.1, -4, M.b, -0.25, 1) # for Mb
    vec.tmp2 <- c(profile.vals2[iboot2], 10.0, 16.0, -1, 12.86, -0.25, 1) # for R0
    
    # dat.tmp[164] <- as.character(paste(vec.tmp, collapse = "\t")) # use for full model restricted data styr1975 (S4)
    # dat.tmp[345] <- as.character(paste(vec.tmp1, collapse = "\t")) # use for P1 or P3; NOTE, this line number found manually in countlines.dat, 337 for mscale, 338 for mb
    # dat.tmp[346] <- as.character(paste(vec.tmp2, collapse = "\t")) # for m.b profiling only
    dat.tmp[348] <- as.character(paste(vec.tmp2, collapse = "\t")) # for profiling steepness
    # dat.tmp[285] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2a
    # dat.tmp[210] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2b
    # dat.tmp[170] <- as.character(paste(vec.tmp, collapse = "\t")) # use for P2c
    # dat.tmp[255] <- as.character(paste(vec.tmp, collapse = "\t")) # use for B1 
    ###############################################################
    ####### Run admb. -ind switch changes the name of the data input file each bootstrap iteration
    write(file = newname.b.dat, dat.tmp)
    bamboot <- paste(process.dir2, "/", basename(filename), "-", as.character(iboot1),"-",as.character(iboot2), "_b", ".exe", sep = "")
    file.copy(bamsource, bamboot, overwrite = TRUE)
    bamrun <- paste(basename(filename), "-", as.character(iboot1),"-",as.character(iboot2), "_b", ".exe", sep = "")
    run.command <- paste(bamrun, admb.switch, "-ind", newname.b.dat, sep = " ")
    shell(run.command)
    
    ########## Copy data files to boot.folder
    file.copy(from = newname.b.dat, to = bootout, overwrite = TRUE) 
    file.copy(from = newname.b.rdat, to = bootout, overwrite = TRUE)  
    
    ####### Remove individual processing folders
    setwd(process.dir)
    unlink(process.dir2, recursive = T)
  } # end nboot/inner loop
 
  ####### Record profile summaries for inner loop parameter
  source("C:/Users/mddamiano/Desktop/NOAA projects/Red Snapper/RS M sim_Damiano mod_LorenzenM/Profiling/01compare_runs_iter_mslope.r")
  ll.tot.b[iboot1,] <- lk.tot.b
  ll.dat.b[iboot1,] <- lk.dat.b
  ll.sr.b[iboot1,] <- lk.SR.b
  ll.cpue.b[iboot1,] <- lk.cpue.b
  ll.allagec.b[iboot1,] <- lk.age.b
  ll.sagec.b[iboot1,] <- lk.age.survey1.b
  ll.fl1agec.b[iboot1,] <- lk.age.fleet1.b
  ll.fl2agec.b[iboot1,] <- lk.age.fleet2.b
  ll.fl3lenc.b[iboot1,] <- lk.len.b
  
  ####### Remove individual processing folders
  setwd(OPdir)
  unlink(process.dir, recursive = T)
} # end outer loop
stopCluster(cl)

# *********************************************************************************************
# Create Windows batch files that run a command prompt executable with varying input parameters 
# (input parameters are drawn from model posterior distributions) - Henna Fabritius
# *****************************************************************************************************************

library(coda) # read.openbugs (read posterior distributions from OpenBUGS from coda files)
library(VGAM) # cloglog (model-related data transformation)
library(LaplacesDemon) # rbern (Bernoulli distribution)
library(boot) # inv.logit (model-related data transformation)

samplesize <-1000
model_nr   <-5  # selects which of the model alternatives the simulation will use
hostprop   <-0.25
landscape_sample<-sample(seq(from=0,to=99,by=1),samplesize,replace=TRUE) # randomization of landscape variants used in the simulation

# Read in CODA files for parameter distributions
setwd(choose.dir()) # set working folder
CODA_files<-read.openbugs(stem="model_name_", quiet=T)
Model<-rbind(as.matrix(CODA_files[1]),as.matrix(CODA_files[2])) # bind to a matrix
if(TRUE) Model<-Model[sample.int(nrow(Model),samplesize,replace=FALSE),] # take random draws from the joint posterior distribution

# Read in prior distributions for params that have not yet been estimated
setwd(choose.dir())
priors<-read.csv2("prior_distributions.csv",header=T)

# OPTION A: SIMULATE USING RANDOM DRAWS FROM THE POSTERIOR DISTRIBUTION
# ------------------------------------------------------------------------------------------------------------------
# Read in a existing initialization file that will be modified
text = readLines("input_file.dat",-1) # not a table, so needs to be read without interpreting content
for(i in 1:samplesize){
  zeros<-ifelse(i<10,"000",ifelse(i<100,"00",ifelse(i<1000,"0",""))) # leading zeros in file naming
  # only certain lines of the initialization file need to be modified, note double backslash when writing Windows-style folder paths
  text[1]=paste("C:\\directory\\sub-directory\\folder\\landscape_",landscape_sample[i],"\\    simulation input map (without extensions)",sep="")
  text[6]=paste(model_nr,"\t\t    input file description",sep="")
  text[7]=paste(Model[,"param_01"][i],"            input file description for param. 01, copied from the original file",sep="")
  text[9]=paste(Model[,"param_02"][i],"\t\t    input file description for param. 02, copied from the original file",sep="")
  text[14]=paste(100,"            input file description for param. 03, copied from the original fie",sep="") # we fix this param value
  text[15]=paste(hostprop,"             input file description for param. 02, copied from the original file",sep="")
  text[31]=paste(zeros,i-1,"_param_file.log        log file name, index number will be the same as in the initialization file",sep="")
  writeLines(text,paste("input_file_",zeros,i-1,".dat",sep="")) # initialization files are indexed for the analysis phase
}

# OPTION B: SIMULATE USING POSTERIOR MEDIANS
# ------------------------------------------------------------------------------------------------------------------
# Read in a existing initialization file that will be modified
text = readLines("input_file.dat",-1) # not a table, so needs to be read without interpreting content
for(i in 1:samplesize){
  zeros<-ifelse(i<10,"000",ifelse(i<100,"00",ifelse(i<1000,"0",""))) # leading zeros in file naming
  # only certain lines of the initialization file need to be modified, note double backslash when writing Windows-style folder paths
  text[1]=paste("C:\\directory\\sub-directory\\folder\\landscape_",landscape_sample[i],"\\    simulation input map (without extensions)",sep="")
  text[6]=paste(model_nr,"\t\t    input file description",sep="")
  text[7]=paste(median(Model[,"param_01"]),"            input file description for param. 01, copied from the original file",sep="")
  text[9]=paste(median(Model[,"param_02"]),"\t\t    input file description for param. 02, copied from the original file",sep="")
  text[14]=paste(100,"            input file description for param. 03, copied from the original fie",sep="") # we fix this param value
  text[15]=paste(hostprop,"             input file description for param. 02, copied from the original file",sep="")
  text[31]=paste(zeros,i-1,"_param_file.log        log file name, index number will be the same as in the initialization file",sep="")
  writeLines(text,paste("input_file_",zeros,i-1,".dat",sep="")) # initialization files are indexed for the analysis phase
}

# ----------------------------------------------------------------------------------------------------------------------
# CREATE .BAT FILES TO AUTOMATE THE RUNNING OF SAMPLESIZE<=9999 SIMULATIONS
text<-character(length = 0) # initialize a character file
for(i in 1:samplesize){
  zeros<-ifelse(i<10,"000",ifelse(i<100,"00",ifelse(i<1000,"0",""))) # leading zeros in file naming
  text[length(text)+1] = paste("simulation.exe input_file_",zeros,i,".DAT",sep="")
  text[length(text)+1] = paste("rename output.log ",zeros,i,"_output.log",sep="")
  text[length(text)+1] = "timeout 1" # the file system may be slow, safest to wait for a while after file operations
  text[length(text)+1] = paste("mkdir ",zeros,i,sep="")
  text[length(text)+1] = paste("move ",zeros,i,"_coord2000* ",zeros,i,sep="") # move created simulation files to an indexed folder
  text[length(text)+1] = paste("del ",zeros,i,"*",sep="") # delete indexed simulation output files that are not needed
}
writeLines(text,"run_sim_variants.BAT") # save the created batch file and execute in a command line

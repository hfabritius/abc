# ***********************************************************************************************************************************
# Code for choosing population simulations (created with varying extinction-colonization rates) that produce habitat patch occupancies
# similar to (incompletely surveyd multi-year) empirical data. Rejection algorithm modified form Rubin 1984. December 2014 -Henna Fabritius
# ***********************************************************************************************************************************

library(gtools) # invalid function

# SUMMARY STATISTICS OF EMPIRICAL DATA
# ------------------------------------------------------------------
setwd(choose.dir()) # observed field data, n=459 (West Coast)
data<-read.csv2("Habitat_patch_network_&_field_visits.csv") # read in data with four-year (years by columns) survey results
novisit<-is.na(data) # save observation patterns of observed data for the analysis of incomplete simulated data

# Summary statistic 1: average site occupancy
S_A<-c() # The fraction of occupied patches from the pooled data for all analysis years (summary statistic S_A)
for(i in 1:ncol(data)) S_A<-c(S_A,data[!is.na(data[,i]),i]) # pooled abundance data per year for all monitored sites 
S_A<-length(S_A[S_A>0])/length(S_A) # nr. of occupied sites / nr. of sites in the pooled survey data across years

# The fraction of cases in which a site observed occupied was still observed occupied after 1, 2 or 3 years (summary statistics S_B1,S_B2,S_B3)
# The fraction of cases in which a site observed empty was still observed empty after 1, 2 or 3 years (summary statistics S_C1,S_C2,S_C3)
S_BC<-matrix(0,nrow=2,ncol=ncol(data)-1) # initialize matrix for the summary statistics of occupied and empty sites
S_total<-matrix(0,nrow=2,ncol=ncol(data)-1) # initialize matrix for the total number of sequential observation events
for(i in 1:(ncol(data)-1)){
  for(j in (i+1):ncol(data)){
    data_sub<-data[!is.na(data[,i]+data[,j]),c(i,j)]
    if(!invalid(data_sub[data_sub[,1]>0,])) S_total[1,(j-i)]<-S_total[1,(j-i)]+nrow(as.matrix(data_sub[data_sub[,1]>0,])/2) # Patches initially occupied
    if(!invalid(data_sub[data_sub[,1]==0,])) S_total[2,(j-i)]<-S_total[2,(j-i)]+nrow(as.matrix(data_sub[data_sub[,1]==0,])/2) # Patches initially unoccupied
    if(!invalid(data_sub[(data_sub[,1]>0 & data_sub[,2]==0),])) S_BC[1,(j-i)]<-S_BC[1,(j-i)]+nrow(as.matrix(data_sub[(data_sub[,1]>0 & data_sub[,2]>0),])/2) # Occupied patches stayed occupied
    if(!invalid(data_sub[(data_sub[,1]==0 & data_sub[,2]>0),])) S_BC[2,(j-i)]<-S_BC[2,(j-i)]+nrow(as.matrix(data_sub[(data_sub[,1]==0 & data_sub[,2]==0),])/2) # Unoccupied patches stayed unoccupied
}}

S_BC<-S_BC/S_total # calculate fractions of each type (out of all events)
Observed<-c(S_A,S_BC[1,],S_BC[2,])
names(Observed)<-c("S_A","S_B1","S_B2","S_B3","S_C1","S_C2","S_C3")

# SUMMARY STATISTICS OF SIMULATED DATA SETS
# ------------------------------------------------------------------
Table<-data.frame() # initialize a data frame
extRates<-read.csv2("extRates_10000.csv")[,-1] # read in extinction rates used to create simulated data
colRates<-read.csv2("colRates_10000.csv")[,-1] # read in colonization rates used to create simulated data

setwd(choose.dir()) # choose folder for simulated data
DataFiles=list.files(recursive=TRUE,pattern="SimulationResults_Rep") # read a list of indexed simulation file names

# read in simulation results
for (d in 1:length(DataFiles)) {
    
  data<-read.csv2(DataFiles[d])[,-1] # retrieve survey data
  if(ncol(data)<12) for(n in 1:(12-ncol(data))) data<-cbind(data,0)
  data<-data[,9:12] # reduce data to survey results
  data[novisit]<-NA # simulate the same observation pattern as in empirical data
  
  # Turn a specified proportion of sightings into non-sightings, to account for the expected rate of false negatives in the empirical data
  falseneg_rate<-0.1 # adjust the false negative rate according to the species studied
  falsenegatives<-rep(1,length(data[!is.na(data) & data>0])) # create a list of sites with sightings
  falsenegatives[runif(length(falsenegatives),0,1)<falseneg_rate]<-0 # by chance, some of these sites fall below ~prob. of observation by the false negative rate
  data[!is.na(data) & data>0]<-falsenegatives #...and are thus turned negative sightings (NOTE: we skip the assumption that population abundance predicts detectability)
  
  # The fraction of occupied patches from the pooled data for all four years (summary statistic S_A)
  S_A<-c()
  for(i in 1:ncol(data)) S_A<-c(S_A,data[!is.na(data[,i]),i]) # pooled data
  S_A<-length(S_A[S_A>0])/length(S_A) # nr. of occupied patches / nr. of patches in the pooled survey data
  
  # The fraction of cases in which a patch observed occupied was still observed occupied after 1, 2 or 3 years (summary statistics S_B1,S_B2,S_B3)
  # The fraction of cases in which a patch observed empty was still observed empty after 1, 2 or 3 years (summary statistics S_C1,S_C2,S_C3)
  S_BC<-matrix(0,nrow=2,ncol=ncol(data)-1)  # initialize matrix for the summary statistics of occupied and empty sites
  S_total<-matrix(0,nrow=2,ncol=ncol(data)-1) # initialize matrix for the total number of sequential observation events
  for(i in 1:(ncol(data)-1)){
    for(j in (i+1):ncol(data)){
      data_sub<-data[!is.na(data[,i]+data[,j]),c(i,j)]
      if(!invalid(data_sub[data_sub[,1]>0,])) S_total[1,(j-i)]<-S_total[1,(j-i)]+nrow(as.matrix(data_sub[data_sub[,1]>0,])/2) # Patches initially occupied
      if(!invalid(data_sub[data_sub[,1]==0,])) S_total[2,(j-i)]<-S_total[2,(j-i)]+nrow(as.matrix(data_sub[data_sub[,1]==0,])/2) # Patches initially unoccupied
      if(!invalid(data_sub[(data_sub[,1]>0 & data_sub[,2]==0),])) S_BC[1,(j-i)]<-S_BC[1,(j-i)]+nrow(as.matrix(data_sub[(data_sub[,1]>0 & data_sub[,2]>0),])/2) # Occupied patches stayed occupied
      if(!invalid(data_sub[(data_sub[,1]==0 & data_sub[,2]>0),])) S_BC[2,(j-i)]<-S_BC[2,(j-i)]+nrow(as.matrix(data_sub[(data_sub[,1]==0 & data_sub[,2]==0),])/2) # Unoccupied patches stayed unoccupied
    }}
  
  S_BC<-S_BC/S_total # calculate fractions of each type (out of all events)
  Table<-rbind(Table,c(extRates[d],colRates[d],S_A,S_BC[1,],S_BC[2,])) # add results to a table, together with the param. values that produced the simulation
  colnames(Table)<-c("extRate","colRate","S_A","S_B1","S_B2","S_B3","S_C1","S_C2","S_C3")

}
write.csv2(Table,"ResultsTable.csv") # save for future analysis

# REJECTION ALGORITHM (MODIFIED FROM RUBIN 1984)
# ------------------------------------------------------------------
Table<-Table[,1:8] # Summary statistic S_C3 was NA in real data
Table<-Table[complete.cases(Table),] # remove rows with NAs (i.e. simulations where populations have gone extinct)

for(i in 1:nrow(Table)) Table[i,3:8]<-Table[i,3:8]-Observed # calculate difference of summary statistics to empirical data
threshold<-0.3194 # set a random threshold value that can be iterated to accept a desired number of best-matching simulations
S_A<-sort(abs(Table$S_A),decreasing=FALSE)[round(threshold*nrow(Table))] # find a threshold difference of S_A summary statistics that accepts the desired number of simulations
S_B1<-sort(abs(Table$S_B1),decreasing=FALSE)[round(threshold*nrow(Table))] # same for S_B1 (all summary statistics are set to equal weight here)
S_B2<-sort(abs(Table$S_B2),decreasing=FALSE)[round(threshold*nrow(Table))] # and so on
S_B3<-sort(abs(Table$S_B3),decreasing=FALSE)[round(threshold*nrow(Table))]
S_C1<-sort(abs(Table$S_C1),decreasing=FALSE)[round(threshold*nrow(Table))]
S_C2<-sort(abs(Table$S_C2),decreasing=FALSE)[round(threshold*nrow(Table))]

# subset data into simulations that fulfill the threshold criteria for all summary statistics
Table_subset<-Table[(abs(Table$S_A)<=S_A & abs(Table$S_B1)<=S_B1 & abs(Table$S_B2)<=S_B2 & abs(Table$S_B3)<=S_B3 & abs(Table$S_C1)<=S_C1 & abs(Table$S_C2)<=S_C2),]
dim(Table_subset) # check the number of accepted simulations and iterate the threshold accordingly to retrieve the desired number of accepted simulations

# PLOT THE DISTRIBUTIONS OF PARAMETERS THAT PRODUCED ACCEPTED SIMULATIONS
# ------------------------------------------------------------------
par(mfrow=c(1,2))
plot(density(Table_subset$extRate),xlab="Extinction rate",main="") # distribution of extinction rates
plot(density(Table_subset$colRate),xlab="Colonization rate",main="") # distribution of colonization rates
# statistics of the extinction rate distribution
mean(Table_subset$extRate) # mean
median(Table_subset$extRate) # median
sort(Table_subset$extRate)[round(0.025*length(Table_subset$extRate))] # lower 95 % Cr.I.
sort(Table_subset$extRate)[round(0.975*length(Table_subset$extRate))] # upper 95 % Cr.I.
# statistics of the colonization rate distribution
mean(Table_subset$colRate) # mean
median(Table_subset$colRate) # median
sort(Table_subset$colRate)[round(0.025*length(Table_subset$colRate))] # lower 95 % Cr.I.
sort(Table_subset$colRate)[round(0.975*length(Table_subset$colRate))] # upper 95 % Cr.I.

# plot all vs. accepted parameter value combinations
par(mfrow=c(1,1))
plot(Table$extRate,Table$colRate,xlab="Extinction rate",ylab="Colonization rate")
points(Table_subset$extRate,Table_subset$colRate,col="black",pch=16)
legend(x=0.5,y=0.15,legend=c("Tested","Accepted"),pch=21,pt.bg=c("white","black"))
write.csv2(Table_subset,"Posterior.csv")

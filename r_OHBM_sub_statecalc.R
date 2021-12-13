#Usage: 
# Rscript --vanilla <script name> <clean> <order> <pen> <width> <roiname> <subject> <k> <iteration>

#Arguments:
# <clean> Clean or not
# <order> First or second
# <pen> Penalize or not
# <width> Width
# <roiname> ROI set list name
# <subject> Subject ID
# <k> k for k-clustering
# <iteration> Iteration of the 100 loops

#Set code as working directory.

#Loading start time.
start.time <- Sys.time()

#Load packages.
library(tidyverse)
library(Gmedian)
library(reticulate)

#filter,lag,between,first,last,transpose

#Set virtual environment.
use_virtualenv('../envs/FlexEnv')

#Load python packages.
pd <- import("pandas")

#Catch arguments.
args = commandArgs(trailingOnly=T)
clean <- args[1]
order <- args[2]
pen <- args[3]
width <- args[4]
roiname <- args[5]
subject <- args[6]
k <- args[7]
iteration <- args[8]
# clean <- 'clean'
# order <- 'first'
# pen <- 'penalize'
# width <- '139'
# roiname <- 'whole'
# subject <- '101107'
# k <- '12'
# iteration <- '1'
print(paste0('Doing: ',clean,' ',order,' ',pen,' ',width,' ',roiname,' ',subject,
             ' ',k,' ',iteration))

#Set output folder path. If output folder doesn't exist, create it.
outpath <- paste0('../outputs/r_OHBM_stateflex/sub_calcs/ver_kGmedian/',clean,'/',
                  order,'/',pen,'/',width,'/',roiname,'/',subject,'/',k,'/')
dir.create(outpath,recursive=T)

#Set output file names. If they already exist, quit.
clustfile <- paste0(outpath,'subclust_',iteration,'.csv')
medfile <- paste0(outpath,'submedians_',iteration,'.csv')
if (file.exists(clustfile) & file.exists(medfile)) {
  print('Files exist.')
  quit(save='no')
}

#Read in vectorized labels.
veclabels <- read_csv('atlasRvec.csv',col_names=F) %>%
  separate(col='X1',into=c('conn1','conn2'),sep='-')

#Read in the ROI set.
roifile = paste0('roi_',roiname,'.txt')
roiset <- scan(roifile,character(),quote='') 

#For each set of connection indices which contains at least one region that 
#matches the ROI set for extranetwork clustering or contains both regions 
#that match the ROI set for intranetwork clustering, retain the set index.
clust_type <- 'extra'
roi_index <- c()
for (i in c(1:nrow(veclabels))) {
  if (clust_type=='extra') {
    if ((veclabels[i,1] %in% roiset)|(veclabels[i,2] %in% roiset)) {
      roi_index <- c(roi_index,i)
    }
  } else if (clust_type=='intra') {
    if ((veclabels[i,1] %in% roiset)&(veclabels[i,2] %in% roiset)) {
      roi_index <- c(roi_index,i)
    }
  }
}

#Read in the subject dFC matrix. Select the ROI rows and relabel columns and rows.
dFCstorepath <- paste0('../outputs/r_OHBM_slid_dFC/',clean,'/',order,'/',pen,'/',
                       width,'/',subject,'/')
dFCstore <- pd$HDFStore(paste0(dFCstorepath,'r_OHBM_slid_dFC.h5'),'r')
subkey <- paste0('sub_',subject)                      
dFCmat <- dFCstore$select(subkey)[,roi_index]
dFCstore$close()
colnames(dFCmat) <- c(1:length(roi_index))
rownames(dFCmat) <- c(1:nrow(dFCmat))
# dFCmat <- t(h5read(file=paste0(dFCstorepath,'allsubs.hf'),
#                    name=subkey)$block0_values)
# h5closeAll()
print('Read dFCmat.')

#Set ending time for loading and output.
end.time <- Sys.time()
print(paste0('Loading time: ',as.character(end.time-start.time)))

#Define portion of the 100 loops based on the iteration.
nloops = 100
startn = nloops*(as.numeric(iteration)-1) + 1
endn = nloops*(as.numeric(iteration))
loops = c(startn:endn)

#Set best WCS value as the infinity value in the beginning and do the loops.
bestwcs = Inf
for (i in loops) {
  
  #Set starting time for loop.
  start.time <- Sys.time()
  print(paste0('Doing loop: ',as.character(i)))
  
  #Read in the initial points, select k points, convert into a matrix, and relabel.
  kppstorepath <- paste0('../outputs/r_OHBM_stateflex/sub_kpps/',clean,'/',order,'/',
                         pen,'/',width,'/',roiname,'/',subject,'/')
  kppstore <- pd$HDFStore(paste0(kppstorepath,'kppcents.h5'),'r')
  iterkey <- paste0(subkey,'_',as.character(i))                     
  initpts <- as.matrix(kppstore$select(iterkey)[1:as.numeric(k),])
  kppstore$close()
  colnames(initpts) <- c(1:length(roi_index))
  rownames(initpts) <- c(1:as.numeric(k))
  # kppfile <- paste0('../outputs/r_stateflex/sub_kpps/',roiname,'/',width,'/',
  #                   subject,'/kppcents_',as.character(i),'.csv')
  # initpts <- as.matrix(read_csv(kppfile,col_names=F,n_max=as.numeric(k)))
  print('Read initial points.')
  
  #Conduct k-median clustering.
  kmedian <- kGmedian(dFCmat,ncenters=initpts,nstart=1)
  print('K-clustered.')
  
  #Produce the WCS value.
  wcs <- sum(kmedian$withinsrs)
  
  #If it is the first loop or if the current WCS is less than the previous best
  #WCS, replace the best values.
  if ((i == 1)|(wcs < bestwcs)) {
    bestwcs <- wcs
    bestsep_wcs <- kmedian$withinsrs
    bestclusters <- kmedian$cluster
    bestmedians <- kmedian$centers
  }
  
  #Delete all unnecessary variables before the next iteration.
  rm(kppstorepath,kppstore,iterkey,initpts,kmedian,wcs)
  
  #Set ending time for loop and output.
  end.time <- Sys.time()
  print(paste0('Loop done: ',as.character(end.time-start.time)))
}

#Set starting time for saving.
start.time <- Sys.time()

#Package the best WCS and clusters, and the separate WCS and the centroids.
clusterpack <- t(rbind(bestwcs,bestclusters))
medianpack <- cbind(bestsep_wcs,bestmedians)

#Save the cluster package and median centroids.
write_csv(as_tibble(clusterpack),clustfile,col_names=F)
write_csv(as_tibble(medianpack),medfile,col_names=F)

#Set ending time for saving and output.
end.time <- Sys.time()
print(paste0('Saving time: ',as.character(end.time-start.time)))

#test1 <- read_csv(clustfile,col_names=F)
#test2 <- read_csv(medfile,col_names=F)

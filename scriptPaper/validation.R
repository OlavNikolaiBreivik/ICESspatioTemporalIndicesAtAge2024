library(spatioTemporalIndices)
library(spatioTemporalALK)

##############################################################
#Simulation experiment########################################
#Read data
dat_l = readRDS("data/catch_at_length_data_ex_rus.rds")
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")
dat_alk = dat_alk[which(dat_alk$readability==1),]#NB!!!Remove observations with readability issues

#Configurations length part
conf_l = defConf(years = 2011:2020, # years to use, use all years with data by default
                 maxLength = 75,
                 minLength = 20,
                 spatioTemporal =1 ,
                 spatial =1,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,1),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 minDepth=100,maxDepth=450, #Approximately 80% coverage interval
                 applyALK = 1,
                 cutoff =80)

#Define configurations age part
conf_alk = defConf_alk(maxAge = 10,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =1,
                       rwBeta0 = 1)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 20)

# run model
start_time <- Sys.time()
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk, getReportCovariance = FALSE,skip.delta.method = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

set.seed(123)
nSim = 26
nSimSave = 25
sim = simStudy(run,nsim = nSim)
plotSim(sim,ylim = c(-8,2))
source("scriptPaper/utils.R")

abundanceRealizationsTmp = abundanceSim(run,sim)
toRemove = numeric(0)
for(i in 1:length(sim)){
  if(is.na(sim[[i]][1])){
    toRemove = c(toRemove,i)
  }else if(sim[[i]]$opt$convergence !=0){
    toRemove = c(toRemove,i)

  }
}
if(length(toRemove)>0){
  simConverged = sim[-toRemove]
  abundanceRealizations = abundanceRealizationsTmp[-toRemove]
}else{
  simConverged = sim
  abundanceRealizations = abundanceRealizationsTmp
}
simConverged = simConverged[1:nSimSave]
abundanceRealizations = abundanceRealizations[1:nSimSave]
attributes(simConverged)$run = run
save(simConverged,abundanceRealizations, file = "results/simStudy.RData")


png("figures/simExperiment.png")
plotSim(simConverged)
title(main = "Simulation experiment", cex.main = 2)
dev.off()

#Calculate mean relative error
diff = list()
for(i in 1:length(simConverged)){
  rl = as.list(simConverged[[i]]$rep, what = "Est", report = TRUE)
  rlSd = as.list(simConverged[[i]]$rep, "Std", report = TRUE)
  ageLogIndex = rl$logAgeIndex
  ageLogIndexSd = rlSd$logAgeIndex
  diff[[i]] = (ageLogIndex- abundanceRealizations[[i]])/abundanceRealizations[[i]]
}
diffVector = as.vector(do.call(rbind,diff))
diffVector = diffVector[!is.na(diffVector)]
mean(diffVector)#mean relative error
###################################################################################


#Jitter#############################################################################
library(spatioTemporalIndices)
library(spatioTemporalALK)
load("results/run.Rdata")
jj = jit(run,njit = 25,ncores = 1) #Only implemented with use of one core, run must be fitted with twoStage =FALSE
as.matrix(jj$maxVecAll)
####################################################################################





library(spatioTemporalIndices)
library(spatioTemporalALK)
rm(list=ls())

#Read data
dat_l = readRDS("data/catch_at_length_data_ex_rus.rds")
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")

#Configurations length part
conf_l = defConf(years = 1994:2020, # years to use
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
                 minDepth=100,maxDepth=450,
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
run = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk, twoStage = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

#Save indices with covariance structures to be used by SAM
saveIndex(run,file = "indexJoint.txt", folder = "data/NEAhaddockAssessment/")
save(run,file = "results/run.Rdata")


#No ALK uncertainty
load("results/run.Rdata")
pl = run$pl
map = run$map
map$beta0_alk = as.factor(pl$beta0_alk*NA)
map$log_sigma_beta_alk = as.factor(pl$log_sigma_beta_alk*NA)
map$betaLength_alk = as.factor(pl$betaLength_alk*NA)
map$logSigma_alk = as.factor(pl$logSigma_alk*NA)
map$logKappa_alk = as.factor(pl$logKappa_alk*NA)
map$transRho_alk = as.factor(pl$transRho_alk*NA)
map$xS_alk = as.factor(pl$xS_alk*NA)
map$xST_alk = as.factor(pl$xST_alk*NA)

# run model
start_time <- Sys.time()
dat_l = run$dat_l
dat_alk = run$dat_alk
conf_l = run$conf_l
conf_alk = run$conf_alk
confPred = run$confPred
runMAPalk = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,mapSet = map,parSet = pl,runModel = FALSE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

library(TMB)
start_time <- Sys.time()
rep = sdreport(runMAPalk$obj)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed
runMAPalk$rep = rep
saveIndex(runMAPalk,file = "indexJointMAPalk.txt", folder = "data/NEAhaddockAssessment/")
save(runMAPalk,file = "results/runMAPalk.Rdata")



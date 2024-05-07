library(spatioTemporalIndices)
library(spatioTemporalALK)
################################################################################################################
#Finer discretication in space##########################################################
load("results/run.Rdata")
conf_l = run$conf_l
conf_l$cutoff = 55
confPred = run$confPred

start_time <- Sys.time()
run = fitModel(run$dat_l,conf_l, run$confPred,run$dat_alk,run$conf_alk,twoStage = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed
saveIndex(run,file = "indexJointMeshDetailed.txt", folder = "data/NEAhaddockAssessment/")
save(run,file = "results/runMeshDetailed.Rdata")


################################################################################################################
#Finer discretication in length##########################################################
load("results/run.Rdata")
conf_l = run$conf_l
conf_l$reduceLength = 1
conf_l$lengthGroupsReduced = 1:length(conf_l$lengthGroupsReduced)
confPred = run$confPred

start_time <- Sys.time()
runDL1 = fitModel(run$dat_l,conf_l, run$confPred,run$dat_alk,run$conf_alk,twoStage = TRUE)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed
saveIndex(runDL1,file = "indexJointDL1.txt", folder = "data/NEAhaddockAssessment/")
save(runDL1,file = "results/runDL1.Rdata")

###################################################################################
#No spatio-temporal effects
load("results/run.Rdata")
conf_l = run$conf_l
conf_alk = run$conf_alk
conf_l$spatioTemporal = 0
conf_l$spatial = 1
conf_alk$spatioTemporal = 0
conf_alk$spatial = 1
confPred = run$confPred

start_time <- Sys.time()
runNoSST = fitModel(run$dat_l,conf_l, run$confPred,run$dat_alk,conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed
saveIndex(runNoSST,file = "indexJointNoST.txt", folder = "data/NEAhaddockAssessment/")
save(runNoSST,file = "results/runNoST.Rdata")


################################################################################################################
#Run with no covariates#############################
load("results/run.Rdata")
confNoCovariates = run$conf_l
confNoCovariates$splineDepth = c(6,0)
confNoCovariates$sunAlt = c(1,0)
confPred = run$confPred

start_time <- Sys.time()
runNOCovariates = fitModel(run$dat_l,confNoCovariates, run$confPred,run$dat_alk,run$conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

saveIndex(runNOCovariates,file = "indexJointNoCovariates.txt", folder = "data/NEAhaddockAssessment/")
save(runNOCovariates,file = "results/runNoCovariates.Rdata")

#Run with no random walk in beta0######################
load("results/run.Rdata")
conf_lNoRW = run$conf_l
conf_lNoRW$rwBeta0 = 0
conf_alkNoRW = run$conf_alk
conf_alkNoRW$rwBeta0 = 0

start_time <- Sys.time()
runNORW = fitModel(run$dat_l,conf_lNoRW, run$confPred,run$dat_alk,conf_alkNoRW)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed
saveIndex(runNORW,file = "indexJointNoRW.txt", folder = "data/NEAhaddockAssessment/")
save(runNORW,file = "results/runNoRW.Rdata")






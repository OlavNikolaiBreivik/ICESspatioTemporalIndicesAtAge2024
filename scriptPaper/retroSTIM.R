library(spatioTemporalIndices)
library(spatioTemporalALK)

#Retro##############################################
load("results/run.Rdata")
retro = retroSTIM(run,nyears = 10)
save(retro,file = "results/retro.Rdata")
retroMAP = list()
for(i in 1:length(retro)){
  pl = retro[[i]]$pl

  map = retro[[i]]$map
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
  dat_l = retro[[i]]$dat_l
  dat_alk = retro[[i]]$dat_alk
  conf_l = retro[[i]]$conf_l
  conf_alk = retro[[i]]$conf_alk
  confPred = retro[[i]]$confPred
  retroMAP[[i]] = fitModel(dat_l,conf_l, confPred,dat_alk,conf_alk,mapSet = map,parSet = pl,runModel = FALSE)
  end_time <- Sys.time()
  timeUsed = end_time - start_time
  timeUsed

  library(TMB)
  start_time <- Sys.time()
  rep = sdreport(retroMAP[[i]]$obj)
  end_time <- Sys.time()
  timeUsed = end_time - start_time
  timeUsed
  retroMAP[[i]]$rep = rep
  print("--------------------------------------------------------------------------------")
  print(i)
  print("--------------------------------------------------------------------------------")
}

for(i in 1:length(retro)){
  saveIndex(retro[[i]],file = paste0("indexJointRetro",i,".txt"), folder = "data/NEAhaddockAssessment/")
  saveIndex(retroMAP[[i]],file = paste0("indexJointMAPalkRetro",i,".txt"), folder = "data/NEAhaddockAssessment/")
}
#################################################


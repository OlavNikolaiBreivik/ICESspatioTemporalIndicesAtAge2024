library(spatioTemporalIndices)
library(spatioTemporalALK)
library(xtable)
rm(list=ls())

#Read data
dat_l = readRDS("data/catch_at_length_data_ex_rus.rds")
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")

##############################################################
#Latent spatio-temporal effects on and off
##############################################################

#Configurations length part
conf_l = defConf(years = 1994:2020,
                 maxLength = 75,
                 minLength = 20,
                 spatioTemporal = 0,
                 spatial =0,
                 nugget = 1,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,2),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 minDepth=100,maxDepth=450,
                 applyALK = 0,
                 cutoff =80)


confPred = defConfPred(conf=conf_l,Depth="data",cellsize = 20)
runs = list()
timeUsed = list()

start_time <- Sys.time()
runs[[1]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsed[[1]] = end_time - start_time

start_time <- Sys.time()
conf_l$spatial=1
runs[[2]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsed[[2]] = end_time - start_time

start_time <- Sys.time()
conf_l$spatial=0
conf_l$spatioTemporal=1
runs[[3]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsed[[3]] = end_time - start_time

start_time <- Sys.time()
conf_l$spatial=1
conf_l$spatioTemporal=1
runs[[4]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsed[[4]] = end_time - start_time

latentTable = AIC(runs[[1]],runs[[2]], runs[[3]],
               runs[[4]])
latentTable = round(latentTable)

latentTable = cbind(c("Non-spatial","$\\alpha$", "$\\gamma$", "$\\alpha + \\gamma$"), latentTable)
rownames(latentTable) = NULL
colnames(latentTable) = c("Latent effects","Parameters", "AIC change")

dAIC = rep(0,length(runs))
for(i in 1:length(runs)){
  dAIC[i] = AIC(runs[[i]])-AIC(runs[[1]])
}
latentTable[,3] = dAIC
mat = xtable(latentTable, type = "latex",digits = 0)
print(mat, sanitize.text.function = function(x) {x},
      file = "tables/latentTab.tex",
      floating = FALSE, include.rownames=FALSE)

save(runs,file = "results/runsLatentEffectsLength.Rdata")



##############################################################
#Find covariates included
##############################################################
conf_l = defConf(years = 1994:2020,
                 maxLength = 75,
                 minLength = 20,
                 spatioTemporal = 1,
                 spatial =1,
                 nugget = 1,
                 rwBeta0 = 1,
                 sunAlt = c(1,0),
                 splineDepth = c(6,0),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 minDepth=100,maxDepth=450,
                 applyALK = 0,
                 cutoff =80)


confPred = defConfPred(conf=conf_l,Depth="data",cellsize = 20)

runsCovariates = list()
timeUsedC = list()
start_time <- Sys.time()
runsCovariates[[1]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[1]] = end_time - start_time

conf_l$sunAlt = c(1,1)
start_time <- Sys.time()
runsCovariates[[2]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[2]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(1,2)
runsCovariates[[3]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[3]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(0,0)
conf_l$splineDepth = c(6,1)
runsCovariates[[4]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[4]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(0,0)
conf_l$splineDepth = c(6,2)
runsCovariates[[5]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[5]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(1,1)
conf_l$splineDepth = c(6,1)
runsCovariates[[6]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[6]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(1,2)
conf_l$splineDepth = c(6,1)
runsCovariates[[7]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[7]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(1,1)
conf_l$splineDepth = c(6,2)
runsCovariates[[8]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[8]] = end_time - start_time

start_time <- Sys.time()
conf_l$sunAlt = c(1,2)
conf_l$splineDepth = c(6,2)
runsCovariates[[9]] = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsedC[[9]] = end_time - start_time

library(xtable)
covTable = AIC(runsCovariates[[1]], runsCovariates[[2]],
               runsCovariates[[3]], runsCovariates[[4]],
               runsCovariates[[5]], runsCovariates[[6]],
               runsCovariates[[7]], runsCovariates[[8]],
               runsCovariates[[9]])
tmp = c("No covariates", "Sun effect", "Sun effect*",
                       "Depth effect", "Depth effect*", "Sun effect and depth effect",
                       "Sun effect* and depth effect",
                       "Sun effect and depth effect*",
                       "Sun effect* and depth effect*"
                       )
covTable = round(covTable)
covTable = cbind(tmp, covTable)
rownames(covTable) = NULL
colnames(covTable) = c("Explanatory variables","Parameters", paste0("AIC change"))

dAIC = rep(0,length(runsCovariates))
for(i in 1:length(runsCovariates)){
  dAIC[i] = AIC(runsCovariates[[i]])-AIC(runsCovariates[[1]])
}
covTable[,3] = dAIC
mat = xtable(covTable, type = "latex",digits = 0)
print(mat, sanitize.text.function = function(x) {x},
      file = "tables/covariatesTab.tex",
      floating = FALSE)

save(runsCovariates,file = "results/runsCovariatesLength.Rdata")

#Save run with all covariates as length dependent as separate file, figure of depth effect is in supplementary
run = runsCovariates[[9]]
save(run,file = "results/runLengthDependentDepth.Rdata")
png("figures/depthEffectLengthDepentent.png")
plotResults(run,what = "depth")
dev.off()

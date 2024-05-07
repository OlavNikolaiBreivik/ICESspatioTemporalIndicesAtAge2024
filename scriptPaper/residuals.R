library(spatioTemporalIndices)
rm(list=ls())

#Read data
dat_l = readRDS("data/catch_at_length_data_ex_rus.rds")
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")

#######################################################################################
################################# Residuals without random effects ####################
#Configurations length part
conf_l = defConf(years = 1994:2020,
                 maxLength = 75,
                 minLength = 20,
                 spatioTemporal =0 ,
                 spatial =0,
                 nugget = 0,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,1),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 0,
                 cutoff =80)

confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 100)

# run model
start_time <- Sys.time()
run = fitModel(dat_l,conf_l, confPred)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


#Neet to unload spatioTemporalALK
#OSA residuals
options(future.availableCores.methods = "mc.cores")
options(mc.cores = 10)
nHualsRemove = dim(run$data$fishObsMatrix)[1]-50 #50 last hauls, resulting in 600 residuals
mm = max(run$data$fishObsMatrix[-(1:nHualsRemove),])
res <- oneStepPredict(run$obj, observation.name ="obsVector",
                      data.term.indicator="keep",discrete = TRUE, method = "oneStepGeneric",
                      conditional = 1:(nHualsRemove*length(run$data$lengthGroups)), parallel = TRUE,
                      range = c(0,round(10*mm)),
                      reverse = TRUE)

png("figures/noRandomQQ.png")
qqnorm(res$residual,main = "Quantile-quantile no random effects",
       cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
abline(0,1)
dev.off()

png("figures/noRandomACF.png")
acf(res$residual,main = " ", ylab = " ",
    cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
title(main="Autocorrelation no random effects",cex.main = 1.5, ylab = "Autocorrelation", cex.lab = 1.5)
dev.off()

#######################################################################################
################################# Residuals with random effects #######################
library(spatioTemporalIndices)
dat_l = readRDS("data/catch_at_length_data_ex_rus.rds")

conf_l = defConf(years = 1994:2020,
                 maxLength = 75,
                 minLength = 20,
                 spatioTemporal =1 ,
                 spatial =1,
                 nugget = 1,
                 rwBeta0 = 1,
                 sunAlt = c(1,2),
                 splineDepth = c(6,1),
                 dLength = 5,
                 reduceLength = 3,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 minDepth=150,maxDepth=400,
                 applyALK = 0,
                 cutoff =80)
confPred = defConfPred(conf=conf_l,Depth="DATA",cellsize = 100)

# run model
start_time <- Sys.time()
run = fitModel(dat_l,conf_l, confPred)#,dat_alk,conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


#Neet to unload spatioTemporalALK
#OSA residuals
options(future.availableCores.methods = "mc.cores")
options(mc.cores = 10)
nHualsRemove = dim(run$data$fishObsMatrix)[1]-50
mm = max(run$data$fishObsMatrix[-(1:nHualsRemove),])
res <- oneStepPredict(run$obj, observation.name ="obsVector",
                      data.term.indicator="keep",discrete = TRUE, method = "oneStepGeneric",
                      conditional = 1:(nHualsRemove*length(run$data$lengthGroups)), parallel = TRUE,
                      range = c(0,round(10*mm)),
                      reverse = TRUE)
png("figures/withRandomQQ.png")
qqnorm(res$residual,main = "Quantile-quantile with random effects",
       cex.main = 1.5, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
abline(0,1)
dev.off()
png("figures/withRandomACF.png")
acf(res$residual,main = " ", ylab = " ",
    cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
title(main="Autocorrelation with random effects",cex.main = 1.5, ylab = "Autocorrelation", cex.lab = 1.5)
dev.off()





###########################ALK without random effects#########################################
#Run only ALK model
library(spatioTemporalIndices)
library(spatioTemporalALK)
rm(list=ls())

#Read data
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")

runs_alk= list()
#Needed for using the same mesh
conf_l = defConf(years = 1994:2020,
                 maxLength = 75,
                 minLength = 20,
                 dLength = 5,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 cutoff =80)
#Define configurations age part
conf_alk = defConf_alk(years = 1994:2020,
                       maxAge = 10,
                       minAge = 3,
                       spatioTemporal =0,
                       spatial =0,
                       rwBeta0 = 1,
                       readability = 1,
                       cutoff =80)
conf_alk$meshSimilar = TRUE

#Set up data and run model
data_alk = setUpData_alk(dat_alk,conf_alk,conf_l  = conf_l )
par_alk = defpar_alk(data_alk,conf_alk)
map_alk = setMap_alk(conf_alk,par_alk)
start_time <- Sys.time()
#data_alk$readability = rep(1,length(data_alk$readability))
run = fitALK(data = data_alk,par = par_alk,conf = conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed


#OSA residuals
options(future.availableCores.methods = "mc.cores")
options(mc.cores = 10)
conditional = 1:(length(run$data$age)-600)
ageRange = seq(min(run$data$age), max(run$data$age))
res <- oneStepPredict(run$obj, observation.name ="age",
                      data.term.indicator="keep",
                      discrete = TRUE,
                      discreteSupport = ageRange,
                      method = "cdf",
                      conditional = conditional,
                      trace = 3)

png("figures/withoutRandomQQ_ALK.png")
qqnorm(res$residual,main = "Quantile-quantile with random effects",
       cex.main = 1.5, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
abline(0,1)
dev.off()
png("figures/withoutRandomACF_ALK.png")
acf(res$residual,main = " ", ylab = " ",
    cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
title(main="Autocorrelation with random effects",cex.main = 1.5, ylab = "Autocorrelation", cex.lab = 1.5)
dev.off()




###########################ALK with random effects#########################################
#Run only ALK model
library(spatioTemporalIndices)
library(spatioTemporalALK)
#library(xtable)
rm(list=ls())

#Read data
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")
runs_alk= list()
#Needed for using the same mesh
conf_l = defConf(years = 1994:2020,
                 maxLength = 75,
                 minLength = 20,
                 dLength = 5,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 cutoff =80)
#Define configurations age part
conf_alk = defConf_alk(years = 1994:2020,
                       maxAge = 10,
                       minAge = 3,
                       spatioTemporal = 2,
                       spatial =1,
                       rwBeta0 = 1,
                       readability = 1,
                       cutoff =80 )
conf_alk$meshSimilar = TRUE

#Set up data and run model
data_alk = setUpData_alk(dat_alk,conf_alk,conf_l  = conf_l )
par_alk = defpar_alk(data_alk,conf_alk)
map_alk = setMap_alk(conf_alk,par_alk)
start_time <- Sys.time()
#data_alk$readability = rep(1,length(data_alk$readability))
run = fitALK(data = data_alk,par = par_alk,conf = conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

#OSA residuals
options(future.availableCores.methods = "mc.cores")
options(mc.cores = 5)
conditional = 1:(length(run$data$age)-600)
ageRange = seq(min(run$data$age), max(run$data$age))
res <- oneStepPredict(run$obj, observation.name ="age",
                      data.term.indicator="keep",
                      discrete = TRUE,
                      discreteSupport = ageRange,
                      method = "cdf",
                      conditional = conditional,
                      trace = 3)

png("figures/withRandomQQ_ALK.png")
qqnorm(res$residual,main = "Quantile-quantile with random effects",
       cex.main = 1.5, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
abline(0,1)
dev.off()
png("figures/withRandomACF_ALK.png")
acf(res$residual,main = " ", ylab = " ",
    cex.main = 2, cex.lab = 1.5, cex = 1.5,cex.axis = 1.5)
title(main="Autocorrelation with random effects",cex.main = 1.5, ylab = "Autocorrelation", cex.lab = 1.5)
dev.off()


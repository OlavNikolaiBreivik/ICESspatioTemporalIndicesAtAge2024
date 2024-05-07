#Run only ALK model
library(spatioTemporalIndices)
library(spatioTemporalALK)
library(xtable)
rm(list=ls())

#Read data
dat_alk = readRDS("data/catch_at_age_data_ex_rus.rds")
runs_alk= list()
#Needed for using the same mesh
conf_l = defConf(years = 1994:2020, # years to use, use all years with data by default
                 maxLength = 75, # Numeric = use directly; NULL = use input data to determine
                 minLength = 20, # Numeric = use directly; NULL = use input data to determine
                 dLength = 5,
                 stratasystem = list(dsn="shapefiles/Vintertokt_strata_system", layer = "Vintertoktet_nye_strata"),
                 cutoff =80)
#Define configurations age part
conf_alk = defConf_alk(years = 1994:2020,
                       maxAge = 10,
                       minAge = 3,
                       spatioTemporal = 0,
                       spatial =0,
                       rwBeta0 = 1,
                       cutoff =80, cbound = 130 )
conf_alk$meshSimilar = TRUE

#Set up data and run model
data_alk = setUpData_alk(dat_alk,conf_alk,conf_l  = conf_l )
par_alk = defpar_alk(data_alk,conf_alk)
map_alk = setMap_alk(conf_alk,par_alk)
start_time <- Sys.time()
runs_alk[[1]] = fitALK(data = data_alk,par = par_alk,conf = conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

###Spatial effect###
conf_alk$spatial = 1
data_alk = setUpData_alk(dat_alk,conf_alk,conf_l  = conf_l )
par_alk = defpar_alk(data_alk,conf_alk)
map_alk = setMap_alk(conf_alk,par_alk)
start_time <- Sys.time()
runs_alk[[2]] = fitALK(data = data_alk,par = par_alk,conf = conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

###Spatio-temporal effect###
conf_alk$spatial = 0
conf_alk$spatioTemporal = 2
data_alk = setUpData_alk(dat_alk,conf_alk,conf_l  = conf_l )
par_alk = defpar_alk(data_alk,conf_alk)
map_alk = setMap_alk(conf_alk,par_alk)
start_time <- Sys.time()
runs_alk[[3]] = fitALK(data = data_alk,par = par_alk,conf = conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

###Spatial and spatio-temporal effect##
conf_alk$spatial = 1
conf_alk$spatioTemporal = 2
data_alk = setUpData_alk(dat_alk,conf_alk,conf_l  = conf_l )
par_alk = defpar_alk(data_alk,conf_alk)
map_alk = setMap_alk(conf_alk,par_alk)
start_time <- Sys.time()
runs_alk[[4]] = fitALK(data = data_alk,par = par_alk,conf = conf_alk)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed

latentTable_alk =   round(AIC(runs_alk[[1]],runs_alk[[2]],runs_alk[[3]],runs_alk[[4]]))
latentTable_alk = cbind(c("Non-spatial" ,"$\\alpha^{(alk)}$", "$\\gamma^{(alk)}$",
                          "$\\alpha^{(alk)}$ + $\\gamma^{(alk)}$"), latentTable_alk)

dAIC = rep(0,length(runs_alk))
for(i in 1:length(runs_alk)){
  dAIC[i] = AIC(runs_alk[[i]])-AIC(runs_alk[[1]])
}
latentTable_alk[,3] = dAIC

colnames(latentTable_alk) = c("Latent effects","Parameters", "AIC change")
mat = xtable(latentTable_alk, type = "latex",digits = 0)
print(mat, sanitize.text.function = function(x) {x},
      file = "tables/latentTab_alk.tex",
      floating = FALSE, include.rownames=FALSE)

save(runs_alk,file = "results/runsLatentEffectsAge.Rdata")

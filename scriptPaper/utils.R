plotSim = function(sim,...){
  run = attributes(sim)$run
  pl = run$pl
  plSd = run$plSd
  what = c("","log_sigma", "log_kappa","tan_rho_l")
  if(run$conf_l$spatioTemporal==1){
    what = c(what,"tan_rho_t")
  }
  what = c(what,"betaLength_alk")
  if(run$conf_alk$spatial!=0){
    what = c(what,"logSigma_alk")
    what = c(what,"logKappa_alk")
  }
  if(run$conf_alk$spatioTemporal==1){
    what = c(what,"transRho_alk")
  }
  if(run$conf_l$rwBeta0==0){
    what = c(what,"beta0")
  }
  if(run$conf_alk$rwBeta0==0){
    what = c(what,"beta0_alk")
  }


  namesWhat = NULL;
  for(i in 1:length(what)){
    namesWhat = c(namesWhat, rep(what[i], length(unlist(pl[what[i]]))))
  }
  namesGreek = namesWhat
  namesGreek[which(namesWhat=="log_sigma")] = c(expression(log~sigma[alpha]),expression(log~sigma[gamma]),expression(log~sigma[psi]))
  namesGreek[which(namesWhat=="log_kappa")] = c(expression(log~kappa[alpha]),expression(log~kappa[gamma]))
  namesGreek[which(namesWhat=="tan_rho_l")] = c(expression("internal"~rho[alpha][",l"]),expression("internal"~rho[gamma][",l"]),expression("internal"~rho[psi][",l"]))
  namesGreek[which(namesWhat=="betaLength_alk")] = c(expression(beta[2]^"(a)"),expression(beta[3]^"(a)"),expression(beta[4]^"(a)"),
                                                        expression(beta[5]^"(a)"),expression(beta[6]^"(a)"),expression(beta[7]^"(a)"),
                                                     expression(beta[8]^"(a)"),expression(beta[9]^"(a)"))
  namesGreek[which(namesWhat=="tan_rho_t")] = c(expression("internal"~rho[gamma][",y"]))
  namesGreek[which(namesWhat=="logSigma_alk")] = c(expression(log~sigma[alpha]^"(alk)"),expression(log~sigma[gamma]^"(alk)"))
  namesGreek[which(namesWhat=="logKappa_alk")] = c(expression(log~kappa[alpha]^"(alk)"),expression(log~kappa[gamma]^"(alk)"))


  par(mar = c(9, 5, 3, 2))
  library(plotrix)
  pointEst = unlist(pl[what])
  li =  unlist(pl[what]) - 2*unlist(plSd[what])
  ui = unlist(pl[what]) + 2*unlist(plSd[what])

  pointEst["tan_rho_l1"] =    2/(1 + exp(-2 * pointEst["tan_rho_l1"] )) - 1
  pointEst["tan_rho_l2"] =   2/(1 + exp(-2 * pointEst["tan_rho_l2"] )) - 1
  li["tan_rho_l1"] =   2/(1 + exp(-2 * li["tan_rho_l1"] )) - 1
  li["tan_rho_l2"] =    2/(1 + exp(-2 * li["tan_rho_l2"] )) - 1
  ui["tan_rho_l1"] =   2/(1 + exp(-2 * ui["tan_rho_l1"] )) - 1
  ui["tan_rho_l2"] =   2/(1 + exp(-2 * ui["tan_rho_l2"] )) - 1
  pointEst["tan_rho_l1"] = pointEst["tan_rho_l1"]^(1/run$conf_l$reduceLength) #Accommodates for reduced length dimension
  pointEst["tan_rho_l2"]= pointEst["tan_rho_l2"]^(1/run$conf_l$reduceLength)
  li["tan_rho_l1"] = li["tan_rho_l1"]^(1/run$conf_l$reduceLength)
  li["tan_rho_l2"] = li["tan_rho_l2"]^(1/run$conf_l$reduceLength)
  ui["tan_rho_l1"] = ui["tan_rho_l1"]^(1/run$conf_l$reduceLength)
  ui["tan_rho_l2"] = ui["tan_rho_l2"]^(1/run$conf_l$reduceLength)
  pointEst["tan_rho_l1"]  = -0.5*log(2/(pointEst["tan_rho_l1"]+1) - 1)
  pointEst["tan_rho_l2"]  = -0.5*log(2/(pointEst["tan_rho_l2"]+1) - 1)
  li["tan_rho_l1"]  = -0.5*log(2/(li["tan_rho_l1"]+1) - 1)
  li["tan_rho_l2"]  = -0.5*log(2/(li["tan_rho_l2"]+1) - 1)
  ui["tan_rho_l1"]  = -0.5*log(2/(ui["tan_rho_l1"]+1) - 1)
  ui["tan_rho_l2"]  = -0.5*log(2/(ui["tan_rho_l2"]+1) - 1)

  plotCI(pointEst,
         li =  li,
         ui = ui,
         xlab = "Parameter",
         ylab = "Value",cex.lab = 1.4,
         xaxt="n",...)

  axis(1,at = 1:length(unlist(pl[what])), namesGreek,  las=2,cex.lab= 2,cex.axis = 1.3)
  k = lapply(sim, function(f){
    if(length(f)>1){
      plSim = unlist(f$pl[what])
      plSim["tan_rho_l1"] =    2/(1 + exp(-2 * plSim["tan_rho_l1"] )) - 1
      plSim["tan_rho_l2"] =   2/(1 + exp(-2 * plSim["tan_rho_l2"] )) - 1
      plSim["tan_rho_l1"] = plSim["tan_rho_l1"]^(1/run$conf_l$reduceLength) #Accommodates for reduced length dimension
      plSim["tan_rho_l2"]= plSim["tan_rho_l2"]^(1/run$conf_l$reduceLength)
      plSim["tan_rho_l1"]  = -0.5*log(2/(plSim["tan_rho_l1"]+1) - 1)
      plSim["tan_rho_l2"]  = -0.5*log(2/(plSim["tan_rho_l2"]+1) - 1)

      points(plSim, col = 'red', pch = 2, cex = 0.5)
    }
  })
}

abundanceSim = function(run,sim){
  par = run$pl
  simData = attributes(sim)$simData
  abundanceRealizations = list()
  for(i in 1:length(sim)){
    par$xS = simData[[i]]$xS
    par$xST = simData[[i]]$xST
    par$xS_alk = simData[[i]]$xS_alk
    par$xST_alk = simData[[i]]$xST_alk
    simObj =  fitModel(run$dat_l,run$conf_l, run$confPred,run$dat_alk,run$conf_alk ,parSet = par, runModel = FALSE)$obj
#    rep =simObj$obj$report()
    sdrep = sdreport(simObj,ignore.parm.uncertainty = TRUE,skip.delta.method=TRUE)
    rl = as.list(sdrep,what = "Est", report = TRUE)
    abundanceRealizations[[i]] = rl$logAgeIndex
  }
  return(abundanceRealizations)
}


partable<-function(run){
  pl = as.list(run$rep,"Est")
  plsd = as.list(run$rep,"Std. Error")

  digits = 3

  rho_l =  2/(1 + exp(-2 * pl$tan_rho_l)) - 1
  rho_l_U =  2/(1 + exp(-2 * (pl$tan_rho_l + 1.96*plsd$tan_rho_l))) - 1
  rho_l_L =  2/(1 + exp(-2 * (pl$tan_rho_l - 1.96*plsd$tan_rho_l))) - 1
  rho_t =  2/(1 + exp(-2 * pl$tan_rho_t)) - 1
  rho_t_U =  2/(1 + exp(-2 * (pl$tan_rho_t+ 1.96*plsd$tan_rho_t))) - 1
  rho_t_L =  2/(1 + exp(-2 * (pl$tan_rho_t- 1.96*plsd$tan_rho_t))) - 1

  rho_l[1:2] = rho_l[1:2]^(1/run$conf_l$reduceLength) #Accommodates for reduced length dimension
  rho_l_U[1:2] = rho_l_U[1:2]^(1/run$conf_l$reduceLength)
  rho_l_L[1:2] = rho_l_L[1:2]^(1/run$conf_l$reduceLength)

  rho_l = round(rho_l,digits)
  rho_l_U = round(rho_l_U,digits)
  rho_l_L = round(rho_l_L,digits)
  rho_t = round(rho_t,digits)
  rho_t_U = round(rho_t_U,digits)
  rho_t_L = round(rho_t_L,digits)


  sigma = round(exp(pl$log_sigma),digits)
  sigmaU = round(exp(pl$log_sigma + 1.96*plsd$log_sigma),digits)
  sigmaL = round(exp(pl$log_sigma - 1.96*plsd$log_sigma),digits)

  kappa = round(exp(pl$log_kappa),4)
  kappaU = round(exp(pl$log_kappa +  1.96*plsd$log_kappa),digits)
  kappaL = round(exp(pl$log_kappa -  1.96*plsd$log_kappa),digits)

  sigmaRW = c(exp(pl$log_sigma_beta0), exp(pl$log_sigma_beta0 + 1.96*plsd$log_sigma_beta0*c(-1,1)))
  sigmaRW_alk = c(exp(pl$log_sigma_beta_alk[1]), exp(pl$log_sigma_beta_alk[1] + 1.96*plsd$log_sigma_beta_alk[1]*c(-1,1)))

  sigma_alkS = c(exp(pl$logSigma_alk)[1],exp(pl$logSigma_alk[1] + 1.96*plsd$logSigma_alk[1]*c(-1,1)))
  sigma_alkST = c(exp(pl$logSigma_alk)[2],exp(pl$logSigma_alk[2] + 1.96*plsd$logSigma_alk[2]*c(-1,1)))

  kappa_alkS = c(exp(pl$logKappa_alk)[1],exp(pl$logKappa_alk[1] + 1.96*plsd$logKappa_alk[1]*c(-1,1)))
  kappa_alkST = c(exp(pl$logKappa_alk)[2],exp(pl$logKappa_alk[2] + 1.96*plsd$logKappa_alk[2]*c(-1,1)))

  sigmaRW = round(sigmaRW,digits);
  sigmaRW_alk = round(sigmaRW_alk,digits);
  sigma_alkS = round(sigma_alkS,digits);
  sigma_alkST = round(sigma_alkST,digits);
  kappa_alkS = round(kappa_alkS,digits);
  kappa_alkST = round(kappa_alkST,digits);

  if(run$conf_l$spatial==0){
    sigma[1] = "-";sigmaL[1] = "-";sigmaU[1] = "-"
    kappa[1] = "-";kappaL[1] = "-";kappaU[1] = "-"
    rho_l[1] = "-";rho_l_L[1] = "-";rho_l_U[1] = "-"
  }
  if(run$conf_l$spatioTemporal==0){
    sigma[2] = "-";sigmaL[2] = "-";sigmaU[2] = "-"
    kappa[2] = "-";kappaL[2] = "-";kappaU[2] = "-"
    rho_l[2] = "-";rho_l_L[2] = "-";rho_l_U[2] = "-"
    rho_t = "-";rho_t_L = "-";rho_t_U = "-"
  }
  if(run$conf_alk$spatial==0){
    sigma_alkS = rep("-",3);
    kappa_alkS = rep("-",3);
  }
  if(run$conf_alk$spatioTemporal==0){
    sigma_alkST = rep("-",3);
    kappa_alkST = rep("-",3);
  }

  log_lambda = c(pl$log_lambda[1],pl$log_lambda[1] + 1.96*plsd$log_lambda[1]*c(-1,1))


  par = matrix(0,16,3)
  par[1,1] = rho_l[1]; par[1,2:3] = c(rho_l_L[1],rho_l_U[1])
  par[2,1] = rho_l[2]; par[2,2:3] = c(rho_l_L[2],rho_l_U[2])
  par[3,1] = rho_l[3]; par[3,2:3] = c(rho_l_L[3],rho_l_U[3])
  par[4,1] = rho_t; par[4,2:3] = c(rho_t_L,rho_t_U)

  par[5,1] = sigma[1]; par[5,2:3] = c(sigmaL[1],sigmaU[1])
  par[6,1] = sigma[2]; par[6,2:3] = c(sigmaL[2],sigmaU[2])
  par[7,1] = sigma[3]; par[7,2:3] = c(sigmaL[3],sigmaU[3])
  par[8,1] = kappa[1]; par[8,2:3] = c(kappaL[1],kappaU[1])
  par[9,1] = kappa[2]; par[9,2:3] = c(kappaL[2],kappaU[2])

  par[10,] = sigmaRW
  par[11,] = sigmaRW_alk
  par[12,] = sigma_alkS
  par[13,] = sigma_alkST
  par[14,] = kappa_alkS
  par[15,] = kappa_alkST
  par[16,] = log_lambda

  rownames(par) = c("rho_l_S","rho_l_ST","rho_l_nugg","rho_t","sigmaS", "sigmaST","sigmaNugget", "kappaS", "kappaST", "sigma_beta0_RW",
                    "sigmaRW_alk","sigmaS_alk","sigmaST_alk", "kappaS_alk", "kappaST_alk", "log_lambda")
  colnames(par) = c("MLE", "0.025Percentile", "0.975Percentile")
  return(par)
}



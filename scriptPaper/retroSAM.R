#First run runAssessment.R
fitSepCovRet = list()
fitJointCovRet = list()
fitJointRet = list()
fitRet = retro(fit, year = 10, ncores = 1)

fitJointAtt = fitJoint
fitJointCovAtt = fitJointCov
fitSepCovAtt = fitSepCov

for(y in 1:length(fitRet)){
  cn<-read.ices("data/NEAhaddockAssessment/cn.dat")
  cn = cn[1:(dim(cn)[1]-y),]
  cw<-read.ices("data/NEAhaddockAssessment/cw.dat")
  cw = cw[1:(dim(cw)[1]-y),]
  dw<-read.ices("data/NEAhaddockAssessment/dw.dat")
  dw = dw[1:(dim(dw)[1]-y),]
  lf<-read.ices("data/NEAhaddockAssessment/lf.dat")
  lf = lf[1:(dim(lf)[1]-y),]
  lw<-read.ices("data/NEAhaddockAssessment/lw.dat")
  lw = lw[1:(dim(lw)[1]-y),]
  mo<-read.ices("data/NEAhaddockAssessment/mo.dat")
  mo = mo[1:(dim(mo)[1]-y),]
  nm<-read.ices("data/NEAhaddockAssessment/nm.dat")
  nm = nm[1:(dim(nm)[1]-y),]
  pf<-read.ices("data/NEAhaddockAssessment/pf.dat")
  pf = pf[1:(dim(pf)[1]-y),]
  pm<-read.ices("data/NEAhaddockAssessment/pm.dat")
  pm = pm[1:(dim(pm)[1]-y),]
  sw<-read.ices("data/NEAhaddockAssessment/sw.dat")
  sw = sw[1:(dim(sw)[1]-y),]
  surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
  time = attributes(surveys[[1]])$time

  #Read joint index
  indexJoint= as.matrix(read.table(paste0("data/NEAhaddockAssessment/indexJointRetro",y,".txt", sep = " ")))[,-1]
  load(paste0("data/NEAhaddockAssessment/cov_indexJointRetro",y,".Rda"))
  covJoint = lapply(covYears,function(f){
    cov = f[-1,]
    cov = cov[,-1]
    cov
  } )

  #Read mapALK index
  indexJoint= as.matrix(read.table(paste0("data/NEAhaddockAssessment/indexJointMAPalkRetro",y,".txt", sep = " ")))[,-1]
  load(paste0("data/NEAhaddockAssessment/cov_indexJointMAPalkRetro",y,".Rda"))
  covSep = lapply(covYears,function(f){
    cov = f[-1,]
    cov = cov[,-1]
    cov
  } )
  #Fit joint
  surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
  att = attributes(surveys[[1]])
  surveys$`BS-NoRu-Q1(BTr)` = indexJoint
  attributes(surveys[[1]])$time = att$time
  attributes(surveys[[1]])$twofirst = att$twofirst
  attributes(surveys[[1]])$dimnames[[2]] = att$dimnames[[2]]
  dat<-setup.sam.data(surveys=surveys,
                      residual.fleet=cn,
                      prop.mature=mo,
                      stock.mean.weight=sw,
                      catch.mean.weight=cw,
                      dis.mean.weight=dw,
                      land.mean.weight=lw,
                      prop.f=pf,
                      prop.m=pm,
                      natural.mortality=nm,
                      land.frac=lf)
  conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
  fitJoint<-sam.fit(dat,conf,par)

  #Fit joint with provided covariance structure
  surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
  att = attributes(surveys[[1]])
  surveys$`BS-NoRu-Q1(BTr)` = indexJoint
  attributes(surveys[[1]])$time = att$time
  attributes(surveys[[1]])$twofirst = att$twofirst
  attributes(surveys[[1]])$dimnames[[2]] = att$dimnames[[2]]
  attr(surveys[[1]], "cov-weight") <- covJoint
  dat<-setup.sam.data(surveys=surveys,
                      residual.fleet=cn,
                      prop.mature=mo,
                      stock.mean.weight=sw,
                      catch.mean.weight=cw,
                      dis.mean.weight=dw,
                      land.mean.weight=lw,
                      prop.f=pf,
                      prop.m=pm,
                      natural.mortality=nm,
                      land.frac=lf)
  conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
  par<-defpar(dat,conf)
  fitJointCov<-sam.fit(dat,conf,par)

  #Fit separate with provided covariance
  surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
  surveys$`BS-NoRu-Q1(BTr)` = indexJoint
  attributes(surveys[[1]])$time = att$time
  attributes(surveys[[1]])$twofirst = att$twofirst
  attributes(surveys[[1]])$dimnames[[2]] = att$dimnames[[2]]
  attr(surveys[[1]], "cov-weight") <- covSep
  dat<-setup.sam.data(surveys=surveys,
                      residual.fleet=cn,
                      prop.mature=mo,
                      stock.mean.weight=sw,
                      catch.mean.weight=cw,
                      dis.mean.weight=dw,
                      land.mean.weight=lw,
                      prop.f=pf,
                      prop.m=pm,
                      natural.mortality=nm,
                      land.frac=lf)

  par<-defpar(dat,conf)
  fitSepCov<-sam.fit(dat,conf,par)

  fitJointRet[[y]] = fitJoint
  fitJointCovRet[[y]] =fitJointCov
  fitSepCovRet[[y]] =fitSepCov
}

attr(fitJointRet, "fit")= fitJointAtt
attr(fitJointCovRet, "fit")= fitJointCovAtt
attr(fitSepCovRet, "fit")= fitSepCovAtt

class(fitJointRet)<-"samset"
class(fitJointCovRet)<-"samset"
class(fitSepCovRet)<-"samset"


library(scales)
png("figures/ssbSepCovRet.png")
ssbplot(fitWeb, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Proposed~index~with~bold(R)[y]^fixedALK ),cex.main = 2,cex.lab = 1.4, ylim =c(0,650000))
ssbplot(fitSepCovRet,xlim = c(2000,2020), add = TRUE)
dev.off()

png("figures/ssbJointCovRet.png")
ssbplot(fitWeb, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Proposed~index~with~bold(R)[y] ),cex.main = 2,cex.lab = 1.4, ylim =c(0,650000))
ssbplot(fitJointCovRet,xlim = c(2000,2020), add = TRUE)
dev.off()

png("figures/ssbJointRet.png")
ssbplot(fitWeb, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Proposed~index~with~no~covariance~data),cex.main = 2,cex.lab = 1.4, ylim =c(0,650000))
ssbplot(fitJointRet,xlim = c(2000,2020), add = TRUE)
dev.off()

png("figures/ssbRet.png")
ssbplot(fitWeb, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Official~index~currently~used),cex.main = 2,cex.lab = 1.4, ylim =c(0,650000))
ssbplot(fitRet,xlim = c(2000,2020), add = TRUE)
dev.off()

png("figures/fRet.png")
fbarplot(fitWeb,partial = FALSE, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Official~index~currently~used),cex.main = 2,cex.lab = 1.2, ylim =c(0,0.6))
fbarplot(fitRet, add = TRUE)
dev.off()

png("figures/fSepCovRet.png")
fbarplot(fitWeb,partial = FALSE, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Proposed~index~with~bold(R)[y]^fixedALK ),cex.main = 2,cex.lab = 1.2, ylim =c(0,0.6))
fbarplot(fitSepCovRet, add = TRUE)
dev.off()

png("figures/fJointCovRet.png")
fbarplot(fitWeb,partial = FALSE, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Proposed~index~with~bold(R)[y] ),cex.main = 2,cex.lab = 1.2, ylim =c(0,0.6))
fbarplot(fitJointCovRet, add = TRUE)
dev.off()

png("figures/fJointRet.png")
fbarplot(fitWeb,partial = FALSE, col = 'green',cicol = alpha("green",0.1),xlim = c(2000,2022), main = expression(Proposed~index),cex.main = 2, ylim =c(0,0.6))
fbarplot(fitJointRet, add = TRUE)
dev.off()


#Joint desnity of SSB_Y and F_{Y-1}
xl <- expression(log(SSB[T]))
yl <- expression(log(bar(F)[4-7~','][T-1]))
mtext(yl,side=2, line =2, cex = 1.2)  # "line" sets where the text appears relative to the axis.

library(TMB)
setEPS()
postscript("figures/retroSFF.eps")
mfrow = c(3,4)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
sdrepWeb = sdreport(fitWeb$obj,getReportCovariance = TRUE)
for(i in 0:length(fitSepCovRet)){
  yy = which(fitWeb$data$year==(2020-i))
  ssbIndexFull = which(names(fitWeb$sdrep$value)=="logssb")[yy]
  fbarIndexFull = which(names(fitWeb$sdrep$value)=="logfbar")[yy-1]
  if(i==0){
    ssbIndex = rev(which(names( attr(fitJointCovRet, "fit")$sdrep$value)=="logssb"))[1]
    fbarIndex = rev(which(names( attr(fitJointCovRet, "fit")$sdrep$value)=="logfbar"))[2]

    sdrepJointCov = sdreport(attr(fitJointCovRet, "fit")$obj,getReportCovariance = TRUE)
    sdrepJoint = sdreport(attr(fitJointRet, "fit")$obj,getReportCovariance = TRUE)
    sdrepFit = sdreport(attr(fitRet, "fit")$obj,getReportCovariance = TRUE)
    sdrepSepCov = sdreport(attr(fitSepCovRet, "fit")$obj,getReportCovariance = TRUE)
  }else{
    ssbIndex = rev(which(names(fitSepCovRet[[i]]$sdrep$value)=="logssb"))[1]
    fbarIndex = rev(which(names(fitSepCovRet[[i]]$sdrep$value)=="logfbar"))[2]

    sdrepJointCov = sdreport(fitJointCovRet[[i]]$obj,getReportCovariance = TRUE)
    sdrepJoint = sdreport(fitJointRet[[i]]$obj,getReportCovariance = TRUE)
    sdrepFit = sdreport(fitRet[[i]]$obj,getReportCovariance = TRUE)
    sdrepSepCov = sdreport(fitSepCovRet[[i]]$obj,getReportCovariance = TRUE)
  }

  muWeb = sdrepWeb$value[c(ssbIndexFull,fbarIndexFull)]
  covWeb = sdrepWeb$cov[c(ssbIndexFull,fbarIndexFull),c(ssbIndexFull,fbarIndexFull)]
  muJointCov = sdrepJointCov$value[c(ssbIndex,fbarIndex)]
  covJointCov = sdrepJointCov$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]
  muJoint = sdrepJoint$value[c(ssbIndex,fbarIndex)]
  covJoint = sdrepJoint$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]
  muSepCov = sdrepSepCov$value[c(ssbIndex,fbarIndex)]
  covSepCov = sdrepSepCov$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]

  muStox = sdrepFit$value[c(ssbIndex,fbarIndex)]
  covStox = sdrepFit$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]
  if(i ==0 | i==4){
    plot(-99,-99,xlim = c(11.8,13.8),ylim = c(-2.3,-0.40),
         xaxt="n",
         main = "")
  }else if(i==8){
    plot(-99,-99,xlim = c(11.8,13.8),ylim = c(-2.3,-0.40),
         main = "")
  }else if(i>8){
    plot(-99,-99,xlim = c(11.8,13.8),ylim = c(-2.3,-0.40),
         yaxt="n",
         main = "")
  }else{
    plot(-99,-99,xlim = c(11.8,13.8),ylim = c(-2.3,-0.40),
         yaxt="n",xaxt="n",
         main = "")
  }
  mtext(paste0("T = ", 2020-i),side=3, line =-2.5, cex = 1.5)  # "line" sets where the text appears relative to the axis.

  if(i==4)mtext(yl,side=2, line =2, cex = 1.2)  # "line" sets where the text appears relative to the axis.
  if(i==9)mtext(xl,side=1, line =2.8, cex = 1.2)  # "line" sets where the text appears relative to the axis.

  lines( ellipse::ellipse( covWeb, centre = muWeb,level = 0.95) ,type = 'l',lt =2,lw = 1)
  points(muWeb[1],muWeb[2],cex = 1,pch = 19)
  lines( ellipse::ellipse( covStox, centre = muStox,level = 0.95) , col="blue",type = 'l',lt = 2,lw = 1)
  points(muStox[1],muStox[2],cex = 2,pch = 18, col = 'blue')
  lines( ellipse::ellipse( covJoint, centre = muJoint,level = 0.95) , col="purple",type = 'l',lt = 2,lw = 1)
  points(muJoint[1],muJoint[2],cex = 2,pch = 20, col = 'purple')
  lines( ellipse::ellipse( covJointCov, centre = muJointCov,level = 0.95) , col="#66C2A5",type = 'l',lw = 1)
  points(muJointCov[1],muJointCov[2],cex = 1,pch = 15, col = "#66C2A5")
  lines( ellipse::ellipse( covSepCov, centre = muSepCov,level = 0.95) , col="red",type = 'l',lw = 1)
  points(muSepCov[1],muSepCov[2],cex = 1,pch = 17, col = 'red')

  if(i==0){
    legend("bottomright", inset=.01,
           c("Assessment in year 2023", expression(bold(I)[y]~no~covariance~data), expression(bold(I)[y]~with~R[y]), expression(bold(I)[y]~with~R[y]^fixedALK ),"Official index"),
           pch = c(19,20,15,17,18), col = c("black", "purple","#66C2A5","red","blue"), cex = 0.75, lwd = 1.6, lty = c(2,2,1,1,2))
  }
}
dev.off()


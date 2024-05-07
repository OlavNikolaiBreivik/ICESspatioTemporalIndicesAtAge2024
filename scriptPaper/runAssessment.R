library(stockassessment)
load("data/NEAhaddockAssessment/NEA_haddock_2023_officiall.Rds")#Load officiall assessment model called fitWeb
fitWeb<-stockassessment:::refit(fitWeb) #need local version

#Read data
cn<-read.ices("data/NEAhaddockAssessment/cn.dat")
cw<-read.ices("data/NEAhaddockAssessment/cw.dat")
dw<-read.ices("data/NEAhaddockAssessment/dw.dat")
lf<-read.ices("data/NEAhaddockAssessment/lf.dat")
lw<-read.ices("data/NEAhaddockAssessment/lw.dat")
mo<-read.ices("data/NEAhaddockAssessment/mo.dat")
nm<-read.ices("data/NEAhaddockAssessment/nm.dat")
pf<-read.ices("data/NEAhaddockAssessment/pf.dat")
pm<-read.ices("data/NEAhaddockAssessment/pm.dat")
sw<-read.ices("data/NEAhaddockAssessment/sw.dat")
surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")

#Read joint index
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJoint.txt", sep = " "))[,-1]
varJoint=  as.matrix(read.table("data/NEAhaddockAssessment/sdindexJoint.txt", sep = " ",header = T))[,-1]^2
load("data/NEAhaddockAssessment/cov_indexJoint.Rda")

covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )

indexSep=  as.matrix(read.table("data/NEAhaddockAssessment/indexJointMAPalk.txt", sep = " "))[,-1]
varSep=  as.matrix(read.table("data/NEAhaddockAssessment/sdindexJointMAPalk.txt", sep = " ",header = T))[,-1]^2
load("data/NEAhaddockAssessment/cov_indexJointMAPalk.Rda")
covSep = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )

dim = dim(surveys$`BS-NoRu-Q1(BTr)`)
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
fit<-sam.fit(dat,conf,par)

#############Fit joint
surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
surveys$`BS-NoRu-Q1(BTr)`[1:dim[1],1:dim[2]] = indexJoint
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
fitJoint<-sam.fit(dat,conf,par)

#Fit joint with provided covariance structure
surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
surveys$`BS-NoRu-Q1(BTr)`[1:dim[1],1:dim[2]] = indexJoint
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

par<-defpar(dat,conf)
fitJointCov<-sam.fit(dat,conf,par)

#Fit separate with provided covariance
surveys<-read.ices("data/NEAhaddockAssessment/survey.dat")
surveys$`BS-NoRu-Q1(BTr)`[1:dim[1],1:dim[2]] = indexSep
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

#Joint desnity of SSB_Y and F_{Y-1}
yy = which(rev(fitWeb$data$year==2020))
ssbIndexFull = rev(which(names(fitWeb$sdrep$value)=="logssb"))[yy]
fbarIndexFull = rev(which(names(fitWeb$sdrep$value)=="logfbar"))[yy-1]
ssbIndex = rev(which(names(fit$sdrep$value)=="logssb"))[1]
fbarIndex = rev(which(names(fit$sdrep$value)=="logfbar"))[2]

library(TMB)
sdrepStox = sdreport(fit$obj,getReportCovariance = TRUE)
sdrepWeb = sdreport(fitWeb$obj,getReportCovariance = TRUE)
sdrepJoint = sdreport(fitJoint$obj,getReportCovariance = TRUE)
sdrepJointCov = sdreport(fitJointCov$obj,getReportCovariance = TRUE)
sdrepSepCov = sdreport(fitSepCov$obj,getReportCovariance = TRUE)

muWeb = sdrepWeb$value[c(ssbIndexFull,fbarIndexFull)]
covWeb = sdrepWeb$cov[c(ssbIndexFull,fbarIndexFull),c(ssbIndexFull,fbarIndexFull)]
muStox = sdrepStox$value[c(ssbIndex,fbarIndex)]
covStox = sdrepStox$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]
muJoint = sdrepJoint$value[c(ssbIndex,fbarIndex)]
covJoint = sdrepJoint$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]
muJointCov = sdrepJointCov$value[c(ssbIndex,fbarIndex)]
covJointCov = sdrepJointCov$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]
muSepCov = sdrepSepCov$value[c(ssbIndex,fbarIndex)]
covSepCov = sdrepSepCov$cov[c(ssbIndex,fbarIndex),c(ssbIndex,fbarIndex)]

setEPS()
postscript("figures/ssbAndF.eps")
xl <- expression(log(SSB[2020]))
yl <- expression(log(bar(F)[4-7~","~2019]))
plot(-99,-99,xlim = c(11.95,12.9),ylim = c(-1.5,-0.55),
     xlab ="",
     ylab = "",
     main = expression(Estimated~log(SSB[2020])~and~log(bar(F)[4-7~","~2019]) ),
     cex.main = 2)

mtext(xl,side=1, line =2.8, cex = 1.5)  # "line" sets where the text appears relative to the axis.
mtext(yl,side=2, line =2, cex = 1.5)  # "line" sets where the text appears relative to the axis.

lines( ellipse::ellipse( covWeb, centre = muWeb,level = 0.95) ,type = 'l',lt =2,lw = 4)
points(muWeb[1],muWeb[2],cex = 2,pch = 19)
lines( ellipse::ellipse( covStox, centre = muStox,level = 0.95) , col="blue",type = 'l',lt = 2,lw = 4)
points(muStox[1],muStox[2],cex = 2,pch = 18, col = 'blue')
lines( ellipse::ellipse( covJoint, centre = muJoint,level = 0.95) , col="brown",type = 'l',lt = 3,lw = 4)
points(muJoint[1],muJoint[2],cex = 2,pch = 20, col = 'brown')
lines( ellipse::ellipse( covJointCov, centre = muJointCov,level = 0.95) , col="green",type = 'l',lw = 4)
points(muJointCov[1],muJointCov[2],cex = 2,pch = 15, col = 'green')
lines( ellipse::ellipse( covSepCov, centre = muSepCov,level = 0.95) , col="red",type = 'l',lw = 4)
points(muSepCov[1],muSepCov[2],cex = 2,pch = 17, col = 'red')

legend("topright", inset=.01,
#       c("Assessment in year 2023", "Proposed index, no data on covariance", expression(Proposed~index~with~ Psi ), expression(Proposed~index~with~ Psi^fixedALK ),"Official index"),
       c("Assessment in year 2023", "Proposed index, no data on covariance", expression(Proposed~index~with~cov(bold(I)[y]) ), expression(Proposed~index~with~cov(bold(I)[y]^fixedALK) ),"Official index"),
       pch = c(19,20,15,17,18), col = c("black", "brown","green","red","blue"), cex = 0.95, lwd = 2, lty = c(2,3,1,1,2))
dev.off()

##########################################################################
#Make table with variance and scaling parameter###########################
sdIndex = unique(fit$conf$keyVarObs[2,])+1
sdIndex = sdIndex[sdIndex>0]
logSdStox = c(fit$pl$logSdLogObs[sdIndex],fit$pl$logSdLogObs[sdIndex]-2*fit$plsd$logSdLogObs[sdIndex], fit$pl$logSdLogObs[sdIndex]+2*fit$plsd$logSdLogObs[sdIndex])
logSdJoint = c(fitJoint$pl$logSdLogObs[sdIndex],fitJoint$pl$logSdLogObs[sdIndex]-2*fitJoint$plsd$logSdLogObs[sdIndex], fitJoint$pl$logSdLogObs[sdIndex]+2*fitJoint$plsd$logSdLogObs[sdIndex])

varTab = exp(logSdStox)
varTab = rbind(varTab,exp(logSdJoint))
rownames(varTab) = c("Official index","Proposed index")
colnames(varTab) = c("$\\sigma_{SAM}$", "Q025", "Q0975")
varTab = round(varTab,2)
varTab

logJointCovScale = c(fitJointCov$pl$logSdLogObs[sdIndex],fitJointCov$pl$logSdLogObs[sdIndex]-2*fitJointCov$plsd$logSdLogObs[sdIndex], fitJointCov$pl$logSdLogObs[sdIndex]+2*fitJointCov$plsd$logSdLogObs[sdIndex])
logSepCovScale = c(fitSepCov$pl$logSdLogObs[sdIndex],fitSepCov$pl$logSdLogObs[sdIndex]-2*fitSepCov$plsd$logSdLogObs[sdIndex], fitSepCov$pl$logSdLogObs[sdIndex]+2*fitSepCov$plsd$logSdLogObs[sdIndex])

varScalingTab = exp(logJointCovScale)
varScalingTab = rbind(varScalingTab,exp(logSepCovScale))
rownames(varScalingTab) = c("Proposed index","Proposed index with fixed ALK")
colnames(varScalingTab) = c("$k$", "Q025", "Q0975")
varScalingTab

combinedTab = as.data.frame(matrix(0,4,3))
colnames(combinedTab) = c( "$\\sigma_{SAM}$", "$k_{SAM}$", "$k_{SAM}*mean(\\sqrt{diag(\\Phi)})$")

rownames(combinedTab) = c("Official index", "$\\bf{I}_y$",
                    "$\\bf{I}_y$ with cov$({\\bf I}_y)$",
                    "$\\bf{I}_y$ with cov$({\\bf I}^{fixedALK}_y)$")

varScalingTab = round(varScalingTab,2)
combinedTab[1,1] = paste0( varTab[1,1] , " (",varTab[1,2],",",varTab[1,3], ")")
combinedTab[2,1] = paste0( varTab[2,1] , " (",varTab[2,2],",",varTab[2,3], ")")
combinedTab[1,2] = "-"
combinedTab[2,2] = "-"
combinedTab[1,3] = "-"
combinedTab[2,3] = "-"
combinedTab[3,2] = paste0( varScalingTab[1,1] , " (",varScalingTab[1,2],",",varScalingTab[1,3], ")")
combinedTab[4,2] = paste0( varScalingTab[2,1] , " (",varScalingTab[2,2],",",varScalingTab[2,3], ")")
combinedTab[3,3] = round(mean(varScalingTab[1,1]*sqrt(varJoint)),2)
combinedTab[4,3] = round(mean(varScalingTab[2,1]*sqrt(varSep)),2)
combinedTab[3,1] = "-"
combinedTab[4,1] = "-"

matQ = xtable::xtable(combinedTab, type = "latex")
print(matQ, sanitize.text.function = function(x) {x},
      file = "tables/varSAMTab.tex",
      floating = FALSE)



################################################################################
#Plot relative difference in variance estimates.
pdf("figures/relativeSdChange.pdf", width = 7,height = 5)
varChangeRelAll = (sqrt(varJoint)-sqrt(varSep)) /sqrt(varSep)
varChangeRel = varChangeRelAll

matplot(1+varChangeRel, xaxt = "n", main = "Relative change in standard deviation",ylab = "Scaling factor",xlab = "Year",
        cex.main = 1.8,cex.lab = 1.5, lwd=1.8, type = "l", las=1)
abline(h = 1)
text(0.75,1+varChangeRel[1,-ncol(varChangeRel)],  c(3:9), , col=1:6, cex = 1.2)
text(0.75,1+varChangeRel[1,ncol(varChangeRel)], expression(10^"+"), bg="white", col=2, cex = 1.2)
axis(side =1, at=1:27, labels=1994:2020)
#abline(h =exp(fitSepCov$pl$logSdLogObs[4]),lwd = 2,lty = 2)
dev.off()


#Plot of raw covariances, in supplementary.
load("data/NEAhaddockAssessment/cov_indexJoint.Rda")
library(ellipse)
ccolors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white","#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")
#pdf("figures/correlation.pdf", width = 5,height = 7)
pdf("figures/correlationMGWG.pdf", width = 6,height = 5)
mfrow = c(4,7)
par(oma=c(1,1,1,1),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
for(i in 1:length(covYears)){
  covYearsTmp = covYears[[i]][-1,-1]
  xx = cov2cor(covYearsTmp)
  plotcorr(xx, col = ccolors[5*xx+6],main = paste0("Year ", 1993 + i), mar = c(0,0,0.7,0.7))
}
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.3, y = seq(0,1,l=5), labels = round(seq(-1,1,l=5),2))
rasterImage(rev(ccolors), 0, 0, 1,1)
dev.off()



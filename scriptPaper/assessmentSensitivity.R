library(stockassessment)
library(scales)
load("data/NEAhaddockAssessment/NEA_haddock_2023_officiall.Rds")
fitWeb<-stockassessment:::refit(fitWeb) #need local version

####################################################Assessment in paper######################################
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
load("data/NEAhaddockAssessment/cov_indexJoint.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
dim = dim(surveys$`BS-NoRu-Q1(BTr)`)
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

conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
par<-defpar(dat,conf)
fitJointCov<-sam.fit(dat,conf,par)

####################################################Finer spatial resolution######################################
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJointMeshDetailed.txt", sep = " "))[,-1]
load("data/NEAhaddockAssessment/cov_indexJointMeshDetailed.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
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
fitJointCovMeshDetailed<-sam.fit(dat,conf,par)

###No discretication in length in latent space
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJointDL1.txt", sep = " "))[,-1]
load("data/NEAhaddockAssessment/cov_indexJointDL1.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
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
fitJointCovDL1<-sam.fit(dat,conf,par)

###################################################No covariates############################################
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJointNoCovariates.txt", sep = " "))[,-1]
load("data/NEAhaddockAssessment/cov_indexJointNoCovariates.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
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

conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
par<-defpar(dat,conf)
fitJointCovNoCovariates<-sam.fit(dat,conf,par)


###################################################No random walk############################################
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJointNoRW.txt", sep = " "))[,-1]
load("data/NEAhaddockAssessment/cov_indexJointNoRW.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
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

conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
par<-defpar(dat,conf)
fitJointCovNoRW<-sam.fit(dat,conf,par)


###################################################No spatio-temporal structure#############################
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJointNoST.txt", sep = " "))[,-1]
load("data/NEAhaddockAssessment/cov_indexJointNoST.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
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

conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
par<-defpar(dat,conf)
fitJointCovNoST<-sam.fit(dat,conf,par)


###################################################Finer integration resolution#############################
indexJoint= as.matrix(read.table("data/NEAhaddockAssessment/indexJointMoreIntPoints.txt", sep = " "))[,-1]
load("data/NEAhaddockAssessment/cov_indexJointMoreIntPoints.Rda")
covJoint = lapply(covYears,function(f){
  cov = f[-1,]
  cov = cov[,-1]
  cov
} )
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

conf = loadConf(dat, "data/NEAhaddockAssessment/conf.cfg")
par<-defpar(dat,conf)
fitJointCovMoreIntPoints<-sam.fit(dat,conf,par)


#################################################################################################
#Plot them all
#################################################################################################
colSet <- c("#332288"  , "#88CCEE"  , "#44AA99"  , "#117733"  , "#999933"  , "#DDCC77"  , "#661100"  , "#CC6677"  , "#882255"  , "#AA4499")

#No ST
#png("figures/sensitivityNoST.png",height = 1000,width = 900)
pdf("figures/sensitivityNoST.pdf",height = 8,width = 15)
par(mfrow = c(1,2),oma = c(1,1,1,1))
leg = c("In paper","No spatio-temporal structures")
cex = 1;cex.main = 2
ssbplot(c(fitJointCov, fitJointCovNoST),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "SSB",cex.main =cex.main,)
legend(x = 1995, y = 60000,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]),  bty="n",cex = cex)
legend(x = 1995, y = 60000,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted",  bty="n",cex = cex)

fbarplot(c(fitJointCov,fitJointCovNoST),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "F-bar",cex.main =cex.main)
legend(x = 1995, y = 0.1,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 0.1,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)
dev.off()

#No RW
pdf("figures/sensitivityNoRW.pdf",height = 8,width = 15)
par(mfrow = c(1,2))
leg = c("In paper","No random walk")
cex = 1.5;cex.main = 2
ssbplot(c(fitJointCov, fitJointCovNoRW),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "SSB",cex.main =cex.main)
legend(x = 1995, y = 60000,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 60000,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)

fbarplot(c(fitJointCov,fitJointCovNoRW),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "F-bar",cex.main =cex.main)
legend(x = 1995, y = 0.1,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 0.1,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)
dev.off()


#No covariates
pdf("figures/sensitivityNoCovariates.pdf",height = 8,width = 15)
par(mfrow = c(1,2))
leg = c("In paper","No covariates")
cex = 1.5;cex.main = 2
ssbplot(c(fitJointCov, fitJointCovNoCovariates),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "SSB",cex.main =cex.main)
legend(x = 1995, y = 60000,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 60000,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)

fbarplot(c(fitJointCov,fitJointCovNoCovariates),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "F-bar",cex.main =cex.main)
legend(x = 1995, y = 0.1,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 0.1,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)
dev.off()


#More detailed latent resolution
pdf("figures/sensitivityDetailedMesh.pdf",height = 8,width = 15)
par(mfrow = c(1,2))
leg = c("In paper","Detailed spatial mesh")
cex = 1.5;cex.main = 2
ssbplot(c(fitJointCov, fitJointCovMeshDetailed),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "SSB",cex.main =cex.main)
legend(x = 1995, y = 60000,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 60000,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)

fbarplot(c(fitJointCov,fitJointCovMeshDetailed),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "F-bar",cex.main =cex.main)
legend(x = 1995, y = 0.1,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 0.1,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)
dev.off()

pdf("figures/sensitivityDL1.pdf",height = 8,width = 15)
par(mfrow = c(1,2))
leg = c("In paper","5cm resolution in latent length space")
cex = 1.5;cex.main = 2
ssbplot(c(fitJointCov, fitJointCovDL1),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "SSB",cex.main =cex.main)
legend(x = 1995, y = 60000,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 60000,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)

fbarplot(c(fitJointCov,fitJointCovDL1),addCI = rep(TRUE,2),xlim = c(1994,2020), cex.lab =1.3, cex = 2, main = "F-bar",cex.main =cex.main)
legend(x = 1995, y = 0.1,legend=leg, lwd=3, col=c(par("col"),colSet[((2:2)-1)%%length(colSet)+1]), bty="n",cex = cex)
legend(x = 1995, y = 0.1,legend=leg, lwd=2, col=rep("black", length(leg)), lty="dotted", bty="n",cex = cex)
dev.off()




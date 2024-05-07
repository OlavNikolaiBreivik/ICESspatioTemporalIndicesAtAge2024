library(spatioTemporalALK)
library(spatioTemporalIndices)
library(RColorBrewer)
library(xtable)
#Make figures
load("results/run.Rdata") #File created in run.R script.

############################################################################################
###################Plot fixed effects#######################################################
pdf("figures/sunEffect.pdf")
plotSunAlt(run)
dev.off()

pdf("figures/depthEffect.pdf")
plotResults(run,what = "depth")
dev.off()



############################################################################################
#################mesh and survey map####################################################
newmap <- maps::map("world", c("Norway","Sweden","Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
mapTmp = data.frame(newmap$x,newmap$y)
mapTmp[which(is.na(mapTmp[,1])),] = 3.141592 #Need something different from NA
names(mapTmp) = c("X","Y")
attr(mapTmp, "projection") = "LL"
attr(mapTmp, "zone") = run$conf_l$zone
ddpcr::quiet(mapTmp <- PBSmapping::convUL(mapTmp))
colnames(mapTmp) = c("UTMX", "UTMY")
mapTmp[which(is.na(newmap$x)),] = NA

pdf("figures/spatialMeshUTM.pdf")
plot(-1,-1, xlim = c(-250,1300), ylim = c(7500,9000),
     xlab = "East (km)", ylab = "North (km)",
     cex.lab = 1.5,cex.main = 2)
polygon(mapTmp,col = 'lightgrey')
plot(attributes(run$data)$meshS, add = TRUE)
title("Spatial mesh",cex.main  = 2)
dev.off()

pdf("figures/surveyMap.pdf")
plot(-10,-10,ylim = c(66,80.3),xlim = c(8,58),
     cex.lab = 1.5,cex.axis = 1.8, cex.main=2, cex.sub= 1.8,
     main = "Stations in year 2015 and 2016",
     xlab = "Longitude",
     ylab = "Latitude")
plot(run$conf_l$strata,add = TRUE,lwd = 2,col = 0)
points(newmap,type = 'l')
points(attributes(run$data)$locObsLatLon[which(attributes(run$data)$year == 2015),],cex = 1.5, col = 'blue')
points(attributes(run$data)$locObsLatLon[which(attributes(run$data)$year == 2016),],cex = 1.5,col = 'red', pch = 5,)
polygon(newmap,col = 'lightgrey')
dev.off()



############################################################################################
#############ALK plot in 2015 and 2016######################################################
col = brewer.pal(n = 10, name = 'Set1')
col2 = brewer.pal(n = 10, name = 'Paired')
col[6] = col2[3]

yearsPlot = c(2015,2016)
lab = FALSE
lengthInt = seq(min(run$data$length), max(run$data$length), length.out = 200)
xlim = c(min(lengthInt)-3, max(lengthInt))
for(year in yearsPlot){
  pdf(paste0("figures/alk", year,".pdf"))
  yearIndex = which(run$conf_alk$years == year)-1
  nAges = run$data$ageRange[2]- run$data$ageRange[1] + 1
  linPredMatrix = matrix(0,length(lengthInt), nAges-1)
  for(a in 1:(nAges-1)){
    linPredMatrix[,a] = run$pl$beta0_alk[yearIndex+1, a] + run$pl$betaLength_alk[a]*lengthInt
  }
  ALK = matrix(0,length(lengthInt), nAges)
  for(a in 1:nAges){
    probLess = rep(0,dim(ALK)[1])
    if(a>1){
      for(b in 1:(a-1)){
        tmp = ALK[,b];
        probLess = probLess + tmp;
      }
    }
    if(a <(nAges)){
      tmp2 = linPredMatrix[,a];
      ALK[,a] = plogis(tmp2)*(1-probLess);
    }else{
      ALK[,nAges] = (1-probLess);
    }
  }
  plot(lengthInt,ALK[,1], type = "l", ylab = "Probability", xlab = "Length",
       ylim = c(0,1), main = paste0("ALK year ",year),
       col = col[1],lwd = 3, cex.main = 2,cex.lab = 1.5, xlim = xlim)
  for(l in 2:nAges){
      lines(lengthInt,ALK[,l],col = col[l],lwd = 3)
  }
  #Include labels
  if(lab){
    ages = run$conf_alk$minAge:run$conf_alk$maxAge
    nAges =  length(ages)
    for(i in 0:nAges){
      text(7, i/(nAges+1)+0.05, paste("a =", i + ages[1]-1), col = col[i+1],cex = 1.3)
    }
  }
  lab = TRUE
  dev.off()
}

############################################################################################
#Contour plot length 40 in years 2015 and 2016##############################################
ll = run$obj$report()$lengthIndexDetailed
aa = run$obj$report()$ageIndexDetailed
length = 40
age = 5
col = brewer.pal( 9,"YlOrRd")
nCol = 20
col = colorRampPalette(col)(nCol)
xlim = c(-200,1170)
ylim = c(7400,9050)
l = which(run$conf_l$lengthGroups== length)
years = which(run$conf_l$years %in% yearsPlot)
zlim = c(0,log(max(ll[years,l,])))
breaks = c(-1,seq(0,zlim[2],by  = 0.1))
lab = FALSE
for(year in years){
  pdf(paste0("figures/length",length,"year", 1993 +year,".pdf"))
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))
  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
        z[xC,yC] = log(ll[year,l,which(run$data$xInt ==xx & run$data$yInt ==yy)])
      }
      yC = yC+1
    }
    xC = xC +1
  }
  image(x,y, z,col =  col,
             legend.width = 2, cex = 1.6,zlim = zlim, xlim = xlim,ylim = ylim,
             axis.args = list(cex.axis = 1.3),xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
             main = paste0("Log CPUE at length ", length, " cm in year ",1993 + year ),
             cex.lab = 1.5,cex.main = 2)


  polygon(mapTmp,col = 'lightgrey')
  cc = log(run$data$fishObsMatrix[which(attributes(run$data)$year == (year +1993)),l] / run$data$dist[which(attributes(run$data)$year == (year +1993))])
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),],cex = cc)
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),][which(cc< -9999),],pch = 4)

  if(year ==years[2]){
    image.plot(zlim = zlim, col = col,#colorRampPalette(c("white", "yellow", "red"))(length(breaks)-1),
               legend.only = TRUE, horizontal = FALSE, axis.args = list(cex.axis = 2),
               main = "Log CPUE Scale", legend.mar = 6,legend.width = 2)
  }

  dev.off()
}

############################################################################################
#Contour plot age 5 in years 2015 and 2016##################################################
Ahaul = run$obj$report()$ALK_haul
zlim = c(0,log(max(aa[years,age-1,])))
breaks = c(-1,seq(0,zlim[2],by  = 0.1))
col = brewer.pal( 9,"YlOrRd")
nCol = 20
col = colorRampPalette(col)(nCol)
xlim = c(-200,1170)
ylim = c(7400,9050)
for(year in years){
  pdf(paste0("figures/age",age,"year",1993 + year,".pdf"))
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))
  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
        z[xC,yC] = log(aa[year,age-1,which(run$data$xInt ==xx & run$data$yInt ==yy)])
      }
      yC = yC+1
    }
    xC = xC +1
  }
  image(x,y, z,col =  col,
             legend.width = 2, cex = 1.6,zlim = zlim,xlim = xlim,ylim = ylim,
              xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
             main = paste0("Log CPUE at age ", age, " in year ",1993 + year ),
             cex.lab = 1.5,cex.main = 2)

  polygon(mapTmp,col = 'lightgrey')
  cc = run$data$fishObsMatrix[which(attributes(run$data)$year == (year +1993)),]
  AA = matrix(NA,dim(cc)[1],dim(Ahaul)[2])
  for(sss in 1:dim(cc)[1]){
    for(aaa in 1:dim(Ahaul)[2]){
      AA[sss,]= cc[sss,]%*%Ahaul[,,sss,year]
    }
  }
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),],cex = log(AA[,age-1]/ run$data$dist[which(attributes(run$data)$year == (year +1993))] ))
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),][which(AA[,age-1]<1),],pch = 4)

  if(year ==years[2]){
    image.plot(zlim = zlim, col = col,
               legend.only = TRUE, horizontal = FALSE, axis.args = list(cex.axis = 2),
               main = "Log CPUE Scale", legend.mar = 6,legend.width = 2)
  }

  dev.off()
}

############################################################################################
#Abundance boble plot#######################################################################
minAge = 3
maxAgeInternal = run$conf_alk$maxAge - minAge + 2
nYears = length(run$conf_l$years)
ageLogIndex = run$rl$logAgeIndex[,2:maxAgeInternal]
ageLogIndexSd =  run$rlSd$logAgeIndex[,2:maxAgeInternal]
reduce = 1.5
pdf("figures/cohort.pdf",height = 12,width = 14)
tmp = apply(ageLogIndex,2,max)
tmp2 =  apply(ageLogIndex,2,min)
increase = 10
indexSdTmp = ageLogIndexSd
maxVar = max(indexSdTmp) + 0.01
nCol = 100
col = rev(brewer.pal( 9,"YlOrRd"))
col = colorRampPalette(col)(nCol)

par(oma=c(4,4,3,3),mar=c(0.5,0.5,0.5,0.0))
plot(rep(1, dim(ageLogIndex)[2]),2:run$data$ageRange[2],cex = (ageLogIndex[1,]-0.99*tmp2)/(tmp-tmp2) *increase +0.5, #Add a small value to make the smallest indices visible.
     xlim =c(1,dim(ageLogIndex)[1]+2),pch = 16,
     yaxt="n",xaxt = "n",xlab = "Year",ylab = "Age",cex.lab = 1.5,
     col = col[indexSdTmp[1,]/maxVar*nCol])
for(i in 2:nYears){
  points(rep(i, dim(ageLogIndex)[2]),2:run$data$ageRange[2],
         cex = (ageLogIndex[i,] -0.99*tmp2)/(tmp-tmp2) *increase+0.5,
         pch = 16,col = col[indexSdTmp[i,]/maxVar*nCol])
}
axis(1,1:length(run$conf_l$years),run$conf_l$years,cex.axis = 1.5)
axis(2,2:10,3:11,cex.axis = 1.5)
x = sort(unique(run$data$xInt))
y = sort(unique(run$data$yInt))
z = matrix(NA,nrow = length(x), ncol = length(y))
image.plot(x,y = 1, z,col =  col[1:nCol],
           legend.width = 2, cex = 1.2,zlim = c(min(indexSdTmp),max(indexSdTmp^2)),
           axis.args = list(cex.axis = 1.5),
           main = "Var" ,legend.cex = 1,legend.only = TRUE,add = TRUE,
           smallplot = c(0.95, 1, .05, .95))
mtext(text="Year",cex=2.5,side=1,line=2.5,outer=TRUE)
mtext(text="Age",cex=2.5,side=2,line=2,outer=TRUE)
mtext("Age index plot", outer=TRUE,  cex=3.5, line=0.6)
dev.off()
############################################################################################

############################################################################################
#Plot mean age at length 40 cm in 2015 and 2016##############################################
alk = run$obj$report()$ALK_int
length = 40
zlim = c(4,5)
l = which(run$conf_l$lengthGroups== length)
yearsPlot = c(2015,2016)
yIndex = which(run$conf_l$years %in% yearsPlot)
for(year in yIndex){
  pdf(paste0("figures/meanAgeAtLength", length,"Year",1993 + year,".pdf"))
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))
  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      ii = which(run$data$xInt ==xx & run$data$yInt ==yy)
      if(length(ii)>0){
        z[xC,yC] =  sum(alk[l,,ii,year]*(2:10))
        #          log(aa[year,age-1,which(run$data$xInt ==xx & run$data$yInt ==yy)])
      }
      yC = yC+1
    }
    xC = xC +1
  }
  breaks = c(-1,seq(zlim[1],zlim[2],by  = 0.01))
  col = brewer.pal(9, "YlOrRd")
  image(x,y, z,col =  colorRampPalette(col)(length(breaks)-1), cex = 1.6,zlim = zlim,
             xlab = "East (km)", ylab = "North (km)",
             main = "",
             cex.lab = 1.4,cex.main = 1.5)
  title(paste0("Mean age at ",length,"cm in year ", 1993 + year),cex.main=2)
  polygon(mapTmp,col = 'lightgrey')
  cc = log(run$data$fishObsMatrix[which(attributes(run$data)$year == (year +1993)),l]/run$data$dist[which(attributes(run$data)$year == (year +1993))])
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),],cex = cc)
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),][which(cc< -9999),],pch = 4)
  dev.off()
}

#Plot average mean age at length 40 in all years
pdf("figures/meanAgeAtLengthSpatial.pdf")
x = sort(unique(run$data$xInt))
y = sort(unique(run$data$yInt))
z = matrix(NA,nrow = length(x), ncol = length(y))
xC = 1
for(xx in x){
  yC = 1
  for(yy in y){
    ii = which(run$data$xInt ==xx & run$data$yInt ==yy)
    if(length(ii)>0){
      z[xC,yC] = 0
      for(year in  1:length(run$conf_l$years)){
        z[xC,yC] = z[xC,yC]+  sum(alk[l,,ii,year]*(2:10))/length(run$conf_l$years)

      }
    }
    yC = yC+1
  }
  xC = xC +1
}
breaks = c(-1,seq(zlim[1],zlim[2],by  = 0.01))
col = brewer.pal(9, "YlOrRd")
image.plot(x,y, z,col = colorRampPalette(col)(length(breaks)-1),
           legend.width = 2, cex = 1.6,zlim = zlim,
           xlab = "East (km)", ylab = "North (km)",
           main = "",
           cex.lab = 1.4,cex.main = 1.5)

title(paste0("Mean age at ",length,"cm in years 1994-2020"),cex.main=2)
polygon(mapTmp,col = 'lightgrey')
dev.off()


############################################################################################
##############################SUPPLEMENTARY#################################################
############################################################################################

#Parameter table
source("scriptPaper/utils.R")
parameter = partable(run)
colnames(parameter) = c("Estimate","$Q_{0.025}$","$Q_{0.975}$")
rownames(parameter)[1] = "$\\rho_{\\alpha,l}$"
rownames(parameter)[2] = "$\\rho_{\\gamma,l}$"
rownames(parameter)[3] = "$\\rho_{\\psi,l}$"
rownames(parameter)[4] = "$\\rho_{\\gamma,y}$"
rownames(parameter)[5] = "$\\sigma_{\\alpha}$"
rownames(parameter)[6] = "$\\sigma_{\\gamma}$"
rownames(parameter)[7] = "$\\sigma_{\\psi}$"
rownames(parameter)[8] = "$\\log\\kappa_{\\alpha}$"
rownames(parameter)[9] = "$\\log\\kappa_{\\gamma}$"
parameter[8,] = log(parameter[8,])
parameter[9,] = log(parameter[9,])
rownames(parameter)[10] = "$\\sigma_{\\beta}$"
rownames(parameter)[11] = "$\\sigma_{\\beta^{(alk)}}$"
rownames(parameter)[12] = "$\\sigma_{\\alpha}^{(alk)}$"
rownames(parameter)[13] = "$\\sigma_{\\gamma}^{(alk)}$"
rownames(parameter)[14] = "$\\log\\kappa_{\\alpha}^{(alk)}$"
rownames(parameter)[15] = "$\\log\\kappa_{\\gamma}^{(alk)}$"
parameter[14,] = log(parameter[14,])
parameter[15,] = log(parameter[15,])
rownames(parameter)[16] = "$\\log\\lambda$"

for(i in 1:dim(parameter)[1]){
  parameter[i,1] = paste0(round(as.numeric(parameter[i,1]),3), " (",round(as.numeric(parameter[i,2]),3), ",",round(as.numeric(parameter[i,3]),3),")")
}
parameter = parameter[,-2]
parameter = data.frame(parameter[,-2])
colnames(parameter) = "Estimate"

parameter$Parameter = rownames(parameter)

parPrint = data.frame(rownames(parameter)[1:8])
parPrint$Estimate = parameter$Estimate[1:8]
parPrint$par2 = c(rownames(parameter)[9:16])
parPrint$Estimate2 = c(parameter$Estimate[9:16])

names = c("Parameter","Estimate","Parameter","Estimate")
parPrint = rbind(names,parPrint)


pp = xtable(parPrint, type = "latex",digits = 3,lable = NULL)
print(pp, sanitize.text.function = function(x) {x},
      file = "tables/partab.tex",
      floating = FALSE,include.rownames=FALSE, include.colnames = FALSE)




#Plot mean age at age 40 in all years#######################################################
alk = run$obj$report()$ALK_int
length = 40
l = which(run$conf_l$lengthGroups== length)
pdf(paste0("figures/meanAgeAllYearsLength", length,".pdf"),height = 18,width = 14)
mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
yaxt = rep(c("s",rep("n",4)),6)
xaxt = rep("n",27)
xaxt[1:5] = "s"
xaxt= rev(xaxt)
col = brewer.pal(n = 10, name = 'Paired')
zlim = c(3.5,5.9)
for(year in 1:27){
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))
  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      ii = which(run$data$xInt ==xx & run$data$yInt ==yy)
      if(length(ii)>0){
        z[xC,yC] =  sum(alk[l,,ii,year]*(2:10))
      }
      yC = yC+1
    }
    xC = xC +1
  }
  breaks = c(-1,seq(zlim[1],zlim[2],by  = 0.1))
  col = brewer.pal(9, "YlOrRd")

  image(x,y, z,col =  colorRampPalette(col)(length(breaks)-1), cex = 1.6,zlim = zlim,
        xlab = "Eastern direction (km)", ylab = "Northern direction (km)",
        axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
        main = "",
        cex.lab = 1.4,cex.main = 1.5)
  title(1993 + year,line = -2,cex.main=2)
  polygon(mapTmp,col = 'lightgrey')
  cc = log(run$data$fishObsMatrix[which(attributes(run$data)$year == (year +1993)),l])
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),],cex = cc)
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),][which(cc< -9999),],pch = 4)
  if(year==27){
    plot(1,1,axes=F,ylim = c(-99,-98))
    image.plot(x,y, z,col =  colorRampPalette(col)(length(breaks)-1),
               legend.width = 2, cex = 1.6,zlim = zlim,
               axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
               main = "" ,legend.cex = 2,legend.only = TRUE,
               smallplot = c(0.1, 0.15, .1, .9))

  }
}
mtext(text="Eastern direction (km)",cex=2,side=1,line=2,outer=TRUE)
mtext(text="Northern direction (km)",cex=2,side=2,line=2.4,outer=TRUE)
mtext(paste0("Mean age at length = ",length," in each year "  ), outer=TRUE,  cex=2, line=0.5)
dev.off()

############################################################################################
#############ALK all years##########################################
col = brewer.pal(n = 10, name = 'Set1')
col2 = brewer.pal(n = 10, name = 'Paired')
col[6] = col2[3]
yearsPlot = 1994:2020
lab = FALSE
lengthInt = seq(min(run$data$length), max(run$data$length), length.out = 200)
xlim = c(min(lengthInt)-3, max(lengthInt))
yaxt = rep(c("s",rep("n",4)),6)
xaxt = rep("n",27)
xaxt[1:5] = "s"
xaxt= rev(xaxt)
pdf(paste0("figures/alkAllYears.pdf"))
mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
for(year in yearsPlot){
  yearIndex = which(run$conf_alk$years == year)-1
  nAges = run$data$ageRange[2]- run$data$ageRange[1] + 1
  linPredMatrix = matrix(0,length(lengthInt), nAges-1)
  for(a in 1:(nAges-1)){
    linPredMatrix[,a] = run$pl$beta0_alk[yearIndex+1, a] + run$pl$betaLength_alk[a]*lengthInt
  }
  ALK = matrix(0,length(lengthInt), nAges)
  for(a in 1:nAges){
    probLess = rep(0,dim(ALK)[1])
    if(a>1){
      for(b in 1:(a-1)){
        tmp = ALK[,b];
        probLess = probLess + tmp;
      }
    }
    if(a <(nAges)){
      tmp2 = linPredMatrix[,a];
      ALK[,a] = plogis(tmp2)*(1-probLess);
    }else{
      ALK[,nAges] = (1-probLess);
    }
  }
  plot(lengthInt,ALK[,1], type = "l", ylab = "Probability", xlab = "Length",
       ylim = c(0,1), main = "",yaxt = yaxt[year-1993],xaxt = xaxt[year-1993],
       col = col[1],lwd = 1, cex.main = 1,cex.lab = 1, xlim = xlim)
  title(year,line = -1,cex.main=1)
  for(l in 2:nAges){
    lines(lengthInt,ALK[,l],col = col[l],lwd = 1)
  }#Include labels
  if(lab){
    ages = run$conf_alk$minAge:run$conf_alk$maxAge
    nAges =  length(ages)
    for(i in 0:nAges){
      text(12.5, i/(nAges+1)+0.05, paste("a =", i + ages[1]-1), col = col[i+1],cex = 0.7)
    }
  }
  if(year==2019)lab = TRUE
}
dev.off()

############################################################################################
#############Contour plot CPUE length 40 all years##########################################
ll = run$obj$report()$lengthIndexDetailed
length = 40
l = which(run$conf_l$lengthGroups== length)
pdf(paste0("figures/densityLength", length,".pdf"),height = 18,width = 14)
zlim = c(0,log(max(ll[,l,])))
mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
yaxt = rep(c("s",rep("n",4)),6)
xaxt = rep("n",27)
xaxt[1:5] = "s"
xaxt= rev(xaxt)
col = brewer.pal( 9,"YlOrRd")
nCol = 20
col = colorRampPalette(col)(nCol)
for(year in 1:27){
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))
  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
        z[xC,yC] = log(ll[year,l,which(run$data$xInt ==xx & run$data$yInt ==yy)])
      }
      yC = yC+1
    }
    xC = xC +1
  }
  breaks = c(-1,seq(0,zlim[2],by  = 1))
  image(x,y, z,col =  col,
        legend.width = 2, cex = 1.6,zlim = zlim,
        axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
        main = "")

  title(1993 + year,line = -2,cex.main=2)
  polygon(mapTmp,col = 'lightgrey')
  cc = log(run$data$fishObsMatrix[which(attributes(run$data)$year == (year +1993)),l])
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),],cex = cc)
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),][which(cc< -9999),],pch = 4)
  if(year==27){
    plot(1,1,axes=F,ylim = c(-99,-98))
    image.plot(x,y, z,col =  col,
               legend.width = 2, cex = 1.6,zlim = zlim,
               axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
               main = "" ,legend.cex = 2,legend.only = TRUE,
               smallplot = c(0.1, 0.15, .1, .9))
  }
}
mtext(text="Eastern direction (km)",cex=2,side=1,line=2,outer=TRUE)
mtext(text="Northern direction (km)",cex=2,side=2,line=2.4,outer=TRUE)
mtext(paste0("Predicted log CPUE of length ", length, " cm NEA haddock"  ), outer=TRUE,  cex=2, line=0.5)
dev.off()

############################################################################################
#Contour plot age all years
aa = run$obj$report()$ageIndexDetailed
age = 5
pdf(paste0("figures/densityAge", age,".pdf"),height = 18,width = 14)
zlim = c(0,log(max(aa[,age-1,]))/1)
mfrow = c(6,5)
par(oma=c(5,5,4,0.5),mar=c(0.5,0.5,0.0,0.0),mfrow = mfrow)
yaxt = rep(c("s",rep("n",4)),6)
xaxt = rep("n",27)
xaxt[1:5] = "s"
xaxt= rev(xaxt)
col = brewer.pal( 9,"YlOrRd")
nCol = 20
col = colorRampPalette(col)(nCol)
for(year in 1:27){
  x = sort(unique(run$data$xInt))
  y = sort(unique(run$data$yInt))
  z = matrix(NA,nrow = length(x), ncol = length(y))
  xC = 1
  for(xx in x){
    yC = 1
    for(yy in y){
      if(length(which(run$data$xInt ==xx & run$data$yInt ==yy))>0){
        z[xC,yC] = log(aa[year,age-1,which(run$data$xInt ==xx & run$data$yInt ==yy)])
      }
      yC = yC+1
    }
    xC = xC +1
  }
  breaks = c(-1,seq(0,zlim[2],by  = 1))
  image(x,y, z,col =  col,
        legend.width = 2, cex = 1.6,zlim = zlim,
        axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
        main = "")
  title(1993 + year,line = -2,cex.main=2)
  polygon(mapTmp,col = 'lightgrey')
  cc = log(run$data$fishObsMatrix[which(attributes(run$data)$year == (year +1993)),l])
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),],cex = cc)
  points(attributes(run$data)$locObs[which(attributes(run$data)$year == (year +1993)),][which(cc< -9999),],pch = 3)
  if(year==27){
    plot(1,1,axes=F,ylim = c(-99,-98))
    image.plot(x,y, z,col =  col,
               legend.width = 2, cex = 1.6,zlim = zlim,
               axis.args = list(cex.axis = 1.3),yaxt = yaxt[year],xaxt = xaxt[year],
               main = "" ,legend.cex = 2,legend.only = TRUE,
               smallplot = c(0.1, 0.15, .1, .9))
  }
}
mtext(text="Eastern direction (km)",cex=2,side=1,line=2,outer=TRUE)
mtext(text="Northern direction (km)",cex=2,side=2,line=2.4,outer=TRUE)
mtext(paste0("Predicted log CPUE of age ", age  , " NEA haddock"), outer=TRUE,  cex=2, line=0.5)
dev.off()




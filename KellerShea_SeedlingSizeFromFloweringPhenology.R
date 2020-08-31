#Supporting Information. Keller, J.A., and K. Shea. 2020. Warming and shifting phenology accelerate an invasive plant life cycle. Ecology

#Written using R version 3.6.3 (2020-02-29)

#Create a dataframe containing the size (z) and date of seed release (PlDate) of all newly recruited individuals
#Here we use simulated data for illustration
d.seedlings<-data.frame(z=NA, PlDate=rep(c( 180, 190, 200, 210, 220, 230, 240, 250), each=50))
head(d.seedlings)


for(i in 1:length(d.seedlings[,1])){
  d.seedlings$z[i]<--0.04*d.seedlings$PlDate[i]+10+rnorm(1,mean=0,sd=0.2)
}
plot(z~PlDate,data=d.seedlings)

sizes<-d.seedlings$z
dates<-d.seedlings$PlDate
head(sizes)
head(dates)
dateband<-10 #Set the bandwidth for the 2-d kernel fitting relating date of seed release and seedling size
sizeband<-0.5 #Set the bandwidth for the 2-d kernel fitting relating date of seed release and seedling size
m<-300 #The number of meshpoints used in the integral projection model
L<--1 #The lower size limit of the integral projection model. Here, size is ln(longest leaf length in centimeters)
U<-4.5 #The upper size limit of the integral projection model.Here, size is ln(longest leaf length in centimeters)

#Define the kernel relating date of flowering to date of seed release
#On each day following flowering, what is the probability that seed from that flower will be released?
days<-0:30
p.release<-dnorm(days, mean=15, sd=4)
p.release<-p.release/sum(p.release)
plot(p.release)
sum(p.release)

#load libraries
library(KernSmooth) #version 2.23-16
library(RSAGA) #version 1.3.0
library(ggplot2) #version 3.3.0

#Fit the 2-d kernel relating date of seed release and seedling sizes. A similar result could also be acheived using regression
kskern<-bkde2D(cbind.data.frame(dates,sizes),bandwidth = c(dateband,sizeband),gridsize = c(max(dates)-min(dates)+1,m), range.x=list(c(min(dates),max(dates)),c(L,U)))

#define the meshpoints of the integral projection model
h <- (U - L)/m
meshpts <- L + ((1:m) - 1/2) * h

#Convert the fitted 2-d kernel to a more easily worked with format
kskern.xyz<-grid.to.xyz(kskern$fhat,varname="density",colnames=c("z","PlDate","density"))
head(kskern.xyz)
kskern.xyz$z<-rep(meshpts,max(dates)-min(dates)+1)
kskern.xyz$PlDate<-rep(seq(min(dates),max(dates),by=1),each=m)

#standardize each date's probability density to sum to 1
for(i in min(kskern.xyz$PlDate):max(kskern.xyz$PlDate)){
  kskern.xyz$density[kskern.xyz$PlDate==i]<-kskern.xyz$density[kskern.xyz$PlDate==i]/sum(kskern.xyz$density[kskern.xyz$PlDate==i])
}

#And plot
p<-ggplot(kskern.xyz,aes(x=PlDate, y=z))
p+geom_raster(aes(fill=density))+geom_point(data=d.seedlings)

#If seeds may be released before the first tested planting date, we assume they behave the same as those released at the first tested date
#Build density estimates back to an early first flower date

fflower <- 170 #Set the date of first flowering to day of year 170

if(fflower<min(dates)){
  earlydays<-rep(fflower:min(dates),each=m)
  earlysizes<-rep(meshpts,min(dates)-fflower+1)
  earlydens<-rep(kskern.xyz$density[kskern.xyz$PlDate==min(dates)],min(dates)-fflower+1)
  early<-cbind.data.frame(earlydays,earlysizes,earlydens)
  colnames(early)<-c("PlDate","z","density")
  kskern.xyz<-rbind(kskern.xyz,early)
}
p<-ggplot(kskern.xyz,aes(x=PlDate, y=z))
p+geom_raster(aes(fill=density))+geom_point(data=d.seedlings)

#Weight densities for each date by the probability that seeds would be released on that date
release<-cbind.data.frame(seq(from=fflower,to=fflower+length(p.release)-1,by=1),p.release)
colnames(release)<-c("doy","p.release")
release
plot(p.release~doy,release)

if(max(dates)>max(release$doy)){
  fillreleasehigh<-cbind.data.frame(seq(from=max(release$doy), to=max(dates), by=1),rep(0, max(dates)-max(release$doy)+1))
  } else {
  fillreleasehigh<-cbind.data.frame(NA, NA)  
  }
colnames(fillreleasehigh)<-c("doy","p.release")

if(min(release$doy)>min(dates)){
  fillreleaselow<-cbind.data.frame(seq(from=min(dates), to=min(release$doy), by=1),rep(0, min(release$doy)-min(dates)+1))
  } else {
  fillreleaselow<-cbind.data.frame(NA, NA)  
  }
colnames(fillreleaselow)<-c("doy","p.release")

release<-rbind.data.frame(release,fillreleasehigh, fillreleaselow)
release<-na.omit(release)
plot(p.release~doy,release)

for(i in min(release$doy):max(release$doy)){
  kskern.xyz$density[kskern.xyz$PlDate==i]<-kskern.xyz$density[kskern.xyz$PlDate==i]*release$p.release[release$doy==i]
}

p<-ggplot(kskern.xyz,aes(x=PlDate, y=z))
p+geom_raster(aes(fill=density))

#Sum the probabilities for each size across dates and standardize
p.size<-rep(NA,length(meshpts))
for(i in 1:length(meshpts)){
  p.size[i]<-sum(kskern.xyz$density[kskern.xyz$z==meshpts[i]])
}
p.size<-p.size/(sum(p.size)*h)

plot(p.size~meshpts)
sum(p.size*h)

#Package this process in a function:

size_from_phenology<-function(sizes,dates,dateband,sizeband,m,L,U,fflower,p.release){
  kskern<-bkde2D(cbind.data.frame(dates,sizes),bandwidth = c(dateband,sizeband),gridsize = c(max(dates)-min(dates)+1,m), range.x=list(c(min(dates),max(dates)),c(L,U)))
  
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  kskern.xyz<-grid.to.xyz(kskern$fhat,varname="density",colnames=c("z","PlDate","density"))
  head(kskern.xyz)
  kskern.xyz$z<-rep(meshpts,max(dates)-min(dates)+1)
  kskern.xyz$PlDate<-rep(seq(min(dates),max(dates),by=1),each=m)
  
  for(i in min(kskern.xyz$PlDate):max(kskern.xyz$PlDate)){
    kskern.xyz$density[kskern.xyz$PlDate==i]<-kskern.xyz$density[kskern.xyz$PlDate==i]/sum(kskern.xyz$density[kskern.xyz$PlDate==i])
  }
  
  if(fflower<min(dates)){
    earlydays<-rep(fflower:min(dates),each=m)
    earlysizes<-rep(meshpts,min(dates)-fflower+1)
    earlydens<-rep(kskern.xyz$density[kskern.xyz$PlDate==min(dates)],min(dates)-fflower+1)
    early<-cbind.data.frame(earlydays,earlysizes,earlydens)
    colnames(early)<-c("PlDate","z","density")
    kskern.xyz<-rbind(kskern.xyz,early)
  }
  
  release<-cbind.data.frame(seq(from=fflower,to=fflower+length(p.release)-1,by=1),p.release)
  colnames(release)<-c("doy","p.release")
  
  if(max(dates)>max(release$doy)){
    fillreleasehigh<-cbind.data.frame(seq(from=max(release$doy), to=max(dates), by=1),rep(0, max(dates)-max(release$doy)+1))
  } else {
    fillreleasehigh<-cbind.data.frame(NA, NA)  
  }
  colnames(fillreleasehigh)<-c("doy","p.release")
  
  if(min(release$doy)>min(dates)){
    fillreleaselow<-cbind.data.frame(seq(from=min(dates), to=min(release$doy), by=1),rep(0, min(release$doy)-min(dates)+1))
  } else {
    fillreleaselow<-cbind.data.frame(NA, NA)  
  }
  colnames(fillreleaselow)<-c("doy","p.release")
  
  release<-rbind.data.frame(release,fillreleasehigh, fillreleaselow)
  release<-na.omit(release)
  
  for(i in min(release$doy):max(release$doy)){
    kskern.xyz$density[kskern.xyz$PlDate==i]<-kskern.xyz$density[kskern.xyz$PlDate==i]*release$p.release[release$doy==i]
  }
  p.size<-rep(NA,length(meshpts))
  for(i in 1:length(meshpts)){
    p.size[i]<-sum(kskern.xyz$density[kskern.xyz$z==meshpts[i]])
  }
  p.size<-p.size/(sum(p.size)*h)
  return(p.size)
}

#Plot for various dates of first flowering
par(mfrow=c(3,1))
plot(size_from_phenology(sizes=sizes,dates=dates,m=300,L=-1,U=4.5,sizeband=0.3,dateband=10,p.release=p.release, fflower=180) ~ meshpts, ylim=c(0,0.9), xlab="ln(longest leaf length)", ylab="probability density", main="180")
abline(v=-0.04*(180+15)+10) #expected value based on mean value of simulated data at the date of flowering plus 15 days (which is the most probable number of days in between flowering and seed release)
plot(size_from_phenology(sizes=sizes,dates=dates,m=300,L=-1,U=4.5,sizeband=0.3,dateband=10,p.release=p.release, fflower=200) ~ meshpts, ylim=c(0,0.9), xlab="ln(longest leaf length)", ylab="probability density", main="200")
abline(v=-0.04*(200+15)+10)
plot(size_from_phenology(sizes=sizes,dates=dates,m=300,L=-1,U=4.5,sizeband=0.3,dateband=10,p.release=p.release, fflower=220) ~ meshpts, ylim=c(0,0.9), xlab="ln(longest leaf length)", ylab="probability density", main="220")
abline(v=-0.04*(220+15)+10)


library(ggplot)
library(randomForest)
library(spam)
library(viridis)
library(raster)

################################################################
#  Set up
################################################################

# Coordinates of location
lng<-153.3150
lat<--28.7270


# Files with day temperature, night temperature and precipitation
flDT<-"NSW_Daytime_temp.csv"
flNT<-"NSW_Nighttime_temp.csv"
flPR<-"NSW_Precipitation.csv"


# Location of Landsat tiles
dirL8<-"NSW_L8"
# Prefix of RS data
prefL8<-"NSW_L8_Bands_"

# Divide data into RF training and forecasting
train.date<-as.Date("2014-04-11")
test.date<-as.Date("2015-11-01")

################################################################
#  Functions
################################################################

smooth_Whittaker<-function(y, lambda=10^1.3){
  m = length(y)
  w<-rep(1, m)
  D = diff(diag.spam(m), diff = 2)
  W = diag.spam(w)
  z = solve(W + lambda * t(D) %*% D, w * y)
  return(z)
}

kelvin_to_celsius<-function(K){
  K- 273.15  
}

################################################################
#  Temperature 
################################################################

tmpD<-read.csv(fileDT, header = TRUE)
tmpD$date<-as.Date(tmpD$date)
tmpD$SmoothTemp<-smooth_Whittaker(tmpD$temp, lambda=100)

ggplot()+
  geom_path(data=tmpD, aes(x=date, y=temperature), col="green")+
  geom_path(data=tmpD, aes(x=date, y=SmoothTemp), col="red")+
  theme_bw()+ xlab("Date")+ ylab("Temperature")

tmpN<-read.csv(flNT, header = TRUE)
tmpN$date<-as.Date(tmpN$date)
tmpN$SmoothTemp<-smooth_Whittaker(tmpN$temp, lambda=100)

ggplot()+
  geom_path(data=tmpN, aes(x=date, y=temperature), col="green")+
  geom_path(data=tmpN, aes(x=date, y=SmoothTemp), col="red")+
  theme_bw()+ xlab("Date")+ ylab("Temperature")

################################################################
#  Precipitation
################################################################

precip<-read.csv(flPR, header = TRUE)
precip$date<-as.Date(precip$date)

ggplot()+
  geom_line(data=precip, aes(x=date, y=precip), col="green")+
  theme_bw()+ xlab("Date")+ ylab("Precipitation")

################################################################
#  Vegetation Indices
################################################################

lf<-list.files(dirL8)
dates<-c()
for(i in 1:length(lf)){
  x<-strsplit(lf[i],'\\.')
  x<-x[[1]][1]
  x<-strsplit(x,prefL8)
  if(length(x[[1]])==2){
    x<-x[[1]][2]
    dates<-c(dates, x)
  }
}
datesAUS<-sort(dates)
range(datesAUS)

dataAUS<-data.frame(date=c(), SR_B1=c(), SR_B2=c(), SR_B3=c(),
                    SR_B4=c(), SR_B5=c(), SR_B6=c(), SR_B7=c())
for(ii in 1:length(datesAUS)){
  fln<-paste0(dirL8, prefL8, datesAUS[ii],".tif")
  l8bds <- rast(x=fln)
  rasValue <- extract(l8bds, data.frame(lng, lat))
  dataAUS<-rbind(dataAUS, data.frame(date=datesAUS[ii], 
                                     SR_B1=rasValue$SR_B1, SR_B2=rasValue$SR_B2, SR_B3=rasValue$SR_B3,
                                     SR_B4=rasValue$SR_B4, SR_B5=rasValue$SR_B5, SR_B6=rasValue$SR_B6,
                                     SR_B7=rasValue$SR_B7))
}

smoothing_coeff=1
dts<-seq(min(as.Date(datesAUS)), max(as.Date(datesAUS)), by=10)
clmnm<-colnames(dataAUS)

for(jj in 2:ncol(dataAUS)){
    xx<-dataAUS[, c(2,jj)]
    xx<-xx[-which(is.na(xx[,2])),]
    xx<-xx[order(xx[,1]),]
    yy<-smooth_Whittaker(xx[,2], lambda=smoothing_coeff)
    ff<-approx(x=as.Date(xx$date), y = yy, xout=dts, method = "linear")
    if(clmnm[jj]=="SR_B1") SR_B1<-ff$y
    if(clmnm[jj]=="SR_B2") SR_B2<-ff$y
    if(clmnm[jj]=="SR_B3") SR_B3<-ff$y
    if(clmnm[jj]=="SR_B4") SR_B4<-ff$y
    if(clmnm[jj]=="SR_B5") SR_B5<-ff$y
    if(clmnm[jj]=="SR_B6") SR_B6<-ff$y
    if(clmnm[jj]=="SR_B7") SR_B7<-ff$y
  }
  dataAUSsmooth<-data.frame(date=dts, SR_B1=SR_B1, 
                 SR_B2=SR_B2, SR_B3=SR_B3, SR_B4=SR_B4, SR_B5=SR_B5, SR_B6=SR_B6, SR_B7=SR_B7)

# Ratio vegetation index
dataAUSsmooth$RVI <- dataAUSsmooth$SR_B5/dataAUSsmooth$SR_B4
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=RVI))+
  geom_path(aes(group=ID))

# Difference vegetation index
dataAUSsmooth$DVI <- dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=DVI))+
  geom_path(aes(group=ID))

#  Normalized difference vegetation index
dataAUSsmooth$NDVI<-(dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4)/(dataAUSsmooth$SR_B5+dataAUSsmooth$SR_B4)
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=NDVI))+
  geom_path(aes(group=ID))+xlab("Date")

# Enhanced Vegetation Index
G<-2.5; C1<-6; C2<-7; L<-1
dataAUSsmooth$EVI<-G*(dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4)/
  (dataAUSsmooth$SR_B5+C1*dataAUSsmooth$SR_B4+
     C2*dataAUSsmooth$SR_B2 +L)
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=EVI))+
  geom_path(aes(group=ID))

# Soil-adjusted vegetation index
L=0.5
dataAUSsmooth$SAVI<-(1+L)*(dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4)/
  (dataAUSsmooth$SR_B5+dataAUSsmooth$SR_B4+L)
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=SAVI))+
  geom_path(aes(group=ID))

# Modified soil-adjusted vegetation index
dataAUSsmooth$MSAVI<-0.5*((2*dataAUSsmooth$SR_B5+1)-
                            sqrt((2*dataAUSsmooth$SR_B5+1)^2)+8*(dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4))
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=MSAVI))+
  geom_path(aes(group=ID)) +xlab("Date") +
  xlim(as.Date("2013-01-01"), as.Date("2018-01-01"))

# Optimized soil-adjusted vegetation index
X<-0.16
dataAUSsmooth$OSAVI<-(dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4)/
  (dataAUSsmooth$SR_B5+dataAUSsmooth$SR_B4+X)
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=OSAVI))+
  geom_path(aes(group=ID))

# Normalized difference phenology index
dataAUSsmooth$NDPI<-(dataAUSsmooth$SR_B5-(0.74*dataAUSsmooth$SR_B4+0.26*dataAUSsmooth$SR_B6))/
  (dataAUSsmooth$SR_B5+(0.74*dataAUSsmooth$SR_B4+0.26*dataAUSsmooth$SR_B6))
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=NDPI))+
  geom_path(aes(group=ID))

# kernel NDVI
dataAUSsmooth$kNDVI<-tanh(((dataAUSsmooth$SR_B5-dataAUSsmooth$SR_B4)^2)/
                            (2*0.5*dataAUSsmooth$SR_B5+dataAUSsmooth$SR_B4))
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=kNDVI))+
  geom_path(aes(group=ID))

# NIRv
dataAUSsmooth$NIRv<-dataAUSsmooth$SR_B5*dataAUSsmooth$NDVI
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=NIRv))+
  geom_path(aes(group=ID))

# Global Environment Monitoring Index
eta<-(2*(dataAUSsmooth$SR_B5^2 -dataAUSsmooth$SR_B4^2)+1.5*dataAUSsmooth$SR_B5+
        0.5*dataAUSsmooth$SR_B4)/(dataAUSsmooth$SR_B5+dataAUSsmooth$SR_B4+0.5)
dataAUSsmooth$GEMI<-eta*(1-0.25*eta)-
  (dataAUSsmooth$SR_B4-0.0125)/(1-dataAUSsmooth$SR_B4)
ggplot(dataAUSsmooth, aes(x=as.Date(date), y=GEMI))+
  geom_path(aes(group=ID))

################################################################
#  Train RF model
################################################################

KK<-100
ntree <- 100

frcst<-data.frame(date=c(), obs=c(), pred=c(), diff=c(), VI=c())
for(nm.vi in colnames(dataAUSsmooth)[10:ncol(dataAUSsmooth)]){
  cat("\n");print(nm.vi)
    xx<-dataAUSsmooth[, c('date', nm.vi)]
    d0<-min(xx$date)
    xx$day<-as.numeric(difftime(xx$date, d0, units = c("days")))
    for(ii in 1:nrow(xx)){
      yy<-tmpD
      yy$day<-as.numeric(difftime(yy$date, d0, units = c("days")))
      y.d<-approx(yy$day, yy$temp,
                  seq(xx$day[ii]-KK, xx$day[ii]-1, by= 1))
      y.d<-y.d$y
      yy<-tmpN
      yy$day<-as.numeric(difftime(yy$date, d0, units = c("days")))
      y.n<-approx(yy$day[yy$Type=="Night time"], yy$temp[yy$Type=="Night time"],
                  seq(xx$day[ii]-KK, xx$day[ii]-1, by= 1))
      y.n<-y.n$y
      yy<-precip
      yy$day<-as.numeric(difftime(yy$date, d0, units = c("days")))
      y.r<-approx(yy$day, yy$prec,
                  seq(xx$day[ii]-KK, xx$day[ii]-30, by= 30))
      y.r<-y.r$y
      y<-c(xx[ii, nm.vi], y.d, y.n, y.r)
      y<-as.data.frame(t(y))
      colnames(y)<-c("X", sapply(1:length(y.d), function(a) paste0("D", a)),
                     sapply(1:length(y.n), function(a) paste0("N", a)),
                     sapply(1:length(y.r), function(a) paste0("R", a)))
      if(ii==1){
        dat<-y
      } else {
        dat<-rbind(dat, y)
      }
    }
    wh.tr<-which(xx$date<=train.date)
    wh.tst<-which(xx$date>train.date & xx$date<=test.date)
    wh.frc<-which(xx$date>test.date)
    
    rf.fit <- randomForest(X ~ ., data=dat[c(wh.tr, wh.tst),], ntree=ntree)     
    prediction <- predict(rf.fit, dat[,-1])
    frcst<-rbind(frcst, data.frame(date=xx$date, 
                                   obs=dat$X, pred=prediction,
                                   diff=prediction-dat$X, VI=nm.vi))
    
  }

################################################################
#  Check for VI anomaly
################################################################

ggplot(frcst, aes(x=as.Date(date), y=diff)) + 
  geom_path(aes(col=diff)) +
  scale_colour_gradient2(name="", 
                         low = "#a50026", mid = "#fee08b", high = "black", 
                         midpoint = 0)+
  xlab("Date") + theme_minimal() + ylab(paste0("\u0394", "NDVI"))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") 




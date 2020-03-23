### storms are defined as:
## 1 - windspeed threshold 10 m/s
## 2 - winds below threshold for at least 24 hours before event
## 3 - at least 12 hours of winds greater than or equal to threshold

### calms are defined as:
## 1 - winds greater than or equal to threshold for minimum 12 hours before event
## 2 - at least 10 days of winds below threshold

library(lubridate)
library(R.matlab)

# ### extract the raw winds
# setwd("~/Desktop/Professional/projects/dissertation/data/storms")
# ### this is a matlab file format, I originially downloaded and concatenated the data with matlab
# data <- readMat('20180601_cfsr_full.mat')
# u <- data$u
# v <- data$v
# time <- data$t.dv
# lon <- data$lon.a-360
# lat <- rev(data$lat.a)
# rm(data)
# uv <- array(NA,c(dim(u)[1],dim(u)[2],dim(u)[3]))
# for(i in 1:dim(u)[2]){
#   for(j in 1:dim(u)[3]){
#     uv[,i,j] <- sqrt(u[,i,j]^2 + v[,i,j]^2)
#   }
# }
# image(uv[1,,])
# rm(u,v,i,j)

#### previous code makes the data being loaded below
load("~/Desktop/professional/projects/dissertation/data/storms/cfsr_2week_part_a.RData")
### make a land mask
mask <- uv[1,,]
mask[!is.na(mask)] <- 1
mask[is.na(mask)] <- 0

### hake yolk sack period in days
ys <- 10
### windspeed threshold (meters per second)
threshold <- 10
storm <- matrix(NA,1150450,15)
calm <- matrix(NA,1150450,15)
a <- 1
b <- 0
c <- 1
d <- 0
e <- 1
f <- 0
pb <- txtProgressBar(min = 0, max = dim(uv)[2], style = 3)
for(i in 1:dim(uv)[2]){
  for(j in 1:dim(uv)[3]){
    windsp <- uv[,i,j]
    if (sum(is.na(windsp))==0){
      ### new storms
      m <- 1
      n <- 1
      io <- 0
      start <- NULL
      stop <- NULL
      scount <- NULL
      csum <- NULL
      mws <- NULL
      sdws <- NULL
      count <- 0
      ## 24 hours of calm before
      r <- 4
      ## 12 hours of storms
      s <- 1
      for(k in r:(length(windsp)-r)){
        if(io==0 & mean(windsp[k:(k+s)]>=threshold)==1 & mean(windsp[(k-r):(k-1)]<threshold)==1){
          ## identifies the start of a storm event
          start[m] <- k
          m <- m + 1
          io <- 1
        } 
        if(io==1 & mean(windsp[(k-s):k]>=threshold)==1 & mean(windsp[(k+1)]<threshold)==1){
          ## identifies the stop of a storm event
          stop[n] <- k
          csum[m-1] <- sum(unlist(lapply(windsp[start[m-1]:(stop[n]-1)],function(x){x^3})))
          mws[m-1] <- mean(windsp[start[m-1]:(stop[n]-1)])
          sdws[m-1] <- sd(windsp[start[m-1]:(stop[n]-1)])
          n <- n + 1
          scount[m-1] <- count
          io <- 0
          count <- 0
        }
        if(io==1){
          count <- count + 1
        }
      }
      if(length(start)>0){
        if(length(scount)==(length(start)-1)){
          scount <- c(scount,NA)
          stop <- c(stop,NA)
          csum <- c(csum,NA)
          mws <- c(mws,NA)
          sdws <- c(sdws,NA)
        }
        d <- d + length(start)
        storm[c:d,1] <- start # start index
        storm[c:d,2] <- lat[i]
        storm[c:d,3] <- lon[j]
        storm[c:d,4:9] <- time[start,] # year, month, day, hour, minute, second
        storm[c:d,10] <- scount # length of event
        storm[c:d,11] <- stop # stop index
        storm[c:d,12] <- csum # cummulative windspeed cubed, index of strength
        storm[c:d,13] <- mws # mean winds speed during event
        storm[c:d,14] <- sdws # standard deviation of wind speed during event
        storm[(c+1):d,15] <- start[2:length(start)]-stop[1:(length(stop)-1)] # duration between events
        c <- d + 1
      }
      ### calm periods
      m <- 1
      n <- 1
      io <- 0
      start <- NULL
      stop <- NULL
      scount <- NULL
      csum <- NULL
      mws <- NULL
      sdws <- NULL
      count <- 0
      ## 12 hours of storms before
      r <- 2
      ## ten days of calm
      s <- 4*ys
      for(k in s:(length(windsp)-s)){
        if(io==0 & mean(windsp[k:(k+s)]<threshold)==1 & mean(windsp[(k-r):(k-1)]>=threshold)==1){
          start[m] <- k
          m <- m + 1
          io <- 1
        } 
        if(io==1 & mean(windsp[(k-s):k]<threshold)==1 & mean(windsp[(k+1)]>=threshold)==1){
          stop[n] <- k
          csum[m-1] <- sum(unlist(lapply(windsp[start[m-1]:(stop[n]-1)],function(x){x^3})))
          mws[m-1] <- mean(windsp[start[m-1]:(stop[n]-1)])
          sdws[m-1] <- sd(windsp[start[m-1]:(stop[n]-1)])
          n <- n + 1
          scount[m-1] <- count
          io <- 0
          count <- 0
        }
        if(io==1){
          count <- count + 1
        }
      }
      if(length(start)>0){
        if(length(scount)==(length(start)-1)){
          scount <- c(scount,NA)
          stop <- c(stop,NA)
          csum <- c(csum,NA)
          mws <- c(mws,NA)
          sdws <- c(sdws,NA)
        }
        f <- f + length(start)
        calm[e:f,1] <- start
        calm[e:f,2] <- lat[i]
        calm[e:f,3] <- lon[j]
        calm[e:f,4:9] <- time[start,]
        calm[e:f,10] <- scount
        calm[e:f,11] <- stop
        calm[e:f,12] <- csum
        calm[e:f,13] <- mws
        calm[e:f,14] <- sdws
        calm[(e+1):f,15] <- start[2:length(start)]-stop[1:(length(stop)-1)]
        e <- f + 1
      }
    }####
  }
  setTxtProgressBar(pb, i)
}

### remove unused rows to save memory
storm <- storm[-which(is.na(storm[,1])),]
calm <- calm[-which(is.na(calm[,1])),]

### monthly aggregates
lon<-as.vector(lon)
### storms per year and mean storm strength and length
stormsmth <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
stormmth.c <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
stormmth.s <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
stormmth.l <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
stormmth.b <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
m <- 1
n <- 0
pb <- txtProgressBar(min = 0, max = length(lat), style = 3)
for(i in 1:length(lat)){
  for(j in 1:length(lon)){
    x <- as.data.frame(storm[storm[,2]==lat[i] & storm[,3]==lon[j],])
    for(k in 1979:2015){
      storms <- rep(NA,length(1:12))
      storm.c <- rep(NA,length(1:12))
      storm.s <- rep(NA,length(1:12))
      storm.l <- rep(NA,length(1:12))
      storm.b <- rep(NA,length(1:12))
      for(l in 1:12){
        y <- x[which(x[,4]==k & x[,5]==l),]
        storms[l] <- nrow(y)
        storm.c[l] <- mean(y[,12],na.rm=T)
        storm.s[l] <- mean(y[,13],na.rm=T)
        storm.l[l] <- mean(y[,10],na.rm=T)
        storm.b[l] <- mean(y[,15],na.rm=T)
      }
      n <- n + 12
      stormsmth[i,j,m:n] <- storms
      stormmth.c[i,j,m:n] <- storm.c
      stormmth.s[i,j,m:n] <- storm.s
      stormmth.l[i,j,m:n] <- storm.l
      stormmth.b[i,j,m:n] <- storm.b
      m <- n + 1
    }
    m <- 1
    n <- 0
  }
  setTxtProgressBar(pb, i)
}
stormsmth[which(is.na(stormsmth))] <- 0
stormmth.c[which(is.na(stormmth.c))] <- 0
stormmth.s[which(is.na(stormmth.s))] <- 0
stormmth.l[which(is.na(stormmth.l))] <- 0
stormmth.b[which(is.na(stormmth.b))] <- 0

stormsmth[mask==0] <- NA
stormmth.c[mask==0] <- NA
stormmth.s[mask==0] <- NA
stormmth.l[mask==0] <- NA
stormmth.b[mask==0] <- NA

### calms per year and mean storm strength and length
calmsmth <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
calmmth.c <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
calmmth.s <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
calmmth.l <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
calmmth.b <- array(NA,c(length(lat),length(lon),length(1:12)*length(1979:2015)))
m <- 1
n <- 0
pb <- txtProgressBar(min = 0, max = length(lat), style = 3)
for(i in 1:length(lat)){
  for(j in 1:length(lon)){
    x <- as.data.frame(calm[calm[,2]==lat[i] & calm[,3]==lon[j],])
    for(k in 1979:2015){
      calms <- rep(NA,length(1:12))
      calm.c <- rep(NA,length(1:12))
      calm.s <- rep(NA,length(1:12))
      calm.l <- rep(NA,length(1:12))
      calm.b <- rep(NA,length(1:12))
      for(l in 1:12){
        y <- x[which(x[,4]==k & x[,5]==l),]
        calms[l] <- nrow(y)
        calm.c[l] <- mean(y[,12],na.rm=T)
        calm.s[l] <- mean(y[,13],na.rm=T)
        calm.l[l] <- mean(y[,10],na.rm=T)
        calm.b[l] <- mean(y[,15],na.rm=T)
      }
      n <- n + 12
      calmsmth[i,j,m:n] <- calms
      calmmth.c[i,j,m:n] <- calm.c
      calmmth.s[i,j,m:n] <- calm.s
      calmmth.l[i,j,m:n] <- calm.l
      calmmth.b[i,j,m:n] <- calm.b
      m <- n + 1
    }
    m <- 1
    n <- 0
  }
  setTxtProgressBar(pb, i)
}
calmsmth[which(is.na(calmsmth))] <- 0
calmmth.c[which(is.na(calmmth.c))] <- 0
calmmth.s[which(is.na(calmmth.s))] <- 0
calmmth.l[which(is.na(calmmth.l))] <- 0
calmmth.b[which(is.na(calmmth.b))] <- 0

calmsmth[mask==0] <- NA
calmmth.c[mask==0] <- NA
calmmth.s[mask==0] <- NA
calmmth.l[mask==0] <- NA
calmmth.b[mask==0] <- NA

time.mmyyyy <- expand.grid(month=1:12,year=1979:2015)


### creating the Netcdf file
library(ncdf4)
dimlon <- ncdim_def('Lon','degreesE',seq(lon[1],lon[21],.5))
dimlat <- ncdim_def('Lat','degreesN',seq(lat[17],lat[1],-.5))
dimtime <- ncdim_def('Time','decimal years',seq(1979,2015.999,1/12))
calms_var <- ncvar_def('calms','count',list(dimlat,dimlon,dimtime),-1,"Number of calms per month",prec='integer')
storms_var <- ncvar_def('storms','count',list(dimlat,dimlon,dimtime),-1,"Number of storms per month",prec='integer')

wcalms_var <- ncvar_def('calms_ws','meters per second',list(dimlat,dimlon,dimtime),-1,"Mean wind speed during storms",prec='double')
wstorms_var <- ncvar_def('storms_ws','meters per second',list(dimlat,dimlon,dimtime),-1,"Mean wind speed during storms",prec='double')
lcalms_var <- ncvar_def('calms_len','days',list(dimlat,dimlon,dimtime),-1,"Duration of calm events per month",prec='double')
lstorms_var <- ncvar_def('storms_len','days',list(dimlat,dimlon,dimtime),-1,"Duration of storm events per month",prec='double')
bcalms_var <- ncvar_def('calms_bt','days',list(dimlat,dimlon,dimtime),-1,"Duration between calms per month",prec='double')
bstorms_var <- ncvar_def('storms_bt','days',list(dimlat,dimlon,dimtime),-1,"Duration between storms per month",prec='double')
events <- nc_create('events.nc',list(calms_var,wcalms_var,lcalms_var,bcalms_var,
                                     storms_var,wstorms_var,lstorms_var,bstorms_var))
ncvar_put(events,calms_var,calmsmth)
ncvar_put(events,wcalms_var,calmmth.s)
ncvar_put(events,lcalms_var,calmmth.l/4)
ncvar_put(events,bcalms_var,calmmth.b/4)

ncvar_put(events,storms_var,stormsmth)
ncvar_put(events,wstorms_var,stormmth.s)
ncvar_put(events,lstorms_var,stormmth.l/4)
ncvar_put(events,bstorms_var,stormmth.b/4)
# add global attributes
title <- "Summation of storm and calms events in Southern California Current per month per year relevant to larval Pacific hake (Merluccius productus)"
ncatt_put(events,0,"title",title)
ncatt_put(events,0,"wind input used",'NOAA CFSR, 6-hourly')
# ncatt_put(events,0,"institution",institution$value)
source <- "https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2#CFS%20Reanalysis%20(CFSR)"
ncatt_put(events,0,"source",source)
refs <- 'Turley & Rykaczewski, 2019, Can. J. Fish, 76(12): 2418-2432, https://doi.org/10.1139/cjfas-2018-0458'
ncatt_put(events,0,"references",refs)
history <- paste("B.D. Turley", date(), sep=", ")
ncatt_put(events,0,"history",history)
nc_close(events)

          
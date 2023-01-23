library(zoo)
library(dplyr)
library(hydroGOF)
library(trend)
library(sp)
library(raster)
library(gplots)
library(stats)
library(car)
library(reshape2)
library(ggplot2)





##### GAP FILLING DATA SELECTION #####
### Reading data ###
meta <- read.csv("Meta_gap0.csv", header=T, sep=";")
data.nogap <- read.csv("Data_nogap.csv", header=T, sep=";")
data.gap <- read.csv("Data_gap0.csv", header=T, sep=";")



### Fraction of data gain thks to gap filling ###
data.nogap.v2 <- data.nogap[which(data.nogap$Month == 10 | data.nogap$Month == 11 |
                                  data.nogap$Month == 12 | data.nogap$Month == 1 |
                                  data.nogap$Month == 2 | data.nogap$Month == 3 |
                                  data.nogap$Month == 4),]
data.nogap.v2 <- data.nogap.v2[-c(1:4),]
rownames(data.nogap.v2) <- NULL
N.data.nogap <- colSums(!is.na(data.nogap.v2[,-c(1:2)]))/length(data.nogap.v2$Year)*100

data.gap.v2 <- data.gap[which(data.gap$Month == 10 | data.gap$Month == 11 |
                              data.gap$Month == 12 | data.gap$Month == 1 |
                              data.gap$Month == 2 | data.gap$Month == 3 |
                              data.gap$Month == 4),]
data.gap.v2 <- data.gap.v2[-c(1:4),]
rownames(data.gap.v2) <- NULL
N.data.gap <- colSums(!is.na(data.gap.v2[,-c(1:2)]))/length(data.gap.v2$Year)*100

delta.N <- N.data.gap - N.data.nogap



### Data selection ###
tt <- c(10,20,30)

for (kk in 1:length(tt)) {
  
  meta.x <- meta[which(delta.N <= tt[kk]),]
  data.x <- data.nogap[,c(1:2,as.numeric(rownames(meta.x))+2)]
  
  N <- length(meta.x[,1]) # stations in data.x and meta.x need to be in the same order
  mat_filled <- data.x # matrix were the filled series will be saved
  common_years <- 10 # minimum number of years with common data between test and reference series
  min.weight <- 0.00005 # minimum weight for reference station selection
  
  ### Weigths matrix ###
  # Computing station weights based on distance and elevation difference between station pairs #
  Rt <- 6371
  conv <- pi/180
  dist_r <- matrix(nrow=length(meta.x[,1]), ncol=N)
  delta_z <- matrix(nrow=length(meta.x[,1]), ncol=N)
  
  for (ii in 1:N) {
    angle <- sin(conv*meta.x$Latitude) * sin(conv*meta.x$Latitude[ii]) + cos(conv*meta.x$Latitude) * cos(conv*meta.x$Latitude[ii]) * cos(abs(conv*meta.x$Longitude-conv*meta.x$Longitude[ii])) 
    angle[which(angle > 1)] <- 1
    dist_r[,ii] <- Rt*acos(angle)
    delta_z[,ii] <- abs(meta.x$Elevation-meta.x$Elevation[ii])
  }
  
  # Halving coefficients for weights: 20km for distance and 200m for elevation difference #
  halv_dist <- 20
  halv_elev <- 200
  c_dist <- (halv_dist^2)/log(2)
  c_h <- (halv_elev^2)/log(2)
  
  wd <- exp(-(dist_r^2)/c_dist)
  
  wh <- matrix(nrow=length(meta.x[,1]), ncol=N)
  for (t in 1:N) {
    wh[,t] <- exp(-(delta_z[,t]^2)/c_h)
  }
  
  pesi <- wd*wh
  diag(pesi) <- NA 
  
  ### Gap filling ###
  for (i in 1:N) { # iter on stations
    print(paste("station",meta.x$Name[i],i))
    NMAX <- 10 # max number of reference stations used for the reconstruction 
    missing <- which(is.na(data.x[,i+2]))
    
    
    if (length(missing) > 0) {
      d <- which(pesi[,i] > min.weight) # only stations with weights above the threshold are considered 
      col_ref <- length(d)
      
      if (col_ref > 2) { # stations are sorted by decreasing weights (the reconstruction is performed only if the selected stations are more than two)
        l <- sort.int(as.numeric(pesi[d,i]), decreasing=TRUE, index.return=TRUE)
        d <- d[l$ix]
        
        if (col_ref < NMAX) {NMAX <- col_ref} # if selected stations are less than NMAX, NMAX is set to the number of stations with relevant weight
        else {NMAX <- 10} 
        
        
        for (k in 1:12) { # iter on months
          mes <- which(data.x[,2] == k) # select data for month k
          vet_test <- data.x[mes,i+2] # extract the values of the test series
          vet_filled <- vet_test # create a duplicate of reference series for the month k where the filled gaps will be added
          gaps <- which(is.na(vet_test)) # select the years with missing data
          
          if (length(gaps) > 0) {
            mat_ref <- data.x[mes,d+2] # matrix of the sorted reference series for the month k
            stat <- matrix(nrow=length(gaps), ncol=NMAX) # matrix were reconstruction statistics will be stored
            filled <- numeric(length(gaps)) # vector for storing only reconstructed gaps of test series
            filled[] <- NA
            
            for (k1 in 1:length(gaps)) { # iter on gaps
              NN <- 1 
              j <- 1
              
              while (NN<=NMAX && j<=col_ref) { # cycle on reference stations until NMAX is reached or until the total number of available reference stations is exploited 
                if (!is.na(mat_ref[gaps[k1],j])) { # not missing value in the reference 
                  sel <- which(!is.na(vet_test) & !is.na(mat_ref[,j])) # years with common data
                  
                  if (length(sel) >= common_years) {
                    sum_test <- sum(vet_test[sel])
                    sum_ref <- sum(mat_ref[sel,j])
                    ratio_sum <- sum_test/sum_ref
                    
                    if (sum_ref == 0) {ratio_sum <- 1}
                    
                    stat[k1,NN] <- round((mat_ref[gaps[k1],j]*ratio_sum),4) # reconstructed value is rescaled by the ratio
                    NN <- NN+1    
                  }
                }
                j <- j+1 
              }
            }
            
            for (k2 in 1:length(gaps)) {
              filled[k2] <- mean(stat[k2,],na.rm=TRUE) # filled gap as the average of the NMAX simulations - alternatively, median/weighted mean can be considered  
            } 
            
            vet_filled[gaps] <- filled  
            mat_filled[mes,i+2] <- vet_filled # store the filled series for the month k in the data matrix
            
          }
        }
      }
    }        
  }        
  
  mat_filled2 <- mat_filled
  for (a in 3:ncol(mat_filled2)) {
    mat_filled2[,a][which(is.na(mat_filled2[,a]))] <- NA
    if (is.numeric(mat_filled2[,a]) == T) {mat_filled2[,a] <- round(mat_filled2[,a], digits=0)}
  }
  
  ### Saving the data ###
  write.table(meta.x, paste0("Meta_gap",tt[kk],".csv"), sep=";", row.names=F)
  write.table(mat_filled2, paste0("Data_gap",tt[kk],".csv"), sep=";", row.names=F)
  
}





##### SEASONAL DATA #####
t <- c(10,20,30)

for (k in 1:length(t)) {
  
  Meta <- read.csv(paste0("Meta_gap",t[k],".csv"), header=T, sep=";")
  HN <- read.csv(paste0("Data_gap",t[k],".csv"), header=T, sep=";")
  dt <- seq(7,280,7)
  
  HN_TOT <- as.data.frame(matrix(nrow=40, ncol=length(Meta$Name)+1))
  colnames(HN_TOT) <- c("Date", Meta$Name)
  HN_TOT$Date <- paste(seq(min(HN$Year),max(HN$Year)-1),seq(min(HN$Year)+1,max(HN$Year)),sep="-")
  
  for (i in 1:length(Meta$Name)) {
    HN_X <- HN[,c(1:2,i+2)]
    HN_X <- HN_X[which(HN_X$Month == 10 | HN_X$Month == 11 | HN_X$Month == 12 |
                         HN_X$Month == 1 | HN_X$Month == 2 |
                         HN_X$Month == 3 | HN_X$Month == 4),]
    HN_X <- HN_X[-c(1:4),]
    rownames(HN_X) <- NULL
    
    for (j in 1:length(dt)) {
      HN_X2 <- HN_X[,3][(dt[j]-6):dt[j]]
      NA_X <- (sum(is.na(HN_X2))/length(HN_X2))*100
      if (NA_X <= 30) {
        HN_TOT[j,i+1] <- sum(HN_X2, na.rm=T)
      }
    }
  }
  
  write.table(HN_TOT, paste0("Data_gap",t[k],"_S.csv"), sep=";", row.names=F)
  
}





#### TREND #####
### Data ###
t <- c(0,10,20,30)
m <- c(10:12,1:4)
ms <- c(month.abb[m],"Seas")
ar <- c(0,1000,2000,3000)
a <- 0.05

TTab_tot <- vector("list", length(t))
for (ii in 1:length(t)) {
  TTab_tot[[ii]] <- vector("list", 8)
  names(TTab_tot[[ii]]) <- ms
}
names(TTab_tot) <- t

TAA <- shapefile("TopographicData/TAA.shp")
TAA <- spTransform(TAA, crs("+proj=longlat +datum=WGS84"))
DEM <- raster("TopographicData/DEM_TAA_100m.tif")
Hillshade <- hillShade(terrain(DEM, "slope", "radians", 8),
                       terrain(DEM, "aspect", "radians", 8),
                       angle=45, direction=315, normalize=T)



### Trend computaion & plotting ###
for (k in 1:length(t)) {
  
  # Data reading #
  Meta <- read.csv(paste0("Meta_gap",t[k],".csv"), header=T, sep=";")
  Data <- read.csv(paste0("Data_gap",t[k],".csv"), header=T, sep=";")
  Data_S <- read.csv(paste0("Data_gap",t[k],"_S.csv"), header=T, sep=";")
  
  # Trend calculation #
  TS_s <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=(length(ms)+1)))
  TS_s[,1] <- Meta$Name
  colnames(TS_s) <- c("Name",ms)
  TS_lci <- TS_s
  TS_uci <- TS_s
  MK_p <- TS_s
  Pt <- TS_s
  TS_s_r <- TS_s
  HN_avg <- TS_s
  
  for (j in 1:length(ms)) {
    
    for (i in 1:length(Meta$Name)) {
      
      if (j == 8) {
        
        Data_X <- Data_S[,c(1,i+1)]
        Data_X2 <- Data_X[,2]
        
      } else {
        
        Data_X <- Data[,c(1:2,i+2)]
        Data_X <- Data_X[-c(1:9),]
        rownames(Data_X) <- NULL
        Data_X2 <- Data_X[which(Data_X$Month == m[j]),]
        Data_X2 <- Data_X2[,3]
        rownames(Data_X2) <- NULL
        
      }
      
      
      if (sum(is.na(Data_X2)) == 0 && sum(Data_X2 != 0) >= 5) {
        
        HN_avg[i,j+1] <- round(mean(Data_X2, na.rm=T), digits=3)
        TS_s[i,j+1] <- 10*round(sens.slope(Data_X2, conf.level=0.95)[[1]][[1]], digits=3)
        TS_lci[i,j+1] <- 10*round(sens.slope(Data_X2, conf.level=0.95)[["conf.int"]][[1]], digits=3)
        TS_uci[i,j+1] <- 10*round(sens.slope(Data_X2, conf.level=0.95)[["conf.int"]][[2]], digits=3)
        MK_p[i,j+1] <- round(mk.test(Data_X2)[[2]], digits=5)
        TS_s_r[i,j+1] <- round(TS_s[i,j+1]/mean(Data_X2)*100, digits=0)
        
        if (MK_p[i,j+1] <= a) {Pt[i,j+1] <- 19} else {Pt[i,j+1] <- 1}
        
      } else {
        
        HN_avg[i,j+1] <- NA
        TS_s[i,j+1] <- NA
        TS_lci[i,j+1] <- NA
        TS_uci[i,j+1] <- NA
        MK_p[i,j+1] <- NA
        Pt[i,j+1] <- NA
        TS_s_r[i,j+1] <- NA
        
      }
      
    }
    
  }
  
  # Saving trend data #
  write.table(TS_s, paste0("Trend_gap",t[k],".csv"), sep=";", row.names=F)

  # Trend plot v1 #
  I <- c(-30,-10,-1,1,10,30)
  PtType <- c(6,6,0,2,2)
  CPal <- c("#FF0000","#ff8000","#666666","#2395ff","#0000FF")

  png(paste0("Plots/TrendHN_gap",t[k],"_v1.png"), width=2400, height=1200, res=100)
  par(mar=c(0,0,5,0), mgp=c(3,0,0), mfrow=c(2,4))

  for (j in 1:(length(ms))) {

    IClass_X <- cut(TS_s[,j+1], breaks=I)
    Col_X <- CPal[as.numeric(IClass_X)]
    Pt_X <- PtType[as.numeric(IClass_X)]
    Pt_X[which(TS_s[,j+1] <= I[3] & MK_p[,j+1] <= a)] <- 25
    Pt_X[which(TS_s[,j+1] > I[4] & MK_p[,j+1] <= a)] <- 24
    Pt_X[which(TS_s[,j+1] > I[3] & TS_s[,j+1] <= I[4] & MK_p[,j+1] <= a)] <- 22

    plot(Hillshade, main=ms[j], cex.main=3, col=gray(0:255/255), legend=F, axes=F, box=T, alpha=0.6)
    plot(TAA, lwd=1, border="black", add=T)
    points(Meta$Longitude, Meta$Latitude, type="p", pch=Pt_X, col=Col_X, bg=Col_X, cex=3, lwd=3)
    legend("bottomright", legend=rev(levels(IClass_X)), pch=rev(PtType), col=rev(CPal), pt.lwd=2, cex=2, bty="n")

  }

  dev.off()

  # Trend plot v2 #
  TS_s_TOT <- as.data.frame(matrix(nrow=8*length(Meta$Name), ncol=4))
  TS_s_TOT[,1] <- factor(rep(ms, each=length(Meta$Name)), levels=ms)
  TS_s_TOT[,2] <- rep(Meta$Name, times=8)
  TS_s_TOT[,3] <- rep(Meta$Elevation, times=8)
  colnames(TS_s_TOT) <- c("Period","Name","Elevation","Value")
  TS_lci_TOT <- TS_s_TOT
  TS_uci_TOT <- TS_s_TOT
  Pt_TOT <- TS_s_TOT

  dj <- c(length(Meta$Name)*c(1:length(ms)))

  for (j in 1:length(ms)) {

    TS_s_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],4] <- TS_s[,j+1]
    TS_lci_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],4] <- TS_lci[,j+1]
    TS_uci_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],4] <- TS_uci[,j+1]
    Pt_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],4] <- Pt[,j+1]

  }

  TS_s_TOT <- na.omit(TS_s_TOT)
  rownames(TS_s_TOT) <- NULL
  TS_lci_TOT <- na.omit(TS_lci_TOT)
  rownames(TS_lci_TOT) <- NULL
  TS_uci_TOT <- na.omit(TS_uci_TOT)
  rownames(TS_uci_TOT) <- NULL
  Pt_TOT <- na.omit(Pt_TOT)
  rownames(Pt_TOT) <- NULL

  TS_s_TOT2 <- cbind(TS_s_TOT, TS_lci_TOT$Value, TS_uci_TOT$Value, Pt_TOT$Value)
  colnames(TS_s_TOT2)[5:7] <- c("Lci","Uci","Pt")

  png(paste0("Plots/TrendHN_gap",t[k],"_v2.png"), width=2000, height=1000, res=100)

  print(ggplot(TS_s_TOT2, aes(Elevation, Value)) +
          geom_hline(yintercept=0, size=1, linetype="dashed") +
          geom_linerange(aes(ymin=Lci, ymax=Uci), alpha=0.25, size=0.5, colour="#000000") +
          geom_point(size=3.5, colour="#000000", shape=TS_s_TOT2$Pt) +
          facet_wrap(~ Period, nrow=2, scales="free_y") +
          geom_point(aes(x=2000, y=40), alpha=0) +
          geom_point(aes(x=1000, y=-40), alpha=0) +
          scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
          labs(x="Elevation [m]", y="Trend [cm/decade]") +
          theme_bw() +
          theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
                axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
                strip.text.x = element_text(face="bold", size=15, vjust=1),
                strip.text.y = element_text(face="bold", size=15, vjust=1),
                legend.position = "None", plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines"),
                panel.grid.minor = element_blank()))

  dev.off()

  # Trend table #
  TTab_X <- as.data.frame(matrix(nrow=3, ncol=6))
  TTab_X[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
  colnames(TTab_X) <- c("Altitude [m]","#","AT [cm/decade]","Lci [cm/decade]","Uci [cm/decade]","RT [%]")

  for (j in 1:length(ms)) {

    for (kk in 1:3) {

      TTab_X[kk,2] <- sum(!is.na(TS_s[,j+1][which(Meta$Elevation > ar[kk] & Meta$Elevation <= ar[kk+1])]))
      TTab_X[kk,3] <- round(mean(TS_s[,j+1][which(Meta$Elevation > ar[kk] & Meta$Elevation <= ar[kk+1])], na.rm=T), 2)
      TTab_X[kk,4] <- round(mean(TS_lci[,j+1][which(Meta$Elevation > ar[kk] & Meta$Elevation <= ar[kk+1])], na.rm=T), 2)
      TTab_X[kk,5] <- round(mean(TS_uci[,j+1][which(Meta$Elevation > ar[kk] & Meta$Elevation <= ar[kk+1])], na.rm=T), 2)
      TTab_X[kk,6] <- round(mean(TS_s_r[,j+1][which(Meta$Elevation > ar[kk] & Meta$Elevation <= ar[kk+1])], na.rm=T), 2)

    }

    TTab_tot[[k]][[j]] <- TTab_X

  }
  
}





##### TREND COMPARISON #####
### Data ###
meta <- read.csv("Meta_gap0.csv", header=T, sep=";")
t <- c(0,10,20,30)
m <- c(10:12,1:4)
ms <- c(month.abb[m],"Seas")

trend <- vector("list", length(t))
for (k in 1:length(t)) {
  trend[[k]] <- read.csv(paste0("Trend_gap",t[k],".csv"), header=T, sep=";")
}
names(trend) <- t

stats <- vector("list", 2)
names(stats) <- c("PBIAS","MAE")
stats.x <- as.data.frame(matrix(nrow=3, ncol=9))
stats.x[,1] <- t[-1]
colnames(stats.x) <- c("Threshold",ms)



### Trend comparison ###
trend.ref <- trend[[1]] 
trend.NA <- trend.ref
trend.NA[,-1] <- NA

for (i in 1:2) {
  
  for (k in 1:(length(t)-1)) {
    
    trend.x <- left_join(trend.NA, trend[[k+1]], by="Name")[,c(1,10:17)]
    colnames(trend.x) <- colnames(trend.ref)
    
    if (i == 1) {
      for (j in 1:length(ms)) {
        stats.x[k,j+1] <- pbias(trend.x[,j+1], trend.ref[,j+1], na.rm=T)
      }
    } else {
      for (j in 1:length(ms)) {
        stats.x[k,j+1] <- round(mae(trend.x[,j+1], trend.ref[,j+1], na.rm=T), 5)
      }
    }
    
  }
  
  stats[[i]] <- stats.x
  
}

# write.table(stats[[1]], "PBIAS.csv", sep=";", row.names=F)
# write.table(stats[[2]], "MAE.csv", sep=";", row.names=F)





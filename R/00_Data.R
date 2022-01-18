##### ----- 01 ELABORATED DATA ----- #####
rm(list = ls(all.names = TRUE))

### MTAA DATA ###
Meta_MTAA <- read.csv("00_Data/00_OriginalData/MTAA/MTAA-meta.csv", header=T, sep=";")
Data_MTAA <- read.csv("00_Data/00_OriginalData/MTAA/MTAA-data.csv", header=T, sep=";")

Y1 <- as.numeric(min(Data_MTAA$Year))
Y2 <- as.numeric(max(Data_MTAA$Year))
Year <- rep(Y1:Y2, each=12)
Month <- c(1:12)
Month <- rep(Month, times=(Y2-Y1)+1)

Data_MTAA2 <- as.data.frame(cbind(Year,Month))
for (i in 1:length(Meta_MTAA$Name)) {
  S_X <- Data_MTAA[which(Meta_MTAA$Name[i] == Data_MTAA$Name),]
  S_X2 <- as.data.frame(matrix(NA, nrow=nrow(Data_MTAA2), ncol=1))
  colnames(S_X2) <- Meta_MTAA$Name[i]
  
  for (j in 1:nrow(S_X2)) {
    for (k in 1:nrow(S_X)) {
      if (Data_MTAA2$Year[j] == S_X$Year[k] && Data_MTAA2$Month[j] == S_X$Month[k]) {
        S_X2[j,]=S_X$HN[k]
        k=nrow(S_X)
      }
    }
  }
  
  Data_MTAA2 <- cbind(Data_MTAA2, S_X2)
}

# write.table(Meta_MTAA, "00_Data/01_ElaboratedData/MTAA_MetaData.csv", sep=";", row.names=F)
# write.table(Data_MTAA2, "00_Data/01_ElaboratedData/MTAA_Data_HN.csv", sep=";", row.names=F)



### EURAC DATA ###
Meta_EURAC <- read.csv("00_Data/00_OriginalData/EURAC/TNBZ-meta.csv", header=T, sep=",")
Meta_EURAC <- Meta_EURAC[order(Meta_EURAC$Name),]
rownames(Meta_EURAC) <- NULL
Data_EURAC <- read.csv("00_Data/00_OriginalData/EURAC/TNBZ-data-HNHS.csv", header=T, sep=",")[,c(1:3)]
Data_EURAC <- Data_EURAC[order(Data_EURAC$Name),]
rownames(Data_EURAC) <- NULL
Data_EURAC[,3] <- round(Data_EURAC[,3], digits=0)

# Monthly aggregation #
Data_EURAC <- cbind(sapply(strsplit(Data_EURAC[,2],"-"),"[",1),
                    sapply(strsplit(Data_EURAC[,2],"-"),"[",2),
                    sapply(strsplit(Data_EURAC[,2],"-"),"[",3),
                    Data_EURAC)
colnames(Data_EURAC)[(1:3)] <- c("Year","Month","Day")
Data_EURAC <- Data_EURAC[,c(4,1:3,5:6)]
Data_EURAC_NA <- na.omit(Data_EURAC)

Data_EURACM_2 <- data.frame()
for (i in 1:length(Meta_EURAC$Name)) {
  Data_EURAC_X <- Data_EURAC[which(Meta_EURAC$Name[i] == Data_EURAC$Name),]
  Data_EURAC_X$Date2 <- paste(Data_EURAC_X$Year, Data_EURAC_X$Month, sep="-")
  Data_EURAC_X <- Data_EURAC_X[,c(1:5,7,6)]
  
  y1 <- as.numeric(min(Data_EURAC_X$Year))
  y2 <- as.numeric(max(Data_EURAC_X$Year))
  y <- as.character(rep(y1:y2, each=12))
  m <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  m <- rep(m, times=(y2-y1)+1)
  date <- cbind(y,m)
  date <- paste(date[,1], date[,2], sep="-")
  
  Data_EURACM <- as.data.frame(matrix(nrow=length(date), ncol=3))
  colnames(Data_EURACM) <- c("Name","Date","HN")
  Data_EURACM$Name <- Meta_EURAC$Name[i]
  Data_EURACM$Date <- date
  
  k=1
  for (j in 1:length(date)) {
    Data_EURAC_XM <- Data_EURAC_X$HN[which(date[j] == Data_EURAC_X$Date2)]
    NA_XM <- (sum(is.na(Data_EURAC_XM))/length(Data_EURAC_XM))*100
    if (NA_XM <= 30) {
      Data_EURACM$HN[k] <- sum(Data_EURAC_XM, na.rm=T)
      k=k+1
    } else {k=k+1}
  }
  
  Data_EURACM_2 <- rbind(Data_EURACM_2, Data_EURACM)
}

# Data setting #
y1 <- as.numeric(min(Data_EURAC_NA$Year))
y2 <- as.numeric(max(Data_EURAC_NA$Year))-1
y <- as.character(rep(y1:y2, each=12))
m <- c("01","02","03","04","05","06","07","08","09","10","11","12")
m <- rep(m, times=(y2-y1)+1)
date <- cbind(y,m)
date <- paste(date[,1], date[,2], sep="-")

Data_EURACM2 <- as.data.frame(date)
colnames(Data_EURACM2) <- "Date"
for (i in 1:length(Meta_EURAC$Name)) {
  Data_EURAC_XM2 <- Data_EURACM_2[which(Meta_EURAC$Name[i] == Data_EURACM_2$Name),]
  Data_EURAC_XM3 <- as.data.frame(matrix(NA, nrow=length(date), ncol=1))
  colnames(Data_EURAC_XM3) <- Meta_EURAC$Name[i]
  
  for (j in 1:length(date)) {
    for (k in 1:length(Data_EURAC_XM2$Date)) {
      if (date[j] == Data_EURAC_XM2$Date[k]) {
        Data_EURAC_XM3[j,] <- Data_EURAC_XM2$HN[k]
        k <- length(Data_EURAC_XM2$Date)}
    }
  }
  
  Data_EURACM2 <- cbind(Data_EURACM2, Data_EURAC_XM3)
}

Data_EURACM2 <- cbind(as.numeric(sapply(strsplit(Data_EURACM2[,1],"-"),"[",1)),
                      as.numeric(sapply(strsplit(Data_EURACM2[,1],"-"),"[",2)),
                      Data_EURACM2)
colnames(Data_EURACM2)[(1:2)] <- c("Year","Month")
Data_EURACM2 <- Data_EURACM2[,-3]

# write.table(Meta_EURAC, "00_Data/01_ElaboratedData/EURAC_MetaData.csv", sep=";", row.names=F)
# write.table(Data_EURACM2, "00_Data/01_ElaboratedData/EURAC_Data_HN.csv", sep=";", row.names=F)





##### ----- 02 FULL DATA ----- #####
rm(list = ls(all.names = TRUE))

Meta1 <- read.csv("00_Data/01_ElaboratedData/EURAC_MetaData.csv", header=T, sep=";")
Meta2 <- read.csv("00_Data/01_ElaboratedData/MTAA_MetaData.csv", header=T, sep=";")
Data1 <- read.csv("00_Data/01_ElaboratedData/EURAC_Data_HN.csv", header=T, sep=";")
Data2 <- read.csv("00_Data/01_ElaboratedData/MTAA_Data_HN.csv", header=T, sep=";")

Meta3 <- rbind(Meta1, Meta2)
Meta3 <- Meta3[order(Meta3$Name),]
rownames(Meta3) <- NULL

Data3 <- merge(Data1, Data2, by=c("Year","Month"))
Data3 <- Data3[order(Data3$Year,Data3$Month),]
rownames(Data3) <- NULL
Data4 <- Data3[,-(1:2)]
Data4 <- Data4[,order(colnames(Data4))]
Data5 <- cbind(Data3[,(1:2)], Data4)

# write.table(Meta3, "00_Data/02_FullData/FULL_MetaData.csv", sep=";", row.names=F)
# write.table(Data5, "00_Data/02_FullData/FULL_Data_HN.csv", sep=";", row.names=F)





##### ----- 03 DATA UP 2019-20 ----- #####
rm(list = ls(all.names = TRUE))

### TN ###
Meta_TN <- read.csv("00_Data/03_DataUp19-20/TN/MetaData.csv", header=T, sep=";")

Data1_TN <- vector("list", nrow(Meta_TN))
names(Data1_TN) <- Meta_TN$Name
for (i in 1:nrow(Meta_TN)) {
  Data1_TN[[i]] <- read.csv(paste0("03_DataUp19-20/TN/",Meta_TN$Name[i],".csv"), header=T, sep=";")
}

Date <- seq(as.Date("2019-06-01"), as.Date("2020-05-31"), by="days")
Date <- as.character(Date)

Data2_TN <- data.frame()
for (s in 1:length(Data1_TN)) {
  
  Data_X <- Data1_TN[[s]][-c(1:3),c(1,13)]
  Data_X <- head(Data_X, -7)
  colnames(Data_X) <- c("Date","HN")
  rownames(Data_X) <- NULL
  
  Data_X$HN[which(Data_X$HN == "///")] <- NA
  Data_X$HN[which(Data_X$HN == "8//")] <- 100
  Data_X$HN[which(Data_X$HN == 800 | Data_X$HN == 899 | Data_X$HN == 999)] <- 0
  if (length(which(Data_X$HN > 800)) != 0)  {
    Data_X$HN[which(Data_X$HN > 800)] <- sub('.', '', Data_X$HN[which(Data_X$HN > 800)])
    Data_X$HN <- as.numeric(Data_X$HN)
  }
  Data_X$HN <- as.numeric(Data_X$HN)
  Data_X$Date <- format(strptime(Data_X$Date,"%d/%m/%Y"), "%Y-%m-%d")
  
  Data_X2 <- as.data.frame(matrix(nrow=length(Date), ncol=3))
  colnames(Data_X2) <- c("Name","Date","HN")
  Data_X2$Name <- Meta_TN$Name[s]
  Data_X2$Date <- Date
  Data_X2$HN <- NA
  for (i in 1:length(Data_X$Date)) {
    Data_X2$HN[which(Data_X$Date[i] == Data_X2$Date)] <- Data_X$HN[i]
  }
  
  Data2_TN <- rbind(Data2_TN, Data_X2)
  
}



### BZ ###
Meta_BZ <- read.csv("00_Data/03_DataUp19-20/BZ/MetaData.csv", header=T, sep=";")

Data1_BZ <- vector("list", nrow(Meta_BZ))
names(Data1_BZ) <- Meta_BZ$Name
for (i in 1:nrow(Meta_BZ)) {
  if (i == 6 | i == 11 | i == 25) {
    Data1_BZ[[i]] <- read.csv(paste0("03_DataUp19-20/BZ/",Meta_BZ$Name[i],".csv"), header=T, skip=12, sep=";")
  } else {
    Data1_BZ[[i]] <- read.csv(paste0("03_DataUp19-20/BZ/",Meta_BZ$Name[i],".csv"), header=T, skip=12, sep="\t")
  }
}

Data2_BZ <- data.frame()
for (s in 1:length(Data1_BZ)) {
  
  Data_X <- Data1_BZ[[s]][,c(1,3)]
  Data_X <- Data_X[c(which(Data_X[,1] == "01/06/2019"):which(Data_X[,1] == "31/05/2020")),] 
  colnames(Data_X) <- c("Date","HN")
  rownames(Data_X) <- NULL
  
  Data_X$Date <- format(strptime(Data_X$Date,"%d/%m/%Y"), "%Y-%m-%d")
  Data_X$HN[which(Data_X$HN == "---")] <- NA
  Data_X$HN <- as.numeric(Data_X$HN)
  Data_X <- cbind(Meta_BZ$Name[s], Data_X)
  colnames(Data_X)[1] <- "Name"
  
  Data2_BZ <- rbind(Data2_BZ, Data_X)
  
}



### TN + BZ ###
Meta <- rbind(Meta_TN, Meta_BZ)
Meta <- Meta[order(Meta$Name),]
rownames(Meta) <- NULL
Data <- rbind(Data2_TN, Data2_BZ)
Data <- Data[order(Data$Name),]
rownames(Data) <- NULL



### MONTHLY AGGREGATION ###
Meta_EURAC <- read.csv("00_Data/01_ElaboratedData/EURAC_MetaData.csv", header=T, sep=";")
Date <- c("2019-06","2019-07","2019-08","2019-09","2019-10","2019-11","2019-12",
          "2020-01","2020-02","2020-03","2020-04","2020-05")

Data <- cbind(sapply(strsplit(Data[,2],"-"),"[",1),
              sapply(strsplit(Data[,2],"-"),"[",2),
              sapply(strsplit(Data[,2],"-"),"[",3),
              Data)
colnames(Data)[(1:3)] <- c("Year","Month","Day")
Data <- Data[,c(4,1:3,5:6)]

Data_M <- as.data.frame(matrix(nrow=length(Date), ncol=length(Meta_EURAC$Name)+2))
colnames(Data_M) <- c("Year","Month",Meta_EURAC$Name)
Data_M$Year <- as.numeric(sapply(strsplit(Date,"-"),"[",1))
Data_M$Month <- as.numeric(sapply(strsplit(Date,"-"),"[",2))

for (i in 1:length(Meta$Name)) {
  Data_X <- Data[which(Meta$Name[i] == Data$Name),]
  Data_X$Date2 <- paste(Data_X$Year, Data_X$Month, sep="-")
  Data_X <- Data_X[,c(1:5,7,6)]
  
  Data_X2 <- as.data.frame(matrix(nrow=length(Date), ncol=3))
  colnames(Data_X2) <- c("Name","Date","HN")
  Data_X2$Name <- Meta$Name[i]
  Data_X2$Date <- Date
  
  k=1
  for (j in 1:length(Date)) {
    Data_X3 <- Data_X$HN[which(Date[j] == Data_X$Date2)]
    Data_X3_NA <- (sum(is.na(Data_X3))/length(Data_X3))*100
    if (Data_X3_NA <= 30) {
      Data_X2$HN[k] <- sum(Data_X3, na.rm=T)
      k=k+1
    } else {
      k=k+1
    }
  }
  
  Data_M[,which(colnames(Data_M) == Meta$Name[i])] <- Data_X2$HN
}

### Adding MTAA Data ###
Data_MTAA <- read.csv("00_Data/01_ElaboratedData/MTAA_Data_HN.csv", header=T, sep=";")
Data_MTAA <- head(Data_MTAA, -7)
Data_MTAA <- tail(Data_MTAA, 12)
rownames(Data_MTAA) <- NULL
Data_M <- cbind(Data_M, Data_MTAA[,-c(1:2)])
Data_M2 <- Data_M[,-(1:2)]
Data_M2 <- Data_M2[,order(colnames(Data_M2))]
Data_M3 <- cbind(Data_M[,(1:2)],Data_M2)



### FULL DATA UP ###
Data_M4 <- read.csv("00_Data/02_FullData/FULL_Data_HN.csv", header=T, sep=";")
Data_M4 <- head(Data_M4, -7)
Data_M4 <- rbind(Data_M4,Data_M3)

# write.table(Data_M4, "00_Data/03_DataUp19-20/FULL_DataUp_HN.csv", sep=";", row.names=F)





##### ----- 04 GAP FILLING ----- #####
rm(list = ls(all.names = TRUE))



### DATA ###
dati <- read.csv("00_Data/03_DataUp19-20/FULL_DataUp_HN.csv", header=T, sep=";") # N+2 columns (year, month, N stations)
meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";") # N rows
N <- length(meta[,1]) # stations in dati and meta need to be in the same order
mat_filled <- dati # matrix were the filled series will be saved
common_years <- 10 # minimum number of years with common data between test and reference series
min.weight <- 0.00005 # minimum weight for reference station selection



### WEIGHTS MATRIX ###
# Computing station weights based on distance and elevation difference between station pairs #
Rt <- 6371
conv <- pi/180
dist_r <- matrix(nrow=length(meta[,1]), ncol=N)
delta_z <- matrix(nrow=length(meta[,1]), ncol=N)
wh <- matrix(nrow=length(meta[,1]), ncol=N)

for (ii in 1:N) {
  angle <- sin(conv*meta$Latitude) * sin(conv*meta$Latitude[ii]) + cos(conv*meta$Latitude) * cos(conv*meta$Latitude[ii]) * cos(abs(conv*meta$Longitude-conv*meta$Longitude[ii])) 
  angle[which(angle > 1)] <- 1
  dist_r[,ii] <- Rt*acos(angle)
  delta_z[,ii] <- abs(meta$Elevation-meta$Elevation[ii])
}

# Halving coefficients for weights: 20km for distance and 200m for elevation difference #
halv_dist <- 20
halv_elev <- 200
c_dist <- (halv_dist)^2/log(2)
c_h <- (halv_elev^2)/log(2)
wd <- exp(-(dist_r^2)/c_dist)

for (t in 1:N) {
  wh[,t] <- exp(-(delta_z[,t]^2)/c_h)
}

pesi <- wd*wh
diag(pesi) <- NA 



### GAP FILLING ###
for (i in 1:N) { # iter on stations
  print(paste("station",meta$Name[i],i))
  NMAX <- 10 # max number of reference stations used for the reconstruction 
  missing <- which(is.na(dati[,i+2]))
  
  
  if (length(missing) > 0) {
    d <- which(pesi[,i] > min.weight) # only stations with weights above the threshold are considered 
    col_ref <- length(d)
    
    if (col_ref > 2) { # stations are sorted by decreasing weights (the reconstruction is performed only if the selected stations are more than two)
      l <- sort.int(as.numeric(pesi[d,i]), decreasing=TRUE, index.return=TRUE)
      d <- d[l$ix]
      
      if (col_ref < NMAX) {NMAX <- col_ref} # if selected stations are less than NMAX, NMAX is set to the number of stations with relevant weight
      else {NMAX <- 10} 
      
      
      for (k in 1:12) { # iter on months
        mes <- which(dati[,2] == k) # select data for month k
        vet_test <- dati[mes,i+2] # extract the values of the test series
        vet_filled <- vet_test # create a duplicate of reference series for the month k where the filled gaps will be added
        gaps <- which(is.na(vet_test)) # select the years with missing data
        
        if (length(gaps) > 0) {
          mat_ref <- dati[mes,d+2] # matrix of the sorted reference series for the month k
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

# write.table(mat_filled2, "00_Data/04_GapFilledData/GapFilled_Data_HN.csv", sep=";", row.names=F)



### NOTE ###
# Code given by Michael Matiu (EURAC)





##### ----- 05 SEASONAL DATA ----- #####
rm(list = ls(all.names = TRUE))

library(zoo)

Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";")
HN <- read.csv("00_Data/04_GapFilledData/GapFilled_Data_HN.csv", header=T, sep=";")
P <- read.csv("00_Data/00_OriginalData/EURAC/TNBZMTAA_Data_P.csv", header=T, sep=";")
TMEAN <- read.csv("00_Data/00_OriginalData/EURAC/TNBZMTAA_Data_TMEAN.csv", header=T, sep=";")

dt <- seq(7,280,7)



### HN ###
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

# write.table(HN_TOT, "00_Data/05_SeasonalData/Data_HN_S.csv", row.names=F, sep=";")



### P ###
P_TOT <- as.data.frame(matrix(nrow=40, ncol=length(Meta$Name)+1))
colnames(P_TOT) <- c("Date", Meta$Name)
P_TOT$Date <- paste(seq(min(P$Year),max(P$Year)-1),seq(min(P$Year)+1,max(P$Year)),sep="-")

for (i in 1:length(Meta$Name)) {
  P_X <- P[,c(1:2,i+2)]
  P_X <- P_X[which(P_X$Month == 10 | P_X$Month == 11 | P_X$Month == 12 |
                     P_X$Month == 1 | P_X$Month == 2 |
                     P_X$Month == 3 | P_X$Month == 4),]
  P_X <- P_X[-c(1:4),]
  rownames(P_X) <- NULL
  
  for (j in 1:length(dt)) {
    P_X2 <- P_X[,3][(dt[j]-6):dt[j]]
    NA_X <- (sum(is.na(P_X2))/length(P_X2))*100
    if (NA_X <= 30) {
      P_TOT[j,i+1] <- round(sum(P_X2, na.rm=T), 0)
    }
  }
}

# write.table(P_TOT, "00_Data/05_SeasonalData/Data_P_S.csv", row.names=F, sep=";")



### TMEAN ###
TMEAN_TOT <- as.data.frame(matrix(nrow=40, ncol=length(Meta$Name)+1))
colnames(TMEAN_TOT) <- c("Date", Meta$Name)
TMEAN_TOT$Date <- paste(seq(min(TMEAN$Year),max(TMEAN$Year)-1),seq(min(TMEAN$Year)+1,max(TMEAN$Year)),sep="-")

for (i in 1:length(Meta$Name)) {
  TMEAN_X <- TMEAN[,c(1:2,i+2)]
  TMEAN_X <- TMEAN_X[which(TMEAN_X$Month == 10 | TMEAN_X$Month == 11 | TMEAN_X$Month == 12 |
                             TMEAN_X$Month == 1 | TMEAN_X$Month == 2 |
                             TMEAN_X$Month == 3 | TMEAN_X$Month == 4),]
  TMEAN_X <- TMEAN_X[-c(1:4),]
  rownames(TMEAN_X) <- NULL
  
  for (j in 1:length(dt)) {
    TMEAN_X2 <- TMEAN_X[,3][(dt[j]-6):dt[j]]
    NA_X <- (sum(is.na(P_X2))/length(TMEAN_X2))*100
    if (NA_X <= 30) {
      TMEAN_TOT[j,i+1] <- round(mean(TMEAN_X2, na.rm=T), 1)
    }
  }
}

# write.table(TMEAN_TOT, "00_Data/05_SeasonalData/Data_TMEAN_S.csv", row.names=F, sep=";")





##### ----- 06 DATA ANALYSIS ----- #####
rm(list = ls(all.names = TRUE))



### PRELIMINARY DATA ###
# Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";")
# Meta$Start_Y <- NA
# Meta$End_Y <- NA
# 
# Data1 <- read.csv("00_Data/00_OriginalData/EURAC/TNBZ-data-HNHS.csv", header=T, sep=",")
# Data1 <- cbind(as.numeric(sapply(strsplit(Data1[,2],"-"),"[",1)),
#                as.numeric(sapply(strsplit(Data1[,2],"-"),"[",2)),
#                as.numeric(sapply(strsplit(Data1[,2],"-"),"[",3)),
#                Data1)
# colnames(Data1)[(1:3)] <- c("Year","Month","Day")
# Data1 <- Data1[,c(4,1:3,6)]
# 
# Data2 <- read.csv("00_Data/00_OriginalData/MTAA/MTAA-data.csv", header=T, sep=";")
# Data2$Day <- 1 
# Data2 <- Data2[,c(1:3,5,4)]
# 
# Data <- rbind(Data1,Data2)
# Data <- Data[order(Data$Name),]
# 
# for (i in 1:length(Meta$Name)) {
#   Data_X <- Data[which(Meta$Name[i] == Data$Name),]
#   Meta$Start_Y[i] <- min(Data_X$Year)
#   Meta$End_Y[i] <- max(Data_X$Year)
# }
# 
# write.table(Meta, "00_Data/06_DataAnalysis/MetaPeriod.csv", sep=";", row.names=F)



### DATA ANALYSIS ###
# Number of time series per altitude #
MetaPeriod <- read.csv("00_Data/06_DataAnalysis/MetaPeriod.csv", header=T, sep=";")

jpeg("00_Data/06_DataAnalysis/Plots/Dataset_v1.jpeg", width=1200, height=1000, res=200)
par(mar=c(3,3,0.5,0.5), mgp=c(1.8,0.4,0))
hist(MetaPeriod$Elevation, breaks=seq(0,3000,250), ylim=c(0,25),
     main=NULL, xlab="Elevation [m]", ylab="Number of stations", col="black", border="white", cex.axis=0.8, cex.lab=1.2)
dev.off()

# Number of time series per year #
TS <- data.frame()
for (i in 1:length(MetaPeriod$Name)) {
  Y1_X <- MetaPeriod$Start_Y[i]
  Y2_X <- MetaPeriod$End_Y[i]
  Y_X <- seq(Y1_X,Y2_X)
  
  TS_X <- as.data.frame(cbind(MetaPeriod$Name[i],Y_X))
  colnames(TS_X) <- c("Name","Year")
  
  TS <- rbind(TS, TS_X)
}

TS <- as.data.frame(table(TS$Year))
colnames(TS) <- c("Year","Freq")

jpeg("00_Data/06_DataAnalysis/Plots/Dataset_v2.jpeg", width=1200, height=1000, res=200)
par(mar=c(3,3,0.5,0.5), mgp=c(1.8,0.4,0))
plot(as.numeric(TS$Year), TS$Freq, xaxt="n", yaxt="n",
     xlab="Year", ylab="Number of stations", type="l", cex.axis=0.8, cex.lab=1.2, lend=1, lwd=1.5)
axis(side=1, at=round(seq(1,nrow(TS),length.out=10),0),
     labels=round(seq(min(MetaPeriod$Start_Y),max(MetaPeriod$End_Y),length.out=10),0),
     las=1, cex.axis=0.8, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(0,max(TS$Freq),length.out=8),0), labels=round(seq(0,max(TS$Freq),length.out=8),0),
     las=2, cex.axis=0.8, tck=-0.01, xaxs="i")
dev.off()



### NOTES ###
# MetaPeriod updated after checking manually data update to 2019/20





##### ----- 07 DATA COMPARISON ----- #####
rm(list = ls(all.names = TRUE))



### DATA BEFORE GAPFILLING ###
Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";")
Data10 <- read.csv("00_Data/07_DataComparison/Data_HN_10.csv", header=T, sep=";")
Data10_S <- read.csv("00_Data/07_DataComparison/Data_HN_10_S.csv", header=T, sep=";")
Data30 <- read.csv("00_Data/07_DataComparison/Data_HN_30.csv", header=T, sep=";")
Data30_S <- read.csv("00_Data/07_DataComparison/Data_HN_30_S.csv", header=T, sep=";")
Data50 <- read.csv("00_Data/07_DataComparison/Data_HN_50.csv", header=T, sep=";")
Data50_S <- read.csv("00_Data/07_DataComparison/Data_HN_50_S.csv", header=T, sep=";")
m <- c(10:12,1:4)

# Data10 #
NumNA10 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
NumNA10[,1] <- Meta$Name
colnames(NumNA10) <- c("Name",month.name[m],"Seas")

for (i in 1:length(Meta$Name)) {
  Data10_X <- Data10[,c(1:2,i+2)]
  Data10_X <- Data10_X[-c(1:9),]
  rownames(Data10_X) <- NULL
  
  for (j in 1:length(m)) {
    Data10_X2 <- Data10_X[which(Data10_X$Month == m[j]),]
    Data10_X2 <- Data10_X2[,3]
    rownames(Data10_X2) <- NULL
    
    NumNA10[i,j+1] <- round((sum(is.na(Data10_X2))/length(Data10_X2))*100, 1)
    
  }
  
  NumNA10[,9] <- round((colSums(is.na(Data10_S[,-1]))/nrow(Data10_S))*100, 1)
  
}

# Data30 #
NumNA30 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
NumNA30[,1] <- Meta$Name
colnames(NumNA30) <- c("Name",month.name[m],"Seas")

for (i in 1:length(Meta$Name)) {
  Data30_X <- Data30[,c(1:2,i+2)]
  Data30_X <- Data30_X[-c(1:9),]
  rownames(Data30_X) <- NULL
  
  for (j in 1:length(m)) {
    Data30_X2 <- Data30_X[which(Data30_X$Month == m[j]),]
    Data30_X2 <- Data30_X2[,3]
    rownames(Data30_X2) <- NULL
    
    NumNA30[i,j+1] <- round((sum(is.na(Data30_X2))/length(Data30_X2))*100, 1)
    
  }
  
  NumNA30[,9] <- round((colSums(is.na(Data30_S[,-1]))/nrow(Data30_S))*100, 1)
  
}

# Data50 #
NumNA50 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
NumNA50[,1] <- Meta$Name
colnames(NumNA50) <- c("Name",month.name[m],"Seas")

for (i in 1:length(Meta$Name)) {
  Data50_X <- Data50[,c(1:2,i+2)]
  Data50_X <- Data50_X[-c(1:9),]
  rownames(Data50_X) <- NULL
  
  for (j in 1:length(m)) {
    Data50_X2 <- Data50_X[which(Data50_X$Month == m[j]),]
    Data50_X2 <- Data50_X2[,3]
    rownames(Data50_X2) <- NULL
    
    NumNA50[i,j+1] <- round((sum(is.na(Data50_X2))/length(Data50_X2))*100, 1)
    
  }
  
  NumNA50[,9] <- round((colSums(is.na(Data50_S[,-1]))/nrow(Data50_S))*100, 1)
  
}

# Plotting #
NumM_NA <- as.data.frame(matrix(nrow=3, ncol=9))
NumM_NA[,1] <- c("NumM_NA10","NumM_NA30","NumM_NA50")
NumM_NA[1,c(2:9)] <- round(colMeans(NumNA10[,-1]), 1)
NumM_NA[2,c(2:9)] <- round(colMeans(NumNA30[,-1]), 1)
NumM_NA[3,c(2:9)] <- round(colMeans(NumNA50[,-1]), 1)
colnames(NumM_NA) <- c(NA,month.name[m],"Seas")

png(paste0("00_Data/00_DataComparison/Plots/DataComparison [BEFORE GAPFILLING].png"), width=1100, height=1000, res=100)
par(mar=c(4,5,0.5,0.5), mgp=c(3.5,1,0))

matplot(seq(1:8), t(NumM_NA[,-1]), xaxt="n", yaxt="n", xlab=NA, ylab="Mean_Num_NA [%]", cex.lab=1.5,
        ylim=c(0,100), xlim=c(0.5,8.5), type="b", pch=c(0,1,2), col=c("black","blue","red"), cex=2.5, lwd=2, lty=c(1,1,1))
axis(side=1, at=seq(1:8), labels=c(month.abb[m],"Seas"), las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(0,100,length.out=11), labels=seq(0,100,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topright", legend=rev(c("Mean_Num_NA [10%]","Mean_Num_NA [30%]","Mean_Num_NA [50%]")),
       pch=rev(c(0,1,2)), col=rev(c("black","blue","red")), x.intersp=1.5, y.intersp=1.5, pt.cex=2.5, pt.lwd=2, cex=1.2, bty="n")

dev.off()



### DATA AFTER GAPFILLING ###
Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";")
Data10 <- read.csv("00_Data/07_DataComparison/GapFilled_Data_HN_10.csv", header=T, sep=";")
Data10_S <- read.csv("00_Data/07_DataComparison/GapFilled_Data_HN_10_S.csv", header=T, sep=";")
Data30 <- read.csv("00_Data/07_DataComparison/GapFilled_Data_HN_30.csv", header=T, sep=";")
Data30_S <- read.csv("00_Data/07_DataComparison/GapFilled_Data_HN_30_S.csv", header=T, sep=";")
Data50 <- read.csv("00_Data/07_DataComparison/GapFilled_Data_HN_50.csv", header=T, sep=";")
Data50_S <- read.csv("00_Data/07_DataComparison/GapFilled_Data_HN_50_S.csv", header=T, sep=";")
m <- c(10:12,1:4)

# Data10 #
NumNA10 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
NumNA10[,1] <- Meta$Name
colnames(NumNA10) <- c("Name",month.name[m],"Seas")

for (i in 1:length(Meta$Name)) {
  Data10_X <- Data10[,c(1:2,i+2)]
  Data10_X <- Data10_X[-c(1:9),]
  rownames(Data10_X) <- NULL
  
  for (j in 1:length(m)) {
    Data10_X2 <- Data10_X[which(Data10_X$Month == m[j]),]
    Data10_X2 <- Data10_X2[,3]
    rownames(Data10_X2) <- NULL
    
    NumNA10[i,j+1] <- round((sum(is.na(Data10_X2))/length(Data10_X2))*100, 1)
    
  }
  
  NumNA10[,9] <- round((colSums(is.na(Data10_S[,-1]))/nrow(Data10_S))*100, 1)
  
}

# Data30 #
NumNA30 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
NumNA30[,1] <- Meta$Name
colnames(NumNA30) <- c("Name",month.name[m],"Seas")

for (i in 1:length(Meta$Name)) {
  Data30_X <- Data30[,c(1:2,i+2)]
  Data30_X <- Data30_X[-c(1:9),]
  rownames(Data30_X) <- NULL
  
  for (j in 1:length(m)) {
    Data30_X2 <- Data30_X[which(Data30_X$Month == m[j]),]
    Data30_X2 <- Data30_X2[,3]
    rownames(Data30_X2) <- NULL
    
    NumNA30[i,j+1] <- round((sum(is.na(Data30_X2))/length(Data30_X2))*100, 1)
    
  }
  
  NumNA30[,9] <- round((colSums(is.na(Data30_S[,-1]))/nrow(Data30_S))*100, 1)
  
}

# Data50 #
NumNA50 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
NumNA50[,1] <- Meta$Name
colnames(NumNA50) <- c("Name",month.name[m],"Seas")

for (i in 1:length(Meta$Name)) {
  Data50_X <- Data50[,c(1:2,i+2)]
  Data50_X <- Data50_X[-c(1:9),]
  rownames(Data50_X) <- NULL
  
  for (j in 1:length(m)) {
    Data50_X2 <- Data50_X[which(Data50_X$Month == m[j]),]
    Data50_X2 <- Data50_X2[,3]
    rownames(Data50_X2) <- NULL
    
    NumNA50[i,j+1] <- round((sum(is.na(Data50_X2))/length(Data50_X2))*100, 1)
    
  }
  
  NumNA50[,9] <- round((colSums(is.na(Data50_S[,-1]))/nrow(Data50_S))*100, 1)
  
}

# Plotting #
NumM_NA <- as.data.frame(matrix(nrow=3, ncol=9))
NumM_NA[,1] <- c("NumM_NA10","NumM_NA30","NumM_NA50")
NumM_NA[1,c(2:9)] <- round(colMeans(NumNA10[,-1]), 1)
NumM_NA[2,c(2:9)] <- round(colMeans(NumNA30[,-1]), 1)
NumM_NA[3,c(2:9)] <- round(colMeans(NumNA50[,-1]), 1)
colnames(NumM_NA) <- c(NA,month.name[m],"Seas")

png(paste0("00_DataComparison/Plots/DataComparison [AFTER GAPFILLING].png"), width=1100, height=1000, res=100)
par(mar=c(4,5,0.5,0.5), mgp=c(3.5,1,0))

matplot(seq(1:8), t(NumM_NA[,-1]), xaxt="n", yaxt="n", xlab=NA, ylab="Mean_Num_NA [%]", cex.lab=1.5,
        ylim=c(0,100), xlim=c(0.5,8.5), type="b", pch=c(0,1,2), col=c("black","blue","red"), cex=2.5, lwd=2, lty=c(1,1,1))
axis(side=1, at=seq(1:8), labels=c(month.abb[m],"Seas"), las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(0,100,length.out=11), labels=seq(0,100,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topright", legend=rev(c("Mean_Num_NA [10%]","Mean_Num_NA [30%]","Mean_Num_NA [50%]")),
       pch=rev(c(0,1,2)), col=rev(c("black","blue","red")), x.intersp=1.5, y.intersp=1.5, pt.cex=2.5, pt.lwd=2, cex=1.2, bty="n")

dev.off()



### NOTES ###
# An analysis of different data aggregation before and after gap filling is conducted. Those data are
# produced using NA thresholds of 10, 30 and 50% --> e.g. if one month has > 10% of days with NA, all
# the month is set to NA. The same for seasonal aggregation. In the end it was decided to use 30%.





##### ----- 08 DATA QUALITY CHECK ----- #####
rm(list = ls(all.names = TRUE))

library(trend)
library(sp)
library(raster)



### DATA RECONSTRUCTION ###
dati <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2.csv", header=T, sep=";") # N+2 columns (year, month, N stations)
meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";") # N rows
N <- length(meta[,1]) # stations in dati and meta need to be in the same order
mat_simulated <- dati
mat_simulated[,-(1:2)] <- NA
common_years <- 10 # minimum number of years with common data between test and reference series
min.weight <- 0.00005 # minimum weight for reference station selection

# WEIGHTS MATRIX #
# Computing station weights based on distance and elevation difference between station pairs #
Rt <- 6371
conv <- pi/180
dist_r <- matrix(nrow=length(meta[,1]), ncol=N)
delta_z <- matrix(nrow=length(meta[,1]), ncol=N)

for (ii in 1:N) {
  angle <- sin(conv*meta$Latitude) * sin(conv*meta$Latitude[ii]) + cos(conv*meta$Latitude) * cos(conv*meta$Latitude[ii]) * cos(abs(conv*meta$Longitude-conv*meta$Longitude[ii])) 
  angle[which(angle > 1)] <- 1
  dist_r[,ii] <- Rt*acos(angle)
  delta_z[,ii] <- abs(meta$Elevation-meta$Elevation[ii])
}

# Halving coefficients for weights: 20km for distance and 200m for elevation difference #
halv_dist <- 20
halv_elev <- 200

c_dist <- (halv_dist^2)/log(2)
wd <- exp(-(dist_r^2)/c_dist)

c_h <- (halv_elev^2)/log(2)
wh <- matrix(nrow=length(meta[,1]), ncol=N)
for (t in 1:N) {
  wh[,t] <- exp(-(delta_z[,t]^2)/c_h)
}

pesi <- wd*wh
diag(pesi) <- NA 

# RECONSTRUCTION #
for (i in 1:N) { # iter on stations
  NMAX <- 10 # max number of reference stations used for the reconstruction 
  print(paste(i,"station",meta$Name[i]))
  present <- which(!is.na(dati[,i+2]))  
  
  if (length(present) > 0) {
    d <- which(pesi[,i] > min.weight) # only stations with weights above the threshold are considered 
    col_ref <- length(d)
    
    if (col_ref > 2) { # stations are sorted by decreasing weights (the reconstruction is performed only if the selected stations are more than two)
      l <- sort.int(as.numeric(pesi[d,i]), decreasing=TRUE, index.return=TRUE)
      d <- d[l$ix]
      
      if (col_ref < NMAX) {NMAX <- col_ref} # if selected stations are less than NMAX, NMAX is set to the number of stations with relevant weight
      else {NMAX <- 10} 
      
      for (k in 1:12) { # iter on months
        mes <- which(dati[,2] == k) # select data for month k
        vet_test <- dati[mes,i+2] # extract the values of the test series
        vet_test2 <- vet_test
        data <- which(!is.na(vet_test))
        
        if (length(data) > 0) {
          mat_ref <- dati[mes,d+2] # matrix of the sorted reference series for the month k
          stat <- matrix(nrow=length(data), ncol=NMAX)
          simulated <- numeric(length(data))
          simulated[] <- NA
          
          for (k1 in 1:length(data)) { # iter on data
            NN <- 1 
            j <- 1
            
            while (NN<=NMAX && j<=col_ref) { # cycle on reference stations until NMAX is reached or until the total number of available reference stations is exploited
              if (!is.na(mat_ref[data[k1],j])) { # not missing value in the reference 
                sel <- which(!is.na(vet_test) & !is.na(mat_ref[,j])) # years with common data
                remove <- which(sel == data[k1])
                
                if (length(remove) > 0) {sel <- sel[-remove]}
                
                if (length(sel) >= common_years) {
                  sum_test <- sum(vet_test[sel])
                  sum_ref <- sum(mat_ref[sel,j])
                  ratio_sum <- sum_test/sum_ref
                  
                  if (sum_ref == 0) {ratio_sum <- 1}
                  
                  stat[k1,NN] <- round((mat_ref[data[k1],j]*ratio_sum),4)
                  NN <- NN+1    
                }
              }
              j <- j+1 
            }
          }
          
          for (k2 in 1:length(data)) {
            simulated[k2] <- mean(stat[k2,],na.rm=TRUE)
          } 
          
          vet_test2[data] <- simulated
          mat_simulated[mes,i+2] <- vet_test2 
          
        }
      }
    }
  }        
}        

mat_simulated2 <- mat_simulated
for (a in 3:ncol(mat_simulated2)) {
  mat_simulated2[,a][which(is.na(mat_simulated2[,a]))] <- NA
  if (is.numeric(mat_simulated2[,a]) == T) {mat_simulated2[,a] <- round(mat_simulated2[,a], digits=0)}
}



### SEASONAL DATA ###
HN <- mat_simulated2
dt <- seq(7,280,7)

HN_TOT <- as.data.frame(matrix(nrow=40, ncol=length(meta$Name)+1))
colnames(HN_TOT) <- c("Date", meta$Name)
HN_TOT$Date <- paste(seq(min(HN$Year),max(HN$Year)-1),seq(min(HN$Year)+1,max(HN$Year)),sep="-")

for (i in 1:length(meta$Name)) {
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

# write.table(HN, "00_Data/08_DataQualityCheck/Data_HN_30_v2_simu.csv", sep=";", row.names=F)
# write.table(HN_TOT, "00_Data/08_DataQualityCheck/Data_HN_30_v2_simu_S.csv", sep=";", row.names=F)



### PLOTTING ###
Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";") 
Data_HN <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2.csv", header=T, sep=";")
Data_HN_s <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_simu.csv", header=T, sep=";")
Data_HN_S <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_S.csv", header=T, sep=";")
Data_HN_s_S <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_simu_S.csv", header=T, sep=";")
m <- c(11:12,1:4)

for (i in 1:length(Meta$Name)) {
  Data_X <- cbind(Meta$Name[i], Data_HN[,c(1:2,i+2)], Data_HN_s[,i+2])
  Data_X <- Data_X[-c(1:9),]
  colnames(Data_X)[c(1,4:5)] <- c("Name","Obs","Simu")
  rownames(Data_X) <- NULL
  
  Data_XS <- cbind(Meta$Name[i], Data_HN_S[,c(1,i+1)], Data_HN_s_S[,i+1])
  colnames(Data_XS)[c(1,3:4)] <- c("Name","Obs","Simu")
  rownames(Data_XS) <- NULL
  
  if (sum(is.na(Data_X$Obs)) == length(Data_X$Obs)) {
    
    print(paste(i,"station",Meta$Name[i],"NA"))
    i <- i+1
    
  } else {
    
    print(paste(i,"station",Meta$Name[i]))
    
    # Obs vs Simu 1 #
    max1 <- max(max(Data_X$Obs,na.rm=T),max(Data_X$Simu,na.rm=T))
    png(paste0("00_Data/08_DataQualityCheck/Plots/Obs_vs_Simu_1_v2/",Meta$Name[i],".png"), width=3000, height=2000, res=100)
    par(mar=c(8,10,5,3), mgp=c(7,2.5,0), mfrow=c(2,3))
    
    for (j in 1:length(m)) {
      Data_X2 <- Data_X[which(Data_X$Month == m[j]),]
      rownames(Data_X2) <- NULL
      
      matplot(Data_X2$Year, cbind(Data_X2$Obs,Data_X2$Simu), type="b", col=c("black","red"), pch=c(19,19),
              xlim=c(1979,2021), ylim=c(0,max1),
              main=month.abb[m[j]], xlab=NA, ylab="HN [cm]", cex.main=4, cex.lab=3, cex.axis=3, cex=4, lwd=3, las=1)
      
    }
    
    dev.off()
    
    # Obs vs Simu 2 #
    max2 <- max(max(Data_XS$Obs,na.rm=T),max(Data_XS$Simu,na.rm=T))
    if (max2 == "-Inf") {max2 <- 100}
    png(paste0("00_Data/08_DataQualityCheck/Plots/Obs_vs_Simu_2_v2/",Meta$Name[i],".png"), width=2000, height=1000, res=100)
    par(mar=c(6,7,5,3), mgp=c(4,1,0), mfrow=c(1,2))
    
    plot(Data_X$Obs, Data_X$Simu, type="p", pch=1, cex.lab=1.5, cex.axis=1.5, cex=2, lwd=2, las=1,
         main="Monthly Obs VS Simu Data", cex.main=2,
         xlim=c(0,max1), xlab="Observations [cm]",
         ylim=c(0,max1), ylab="Simulations [cm]")
    abline(0, 1, col="red", lwd=3)
    
    matplot(seq(1:40), cbind(Data_XS$Obs,Data_XS$Simu), type="b", col=c("black","red"), pch=c(19,19),
            xlim=c(0,42), ylim=c(0,max2), main="Seasonal HN series", xaxt="n", xlab="Season",
            ylab="HN [cm]", cex.main=2, cex.lab=1.5, cex.axis=1.5, cex=1.5, lwd=2, las=1)
    axis(side=1, at=seq(1,41,5), labels=c(Data_XS$Date[seq(1,40,5)],"2020-2021"), las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
    
    dev.off()
    
  }
  
}



### STATS ###
# Data #
Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";") 
Data_HN <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2.csv", header=T, sep=";")
Data_HN_s <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_simu.csv", header=T, sep=";")
Data_HN_S <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_S.csv", header=T, sep=";")
Data_HN_s_S <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_simu_S.csv", header=T, sep=";")

m <- c(12,1,2)

Stats_HN <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=7))
Stats_HN[,1] <- Meta$Name
Stats_HN[,2] <- Meta$Elevation
colnames(Stats_HN) <- c("Name","Elevation","BIAS","BIAS_r","MAE","MAE_r","r^2")

# Calculation #
for (i in 1:length(Meta$Name)) {
  Data_X <- cbind(Meta$Name[i], Data_HN[,c(1:2,i+2)], Data_HN_s[,i+2])
  Data_X <- Data_X[-c(1:9),]
  Data_X <- Data_X[which(Data_X$Month == 12 | Data_X$Month == 1 | Data_X$Month == 2),]
  colnames(Data_X)[c(1,4:5)] <- c("Name","Obs","Simu")
  rownames(Data_X) <- NULL
  
  if (sum(is.na(Data_X$Obs)) == length(Data_X$Obs) || sum(is.na(Data_X$Simu)) == length(Data_X$Simu)) {i <- i+1}
  else {
    Stats_HN[i,3] <- round(mean(Data_X$Simu - Data_X$Obs, na.rm=T), 5)
    Stats_HN[i,4] <- round(Stats_HN[i,3]/mean(Data_X$Obs, na.rm=T), 5)*100
    Stats_HN[i,5] <- round(mean(abs(Data_X$Simu - Data_X$Obs), na.rm=T), 5)
    Stats_HN[i,6] <- round(Stats_HN[i,5]/mean(Data_X$Obs, na.rm=T), 5)*100
    Stats_HN[i,7] <- round(summary(lm(Data_X$Obs ~ Data_X$Simu, Data_X))$r.squared, 5)
  }
}

Stats_Dataset <- as.data.frame(matrix(nrow=1, ncol=6))
Stats_Dataset[1,1] <- "Stats_Dataset"
Stats_Dataset[1,-1] <- round(colMeans(Stats_HN[,-c(1:2)], na.rm=T), 5)
colnames(Stats_Dataset) <- c("Name","BIAS","BIAS_r","MAE","MAE_r","r^2")

# Plotting #
MinMax <- as.data.frame(matrix(nrow=2, ncol=5))
MinMax[1,] <- c(-4.5,-12,0,0,0)
MinMax[2,] <- c(4.5,12,50,100,1)
colnames(MinMax) <- c("BIAS","BIAS_r","MAE","MAE_r","r^2")
rownames(MinMax) <- c("Min","Max")

Unit <- c("[cm]","[%]","[cm]","[%]","")

for (j in 1:5) {
  png(paste0("00_Data/08_DataQualityCheck/Plots/Stats_v2/",colnames(Stats_HN)[j+2],".png"), width=1100, height=1000, res=100)
  par(mar=c(6,6,0.5,0.5), mgp=c(4.5,1,0))
  
  plot(Stats_HN$Elevation, Stats_HN[,j+2], xaxt="n", yaxt="n", ylim=c(MinMax[1,j],MinMax[2,j]), xlim=c(0,3000),
       xlab="Altitude [m]", ylab=paste(colnames(Stats_HN)[j+2], Unit[j], sep=" "), type="p",
       pch=1, cex=2, lwd=2.5, cex.axis=1.5, cex.lab=1.5)
  
  if (j <= 2) {abline(h=0, col="red", lty=2, lwd=3)}
  
  axis(side=1, at=seq(0,3000,500), labels=seq(0,3000,500), las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
  axis(side=2, at=seq(MinMax[1,j],MinMax[2,j],length.out=11), labels=seq(MinMax[1,j],MinMax[2,j],length.out=11),
       las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
  
  dev.off()
}

# write.table(Stats_HN, "00_Data/08_DataQualityCheck/Stats_HN_v2.csv", sep=";", row.names=F)
# write.table(Stats_Dataset, "00_Data/08_DataQualityCheck/Stats_Dataset_v2.csv", sep=";", row.names=F)



### ANALYSIS ###
Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";") 
Data_HN <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2.csv", header=T, sep=";")
Data_HN_s <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_simu.csv", header=T, sep=";")
Stats_HN <- read.csv("00_Data/08_DataQualityCheck/Stats_HN_v2.csv", header=T, sep=";")
colnames(Stats_HN)[7] <- "r^2"

# Data_X <- cbind(Meta$Name[i], Data_HN[,c(1:2,i+2)], Data_HN_s[,i+2])
# Data_X <- Data_X[-c(1:9),]
# Data_X <- Data_X[which(Data_X$Month == 12 | Data_X$Month == 1 | Data_X$Month == 2),]
# colnames(Data_X)[c(1,4:5)] <- c("Name","Obs","Simu")
# rownames(Data_X) <- NULL

T_BIAS_r <- Stats_HN[which(Stats_HN$BIAS_r >= 4),]
rownames(T_BIAS_r) <- NULL
T_MAE_r <- Stats_HN[which(Stats_HN$MAE_r >= 45),]
rownames(T_MAE_r) <- NULL
T_r2 <- Stats_HN[which(Stats_HN$`r^2` <= 0.5),]
rownames(T_r2) <- NULL
T_TOT <- Stats_HN[which(Stats_HN$BIAS_r >= 4 & Stats_HN$MAE_r >= 45 & Stats_HN$`r^2` <= 0.5),]
rownames(T_TOT) <- NULL



### PETTITT TEST ###
# Data #
Meta <- read.csv("00_Data/02_FullData/FULL_MetaData.csv", header=T, sep=";") 
Data_HN <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2.csv", header=T, sep=";")
Data_HN_S <- read.csv("00_Data/08_DataQualityCheck/Data_HN_30_v2_S.csv", header=T, sep=";")

m <- rep(c("01","02","12"), times=40)
y <- rep(1980:2020, each=3)
Date <- paste(m, y, sep="-")
Date <- Date[-c(1:2)]
Date <- head(Date, -1)

CP_HN <- cbind(Meta,NA,NA)
colnames(CP_HN)[c(6:7)] <- c("CP","Delta")

TAA <- shapefile("FinalData/TopographicData/TAA.shp")
TAA <- spTransform(TAA, crs("+proj=longlat +datum=WGS84"))
DEM <- raster("FinalData/TopographicData/DEM_TAA_100m.tif")
Hillshade <- hillShade(terrain(DEM, "slope", "radians", 8),
                       terrain(DEM, "aspect", "radians", 8),
                       angle=45, direction=315, normalize=T)

# Calculation and plotting #
for (i in 1:length(Meta$Name)) {
  
  Data_X <- Data_HN[,c(1:2,i+2)]
  Data_X <- Data_X[-c(1:9),]
  Data_X <- Data_X[which(Data_X$Month == 12 | Data_X$Month == 1 | Data_X$Month == 2),]
  Data_X <- as.data.frame(cbind(Date, Data_X))
  rownames(Data_X) <- NULL
  
  if (sum(is.na(Data_X[,4]))/length(Data_X[,4])*100 > 30) {i <- i+1}
  else {
    
    Data_X1 <- Data_X
    Data_X2 <- na.omit(Data_X1)
    rownames(Data_X2) <- NULL
    
    HN1 <- Data_X1[,4]
    HN2 <- Data_X2[,4]
    
    Ptest <- pettitt.test(HN2)
    CP1 <- as.numeric(which(Data_X1$Date == Data_X2$Date[Ptest$estimate]))
    CP2 <- as.numeric(Ptest$estimate)
    m1 <- round(mean(HN2[1:CP2]), 0)
    m2 <- round(mean(HN2[(CP2+1):length(HN2)]), 0)
    d <- m2-m1
    date <- Data_X2$Date[CP2]
    
    M <- Data_X$Month[CP1]
    Y <- Data_X$Year[CP1]
    if (M < 10) {CP_HN[i,6] <- paste(Y,M,sep="0")} else {CP_HN[i,6] <- paste0(Y,M)}
    CP_HN[i,7] <- d
    
    HN1_plot <- HN1
    HN1_plot[which(is.na(HN1_plot))] <- 0.3
    bCol <- array()
    bCol[which(HN1_plot != 0.3)] <- "#000000"
    bCol[which(HN1_plot == 0.3)] <- "#FF0000"
    
    png(paste0("00_Data/08_DataQualityCheck/Plots/CP/",Meta$Name[i],".png"), width=1100, height=1000, res=100)
    par(mar=c(6,6,1.5,1.5), mgp=c(4,1,0))

    plot(HN1_plot, xaxt="n", xlab="Date", ylab="HN [cm]", type="h", col=bCol, ylim=c(0,max(HN2)+50),
         cex.axis=1.5, cex.lab=1.5, lend=1, lwd=2)
    axis(side=1, at=c(1,32,63,94,120), labels=Date[c(1,32,63,94,120)], las=1, xaxs="i", cex.axis=1.5)
    lines(c(1:CP1), rep(m1,CP1), col="blue", lwd=3)
    segments(x0=CP1, x1=CP1, y0=m1, y1=m2, col="red", lwd=3)
    lines(c(CP1:length(HN1)), rep(m2,(length(HN1)-CP1)+1), col="blue", lwd=3)
    points(CP1, m1, type="p", pch=4, col="black", cex=2.5, lwd=2)
    legend("topleft", col=c("blue","blue","red","black"), pch=c(NA,NA,NA,4), lty=c(1,1,1,NA),
           lwd=c(2,2,2,2), pt.cex=2, cex=1.2, bty="n",
           legend=c(paste0("Mean_1: ",round(m1,0)," cm"),
                    paste0("Mean_2: ",round(m2,0)," cm"),
                    paste0("Delta: ",round(d,0)," cm"),
                    paste0("Changing point: ",date)))

    dev.off()
    
  }
}

CP_HN$CP <- as.numeric(CP_HN$CP)

# write.table(CP_HN, "00_Data/08_DataQualityCheck/CP_HN.csv", sep=";", row.names=F)

# Plotting Altitude CP Distribution #
png(paste0("00_Data/08_DataQualityCheck/Plots/CP/Altitude_CP.png"), width=1100, height=1000, res=100)
par(mar=c(6,6,0.5,0.5), mgp=c(4.5,1,0))

plot(CP_HN$CP, CP_HN$Elevation, xaxt="n", yaxt="n", xlim=c(198011,202003), ylim=c(0,3000),
     xlab="Date", ylab="Elevation [m]", type="p", pch=1, cex=2, lwd=2.5, cex.axis=1.5, cex.lab=1.5)
axis(side=1, at=round(seq(198012,202002,length.out=5), 0), labels=c("12-1980","10-1990","07-2000","04-2010","02-2020"),
     las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(0,3000,length.out=7), labels=seq(0,3000,length.out=7),
     las=2, cex.axis=1.5, tck=-0.01, xaxs="i")

dev.off()

# Plotting Spatial CP Distribution #
I <- c(198012,198501,199501,200501,201501,202002)
CPal <- c("red","orange","green","blue","violet")
I_M <- c("12/1980-01/1985","02/1985-01/1995","02/1995-01/2005","02/2005-01/2015","02/2015-02/2020")

png(paste0("00_Data/08_DataQualityCheck/Plots/CP/Spatial_CP.png"), width=1100, height=1000, res=100)
par(mar=c(3,3,2,0), mgp=c(0,1,0))

IClass <- cut(CP_HN$CP, breaks=I)
Colors <- CPal[as.numeric(IClass)]
Pt <- array()
Pt[which(CP_HN$Delta < 0)] <- 25
Pt[which(CP_HN$Delta > 0)] <- 24

plot(Hillshade, col=gray(0:255/255), legend=F, axes=F, box=T, alpha=0.6)
plot(TAA, lwd=1, border="black", add=T)
points(CP_HN$Longitude, CP_HN$Latitude, type="p", pch=Pt, col=Colors, bg=Colors, cex=2.5, lwd=3)
legend("topleft", legend=rev(I_M), col=rev(CPal), lwd=3, cex=1.5, bty="n")
axis(side=1, at=round(seq(Hillshade@extent@xmin,Hillshade@extent@xmax,length.out=5), 2),
     labels=paste0(round(seq(Hillshade@extent@xmin,Hillshade@extent@xmax,length.out=5), 2), "°"),
     las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(Hillshade@extent@ymin,Hillshade@extent@ymax,length.out=5), 2),
     labels=paste0(round(seq(Hillshade@extent@ymin,Hillshade@extent@ymax,length.out=5), 2), "°"),
     las=3, cex.axis=1.5, tck=-0.01, xaxs="i")

dev.off()





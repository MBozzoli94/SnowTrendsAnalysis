library(trend)
library(sp)
library(raster)
library(gplots)
library(stats)
library(car)
library(reshape2)
library(ggplot2)



##### DATA #####
Meta <- read.csv("FinalData/MetaData.csv", header=T, sep=";")
Data <- read.csv("FinalData/Data_HN.csv", header=T, sep=";")
Data_S <- read.csv("FinalData/Data_HN_S.csv", header=T, sep=";")

m <- c(10:12,1:4)
ms <- c(month.abb[m],"Seas")
a <- 0.05
ar <- c(0,1000,2000,3000)

TAA <- shapefile("FinalData/TopographicData/TAA.shp")
TAA <- spTransform(TAA, crs("+proj=longlat +datum=WGS84"))
DEM <- raster("FinalData/TopographicData/DEM_TAA_100m.tif")
Hillshade <- hillShade(terrain(DEM, "slope", "radians", 8),
                       terrain(DEM, "aspect", "radians", 8),
                       angle=45, direction=315, normalize=T)




##### GENERAL TREND #####
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



##### PLOTTING #####
### Spatial distribution plot ###
I <- c(-30,-10,-1,1,10,30)
PtType <- c(6,6,0,2,2)
CPal <- c("#FF0000","#ff8000","#666666","#2395ff","#0000FF")

png(paste0("01_TrendHN/Plots/TrendHN_1_v2.png"), width=2400, height=1200, res=100)
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

### Scatter plot ###
# Adjusting the dataset #
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

# Plotting #
png(paste0("01_TrendHN/Plots/TrendHN_2.png"), width=2000, height=1000, res=100)

ggplot(TS_s_TOT2, aes(Elevation, Value)) +
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
        panel.grid.minor = element_blank())

dev.off()

### Box plot ###
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png(paste0("01_TrendHN/Plots/TrendHN_3.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-30,30), xlim=c(0.5,12.5), xlab=NA, ylab="Trend [cm/decade]", cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)

for (j in 1:(length(m)+1)) {
  
  p <- j+j*0.5
  
  NS[,2] <- rbind(sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)])),
                  sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)])),
                  sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)])))
  
  boxplot(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
          TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
          TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
          add=T, boxwex=0.25, at=c(p-0.3,p,p+0.3), col=c("#a65f02","#04b604","#01a2ff"), xaxt="n", yaxt="n", plot=F)
  text(c(p-0.3,p,p+0.3), -30, labels=NS[,2], cex=1, font=3)
  
}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-30,30,length.out=11), labels=seq(-30,30,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
       fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -30, labels="N.S.:", cex=1, font=3)

dev.off()

### Box plot (relative) ###
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png(paste0("01_TrendHN/Plots/TrendHN_3r.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-32,30), xlim=c(0.5,12.5), xlab=NA, ylab="Relative Trend [%]", cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)

for (j in 1:(length(m)+1)) {
  
  p <- j+j*0.5
  
  NS[,2] <- rbind(sum(!is.na(TS_s_r[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)])),
                  sum(!is.na(TS_s_r[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)])),
                  sum(!is.na(TS_s_r[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)])))
  
  boxplot(TS_s_r[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
          TS_s_r[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
          TS_s_r[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
          add=T, boxwex=0.25, at=c(p-0.3,p,p+0.3), col=c("#a65f02","#04b604","#01a2ff"), xaxt="n", yaxt="n")
  text(c(p-0.3,p,p+0.3), -32, labels=NS[,2], cex=1, font=3)
  
}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-30,30,length.out=11), labels=seq(-30,30,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
       fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -32, labels="N.S.:", cex=1, font=3)

dev.off()



##### TABLE #####
TTab <- vector("list", 8)
names(TTab) <- ms
TTab_X <- as.data.frame(matrix(nrow=3, ncol=6))
TTab_X[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(TTab_X) <- c("Altitude [m]","#","Median [cm]","Lci [cm]","Uci [cm]","Relative [%]")

for (j in 1:(length(m)+1)) {
  
  BP_X <- boxplot(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
                  TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
                  TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
                  plot=F)
  BP_r_X <- boxplot(TS_s_r[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
                    TS_s_r[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
                    TS_s_r[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
                    plot=F)
  TTab_X[,2] <- BP_X[["n"]]
  TTab_X[,3] <- BP_X[["stats"]][3,]
  TTab_X[,4] <- BP_X[["conf"]][1,]
  TTab_X[,5] <- BP_X[["conf"]][2,]
  TTab_X[,6] <- BP_r_X[["stats"]][3,]
  
  TTab[[j]] <- TTab_X
  
}



##### TABLE V2 #####
TTab_v2 <- vector("list", 8)
names(TTab_v2) <- ms
TTab_X_v2 <- as.data.frame(matrix(nrow=3, ncol=6))
TTab_X_v2[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(TTab_X_v2) <- c("Altitude [m]","#","Mean [cm/decade]","Lci [cm/decade]","Uci [cm/decade]","Relative [%]")

for (j in 1:length(ms)) {
  
  for (k in 1:3) {
    
    TTab_X_v2[k,2] <- sum(!is.na(TS_s[,j+1][which(Meta$Elevation > ar[k] & Meta$Elevation <= ar[k+1])]))
    TTab_X_v2[k,3] <- round(mean(TS_s[,j+1][which(Meta$Elevation > ar[k] & Meta$Elevation <= ar[k+1])], na.rm=T), 2)
    TTab_X_v2[k,4] <- round(mean(TS_lci[,j+1][which(Meta$Elevation > ar[k] & Meta$Elevation <= ar[k+1])], na.rm=T), 2)
    TTab_X_v2[k,5] <- round(mean(TS_uci[,j+1][which(Meta$Elevation > ar[k] & Meta$Elevation <= ar[k+1])], na.rm=T), 2)
    TTab_X_v2[k,6] <- round(100*(mean(TS_s[,j+1][which(Meta$Elevation > ar[k] & Meta$Elevation <= ar[k+1])]/
                                      HN_avg[,j+1][which(Meta$Elevation > ar[k] & Meta$Elevation <= ar[k+1])],
                                 na.rm=T)), 2)
    
  }
  
  TTab_v2[[j]] <- TTab_X_v2
  
}





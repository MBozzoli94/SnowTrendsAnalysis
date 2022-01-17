library(trend)
library(sp)
library(raster)
library(gplots)
library(stats)
library(car)
library(reshape2)
library(ggplot2)



##### ----- DATA ----- #####
Meta <- read.csv("FinalData/MetaData.csv", header=T, sep=";")
Var <- c("HN","P","TMEAN")

Data <- vector("list", 3)
for (i in 1:length(Var)) {
  Data[[i]] <- read.csv(paste0("FinalData/Data_",Var[i],".csv"), header=T, sep=";")
}
names(Data) <- Var

Data_S <- vector("list", 3)
for (i in 1:length(Var)) {
  Data_S[[i]] <- read.csv(paste0("FinalData/Data_",Var[i],"_S.csv"), header=T, sep=";")
}
names(Data_S) <- Var

m <- c(10:12,1:4)
ms <- c(month.abb[m],"Seas")
a <- 0.05

TAA <- shapefile("FinalData/TopographicData/TAA.shp")
TAA <- spTransform(TAA, crs("+proj=longlat +datum=WGS84"))
DEM <- raster("FinalData/TopographicData/DEM_TAA_100m.tif")
Hillshade <- hillShade(terrain(DEM, "slope", "radians", 8),
                       terrain(DEM, "aspect", "radians", 8),
                       angle=45, direction=315, normalize=T)





##### ----- TREND FULL ----- #####
### TREND CALCULATION ###
TS <- vector("list", 8)
names(TS) <- ms
TS_X <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=4))
TS_X[,1] <- Meta$Name
colnames(TS_X) <- c("Name",Var)
TS_min <- TS
TS_min_X <- TS_X
TS_max <- TS
TS_max_X <- TS_X
Pv <- TS
Pv_X <- TS_X
Pt <- TS
Pt_X <- TS_X

for (j in 1:length(ms)) {
  
  for (i in 1:length(Meta$Name)) {
    
    if (j == 8) {
      
      Data_X2 <- cbind(Data_S[["HN"]][c(1,i+1)], Data_S[["P"]][i+1], Data_S[["TMEAN"]][i+1])
      colnames(Data_X2)[2:4] <- Var
      
    } else {
      
      Data_X <- cbind(Data[["HN"]][c(1:2,i+2)], Data[["P"]][i+2], Data[["TMEAN"]][i+2])
      Data_X <- Data_X[-c(1:9),]
      Data_X <- head(Data_X,-1)
      colnames(Data_X)[3:5] <- Var
      rownames(Data_X) <- NULL
      
      Data_X2 <- Data_X[which(Data_X$Month == m[j]),]
      rownames(Data_X2) <- NULL
      
    }
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 0) >= 5) {
      
      for (k in 1:length(Var)) {
        
        if (j == 8) {
          TS_X[i,k+1] <- 10*round(sens.slope(Data_X2[,k+1], conf.level=0.95)[[1]][[1]], digits=3)
          TS_min_X[i,k+1] <- 10*round(sens.slope(Data_X2[,k+1], conf.level=0.95)[["conf.int"]][[1]], digits=3)
          TS_max_X[i,k+1] <-10*round(sens.slope(Data_X2[,k+1], conf.level=0.95)[["conf.int"]][[2]], digits=3)
          Pv_X[i,k+1] <- round(mk.test(Data_X2[,k+1])[[2]], digits=5)
        } else {
          TS_X[i,k+1] <- 10*round(sens.slope(Data_X2[,k+2], conf.level=0.95)[[1]][[1]], digits=3)
          TS_min_X[i,k+1] <- 10*round(sens.slope(Data_X2[,k+2], conf.level=0.95)[["conf.int"]][[1]], digits=3)
          TS_max_X[i,k+1] <-10*round(sens.slope(Data_X2[,k+2], conf.level=0.95)[["conf.int"]][[2]], digits=3)
          Pv_X[i,k+1] <- round(mk.test(Data_X2[,k+2])[[2]], digits=5)
        }
        
        Pt_X[i,1+which(Pv_X[i,-1] <= a)] <- 19
        Pt_X[i,1+which(Pv_X[i,-1] > a)] <- 1
        
      }
      
    } else {
      
      TS_X[i,-1] <- NA
      TS_min_X[i,-1] <- NA
      TS_max_X[i,-1] <- NA
      Pv_X[i,-1] <- NA
      Pt_X[i,-1] <- NA
      
    }
    
  }
  
  TS[[j]] <- TS_X
  TS_min[[j]] <- TS_min_X
  TS_max[[j]] <- TS_max_X
  Pv[[j]] <- Pv_X
  Pt[[j]] <- Pt_X
  
}



### PLOTTING ###
# Plot v1 #
I <- as.data.frame(cbind(c(-30,-10,-1,1,10,30),
                         c(-90,-30,-3,3,30,90),
                         c(-1.5,-0.5,-0.1,0.1,0.5,1.5)))
colnames(I) <- Var
PtType <- c(6,6,0,2,2)
CPal <- c("#FF0000","#ff8000","#666666","#2395ff","#0000FF")

png("05_Attribution/FULL/Attribution_1.png", width=1800, height=2400, res=100) # Attribution_2
par(mar=c(0,5.5,5,0), mgp=c(3,0,0), mfrow=c(4,3))

for (j in 1:4) { # 5:8
  
  for (k in 1:length(Var)) {
    
    Title <- NA
    if (j == 1) {Title <- Var[k]} # 5
    
    YLab <- NA
    if (k == 1) {YLab <- ms[j]}
    
    IClass_X <- cut(TS[[j]][[k+1]], breaks=I[,k])
    Col_X <- CPal[as.numeric(IClass_X)]
    if (k == length(Var)) {Col_X <- rev(CPal)[as.numeric(IClass_X)]}
    p_X <- PtType[as.numeric(IClass_X)]
    p_X[which(TS[[j]][[k+1]] <= I[3,k] & Pv[[j]][[k+1]] <= a)] <- 25
    p_X[which(TS[[j]][[k+1]] > I[4,k] & Pv[[j]][[k+1]] <= a)] <- 24
    p_X[which(TS[[j]][[k+1]] > I[3,k] & TS[[j]][[k+1]] <= I[4,k] & Pv[[j]][[k+1]] <= a)] <- 22
    
    plot(Hillshade, col=gray(0:255/255), legend=F, axes=F, box=T, alpha=0.6, main=Title, cex.main=3,
         ylab=YLab, cex.lab=3, font.lab=2)
    plot(TAA, lwd=1, border="black", add=T)
    points(Meta$Longitude, Meta$Latitude, type="p", pch=p_X, col=Col_X, bg=Col_X, cex=3, lwd=3)
    if (k == length(Var)) {
      legend("bottomright", legend=rev(levels(IClass_X)), pch=rev(PtType), col=CPal, pt.lwd=2, cex=2, bty="n")
    } else {
      legend("bottomright", legend=rev(levels(IClass_X)), pch=rev(PtType), col=rev(CPal), pt.lwd=2, cex=2, bty="n")
    }
    
  }
  
}

dev.off()

# Plot v2 #
# Adjusting the dataset #
TS_TOT <- as.data.frame(matrix(nrow=8*length(Meta$Name), ncol=6))
TS_TOT[,1] <- factor(rep(ms, each=length(Meta$Name)), levels=ms)
TS_TOT[,2] <- rep(Meta$Name, times=8)
TS_TOT[,3] <- rep(Meta$Elevation, times=8)
colnames(TS_TOT) <- c("Period","Name","Elevation",colnames(TS_X)[2:4])
TS_min_TOT <- TS_TOT
TS_max_TOT <- TS_TOT
Pt_TOT <- TS_TOT

dj <- c(length(Meta$Name)*c(1:length(ms)))

for (j in 1:length(ms)) {
  
  TS_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:6)] <- TS[[j]][-1]
  TS_min_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:6)] <- TS_min[[j]][-1]
  TS_max_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:6)] <- TS_max[[j]][-1]
  Pt_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:6)] <- Pt[[j]][-1]
  
}

TS_TOT <- na.omit(TS_TOT)
rownames(TS_TOT) <- NULL
TS_min_TOT <- na.omit(TS_min_TOT)
rownames(TS_min_TOT) <- NULL
TS_max_TOT <- na.omit(TS_max_TOT)
rownames(TS_max_TOT) <- NULL
Pt_TOT <- na.omit(Pt_TOT)
rownames(Pt_TOT) <- NULL

TS_TOT2 <- melt(TS_TOT[,-c(1:2)], id.vars="Elevation")
TS_TOT2 <- cbind(rep(TS_TOT$Period, times=3), rep(TS_TOT$Name, times=3), TS_TOT2)
colnames(TS_TOT2) <- c("Period","Name","Elevation","Var","Value")
levels(TS_TOT2$Var) <- Var
TS_TOT2$Value_min <- melt(TS_min_TOT[,-c(1:2)], id.vars="Elevation")$value
TS_TOT2$Value_max <- melt(TS_max_TOT[,-c(1:2)], id.vars="Elevation")$value
TS_TOT2$Pt <- melt(Pt_TOT[,-c(1:2)], id.vars="Elevation")$value

# Plotting #
png(paste0("05_Attribution/FULL/Attribution_3.png"), width=2000, height=1000, res=100)

ggplot(TS_TOT2, aes(Elevation, Value, colour=Var)) +
  geom_hline(yintercept=0, size=0.8, linetype="dashed") +
  geom_point(size=2, shape=TS_TOT2$Pt) +
  geom_linerange(aes(ymin=Value_min, ymax=Value_max), alpha=0.07, size=1) +
  facet_grid(rows=vars(Var), cols=vars(Period), scales="free") +
  scale_colour_manual(values=c("#000000","#3d77ff","#33D81B")) +
  scale_x_continuous(breaks=seq(0,3000,1000), labels=seq(0,3000,1000), limits=c(0,3000)) +
  labs(x="Elevation [m]", y="Trend") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
        strip.text.x = element_text(face="bold", size=15, vjust=1),
        strip.text.y = element_text(face="bold", size=15, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines"),
        panel.grid.minor = element_blank())

dev.off()

# Plot v3 #
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png("05_Attribution/FULL/Attribution_4.png", width=1800, height=2400, res=100) # Attribution_5
par(mar=c(0,8,5,1), mgp=c(5,1.5,0), mfrow=c(4,3))

for (j in 1:4) { # 5:8

  for (k in 1:length(Var)) {

    Title <- NA
    if (j == 1) {Title <- Var[k]} # 5

    YLab <- NA
    if (k == 1) {YLab <- ms[j]}

    NS[,2] <- rbind(sum(!is.na(TS[[j]][[k+1]][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)])),
                    sum(!is.na(TS[[j]][[k+1]][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)])),
                    sum(!is.na(TS[[j]][[k+1]][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)])))

    plot(NULL, xaxt="n", yaxt="n", ylim=c(min(I[,k]),max(I[,k])), xlim=c(0.5,1.5), main=Title, cex.main=3,
         xlab=NA, ylab=YLab, cex.lab=3, font.lab=2)
    abline(h=0, col="black", lty=2, lwd=3)
    boxplot(TS[[j]][[k+1]][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
            TS[[j]][[k+1]][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
            TS[[j]][[k+1]][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
            add=T, boxwex=0.15, at=c(0.75,1,1.25), col=c("#a65f02","#04b604","#01a2ff"), xaxt="n", yaxt="n")
    axis(side=2, at=seq(min(I[,k]),max(I[,k]),length.out=11), labels=seq(min(I[,k]),max(I[,k]),length.out=11),
         las=2, cex.axis=2, tck=-0.02, xaxs="i")
    legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
           fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=2, bty="n")
    text(0.6, min(I[,k]), labels="N.S.:", cex=2, font=3)
    text(c(0.75,1,1.25), min(I[,k]), labels=NS[,2], cex=2, font=3)

  }

}

dev.off()





##### ----- TREND P ----- #####
### GENERAL TREND ###
TS_s <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=(length(ms)+1)))
TS_s[,1] <- Meta$Name
colnames(TS_s) <- c("Name",ms)
TS_lci <- TS_s
TS_uci <- TS_s
MK_p <- TS_s
Pt <- TS_s
TS_s_r <- TS_s

for (j in 1:length(ms)) {
  
  for (i in 1:length(Meta$Name)) {
    
    if (j == 8) {
      
      Data_X2 <- cbind(Data_S[["HN"]][c(1,i+1)], Data_S[["P"]][i+1], Data_S[["TMEAN"]][i+1])
      colnames(Data_X2)[2:4] <- Var
      
    } else {
      
      Data_X <- cbind(Data[["HN"]][c(1:2,i+2)], Data[["P"]][i+2], Data[["TMEAN"]][i+2])
      Data_X <- Data_X[-c(1:9),]
      Data_X <- head(Data_X,-1)
      colnames(Data_X)[3:5] <- Var
      rownames(Data_X) <- NULL
      
      Data_X2 <- Data_X[which(Data_X$Month == m[j]),]
      rownames(Data_X2) <- NULL
      
    }
    
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 0) >= 5) {
      
      TS_s[i,j+1] <- 10*round(sens.slope(Data_X2$P, conf.level=0.95)[[1]][[1]], digits=3)
      TS_lci[i,j+1] <- 10*round(sens.slope(Data_X2$P, conf.level=0.95)[["conf.int"]][[1]], digits=3)
      TS_uci[i,j+1] <- 10*round(sens.slope(Data_X2$P, conf.level=0.95)[["conf.int"]][[2]], digits=3)
      MK_p[i,j+1] <- round(mk.test(Data_X2$P)[[2]], digits=5)
      TS_s_r[i,j+1] <- round(TS_s[i,j+1]/mean(Data_X2$P)*100, digits=0)
      
      if (MK_p[i,j+1] <= a) {Pt[i,j+1] <- 19} else {Pt[i,j+1] <- 1}
      
    } else {
      
      TS_s[i,j+1] <- NA
      TS_lci[i,j+1] <- NA
      TS_uci[i,j+1] <- NA
      MK_p[i,j+1] <- NA
      Pt[i,j+1] <- NA
      TS_s_r[i,j+1] <- NA
      
    }
    
  }
  
}



### PLOTTING ###
### Spatial distribution plot ###
I <- c(-90,-30,-3,3,30,90)
PtType <- c(6,6,0,2,2)
CPal <- c("#FF0000","#ff8000","#666666","#2395ff","#0000FF")

png(paste0("05_Attribution/P/TrendP_1.png"), width=2400, height=1200, res=100)
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
png(paste0("05_Attribution/P/TrendP_2.png"), width=2000, height=1000, res=100)

ggplot(TS_s_TOT2, aes(Elevation, Value)) +
  geom_hline(yintercept=0, size=1, linetype="dashed") +
  geom_point(size=3.5, colour="#3d77ff", shape=TS_s_TOT2$Pt) +
  geom_linerange(aes(ymin=Lci, ymax=Uci), alpha=0.25, size=0.5, colour="#3d77ff") +
  facet_wrap(~ Period, nrow=2, scales="free_y") +
  geom_point(aes(x=2000, y=80), alpha=0) +
  geom_point(aes(x=1000, y=-80), alpha=0) +
  scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
  labs(x="Elevation [m]", y="Trend [mm/decade]") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
        strip.text.x = element_text(face="bold", size=15, vjust=1),
        strip.text.y = element_text(face="bold", size=15, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,1,1.5,1.5), "lines"),
        panel.grid.minor = element_blank())

dev.off()

### Box plot ###
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png(paste0("05_Attribution/P/TrendP_3.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-30,60), xlim=c(0.5,12.5), xlab=NA, ylab="Trend [mm/decade]", cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)

for (j in 1:(length(m)+1)) {
  
  p <- j+j*0.5
  
  NS[,2] <- rbind(sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)])),
                  sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)])),
                  sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)])))
  
  boxplot(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
          TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
          TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
          add=T, boxwex=0.25, at=c(p-0.3,p,p+0.3), col=c("#a65f02","#04b604","#01a2ff"), xaxt="n", yaxt="n")
  text(c(p-0.3,p,p+0.3), -30, labels=NS[,2], cex=1, font=3)
  
}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-30,60,10), labels=seq(-30,60,10), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
       fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -30, labels="N.S.:", cex=1, font=3)

dev.off()

### Box plot (relative) ###
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png(paste0("05_Attribution/P/TrendP_3r.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-20,30), xlim=c(0.5,12.5), xlab=NA, ylab="Relative Trend [%]", cex.lab=1.5)
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
  text(c(p-0.3,p,p+0.3), -20, labels=NS[,2], cex=1, font=3)
  
}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-20,30,10), labels=seq(-20,30,10), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
       fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -20, labels="N.S.:", cex=1, font=3)

dev.off()





##### ----- TREND TMEAN ----- #####
### GENERAL TREND ###
TS_s <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=(length(ms)+1)))
TS_s[,1] <- Meta$Name
colnames(TS_s) <- c("Name",ms)
TS_lci <- TS_s
TS_uci <- TS_s
MK_p <- TS_s
Pt <- TS_s

for (j in 1:length(ms)) {
  
  for (i in 1:length(Meta$Name)) {
    
    if (j == 8) {
      
      Data_X2 <- cbind(Data_S[["HN"]][c(1,i+1)], Data_S[["P"]][i+1], Data_S[["TMEAN"]][i+1])
      colnames(Data_X2)[2:4] <- Var
      
    } else {
      
      Data_X <- cbind(Data[["HN"]][c(1:2,i+2)], Data[["P"]][i+2], Data[["TMEAN"]][i+2])
      Data_X <- Data_X[-c(1:9),]
      Data_X <- head(Data_X,-1)
      colnames(Data_X)[3:5] <- Var
      rownames(Data_X) <- NULL
      
      Data_X2 <- Data_X[which(Data_X$Month == m[j]),]
      rownames(Data_X2) <- NULL
      
    }
    
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 0) >= 5) {
      
      TS_s[i,j+1] <- 10*round(sens.slope(Data_X2$TMEAN, conf.level=0.95)[[1]][[1]], digits=3)
      TS_lci[i,j+1] <- 10*round(sens.slope(Data_X2$TMEAN, conf.level=0.95)[["conf.int"]][[1]], digits=3)
      TS_uci[i,j+1] <- 10*round(sens.slope(Data_X2$TMEAN, conf.level=0.95)[["conf.int"]][[2]], digits=3)
      MK_p[i,j+1] <- round(mk.test(Data_X2$TMEAN)[[2]], digits=5)
      
      if (MK_p[i,j+1] <= a) {Pt[i,j+1] <- 19} else {Pt[i,j+1] <- 1}
      
    } else {
      
      TS_s[i,j+1] <- NA
      TS_lci[i,j+1] <- NA
      TS_uci[i,j+1] <- NA
      MK_p[i,j+1] <- NA
      Pt[i,j+1] <- NA
      
    }
    
  }
  
}



### PLOTTING ###
### Spatial distribution plot ###
I <- c(-1.5,-0.5,-0.1,0.1,0.5,1.5)
PtType <- c(6,6,0,2,2)
CPal <- c("#FF0000","#ff8000","#666666","#2395ff","#0000FF")

png(paste0("05_Attribution/TMEAN/TrendTMEAN_1.png"), width=2400, height=1200, res=100)
par(mar=c(0,0,5,0), mgp=c(3,0,0), mfrow=c(2,4))

for (j in 1:(length(ms))) {
  
  IClass_X <- cut(TS_s[,j+1], breaks=I)
  Col_X <- rev(CPal)[as.numeric(IClass_X)]
  Pt_X <- PtType[as.numeric(IClass_X)]
  Pt_X[which(TS_s[,j+1] <= I[3] & MK_p[,j+1] <= a)] <- 25
  Pt_X[which(TS_s[,j+1] > I[4] & MK_p[,j+1] <= a)] <- 24
  Pt_X[which(TS_s[,j+1] > I[3] & TS_s[,j+1] <= I[4] & MK_p[,j+1] <= a)] <- 22
  
  plot(Hillshade, main=ms[j], cex.main=3, col=gray(0:255/255), legend=F, axes=F, box=T, alpha=0.6)
  plot(TAA, lwd=1, border="black", add=T)
  points(Meta$Longitude, Meta$Latitude, type="p", pch=Pt_X, col=Col_X, bg=Col_X, cex=3, lwd=3)
  legend("bottomright", legend=rev(levels(IClass_X)), pch=rev(PtType), col=CPal, pt.lwd=2, cex=2, bty="n")
  
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
png(paste0("05_Attribution/TMEAN/TrendTMEAN_2.png"), width=2000, height=1000, res=100)

ggplot(TS_s_TOT2, aes(Elevation, Value)) +
  geom_hline(yintercept=0, size=1, linetype="dashed") +
  geom_point(size=3.5, colour="#33D81B", shape=TS_s_TOT2$Pt) +
  geom_linerange(aes(ymin=Lci, ymax=Uci), alpha=0.25, size=0.5, colour="#33D81B") +
  facet_wrap(~ Period, nrow=2) +
  ylim(c(-1.5,1.5)) +
  scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
  labs(x="Elevation [m]", y="Trend [°C/decade]") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
        strip.text.x = element_text(face="bold", size=15, vjust=1),
        strip.text.y = element_text(face="bold", size=15, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,0.7,1.5,1.5), "lines"),
        panel.grid.minor = element_blank())

dev.off()

### Box plot ###
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png(paste0("05_Attribution/TMEAN/TrendTMEAN_3.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-0.4,1), xlim=c(0.5,12.5), xlab=NA, ylab="Trend [°C/decade]", cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)

for (j in 1:(length(m)+1)) {
  
  p <- j+j*0.5
  
  NS[,2] <- rbind(sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)])),
                  sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)])),
                  sum(!is.na(TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)])))
  
  boxplot(TS_s[,j+1][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
          TS_s[,j+1][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
          TS_s[,j+1][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
          add=T, boxwex=0.25, at=c(p-0.3,p,p+0.3), col=c("#a65f02","#04b604","#01a2ff"), xaxt="n", yaxt="n")
  text(c(p-0.3,p,p+0.3), -0.4, labels=NS[,2], cex=1, font=3)
  
}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-0.4,1,0.2), labels=seq(-0.4,1,0.2), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
       fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -0.4, labels="N.S.:", cex=1, font=3)

dev.off()








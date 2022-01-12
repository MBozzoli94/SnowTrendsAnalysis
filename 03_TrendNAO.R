library(pracma)
library(trend)
library(gplots)



##### DATA #####
### Monthly ###
Data1 <- read.csv("03_TrendNAO/Data.csv", header=T, sep=";")
Data2 <- cbind(read.fwf("03_TrendNAO/Data.csv", header=T, sep=";", widths=c(4,2),
                       colClasses="numeric", col.names=c("Year","Month")), 
               Data1[,2])
colnames(Data2)[3] <- "NAO"

### Seasonal ###
Data3 <- as.data.frame(matrix(nrow=40, ncol=2))
colnames(Data3) <- c("Date","NAO")
Data3$Date <- paste(seq(min(Data2$Year),max(Data2$Year)-1),seq(min(Data2$Year)+1,max(Data2$Year)),sep="-")

DataX <- Data2[which(Data2$Month == 10 | Data2$Month == 11 | Data2$Month == 12 |
                     Data2$Month == 1 | Data2$Month == 2 |
                     Data2$Month == 3 | Data2$Month == 4),]
DataX <- DataX[-c(1:4),]
rownames(DataX) <- NULL

dt <- seq(7,280,7)

for (j in 1:length(dt)) {
  DataX2 <- DataX[,3][(dt[j]-6):dt[j]]
  Data3[j,2] <- round(mean(DataX2), 2)
}



##### PLOTTING DATA #####
### Monthly ###
ax <- seq(1,481,length.out=5)
ax_lab <- c("1980-01","1990-01","2000-01","2010-01","2020-01")
bCol <- array()
bCol[which(Data2$NAO == 0)] <- "black"
bCol[which(Data2$NAO > 0)] <- "red"
bCol[which(Data2$NAO < 0)] <- "blue"

MAVG_10y <- movavg(Data2$NAO, n=120, type="e")

png(paste0("03_TrendNAO/Plots/NAO_8020.png"), width=1100, height=1000, res=100)
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3,1,0))

plot(Data2$NAO, xaxt="n", yaxt="n", ylim=c(-3.5,3.5), xlim=c(0,length(Data2$Month)+1),
     xlab="Date", ylab="NAO", type="h", col=bCol, lwd=1.5, cex.axis=1.6, cex.lab=2)
lines(MAVG_10y, col="black", lwd=3)
axis(side=1, at=ax, labels=ax_lab, las=1, cex.axis=1.6, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(-3.5,3.5,length.out=11),1), labels=round(seq(-3.5,3.5,length.out=11),1),
     las=2, cex.axis=1.6, tck=-0.01, xaxs="i")

dev.off()

### Seasonal ###
ax <- c(seq(1,40,10),40)
ax_lab <- c("1980-81","1990-91","2000-01","2010-11","2019-20")
bCol <- array()
bCol[which(Data3$NAO == 0.00)] <- "black"
bCol[which(Data3$NAO > 0)] <- "red"
bCol[which(Data3$NAO < 0)] <- "blue"

MAVG_10y <- movavg(Data3$NAO, n=10, type="e")

png(paste0("03_TrendNAO/Plots/NAO_8020_S.png"), width=1100, height=1000, res=100)
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3,1,0))

plot(Data3$NAO, xaxt="n", yaxt="n", ylim=c(-1.5,1.5), xlim=c(0,length(Data3$Date)+1),
     xlab="Date", ylab="NAO", type="h", col=bCol, lwd=3, cex.axis=1.5, cex.lab=1.5)
lines(MAVG_10y, col="black", lwd=3)
axis(side=1, at=ax, labels=ax_lab, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(-1.5,1.5,length.out=11),1), labels=round(seq(-1.5,1.5,length.out=11),1),
     las=2, cex.axis=1.5, tck=-0.01, xaxs="i")

dev.off()



##### TREND #####
DataM <- Data2[-c(1:9),]
DataM <- head(DataM,-1)
rownames(DataM) <- NULL
DataS <- Data3
m <- c(10:12,1:4)
a <- 0.05

Trend <- as.data.frame(matrix(nrow=4, ncol=8))
rownames(Trend) <- c("TS_s","TS_lci","TS_uci","MK_p")
colnames(Trend) <- c(month.abb[m],"Seas")

### Monthly ###
for (i in 1:length(m)) {
  DataM_X <- DataM[which(DataM$Month == m[i]),]
  DataM_X <- DataM_X[,3]
  
  Trend[1,i] <- 10*round(sens.slope(DataM_X, conf.level=0.95)[[1]][[1]], digits=3)
  Trend[2,i] <- 10*round(sens.slope(DataM_X, conf.level=0.95)[["conf.int"]][[1]], digits=3)
  Trend[3,i] <- 10*round(sens.slope(DataM_X, conf.level=0.95)[["conf.int"]][[2]], digits=3)
  Trend[4,i] <- round(mk.test(DataM_X)[[2]], digits=5)
}

### Seasonal ###
Trend[1,8] <- 10*round(sens.slope(DataS[,2], conf.level=0.95)[[1]][[1]], digits=3)
Trend[2,8] <- 10*round(sens.slope(DataS[,2], conf.level=0.95)[["conf.int"]][[1]], digits=3)
Trend[3,8] <- 10*round(sens.slope(DataS[,2], conf.level=0.95)[["conf.int"]][[2]], digits=3)
Trend[4,8] <- round(mk.test(DataS[,2])[[2]], digits=5)



##### PLOTTING TREND #####
bCol <- array()
bCol[which(Trend[1,] == 0)] <- "black"
bCol[which(Trend[1,] > 0)] <- "red"
bCol[which(Trend[1,] < 0)] <- "blue"

png(paste0("03_TrendNAO/Plots/TrendNAO.png"), width=1100, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plotCI(seq(1:8), t(Trend[1,]), ui=t(Trend[3,]), li=t(Trend[2,]), gap=0, sfrac=0.01, barcol="grey",
       xaxt="n", yaxt="n", ylim=c(-0.5,0.5), xlim=c(0.5,8.5), xlab=NA, ylab="Trend NAO", type="p",
       pch=1, col=bCol, cex=2.5, lwd=2.5, cex.axis=1.5, cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)
axis(side=1, at=seq(1:8), labels=c(month.abb[m],"Seas"), las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-0.5,0.5,length.out=11), labels=seq(-0.5,0.5,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")

dev.off()





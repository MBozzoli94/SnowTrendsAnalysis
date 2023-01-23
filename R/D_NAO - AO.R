library(pracma)
library(zoo)



##### NAO - DATA #####
### Monthly ###
Data1 <- read.csv("FinalData/Data_NAO_raw.csv", header=T, sep=";")
Data2 <- cbind(read.fwf("FinalData/Data_NAO_raw.csv", header=T, sep=";", widths=c(4,2),
                       colClasses="numeric", col.names=c("Year","Month")),
               Data1[,2])
Data2 <- Data2[which(Data2$Month == 10 | Data2$Month == 11 | Data2$Month == 12 |
                     Data2$Month == 1 | Data2$Month == 2 |
                     Data2$Month == 3 | Data2$Month == 4),]
Data2 <- Data2[-c(1:4),]
rownames(Data2) <- NULL
colnames(Data2)[3] <- "NAO"

### Seasonal ###
Data3 <- as.data.frame(matrix(nrow=50, ncol=2))
colnames(Data3) <- c("Date","NAO")
Data3$Date <- paste(seq(min(Data2$Year),max(Data2$Year)-1),seq(min(Data2$Year)+1,max(Data2$Year)),sep="-")

dt <- seq(7,350,7)

for (j in 1:length(dt)) {
  DataX <- Data2[,3][(dt[j]-6):dt[j]]
  Data3[j,2] <- round(mean(DataX), 2)
}



##### NAO - PLOTTING DATA #####
### Monthly ###
ax <- c(1,71,141,211,280)
ax_lab <- c("1980-10","1990-10","2000-10","2010-10","2020-04")
bCol <- array()
bCol[which(Data2$NAO == 0)] <- "black"
bCol[which(Data2$NAO > 0)] <- "red"
bCol[which(Data2$NAO < 0)] <- "blue"

Data2$Mavg <- rollapply(Data2$NAO, 70, mean, partial=T, align="center")

png(paste0("D_NAO - AO/NAO_8020.png"), width=1100, height=1000, res=100)
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3,1,0))

plot(Data2$NAO[-(1:70)], xaxt="n", yaxt="n", ylim=c(-2.8,2.8), xlim=c(0,length(Data2$Month[-(1:70)])+1),
     xlab="Date", ylab="NAO", type="h", col=bCol[-(1:70)], lwd=1.5, cex.axis=1.6, cex.lab=2)
lines(Data2$Mavg[-(1:70)], col="black", lwd=3)
axis(side=1, at=ax, labels=ax_lab, las=1, cex.axis=1.6, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(-2.8,2.8,length.out=11),1), labels=round(seq(-2.8,2.8,length.out=11),1),
     las=2, cex.axis=1.6, tck=-0.01, xaxs="i")

dev.off()

### Seasonal ###
ax <- c(seq(1,40,10),40)
ax_lab <- c("1980-81","1990-91","2000-01","2010-11","2019-20")
bCol <- array()
bCol[which(Data3$NAO == 0.00)] <- "black"
bCol[which(Data3$NAO > 0)] <- "red"
bCol[which(Data3$NAO < 0)] <- "blue"

Data3$Mavg <- rollapply(Data3$NAO, 10, mean, partial=T, align="center")

png(paste0("D_NAO - AO/NAO_8020_S.png"), width=1100, height=1000, res=100)
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3,1,0))

plot(Data3$NAO[-(1:10)], xaxt="n", yaxt="n", ylim=c(-1.2,1.2), xlim=c(0,length(Data3$Date[-(1:10)])+1),
     xlab="Date", ylab="NAO", type="h", col=bCol[-(1:10)], lwd=3, cex.axis=1.6, cex.lab=2)
lines(Data3$Mavg[-(1:10)], col="black", lwd=3)
axis(side=1, at=ax, labels=ax_lab, las=1, cex.axis=1.6, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(-1.2,1.2,0.2),1), labels=round(seq(-1.2,1.2,0.2),1),
     las=2, cex.axis=1.6, tck=-0.01, xaxs="i")

dev.off()




##### AO - DATA #####
### Monthly ###
Data1 <- read.csv("FinalData/Data_AO_raw.csv", header=T, sep=";")
Data2 <- cbind(read.fwf("FinalData/Data_AO_raw.csv", header=T, sep=";", widths=c(4,2),
                        colClasses="numeric", col.names=c("Year","Month")),
               Data1[,2])
Data2 <- Data2[which(Data2$Month == 10 | Data2$Month == 11 | Data2$Month == 12 |
                       Data2$Month == 1 | Data2$Month == 2 |
                       Data2$Month == 3 | Data2$Month == 4),]
Data2 <- Data2[-c(1:4),]
rownames(Data2) <- NULL
colnames(Data2)[3] <- "AO"

### Seasonal ###
Data3 <- as.data.frame(matrix(nrow=50, ncol=2))
colnames(Data3) <- c("Date","AO")
Data3$Date <- paste(seq(min(Data2$Year),max(Data2$Year)-1),seq(min(Data2$Year)+1,max(Data2$Year)),sep="-")

dt <- seq(7,350,7)

for (j in 1:length(dt)) {
  DataX <- Data2[,3][(dt[j]-6):dt[j]]
  Data3[j,2] <- round(mean(DataX), 2)
}



##### AO - PLOTTING DATA #####
### Monthly ###
ax <- c(1,71,141,211,280)
ax_lab <- c("1980-10","1990-10","2000-10","2010-10","2020-04")
bCol <- array()
bCol[which(Data2$AO == 0)] <- "black"
bCol[which(Data2$AO > 0)] <- "red"
bCol[which(Data2$AO < 0)] <- "blue"

Data2$Mavg <- rollapply(Data2$AO, 70, mean, partial=T, align="center")

png(paste0("D_AO - AO/AO_8020.png"), width=1100, height=1000, res=100)
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3,1,0))

plot(Data2$AO[-(1:70)], xaxt="n", yaxt="n", ylim=c(-2.8,2.8), xlim=c(0,length(Data2$Month[-(1:70)])+1),
     xlab="Date", ylab="AO", type="h", col=bCol[-(1:70)], lwd=1.5, cex.axis=1.6, cex.lab=2)
lines(Data2$Mavg[-(1:70)], col="black", lwd=3)
axis(side=1, at=ax, labels=ax_lab, las=1, cex.axis=1.6, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(-2.8,2.8,length.out=11),1), labels=round(seq(-2.8,2.8,length.out=11),1),
     las=2, cex.axis=1.6, tck=-0.01, xaxs="i")

dev.off()

### Seasonal ###
ax <- c(seq(1,40,10),40)
ax_lab <- c("1980-81","1990-91","2000-01","2010-11","2019-20")
bCol <- array()
bCol[which(Data3$AO == 0.00)] <- "black"
bCol[which(Data3$AO > 0)] <- "red"
bCol[which(Data3$AO < 0)] <- "blue"

Data3$Mavg <- rollapply(Data3$AO, 10, mean, partial=T, align="center")

png(paste0("D_AO - AO/AO_8020_S.png"), width=1100, height=1000, res=100)
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3,1,0))

plot(Data3$AO[-(1:10)], xaxt="n", yaxt="n", ylim=c(-1.2,1.2), xlim=c(0,length(Data3$Date[-(1:10)])+1),
     xlab="Date", ylab="AO", type="h", col=bCol[-(1:10)], lwd=3, cex.axis=1.6, cex.lab=2)
lines(Data3$Mavg[-(1:10)], col="black", lwd=3)
axis(side=1, at=ax, labels=ax_lab, las=1, cex.axis=1.6, tck=-0.01, xaxs="i")
axis(side=2, at=round(seq(-1.2,1.2,0.2),1), labels=round(seq(-1.2,1.2,0.2),1),
     las=2, cex.axis=1.6, tck=-0.01, xaxs="i")

dev.off()





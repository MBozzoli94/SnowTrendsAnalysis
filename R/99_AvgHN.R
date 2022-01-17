library(raster)
library(trend)
library(TeachingDemos)



##### DATA #####
Meta <- read.csv("FinalData/MetaData.csv", header=T, sep=";")
Data <- read.csv("FinalData/Data_HN_S.csv", header=T, sep=";")
Var <- c("HN","P","TMEAN")

TAA <- shapefile("FinalData/TopographicData/TAA.shp")
TAA <- spTransform(TAA, crs("+proj=longlat +datum=WGS84"))
DEM <- raster("FinalData/TopographicData/DEM_TAA_100m.tif")
Hillshade <- hillShade(terrain(DEM, "slope", "radians", 8),
                       terrain(DEM, "aspect", "radians", 8),
                       angle=45, direction=315, normalize=T)



##### AVERAGE SEASONAL DATA #####
Data_Avg <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=6))
Data_Avg[,1:5] <- Meta
colnames(Data_Avg) <- c(colnames(Meta),"HN")

for (i in 1:length(Meta$Name)) {
  
  if (sum(is.na(Data[,i+1]))/length(Data[,i+1])*100 <= 30) {
    Data_Avg[i,6] <- round(mean(Data[,i+1],na.rm=T), 0)
  } else {
    Data_Avg[i,6] <- NA
  }
  
}



##### PLOTTING #####
### Spatial plot ###
Data_Avg3 <- na.omit(Data_Avg)
rownames(Data_Avg3) <- NULL
Data_Avg3$Pos_text <- 3
Data_Avg3$Pos_text[c(8,21,22,45,52,62,66,71)] <- 4
Data_Avg3$Pos_text[c(6,30,43,64,67,69)] <- 2
Data_Avg3$Pos_text[c(13,42)] <- 1

# write.table(Data_Avg3, "99_AvgHN/Data_HN_Avg.csv", sep=";", row.names=F)

# CPal <- c("#D3D3D3","#00eeff","#01c8ff","#2395ff","#3d77ff","#2962ff","#0134ff","#2400da","#5500a4","#7a00b2",
#           "#9a00b2","#a800c2","#c400e2","#e601ff","#ff01ee","#ff01bc","#ff019e","#ff0184")
# CSeq <- c(0,5,10,15,20,30,40,50,60,80,100,120,150,200,250,300,400,500,600)
# CPal <- c("#c7d9f5","#8bbff2","#398cff","#3263df","#214691","#008017","#46c846","#90f76d","#fff701","#dcd42c",
#           "#f18d35","#f94943","#911119","#6a1794","#e601ff","#ff0184")
# CSeq <- c(0,5,10,15,20,30,40,50,60,80,100,120,150,200,300,400,600)
# Data_Avg3$Col <- CPal[as.numeric(cut(Data_Avg3$HN, breaks=CSeq))]

jpeg("99_AvgHN/Plots/AvgHN_v1_2.jpeg", width=1200, height=1000, res=200)
par(mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0))
plot(Hillshade, col=gray(0:255/255), legend=F, axes=F, box=T)
plot(TAA, lwd=1, border="black", add=T)
points(Data_Avg3$Longitude, Data_Avg3$Latitude, type="p", pch=22, col="white", bg=Data_Avg3$Col, cex=0.6, lwd=1.5)
shadowtext(Data_Avg3$Longitude, Data_Avg3$Latitude, labels=as.character(Data_Avg3$HN), col=Data_Avg3$Col, bg="white",
           r=0.2, pos=Data_Avg3$Pos, offset=0.25, cex=0.5)
axis(side=1, at=round(seq(Hillshade@extent@xmin,Hillshade@extent@xmax,length.out=5), 2),
     labels=paste0(round(seq(Hillshade@extent@xmin,Hillshade@extent@xmax,length.out=5), 2), "°"),
     las=1, cex.axis=0.6, tck=-0.02, xaxs="i")
axis(side=2, at=round(seq(Hillshade@extent@ymin,Hillshade@extent@ymax,length.out=5), 2),
     labels=paste0(round(seq(Hillshade@extent@ymin,Hillshade@extent@ymax,length.out=5), 2), "°"),
     las=1, cex.axis=0.6, tck=-0.02, xaxs="i")
dev.off()

### Scatter plot ###
Data_Avg2 <- Data_Avg[order(Data_Avg$Elevation),]
Data_Avg2 <- na.omit(Data_Avg2)
rownames(Data_Avg2) <- NULL
Ptest <- pettitt.test(Data_Avg2$HN)
CP <- as.numeric(Ptest$estimate)

LR_HN <- lm(Data_Avg2$HN ~ Data_Avg2$Elevation)
LR_HN_i <- LR_HN[["coefficients"]][[1]]
LR_HN_s <- 100*LR_HN[["coefficients"]][[2]]
LR_HN_se <- 100*summary(LR_HN)$"coefficients"[2,2]

LR_HN1 <- lm(Data_Avg2$HN[1:CP] ~ Data_Avg2$Elevation[1:CP])
LR_HN_i1 <- LR_HN1[["coefficients"]][[1]]
LR_HN_s1 <- 100*LR_HN1[["coefficients"]][[2]]
LR_HN_se1 <- 100*summary(LR_HN1)$"coefficients"[2,2]
X1 <- c(Data_Avg2$Elevation[1],Data_Avg2$Elevation[CP])
Y1 <- c(LR_HN_s1/100*X1[1]+LR_HN_i1,LR_HN_s1/100*X1[2]+LR_HN_i1)

LR_HN2 <- lm(Data_Avg2$HN[(CP+1):length(Data_Avg2$HN)] ~ Data_Avg2$Elevation[(CP+1):length(Data_Avg2$HN)])
LR_HN_i2 <- LR_HN2[["coefficients"]][[1]]
LR_HN_s2 <- 100*LR_HN2[["coefficients"]][[2]]
LR_HN_se2 <- 100*summary(LR_HN2)$"coefficients"[2,2]
X2 <- c(Data_Avg2$Elevation[CP+1],Data_Avg2$Elevation[length(Data_Avg2$HN)])
Y2 <- c(LR_HN_s2/100*X2[1]+LR_HN_i2,LR_HN_s2/100*X2[2]+LR_HN_i2)

png("99_AvgHN/Plots/AvgHN_v2_2.png", width=1550, height=1200, res=100)
par(mar=c(7.5,7.5,0.5,0.5), mgp=c(5.5,1.3,0))

plot(Data_Avg$Elevation, Data_Avg$HN, xaxt="n", yaxt="n", ylim=c(0,600), xlim=c(0,3000),
     xlab="Elevation [m]", ylab="Average total seasonal HN [cm]", type="p", pch=1, cex=2.5, lwd=2, cex.axis=2, cex.lab=2.4)
# abline(LR_HN, col="red", lty=1, lwd=4)
# lines(X1, Y1, col="blue", lty=1, lwd=4)
# lines(X2, Y2, col="green", lty=1, lwd=4)
abline(v=Data_Avg2$Elevation[CP], col="red", lty=2, lwd=5)
axis(side=1, at=seq(0,3000,500), labels=seq(0,3000,500), las=1, cex.axis=2, tck=-0.01, xaxs="i")
axis(side=2, at=seq(0,600,50), labels=seq(0,600,50), las=2, cex.axis=2, tck=-0.01, xaxs="i")
legend("topleft", legend=paste0("CP: ",Data_Avg2$Elevation[CP]," m"), col="red", lty=2, lwd=4, cex=1.8, bty="n")
# legend("topleft", col=c("blue","green","red"), lty=c(1,1,2), lwd=c(3,3,4), cex=1.5, bty="n", y.intersp=1.2,
#        legend=c(# paste0("y = ",round(LR_HN_s,3),"*x - ",abs(round(LR_HN_i,3))),
#                 paste0("y = ",round(LR_HN_s1,3),"*x - ",abs(round(LR_HN_i1,3))),
#                 paste0("y = ",round(LR_HN_s2,3),"*x - ",abs(round(LR_HN_i2,3))),
#                 paste0("CP: ",Data_Avg2$Elevation[CP]," m")))

dev.off()

### Box Plot ###
NS <- as.data.frame(matrix(nrow=3, ncol=2))
NS[,1] <- c("(0 - 1000]","(1000 - 2000]","(2000 - 3000]")
colnames(NS) <- c("Altitude [m]","N_Stations")

png("99_AvgHN/Plots/AvgHN_v3.png", width=1100, height=1000, res=100)
par(mar=c(0.5,7,0.5,0.5), mgp=c(5,1,0))

NS[,2] <- rbind(sum(!is.na(Data_Avg[,6][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)])),
                sum(!is.na(Data_Avg[,6][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)])),
                sum(!is.na(Data_Avg[,6][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)])))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-50,600), xlim=c(0.5,1.5), xlab=NA, ylab="Average seasonal HN [cm]", cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)
boxplot(Data_Avg[,6][which(Meta$Elevation > 0 & Meta$Elevation <= 1000)],
        Data_Avg[,6][which(Meta$Elevation > 1000 & Meta$Elevation <= 2000)],
        Data_Avg[,6][which(Meta$Elevation > 2000 & Meta$Elevation <= 3000)],
        add=T, boxwex=0.1, at=c(0.75,1,1.25), col=c("#a65f02","#04b604","#01a2ff"), xaxt="n", yaxt="n")
axis(side=2, at=seq(0,600,50), labels=seq(0,600,50), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=rev(c("(0 - 1000] m","(1000 - 2000] m","(2000 - 3000] m")),
       fill=rev(c("#a65f02","#04b604","#01a2ff")), x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.6, -40, labels="N.S.:", cex=1.5, font=3)
text(c(0.75,1,1.25), -40, labels=NS[,2], cex=1.5, font=3)

dev.off()





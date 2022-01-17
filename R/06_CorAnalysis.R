library(stats)
library(raster)
library(car)
library(reshape2)
library(ggplot2)



##### DATA #####
Meta <- read.csv("FinalData/MetaData.csv", header=T, sep=";")
Var <- c("HN","P","TMEAN","NAO","AO")
Var_c <- c("#3d77ff","#33D81B","#FF8F00","#ff2600")

Data <- vector("list", 5)
for (i in 1:length(Var)) {
  Data[[i]] <- read.csv(paste0("FinalData/Data_",Var[i],".csv"), header=T, sep=";")
}
names(Data) <- Var

Data_S <- vector("list", 5)
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



##### CORRELATION CALCULATION #####
Cor <- vector("list", 8)
names(Cor) <- ms
Cor_X <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=5))
Cor_X[,1] <- Meta$Name
colnames(Cor_X) <- c("Name",Var[-1])
Cor2 <- Cor
Cor2_X <- Cor_X
Pv <- Cor
Pv_X <- Cor_X
Pt <- Cor
Pt_X <- Cor_X

for (j in 1:length(ms)) {
  
  for (i in 1:length(Meta$Name)) {
    
    if (j == 8) {
      
      Data_X2 <- cbind(Data_S[["HN"]][c(1,i+1)], Data_S[["P"]][i+1], Data_S[["TMEAN"]][i+1], Data_S[["NAO"]][i+1], Data_S[["AO"]][i+1])
      colnames(Data_X2)[2:6] <- Var
      
    } else {
      
      Data_X <- cbind(Data[["HN"]][c(1:2,i+2)], Data[["P"]][i+2], Data[["TMEAN"]][i+2], Data[["NAO"]][i+2], Data[["AO"]][i+2])
      Data_X <- Data_X[-c(1:9),]
      Data_X <- head(Data_X,-1)
      colnames(Data_X)[3:7] <- Var
      rownames(Data_X) <- NULL
      
      Data_X2 <- Data_X[which(Data_X$Month == m[j]),]
      rownames(Data_X2) <- NULL
      
    }
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 0) >= 5) {
      
      for (k in 2:length(Var)) {
        
        if (j == 8) {
          Cor_X[i,k] <- cor.test(Data_X2[,k+1], Data_X2$HN, method="pearson")$estimate
          Pv_X[i,k] <- cor.test(Data_X2[,k+1], Data_X2$HN, method="pearson")$p.value
          Cor2_X[i,k] <- cor.test(Data_X2[,k+1], Data_X2$HN, method="pearson")$estimate^2
        } else {
          Cor_X[i,k] <- cor.test(Data_X2[,k+2], Data_X2$HN, method="pearson")$estimate
          Pv_X[i,k] <- cor.test(Data_X2[,k+2], Data_X2$HN, method="pearson")$p.value
          Cor2_X[i,k] <- cor.test(Data_X2[,k+2], Data_X2$HN, method="pearson")$estimate^2
        }
        
        Pt_X[i,1+which(Pv_X[i,-1] <= a)] <- 19
        Pt_X[i,1+which(Pv_X[i,-1] > a)] <- 1
        
      }
      
    } else {
      
      Cor_X[i,-1] <- NA
      Cor2_X[i,-1] <- NA
      Pv_X[i,-1] <- NA
      Pt_X[i,-1] <- NA
      
    }
    
  }
  
  Cor[[j]] <- Cor_X
  Cor2[[j]] <- Cor2_X
  Pv[[j]] <- Pv_X
  Pt[[j]] <- Pt_X
  
}

# Cor_Avg <- as.data.frame(matrix(nrow=length(Var[-1]), ncol=(1+length(ms))))
# Cor_Avg[,1] <- Var[-1]
# colnames(Cor_Avg) <- c("Var",ms)
# 
# for (j in 1:length(ms)) {
# 
#   Cor_Avg[,j+1] <- round(colMeans(Cor[[j]][-1], na.rm=T), 3)
# 
# }
# 
# write.table(Cor_Avg, "06_CorAnalysis/Cor_Avg_v2.csv", sep=";", row.names=F)



##### PLOTTING #####
### Adjusting the dataset ###
Cor_TOT <- as.data.frame(matrix(nrow=8*length(Meta$Name), ncol=7))
Cor_TOT[,1] <- factor(rep(ms, each=length(Meta$Name)), levels=ms)
Cor_TOT[,2] <- rep(Meta$Name, times=8)
Cor_TOT[,3] <- rep(Meta$Elevation, times=8)
colnames(Cor_TOT) <- c("Period","Name","Elevation",colnames(Cor_X)[2:5])
Pt_TOT <- Cor_TOT

dj <- c(length(Meta$Name)*c(1:length(ms)))

for (j in 1:length(ms)) {
  
  Cor_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:7)] <- Cor[[j]][-1] # 2: squared values
  Pt_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:7)] <- Pt[[j]][-1]
  
}

Cor_TOT <- na.omit(Cor_TOT)
rownames(Cor_TOT) <- NULL
Pt_TOT <- na.omit(Pt_TOT)
rownames(Pt_TOT) <- NULL

Cor_TOT2 <- melt(Cor_TOT[,-c(1:2)], id.vars="Elevation")
Cor_TOT2 <- cbind(rep(Cor_TOT$Period, times=4), rep(Cor_TOT$Name, times=4), Cor_TOT2)
colnames(Cor_TOT2) <- c("Period","Name","Elevation","Var","Value")
levels(Cor_TOT2$Var) <- Var[-1]
Cor_TOT2$Pt <- melt(Pt_TOT[,-c(1:2)], id.vars="Elevation")$value

### Scatter plot ###
png(paste0("06_CorAnalysis/Plots/CorAnalysis_1.png"), width=1500, height=1000, res=100)

ggplot(Cor_TOT2, aes(Elevation, Value, colour=Var)) +
  geom_hline(yintercept=0, size=1, linetype="dashed") +
  geom_point(size=1.5, shape=Cor_TOT2$Pt) +
  facet_grid(rows=vars(Period), cols=vars(Var)) +
  scale_colour_manual(values=c("#3d77ff","#33D81B","#FF8F00","#ff2600")) +
  scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
  ylim(c(-1,1)) +
  labs(x="Elevation [m]", y=expression(r)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
        strip.text.x = element_text(face="bold", size=15, vjust=1),
        strip.text.y = element_text(face="bold", size=15, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines"),
        panel.grid.minor = element_blank())

dev.off()

### Summary plot ###
NS <- as.data.frame(matrix(nrow=(length(Var)-1), ncol=length(ms)))
colnames(NS) <- ms
rownames(NS) <- Var[-1]

png(paste0("06_CorAnalysis/Plots/CorAnalysis_2.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-1,1), xlim=c(0.5,12.5), xlab=NA, ylab=expression(r^2), cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)

for (j in 1:length(ms)) {

  p <- j+j*0.5

  for (k in 2:length(Var)) {
    NS[k-1,j] <- sum(!is.na(Cor[[j]][[k]])) # 2: squared values
  }

  boxplot(Cor[[j]][[2]], Cor[[j]][[3]], Cor[[j]][[4]], Cor[[j]][[5]],
          add=T, boxwex=0.2, at=c(p-0.45,p-0.15,p+0.15,p+0.45),
          col=c("#3d77ff","#33D81B","#FF8F00","#ff2600"), xaxt="n", yaxt="n") # 2: squared values
  text(p, -1, labels=NS[1,j], cex=1, font=3)

}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-1,1,length.out=11), labels=seq(-1,1,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=c("P","TMEAN","NAO","AO"), fill=c("#3d77ff","#33D81B","#FF8F00","#ff2600"),
       x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -1, labels="N.S.:", cex=1, font=3)

dev.off()

### Most correlated spatial plot ###
Cor_sp <- vector("list", 8)
names(Cor_sp) <- ms

for (j in 1:length(ms)) {
  
  Cor_X <- na.omit(Cor[[j]])
  Pv_X <- na.omit(Pv[[j]])
  Cor2_X <- na.omit(Cor2[[j]])
  
  Cor_sp_X <- as.data.frame(matrix(nrow=length(Cor_X$Name), ncol=10))
  Cor_sp_X[,1:5] <- Meta[rownames(Cor_X),]
  colnames(Cor_sp_X) <- c(colnames(Meta),"r","r^2","Color","Pv","Pt")
  
  Cor_X <- Cor_X[,-1]
  Pv_X <- Pv_X[,-1]
  Cor2_X <- Cor2_X[,-1]
  
  for (i in 1:length(Cor_sp_X$Name)) {
    
    Cor_sp_X[i,6] <- Cor_X[i,max.col(Cor2_X)[i]]
    Cor_sp_X[i,7] <- Cor2_X[i,max.col(Cor2_X)[i]]
    Cor_sp_X[i,8] <- Var_c[max.col(Cor2_X)[i]]
    Cor_sp_X[i,9] <- Pv_X[i,max.col(Cor2_X)[i]]
    
    if (Cor_sp_X$Pv[i] > a && Cor_sp_X$r[i] > 0) {Cor_sp_X$Pt[i] <- 2}
    if (Cor_sp_X$Pv[i] <= a && Cor_sp_X$r[i] > 0) {Cor_sp_X$Pt[i] <- 24}
    if (Cor_sp_X$Pv[i] > a && Cor_sp_X$r[i] < 0) {Cor_sp_X$Pt[i] <- 6}
    if (Cor_sp_X$Pv[i] <= a && Cor_sp_X$r[i] < 0) {Cor_sp_X$Pt[i] <- 25}
    
  }
  
  Cor_sp[[j]] <- Cor_sp_X
  
}

png(paste0("06_CorAnalysis/Plots/CorAnalysis_3.png"), width=2400, height=1200, res=100)
par(mar=c(0,0,5,0), mgp=c(3,0,0), mfrow=c(2,4))

for (j in 1:(length(ms))) {
  
  plot(Hillshade, main=ms[j], cex.main=3, col=gray(0:255/255), legend=F, axes=F, box=T, alpha=0.6)
  plot(TAA, lwd=1, border="black", add=T)
  points(Cor_sp[[j]]$Longitude, Cor_sp[[j]]$Latitude, type="p", pch=Cor_sp[[j]]$Pt, col=Cor_sp[[j]]$Color,
         bg=Cor_sp[[j]]$Color, cex=3, lwd=3)
  legend("bottomright", legend=Var[-1], lty=c(1,1,1,1), lwd= c(3,3,3,3), col=Var_c, cex=2, bty="n")
  
}

dev.off()





library(stats)
library(car)
library(reshape2)
library(broom)
library(grid)
library(ggplot2)
library(sp)
library(raster)
# library(rsq)



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



##### EXPLAINED VARIANCE CALCULATION #####
EV <- vector("list", 8)
names(EV) <- ms
EV_X <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=6))
EV_X[,1] <- Meta$Name
colnames(EV_X) <- c("Name",Var[-1],"Residuals")
Coeff <- EV
Coeff_X <- EV_X[,-6]
Pv <- EV
Pv_X <- EV_X
Pt <- EV
Pt_X <- EV_X
# R2_X <- as.data.frame(matrix(nrow=15, ncol=2))
# R2_X[,1] <- c("HN ~ P+TMEAN+NAO+AO", "HN ~ P+TMEAN+NAO", "HN ~ P+TMEAN+AO", "HN ~ P+NAO+AO", "HN ~ TMEAN+NAO+AO",
#               "HN ~ P+TMEAN", "HN ~ P+NAO", "HN ~ P+AO", "HN ~ NAO+AO", "HN ~ TMEAN+NAO", "HN ~ TMEAN+AO",
#               "HN ~ P", "HN ~ TMEAN", "HN ~ NAO", "HN ~ AO")
# colnames(R2_X) <- c("Model","R2_adj")

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
    
    # Data_X2[Var[-1]] <- lapply(Data_X2[Var[-1]], scale)
    # Data_X2$HN <- Data_X2$HN + 1
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 0) >= 5) { # 1
      
      # Anova model #
      LM_X <- lm(HN ~ P+TMEAN+NAO+AO, data=Data_X2)
      Coeff_X[i,-1] <- LM_X$coefficients[-1]
      AN_X <- anova(LM_X)
      AN_X$EV <- AN_X$`Sum Sq`/sum(AN_X$`Sum Sq`)*100
      EV_X[i,-1] <- round(AN_X$EV, 3)
      Pv_X[i,-1] <- AN_X$`Pr(>F)`
      Pt_X[i,1+which(Pv_X[i,-1] <= a)] <- 19
      Pt_X[i,1+which(Pv_X[i,-1] > a)] <- 1
      
      # # Anova best model #
      # for (k in 1:length(R2_X$Model)) {
      #   
      #   R2_X[k,2] <- summary.lm(lm(R2_X[k,1], data=Data_X2))$adj.r.squared
      #   
      # }
      # 
      # m_best <- R2_X[which(R2_X[,2] == max(R2_X[,2])),1]
      # AN_X <- anova(lm(m_best, data=Data_X2))
      # AN_X$EV <- AN_X$`Sum Sq`/sum(AN_X$`Sum Sq`)*100
      # 
      # k2=1
      # for (k1 in 2:length(EV_X)) {
      #   
      #   if (colnames(EV_X)[k1] == rownames(AN_X)[k2]) {
      #     EV_X[i,k1] <- AN_X[k2,6]
      #     Pv_X[i,k1] <- AN_X[k2,5]
      #     k2 <- k2+1
      #   } else {
      #     EV_X[i,k1] <- NA
      #     Pv_X[i,k1] <- NA
      #   }
      #   
      # }
      # 
      # Pt_X[i,1+which(Pv_X[i,-1] <= a)] <- 19
      # Pt_X[i,1+which(Pv_X[i,-1] > a)] <- 1
      
      # # R2 decomposition #
      # # Model 1: HN ~ P+TMEAN+NAO+AO
      # M1_X <- glm(HN ~ P+TMEAN+NAO+AO, data=Data_X2, family=Gamma("log"))
      # R2_M1_X <- rsq(M1_X, adj=FALSE, type="n")
      # Pv_X[i,c(2:3)] <- summary.glm(M1_X)$coefficients[-1,4]
      # Pt_X[i,1+which(Pv_X[i,-1] <= a)] <- 19
      # Pt_X[i,1+which(Pv_X[i,-1] > a)] <- 1
      # 
      # # Model 2: HN ~ TMEAN+NAO+AO
      # M2_X <- glm(HN ~ TMEAN+NAO+AO, data=Data_X2, family=Gamma("log"))
      # R2_M2_X <- rsq(M2_X, adj=FALSE, type="n")
      # EV_X[i,2] <- round((R2_M1_X - R2_M2_X)/R2_M1_X, 4)*100
      # 
      # # Model 3: HN ~ P+AO+NAO
      # M3_X <- glm(HN ~ P+AO+NAO, data=Data_X2, family=Gamma("log"))
      # R2_M3_X <- rsq(M3_X, adj=FALSE, type="n")
      # EV_X[i,3] <- round((R2_M1_X - R2_M3_X)/R2_M1_X, 4)*100
      # 
      # # Model 4: HN ~ P+TMEAN+AO
      # M4_X <- glm(HN ~ P+TMEAN+AO, data=Data_X2, family=Gamma("log"))
      # R2_M4_X <- rsq(M4_X, adj=FALSE, type="n")
      # EV_X[i,4] <- round((R2_M1_X - R2_M4_X)/R2_M1_X, 4)*100
      # 
      # # Model 5: HN ~ P+TMEAN+NAO
      # M5_X <- glm(HN ~ P+TMEAN+NAO, data=Data_X2, family=Gamma("log"))
      # R2_M5_X <- rsq(M5_X, adj=FALSE, type="n")
      # EV_X[i,5] <- round((R2_M1_X - R2_M5_X)/R2_M1_X, 4)*100
      # 
      # # Residuals
      # EV_X[i,6] <- (100 - EV_X[i,2] - EV_X[i,3] - EV_X[i,4] - EV_X[i,5])
      # if (EV_X[i,6] > 100) {EV_X[i,-1] <- NA; Pv_X[i,-1] <- NA; Pt_X[i,-1] <- NA}
      
    } else {
      
      Coeff_X[i,-1] <- NA
      EV_X[i,-1] <- NA
      Pv_X[i,-1] <- NA
      Pt_X[i,-1] <- NA
      
    }
    
  }
  
  Coeff[[j]] <- Coeff_X
  EV[[j]] <- EV_X
  Pv[[j]] <- Pv_X
  Pt[[j]] <- Pt_X
  
}



##### PLOTTING #####
### Adjusting the dataset ###
EV_TOT <- as.data.frame(matrix(nrow=8*length(Meta$Name), ncol=8))
EV_TOT[,1] <- factor(rep(ms, each=length(Meta$Name)), levels=ms)
EV_TOT[,2] <- rep(Meta$Name, times=8)
EV_TOT[,3] <- rep(Meta$Elevation, times=8)
colnames(EV_TOT) <- c("Period","Name","Elevation",colnames(EV_X)[2:6])
Pt_TOT <- EV_TOT

dj <- c(length(Meta$Name)*c(1:length(ms)))

for (j in 1:length(ms)) {
  
  EV_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:8)] <- EV[[j]][-1]
  Pt_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:8)] <- Pt[[j]][-1]
  
}

EV_TOT <- na.omit(EV_TOT)
rownames(EV_TOT) <- NULL
Pt_TOT <- na.omit(Pt_TOT[,-8])
rownames(Pt_TOT) <- NULL
Pt_TOT$Residuals <- 19

EV_TOT2 <- melt(EV_TOT[,-c(1:2)], id.vars="Elevation")
EV_TOT2 <- cbind(rep(EV_TOT$Period, times=5), rep(EV_TOT$Name, times=5), EV_TOT2)
colnames(EV_TOT2) <- c("Period","Name","Elevation","Var","Value")
levels(EV_TOT2$Var) <- c(Var[-1],"Residuals")
EV_TOT2$Pt <- melt(Pt_TOT[,-c(1:2)], id.vars="Elevation")$value

### Scatter plot ###
png(paste0("09_EVAnalysis/Plots/EVAnalysis_1.png"), width=1500, height=1000, res=100)

ggplot(EV_TOT2, aes(Elevation, Value, colour=Var)) +
  # geom_hline(yintercept=0, size=0.8, linetype="dashed") +
  geom_point(size=1.5, shape=EV_TOT2$Pt) +
  facet_grid(rows=vars(Period), cols=vars(Var)) +
  scale_colour_manual(values=c("#3d77ff","#33D81B","#FF8F00","#ff2600","#D00DDC")) +
  scale_x_continuous(breaks=seq(0,3000,1000), labels=seq(0,3000,1000), limits=c(0,3000)) +
  ylim(c(0,100)) +
  labs(x="Elevation [m]", y="EV [%]") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
        strip.text.x = element_text(face="bold", size=15, vjust=1),
        strip.text.y = element_text(face="bold", size=15, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines"),
        panel.grid.minor = element_blank())

dev.off()

### Summary plot ###
NS <- as.data.frame(matrix(nrow=length(Var), ncol=length(ms)))
colnames(NS) <- ms
rownames(NS) <- c(Var[-1],"Res")

png(paste0("09_EVAnalysis/Plots/EVAnalysis_2.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-5,100), xlim=c(0.5,12.5), xlab=NA, ylab="Explained Variance [%]", cex.lab=1.5)

for (j in 1:length(ms)) {

  p <- j+j*0.5

  for (k in 1:length(Var)) {
    NS[k,j] <- sum(!is.na(EV[[j]][[k+1]]))
  }

  boxplot(EV[[j]][[2]], EV[[j]][[3]], EV[[j]][[4]], EV[[j]][[5]], EV[[j]][[6]],
          add=T, boxwex=0.10, at=c(p-0.30,p-0.15,p,p+0.15,p+0.30),
          col=c("#3d77ff","#33D81B","#FF8F00","#ff2600","#D00DDC"), xaxt="n", yaxt="n")
  text(p, -5, labels=NS[1,j], cex=1, font=3)

}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(0,100,10), labels=seq(0,100,10), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=c("P","TMEAN","NAO","AO","Res"), fill=c("#3d77ff","#33D81B","#FF8F00","#ff2600","#D00DDC"),
       x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -5, labels="N.S.:", cex=1, font=3)

dev.off()

### Seasonal scatter plot ###
EV_SEAS <- na.omit(cbind(Meta[,c(2,5)],EV[[8]][2:3]))
EV_SEAS <- EV_SEAS[order(EV_SEAS$Elevation),]
rownames(EV_SEAS) <- NULL   # 28,63
EV_TOT3 <- EV_TOT2[which(EV_TOT2$Period == "Seas"),]
EV_TOT3 <- EV_TOT3[c(which(EV_TOT3$Var == "P"), which(EV_TOT3$Var == "TMEAN")),]

png(paste0("09_EVAnalysis/Plots/EVAnalysis_seas.png"), width=1500, height=1000, res=100)

ggplot(EV_TOT3, aes(Elevation, Value, colour=Var)) +
  # geom_vline(xintercept=EV_SEAS$Elevation[28], size=1.2, linetype="dashed") +
  # geom_vline(xintercept=EV_SEAS$Elevation[63], size=1.2, linetype="dashed") +
  geom_point(size=3, shape=EV_TOT3$Pt) +
  scale_colour_manual(values=c("#3d77ff","#33D81B")) +
  scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
  ylim(c(0,60)) +
  labs(x="Elevation [m]", y="Explained Variance [%]") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, vjust=-0.7), axis.text.y = element_text(size=15, hjust=-0.3),
        axis.title.x = element_text(size=21, vjust=-3), axis.title.y = element_text(size=21, vjust=5),
        plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines")) +
  theme(legend.title = element_blank(), legend.text = element_text(size=14), legend.spacing.y = unit(0, "mm"),
        legend.background = element_blank(), legend.box.background = element_rect(colour = "black"),
        legend.key = element_blank(), legend.key.size = unit(7, "mm"), 
        legend.justification=c(0.75,0), legend.position=c(0.9795, 0.46))
  # annotation_custom(grid.text("TMEAN", x=0.2,  y=0.9, gp=gpar(col="#33D81B", fontsize=18, fontface="bold"))) +
  # annotation_custom(grid.text("P", x=0.45,  y=0.92, gp=gpar(col="#3d77ff", fontsize=18, fontface="bold"))) +
  # annotation_custom(grid.text("TMEAN", x=0.45,  y=0.88, gp=gpar(col="#33D81B", fontsize=18, fontface="bold"))) +
  # annotation_custom(grid.text("P", x=0.75,  y=0.9, gp=gpar(col="#3d77ff", fontsize=18, fontface="bold")))

dev.off()



##### PLOTTING EV + COEFF MAP #####
### Adjusting the dataset ###
EVC <- vector("list", 8)
names(EVC) <- ms

for (j in 1:length(ms)) {
  
  EV2_X <- na.omit(EV[[j]])
  Pv2_X <- na.omit(Pv[[j]][-6])
  Coeff2_X <- na.omit(Coeff[[j]])
  
  EVC_X <- as.data.frame(matrix(nrow=length(EV2_X$Name), ncol=9))
  EVC_X[,1:5] <- Meta[rownames(EV2_X),]
  colnames(EVC_X) <- c(colnames(Meta),"Color","Pv","Coeff","Pt")
  
  EV2_X <- EV2_X[,-c(1,6)]
  Pv2_X <- Pv2_X[,-1]
  Coeff2_X <- Coeff2_X[,-1]
  
  EVC_X[,6] <- Var_c[max.col(EV2_X)]
  
  for (i in 1:length(EVC_X$Name)) {
    
    EVC_X[i,7] <- Pv2_X[i,max.col(EV2_X)[i]]
    EVC_X[i,8] <- Coeff2_X[i,max.col(EV2_X)[i]]
    
    if (EVC_X$Pv[i] > a && EVC_X$Coeff[i] > 0) {EVC_X$Pt[i] <- 2}
    if (EVC_X$Pv[i] <= a && EVC_X$Coeff[i] > 0) {EVC_X$Pt[i] <- 24}
    if (EVC_X$Pv[i] > a && EVC_X$Coeff[i] < 0) {EVC_X$Pt[i] <- 6}
    if (EVC_X$Pv[i] <= a && EVC_X$Coeff[i] < 0) {EVC_X$Pt[i] <- 25}
    
  }
  
  EVC[[j]] <- EVC_X
  
}

### Plotting ###
png(paste0("09_EVAnalysis/Plots/EVAnalysis_3_v2.png"), width=2400, height=1200, res=100)
par(mar=c(0.5,0.5,5,0), mgp=c(3,0,0), mfrow=c(2,4))

for (j in 1:(length(ms))) {
  
  plot(Hillshade, main=ms[j], cex.main=3, col=gray(0:255/255), legend=F, axes=F, box=T, alpha=0.6)
  plot(TAA, lwd=1, border="black", add=T)
  points(EVC[[j]]$Longitude, EVC[[j]]$Latitude, type="p", pch=EVC[[j]]$Pt, col=EVC[[j]]$Color, bg=EVC[[j]]$Color, cex=3, lwd=3)
  legend("bottomright", legend=Var[-1], lty=c(1,1,1,1), lwd= c(3,3,3,3), col=Var_c, cex=2, bty="n")
  
}

dev.off()





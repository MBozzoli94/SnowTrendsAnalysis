library(stats)
library(car)
library(reshape2)
library(ggplot2)



##### DATA #####
Meta <- read.csv("FinalData/MetaData.csv", header=T, sep=";")
Var <- c("HN","P","TMEAN","NAO","AO")

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



##### REGRESSION COEFFICIENTS CALCULATION #####
Coeff <- vector("list", 8)
names(Coeff) <- ms
Coeff_X <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=5))
Coeff_X[,1] <- Meta$Name
colnames(Coeff_X) <- c("Name",Var[-1])
Coeff_min <- Coeff
Coeff_min_X <- Coeff_X
Coeff_max <- Coeff
Coeff_max_X <- Coeff_X
Pv <- Coeff
Pv_X <- Coeff_X
Pt <- Coeff
Pt_X <- Coeff_X
sd_X <- array()

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
    
    Data_X2[Var[-1]] <- lapply(Data_X2[Var[-1]], scale)
    Data_X2$HN <- Data_X2$HN + 1
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 1) >= 5) {
      
      GLM_X <- glm(HN ~ P+TMEAN+NAO+AO, data=Data_X2, family=Gamma("log"))
      
      for (k in 2:length(Var)) {
        
        Coeff_X[i,k] <- summary.lm(GLM_X)$coefficients[k,1]
        Pv_X[i,k] <- summary.lm(GLM_X)$coefficients[k,4]
        Pt_X[i,1+which(Pv_X[i,-1] <= a)] <- 19
        Pt_X[i,1+which(Pv_X[i,-1] > a)] <- 1
        
      }
      
    } else {
  
      Coeff_X[i,-1] <- NA
      Pv_X[i,-1] <- NA
      Pt_X[i,-1] <- NA
      
    }
    
  }
  
  for (k in 2:length(Var)) {
    
    sd_X[k-1] <- sd(Coeff_X[,k], na.rm=T)
    Coeff_min_X[,k] <- Coeff_X[,k]-1.96*sd_X[k-1]
    Coeff_max_X[,k] <- Coeff_X[,k]+1.96*sd_X[k-1]
    
  }
  
  Coeff[[j]] <- Coeff_X
  Coeff_min[[j]] <- Coeff_min_X
  Coeff_max[[j]] <- Coeff_max_X
  Pv[[j]] <- Pv_X
  Pt[[j]] <- Pt_X
  
}



##### PLOTTING #####
### Adjusting the dataset ###
Coeff_TOT <- as.data.frame(matrix(nrow=8*length(Meta$Name), ncol=7))
Coeff_TOT[,1] <- factor(rep(ms, each=length(Meta$Name)), levels=ms)
Coeff_TOT[,2] <- rep(Meta$Name, times=8)
Coeff_TOT[,3] <- rep(Meta$Elevation, times=8)
colnames(Coeff_TOT) <- c("Period","Name","Elevation",colnames(Coeff_X)[2:5])
Coeff_min_TOT <- Coeff_TOT
Coeff_max_TOT <- Coeff_TOT
Pt_TOT <- Coeff_TOT

dj <- c(length(Meta$Name)*c(1:length(ms)))

for (j in 1:length(ms)) {
  
  Coeff_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:7)] <- Coeff[[j]][-1]
  Coeff_min_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:7)] <- Coeff_min[[j]][-1]
  Coeff_max_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:7)] <- Coeff_max[[j]][-1]
  Pt_TOT[(dj[j]-(length(Meta$Name)-1)):dj[j],c(4:7)] <- Pt[[j]][-1]
  
}

Coeff_TOT <- na.omit(Coeff_TOT)
rownames(Coeff_TOT) <- NULL
Coeff_min_TOT <- na.omit(Coeff_min_TOT)
rownames(Coeff_min_TOT) <- NULL
Coeff_max_TOT <- na.omit(Coeff_max_TOT)
rownames(Coeff_max_TOT) <- NULL
Pt_TOT <- na.omit(Pt_TOT)
rownames(Pt_TOT) <- NULL

Coeff_TOT2 <- melt(Coeff_TOT[,-c(1:2)], id.vars="Elevation")
Coeff_TOT2 <- cbind(rep(Coeff_TOT$Period, times=4), rep(Coeff_TOT$Name, times=4), Coeff_TOT2)
colnames(Coeff_TOT2) <- c("Period","Name","Elevation","Var","Value")
levels(Coeff_TOT2$Var) <- Var[-1]
Coeff_TOT2$Value_min <- melt(Coeff_min_TOT[,-c(1:2)], id.vars="Elevation")$value
Coeff_TOT2$Value_max <- melt(Coeff_max_TOT[,-c(1:2)], id.vars="Elevation")$value
Coeff_TOT2$Pt <- melt(Pt_TOT[,-c(1:2)], id.vars="Elevation")$value

### Scatter plot ###
png(paste0("07_GLMAnalysis/Plots/GLMAnalysis_1.png"), width=1500, height=1000, res=100)

ggplot(Coeff_TOT2, aes(Elevation, Value, colour=Var)) +
  geom_hline(yintercept=0, size=0.8, linetype="dashed") +
  geom_point(size=1.5, shape=Coeff_TOT2$Pt) +
  geom_linerange(aes(ymin=Value_min, ymax=Value_max), alpha=0.2, size=1) +
  facet_grid(rows=vars(Period), cols=vars(Var)) +
  scale_colour_manual(values=c("#3d77ff","#33D81B","#FF8F00","#ff2600")) +
  scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
  ylim(c(-1.6,1.6)) +
  labs(x="Elevation [m]", y=expression(beta)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=18, vjust=-3), axis.title.y = element_text(size=18, vjust=5),
        strip.text.x = element_text(face="bold", size=12, vjust=1),
        strip.text.y = element_text(face="bold", size=12, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines"))

dev.off()

### Summary plot ###
NS <- as.data.frame(matrix(nrow=(length(Var)-1), ncol=length(ms)))
colnames(NS) <- ms
rownames(NS) <- Var[-1]

png(paste0("07_GLMAnalysis/Plots/GLMAnalysis_2.png"), width=1500, height=1000, res=100)
par(mar=c(3,5,0.5,0.5), mgp=c(3.5,1,0))

plot(NULL, xaxt="n", yaxt="n", ylim=c(-1.5,1.5), xlim=c(0.5,12.5), xlab=NA, ylab=expression(beta), cex.lab=1.5)
abline(h=0, col="black", lty=2, lwd=3)

for (j in 1:length(ms)) {

  p <- j+j*0.5

  for (k in 2:length(Var)) {
    NS[k-1,j] <- sum(!is.na(Coeff[[j]][[k]]))
  }

  boxplot(Coeff[[j]][[2]], Coeff[[j]][[3]], Coeff[[j]][[4]], Coeff[[j]][[5]],
          add=T, boxwex=0.2, at=c(p-0.45,p-0.15,p+0.15,p+0.45),
          col=c("#3d77ff","#33D81B","#FF8F00","#ff2600"), xaxt="n", yaxt="n")
  text(p, -1.5, labels=NS[1,j], cex=1, font=3)

}

axis(side=1, at=seq(1.5,12,1.5), labels=ms, las=1, cex.axis=1.5, tck=-0.01, xaxs="i")
axis(side=2, at=seq(-1.5,1.5,length.out=11), labels=seq(-1.5,1.5,length.out=11), las=2, cex.axis=1.5, tck=-0.01, xaxs="i")
legend("topleft", legend=c("P","TMEAN","NAO","AO"), fill=c("#3d77ff","#33D81B","#FF8F00","#ff2600"),
       x.intersp=0.5, y.intersp=1.2, cex=1.5, bty="n")
text(0.5, -1.5, labels="N.S.:", cex=1, font=3)

dev.off()





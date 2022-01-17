library(stats)
library(car)
library(reshape2)
library(ggplot2)
library(rsq)



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



##### R2 CALCULATION #####
R2 <- as.data.frame(matrix(nrow=length(Meta$Name), ncol=9))
R2[,1] <- Meta$Name
colnames(R2) <- c("Name",ms)

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
    
    if (sum(is.na(Data_X2$HN)) == 0 && sum(is.na(Data_X2$P)) == 0 && sum(Data_X2$HN != 1) >= 5) { # 1
      
      # LM_X <- lm(HN ~ P+TMEAN+NAO+AO, data=Data_X2)
      # R2[i,j+1] <- summary.lm(LM_X)$r.squared
      GLM_X <- glm(HN ~ P+TMEAN+NAO+AO, data=Data_X2, family=Gamma("log"))
      R2[i,j+1] <- rsq(GLM_X, adj=FALSE, type="n")
      
    } else {
      
      R2[i,j+1] <- NA
      
    }
    
  }
  
}



##### PLOTTING #####
### Adjusting the dataset ###
R2_TOT <- cbind(Meta$Name, Meta$Elevation, R2[,-1])
colnames(R2_TOT)[1:2] <- c("Name","Elevation")
R2_TOT <- melt(R2_TOT[,-1], id.vars="Elevation")
R2_TOT <- cbind(rep(Meta$Name, times=8), R2_TOT)
R2_TOT <- R2_TOT[,c(3,1,2,4)]
colnames(R2_TOT) <- c("Period","Name","Elevation","Value")
R2_TOT <- na.omit(R2_TOT)
rownames(R2_TOT) <- NULL

### Plotting ###
png(paste0("08_R2Analysis/Plots/R2Analysis_v2.png"), width=2000, height=1000, res=100) # v2

ggplot(R2_TOT, aes(Elevation, Value)) +
  geom_point(size=2, colour="#D00DDC") +
  facet_wrap(~ Period, nrow=2) +
  scale_x_continuous(breaks=seq(0,3000,500), labels=seq(0,3000,500), limits=c(0,3000)) +
  ylim(c(0,1)) +
  labs(x="Elevation [m]", y=expression(R^2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=18, vjust=-3), axis.title.y = element_text(size=18, vjust=5),
        strip.text.x = element_text(face="bold", size=12, vjust=1),
        strip.text.y = element_text(face="bold", size=12, vjust=1),
        legend.position = "None", plot.margin = unit(c(0.5,0.5,1.5,1.5), "lines"))

dev.off()





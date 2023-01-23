library(reshape2)
library(stringr)
library(mblm)
library(ggplot2)



### Reading data ###
meta <- read.csv("MetaData_sel.csv", header=T, sep=";")
meta$Long_name <- paste0(meta$Name,"\n[",meta$Elevation," m a.s.l.]")
meta <- meta[,c(1:2,6,3:5)]
data <- read.csv("Data_HN_30_S_sel.csv", header=T, sep=";")
colnames(data)[2:9] <- meta$Long_name



### Trend calculation ###
trend <- data.frame(meta, T_i=NA, T_s=NA)
s <- 1:40

for (i in 1:length(meta$Long_name)) {
  
  data.x <- data[,which(colnames(data) == meta$Long_name[i])]
  
  if (sum(is.na(data.x)) == 0) {
    trend$T_i[i] <- round(mblm(data.x ~ s, repeated=F)[[1]][[1]], 2)
    trend$T_s[i] <- round(mblm(data.x ~ s, repeated=F)[[1]][[2]], 3)
  }
  
}



### Plotting data ###
data.plot <- melt(data, id.vars="Date")
data.plot$Date_s <- 1:40
data.plot$T_i <- rep(trend$T_i, each=40)
data.plot$T_s <- rep(trend$T_s, each=40)
data.plot <- data.plot[,c(2,4,3,5:6)]
colnames(data.plot)[c(1,3)] <- c("Name","Tot_HN")

png("TS_sel.png", width=2000, height=1000, res=200)
ggplot(data.plot, aes(x=Date_s, y=Tot_HN)) +
  geom_point(size=0.8, shape=19) +
  facet_wrap(~Name, nrow=2, scales="free_y") +
  geom_abline(aes(intercept=T_i, slope=T_s), col="red") +
  scale_x_continuous(breaks=c(1,10,20,30,41), labels=c("1980","1990","2000","2010","2020")) +
  labs(x=element_blank(),
       y="Seasonal fresh snow [cm]") +
  theme_bw() +
  theme(axis.text.x=element_text(size=7, family="serif"),
        axis.text.y=element_text(size=7, family="serif"),
        axis.title.y=element_text(size=10, vjust=3.5, family="serif"),
        strip.text.x=element_text(size=7, family="serif"),
        plot.margin=unit(c(1,0.5,0.5,1.5),"lines"))
dev.off()





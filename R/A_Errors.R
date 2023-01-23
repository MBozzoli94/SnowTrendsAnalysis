library(reshape2)
library(ggplot2)





### Data ###
data <- read.csv("Stats_HN.csv", header=T, sep=";")
colnames(data)[3:4] <- c("BIAS [%]","MAE [%]")
data.plot <- data[,-1]
data.plot <- melt(data.plot, id.vars="Elevation")
data.plot$YI <- rep(c(4,45,0.5), each=122)
colnames(data.plot)[2:3] <- c("Var","Val")



### Plot ###
png("Errors.png", width=2400, height=800, res=200)
ggplot(data.plot, aes(x=Elevation, y=Val)) +
  geom_point(size=1, shape=19) +
  facet_wrap(~Var, scales="free_y") +
  geom_hline(aes(yintercept=YI), col="red", lty=2, lwd=0.8) +
  labs(x=element_blank(), y=element_blank()) +
  theme_bw() +
  theme(axis.text.x=element_text(size=10, family="serif"),
        axis.text.y=element_text(size=10, family="serif"),
        strip.text.x=element_text(size=10, family="serif"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"lines"))
dev.off()





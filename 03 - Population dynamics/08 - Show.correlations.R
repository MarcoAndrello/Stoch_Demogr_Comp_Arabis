# SHOW CORRELATIONS AMONG SLTRE CONTRIBUTIONS
# We show 6 correlations
colnames(data)
par.to.show <- data.frame(i.par = c(1,15,15,15,22,22),
                          j.par = c(6,4,6,7,23,24))
par.to.show$i.name <- c("C[S]^mu","C[S]^CV","C[S]^CV","C[S]^CV","C[S]^rho","C[S]^rho")
par.to.show$j.name <- c("C[F[1]]^mu","C[G^'+']^mu","C[F[1]]^mu","C[F[2]]^mu","C[G^'-']^rho","C[G^'=']^rho")

data1 <- data
data1[,c(22,23,24),] <- data1[,c(22,23,24),]*10 # Increase the values to avoid too many decimal digits on the plot

# Set themes for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.x = element_text(size=7, color="black"),
             axis.title.y = element_text(size=7, color="black"),
             plot.title = element_text(size = 7, color = "black", hjust=0.5),
             legend.position = "none")

plots <- list()
for (i.plot in 1 : dim(par.to.show)[1]) {
  i.par <- par.to.show$i.par[i.plot]
  j.par <- par.to.show$j.par[i.plot]
  expr.text.1 <- par.to.show$i.name[i.plot]
  expr.text.2 <- par.to.show$j.name[i.plot]
  data.stat <- data.frame( med = apply(data1[,c(i.par,j.par),],c(1,2),mean),
                           low = apply(data1[,c(i.par,j.par),],c(1,2),quantile,0.025),
                           high = apply(data1[,c(i.par,j.par),],c(1,2),quantile,0.975),
                           Site = dimnames(data1)[[1]]) 
  
  colnames(data.stat) <- c("med.x","med.y","low.x","low.y","high.x","high.y","Site")
  data.stat$Site <- as.character(data.stat$Site)
  data.stat$Site <- factor(data.stat$Site,levels=v.site)
  plots[[i.plot]] <-
    ggplot(data.stat,aes(x=med.x, xmin=low.x, xmax=high.x,y=med.y, ymin=low.y, ymax=high.y, col=Site)) +
    geom_errorbar(width=0,size=0.5) +
    geom_errorbarh(height=0,size=0.5) +
    geom_point(size=1.5) + 
    scale_colour_brewer(type="div",palette=7) +
    labs(x=parse(text=expr.text.1),y=parse(text=expr.text.2)) #+
  # coord_cartesian(xlim=c(-0.25,0.35), ylim=c(-0.25,0.35))
}

png(filename="Figure_4_NegCorr.png",width=16,height=9,units="cm",res=300)
marrangeGrob(plots,nrow=2,ncol=3,top="")
dev.off()

# Plot legends
png("Figure_4_Legend.png",width=2,height=(9+4.5*4/3),units="cm",res=600)
par(mar=c(0,0,0,0))
plot.new()
legend(-0.2,1,c("Neg","Pos"),fill=c("#ef8a62","#67a9cf"),bty="n",ncol=1)
pal <- brewer.pal(7,"RdYlBu")
legend(-0.2,0.6,v.site,fill=pal,bty="n",ncol=1)
dev.off()

# Calculate and plot number of individuals per site per year
# Marco Andrello
# Last modification: 14/01/2020

library(tidyverse)

load("Demographic_data.RData")

data.n <- as.data.frame(table(data$Site,data$Year))
names(data.n) <- c("Site","Year","N")
data <- data.n; rm(data.n)
data$Site <- as.character(data$Site)
data$Site <- factor(data$Site, levels=c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC"))
data$Year <- as.numeric(as.character(data$Year))


theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.x = element_text(size=7, color="black"),
             axis.title.y = element_text(size=7, color="black"),
             plot.margin = margin(0.5,0.5,0,0, unit="cm"),
             legend.position = "right",
             legend.direction = "vertical",
             legend.text = element_text(size=7, color="black"),
             legend.title = element_text(size=0, color="black"))

png("Figure S2-DEMOGRA.png",width=8.5,height=5,units="cm",res=600)
ggplot(data,aes(x=Year,y=N,col=Site)) +
  geom_point() +
  geom_line() +
  scale_colour_brewer(type="div",palette=7)
dev.off()

# Report Initial N in Table 1
data[data$Year==2008,]

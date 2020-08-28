############################################################################
############################################################## SMALL MAMMALS
############################################################################

####################################################
### Read in intra- and inter-annual variation data
####################################################
mam.inter <- readRDS(file="BAT-change_mammals_abs-abund-inter.rds")
mam.intra <- readRDS(file="BAT-change_mammals_abs-abund-intra.rds")

####################################################
### Initial plots
####################################################
### Get basic geographic info on sites
sites4yrs <- readRDS(file="NEONsites_tempreplic_4yrs.rds") #sites that are sampled at least 4 yrs and at least 3 months within each year
colnames(sites4yrs) <- c(colnames(sites4yrs[,1:4]),"siteID",colnames(sites4yrs[,6:ncol(sites4yrs)]))

### Get mean estimates for intra/inter per site (we may, but don't have to, use them later)
sites <- unique(mam.intra$siteID)
mam.intra.m <- matrix(NA, length(sites), 4)
mam.inter.m <- matrix(NA, length(sites), 4)
mam.intra.m <- as.data.frame(mam.intra.m)
mam.inter.m <- as.data.frame(mam.inter.m)

for (i in 1:length(sites)){
  si <- mam.intra %>%
    filter(siteID == sites[i])
  mam.intra.m[i,1] <- sites[i]
  mam.intra.m[i,2] <- mean(si$intra_jacc_diag, na.rm=TRUE)
  mam.intra.m[i,3] <- mean(si$intra_jacctur_diag, na.rm=TRUE)
  mam.intra.m[i,4] <- mean(si$intra_jaccnes_diag, na.rm=TRUE)
}
colnames(mam.intra.m) <- c("siteID",colnames(mam.intra[,5:ncol(mam.intra)]))
  
for (i in 1:length(sites)){
  si <- mam.inter %>%
    filter(siteID == sites[i])
  mam.inter.m[i,1] <- sites[i]
  mam.inter.m[i,2] <- mean(si$inter_jacc_diag, na.rm=TRUE)
  mam.inter.m[i,3] <- mean(si$inter_jacctur_diag, na.rm=TRUE)
  mam.inter.m[i,4] <- mean(si$inter_jaccnes_diag, na.rm=TRUE)
}
colnames(mam.inter.m) <- c("siteID",colnames(mam.inter[,4:ncol(mam.inter)]))

### Combine with geo info
mam.intra.inf <- mam.intra %>%
  left_join(sites4yrs, by = c("siteID"))
mam.inter.inf <- mam.inter %>%
  left_join(sites4yrs, by = c("siteID"))

mam.intra.m.inf <- mam.intra.m %>%
  left_join(sites4yrs, by = c("siteID"))
mam.inter.m.inf <- mam.inter.m %>%
  left_join(sites4yrs, by = c("siteID"))


######################################## Plot relationships w/ latitude
require(ggplot2)
### Intra
# intra jaccard ~ lat
p1 <- ggplot(mam.intra.inf, aes(Latitude,intra_jacc_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=mam.intra.m.inf, aes(Latitude,intra_jacc_diag), col="red", size=3)+
  geom_smooth(data=mam.intra.m.inf, aes(Latitude,intra_jacc_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p1

p1.m <- ggplot(data=mam.intra.m.inf, aes(Latitude,intra_jacc_diag), )+
  geom_point(col="red", size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p1.m

# intra jaccard tur ~ lat
p2 <- ggplot(mam.intra.inf, aes(Latitude,intra_jacctur_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=mam.intra.m.inf, aes(Latitude,intra_jacctur_diag), col="red", size=3)+
  geom_smooth(data=mam.intra.m.inf, aes(Latitude,intra_jacctur_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p2

p2.m <- ggplot(data=mam.intra.m.inf, aes(Latitude,intra_jacctur_diag))+
  geom_point(col="red", size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p2.m

# intra jaccard nes ~ lat
p3 <- ggplot(mam.intra.inf, aes(Latitude,intra_jaccnes_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=mam.intra.m.inf, aes(Latitude,intra_jaccnes_diag), col="red", size=3)+
  geom_smooth(data=mam.intra.m.inf, aes(Latitude,intra_jaccnes_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3

p3.m <- ggplot(data=mam.intra.m.inf, aes(Latitude,intra_jaccnes_diag))+
  geom_point(col="red",size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3.m

require(ggpubr)
n <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, align = "v")
ggsave(file=paste0("mam-intra-lat.jpeg"), n, width=30,height=10, dpi=300)
n <- ggarrange(p1.m, p2.m, p3.m, ncol = 3, nrow = 1, align = "v")
ggsave(file=paste0("mam-intra-lat_mean.jpeg"), n, width=30,height=10, dpi=300)



### Inter-annual variability
# inter jaccard ~ lat
p1 <- ggplot(mam.inter.inf, aes(Latitude,inter_jacc_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=mam.inter.m.inf, aes(Latitude,inter_jacc_diag), col="red",size=3)+
  geom_smooth(data=mam.inter.m.inf, aes(Latitude,inter_jacc_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p1

p1.m <- ggplot(data=mam.inter.m.inf, aes(Latitude,inter_jacc_diag))+
  geom_point(col="red",size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p1.m

# inter jaccard tur ~ lat
p2 <- ggplot(mam.inter.inf, aes(Latitude,inter_jacctur_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=mam.inter.m.inf, aes(Latitude,inter_jacctur_diag), col="red",size=3)+
  geom_smooth(data=mam.inter.m.inf, aes(Latitude,inter_jacctur_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p2

p2.m <- ggplot(data=mam.inter.m.inf, aes(Latitude,inter_jacctur_diag))+
  geom_point(col="red",size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p2.m

# inter jaccard nes ~ lat
p3 <- ggplot(mam.inter.inf, aes(Latitude,inter_jaccnes_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=mam.inter.m.inf, aes(Latitude,inter_jaccnes_diag), col="red",size=3)+
  geom_smooth(data=mam.inter.m.inf, aes(Latitude,inter_jaccnes_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3

p3.m <- ggplot(data=mam.inter.m.inf, aes(Latitude,inter_jaccnes_diag))+
  geom_point(col="red",size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3.m

require(ggpubr)
n <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, align = "v")
ggsave(file=paste0("mam-inter-lat.jpeg"), n, width=30,height=10, dpi=300)
n <- ggarrange(p1.m, p2.m, p3.m, ncol = 3, nrow = 1, align = "v")
ggsave(file=paste0("mam-inter-lat_mean.jpeg"), n, width=30,height=10, dpi=300)



######################################## Plot relationships w/ one another
### Combine 
out_all <- mam.intra %>%
  left_join(mam.inter.m, by = c("siteID"))

out_all.m <- mam.intra.m %>%
  left_join(mam.inter.m, by = c("siteID"))


### intra + inter
require(ggplot2)
# jacc
p4 <- ggplot(out_all, aes(inter_jacc_diag, intra_jacc_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=out_all.m, aes(inter_jacc_diag, intra_jacc_diag), col="red", size=3)+
  geom_smooth(data=out_all.m, aes(inter_jacc_diag, intra_jacc_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p4

p4.m <- ggplot(data=out_all.m, aes(inter_jacc_diag, intra_jacc_diag))+
  geom_point(col="red", size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p4.m

# jacc tur
p5 <- ggplot(out_all, aes(inter_jacctur_diag, intra_jacctur_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=out_all.m, aes(inter_jacctur_diag, intra_jacctur_diag), col="red", size=3)+
  geom_smooth(data=out_all.m, aes(inter_jacctur_diag, intra_jacctur_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p5

p5.m <- ggplot(data=out_all.m, aes(inter_jacctur_diag, intra_jacctur_diag))+
  geom_point(col="red", size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p5.m

# jacc nes
p6 <- ggplot(out_all, aes(inter_jaccnes_diag, intra_jaccnes_diag))+
  geom_point(size=3)+
  geom_smooth(method="glm")+
  geom_point(data=out_all.m, aes(inter_jaccnes_diag, intra_jaccnes_diag), col="red", size=3)+
  geom_smooth(data=out_all.m, aes(inter_jaccnes_diag, intra_jaccnes_diag), col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p6

p6.m <- ggplot(data=out_all.m, aes(inter_jaccnes_diag, intra_jaccnes_diag))+
  geom_point(col="red", size=3)+
  geom_smooth(col="red",method="glm")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=2.5),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=20),
        axis.text.y = element_text(colour='black',size=20),
        axis.title.x = element_text(colour='black',size=25),
        axis.title.y = element_text(colour='black',size=25),
        axis.ticks = element_line(color = 'black', size=2.5),
        axis.ticks.length=unit(0.5,"cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p6.m

n <- ggarrange(p4, p5, p6, ncol = 3, nrow = 1, align = "v")
ggsave(file=paste0("intra-inter.jpeg"), n, width=30,height=10, dpi=300)
n <- ggarrange(p4.m, p5.m, p6.m, ncol = 3, nrow = 1, align = "v")
ggsave(file=paste0("intra-inter_mean.jpeg"), n, width=30,height=10, dpi=300)


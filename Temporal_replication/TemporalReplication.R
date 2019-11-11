
require(devtools)
require(raster)
require(neonUtilities)
require(rhdf5)
require(data.table)
require(tidyverse)
options(stringsAsFactors=F)

#### read in data to explore
m.plot <- read.delim("/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/mam_perplotnight.csv", sep=",")
m.trap <- read.delim(file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/mam_pertrapnight.csv", sep=",")
m.val <- read.delim("/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/validation.csv", sep=",")
m.var <- read.delim(file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/variables.csv", sep=",")

####exploration of mammal ddata
require(dplyr)
t1 <- t(as.data.frame(str_split(m.trap$endDate, "-")))
m.trap$year <- t1[,1]
m.trap$month <- t1[,2]
m.trap$day <- t1[,3]

saveRDS(m.trap, file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/mam_pertrapnight_DateDecomposed-MJ.rds")
m.trap <- readRDS(file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/mam_pertrapnight_DateDecomposed-MJ.rds")
m.t <- as_tibble(m.trap)

#how many temporal replicates within each site within each year
m1 <-  
m.t %>%
  group_by(siteID, year, month) %>%
  select(siteID, year, month) %>%
  summarise(
    count=n()
  ) 

m1 <- arrange(m1, siteID, year, month)

m2 <- 
m1 %>%
  group_by(siteID, year) %>%
  select(siteID, year) %>%
  summarise(
    count=n()
  )

m3 <- 
  m2 %>%
  group_by(siteID) %>%
  select(siteID) %>%
  summarise(
    count=n()
  )



##with end date
m4 <-  
  m.t %>%
  group_by(siteID, year, endDate) %>%
  select(siteID, year, endDate) %>%
  summarise(
    count=n()
  ) 

m4 <- arrange(m4, siteID, year, endDate)

m5 <- 
  m4 %>%
  group_by(siteID, year) %>%
  select(siteID, year) %>%
  summarise(
    count=n()
  )




##plots
sites <- unique(m2[,1])
saveRDS(sites, file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/NEONsites.rds")
require(ggpubr)

##Create some plots
for (i in 1:46){
  m2i <- filter(m2, siteID == sites[i,1])
  p1 <- ggplot(data=m2i, aes(x=year, y=count)) +
    geom_bar(stat="identity")+
    ggtitle(sites[i,1])+
    ylab("No. of Months")+
    xlab("Year")+
    geom_text(aes(label=count), vjust=1.6, color="white", position = position_dodge(0.9), size=10)+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=2.5),
          plot.title = element_text(size=30, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=30,margin=margin(t=10, r=0, b=20, l=0)),
          axis.text.y = element_text(colour='black',size=30,margin=margin(t=0, r=10, b=0, l=0)),
          axis.title.x = element_text(colour='black',size=35),
          axis.title.y = element_text(colour='black',size=35, margin=margin(t=0, r=20, b=0, l=0)),
          axis.ticks = element_line(color = 'black', size=2.5),
          axis.ticks.length=unit(0.5,"cm"),
          legend.position="none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  m5i <- filter(m5, siteID == sites[i,1])
  p2 <- ggplot(data=m5i, aes(x=year, y=count)) +
    geom_bar(stat="identity")+
    ggtitle(sites[i,1])+
    ylab("No. of Visits")+
    xlab("Year")+
    geom_text(aes(label=count), vjust=1.6, color="white", position = position_dodge(0.9), size=10)+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=2.5),
          plot.title = element_text(size=30, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=30,margin=margin(t=10, r=0, b=20, l=0)),
          axis.text.y = element_text(colour='black',size=30,margin=margin(t=0, r=10, b=0, l=0)),
          axis.title.x = element_text(colour='black',size=35),
          axis.title.y = element_text(colour='black',size=35, margin=margin(t=0, r=20, b=0, l=0)),
          axis.ticks = element_line(color = 'black', size=2.5),
          axis.ticks.length=unit(0.5,"cm"),
          legend.position="none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  n <- ggarrange(p1, p2, ncol = 2, nrow = 1, align = "v")
  ggsave(file=paste0("graphics/replic_",sites[i,1],".jpeg"), n, width=20,height=10, dpi=500)
}

#plot the latitudinal position of the sites along with their replication


##Create some more plots
datadir <- "/Users/jarzyna.1/Documents/RESEARCH_NEON"
setwd(datadir)
require(grid)
require(maptools)
require(rgdal)
require(ggplot2)
globe <- readOGR(".", layer = "continent")
us <- readOGR(".", layer = "ContUSA")
neon <- readOGR(".", layer = "NEON_Field_Sites")
alaska <- readOGR(".", layer="tl_2016_02_cousub")

proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
us.proj <- spTransform(us, proj)

# maps
point.all <- neon@data
require (ggrepel)

##getting sites that have at least 4 years of data with at least 3 replicate-months within each year.
m.4y <- 
  m2 %>%
  filter(count >= 3)

m.4y.2 <- 
  m.4y %>%
  group_by(siteID) %>%
  select(siteID) %>%
  summarise(
    count=n()
  )

m.4y.2sel <-
  m.4y.2 %>%
  filter(count >= 4) #these are the 27 sites

##getting additiona sites that have 3 years of data with at least 3 replicate-months within each year.
m.3y.2sel <-
  m.4y.2 %>%
  filter(count == 3) #these are the 3 sites

sites4y <- point.all[point.all$SiteID %in% m.4y.2sel$siteID,]
sites3y <- point.all[point.all$SiteID %in% m.3y.2sel$siteID,]
saveRDS(sites4y, file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/NEONsites_4yrs.rds")
saveRDS(sites3y, file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/NEONsites_3yrs.rds")

point.all.2 <- point.all[point.all$State != "AK" & point.all$State != "HI" & point.all$State != "PR", ]

p.map2 <- ggplot(point.all.2, aes(Longitude,Latitude))+
  geom_polygon(data=us.proj, aes(x=long, y=lat, group=group), fill="grey95", colour="black", alpha=1,size=0.5)+ 
  geom_point(data=point.all.2,aes(Longitude,Latitude), colour="black", fill="lightblue", pch=21, cex=5)+
  geom_point(data=sites4y,aes(Longitude,Latitude), colour="black", fill="yellow", pch=21, cex=5)+
  geom_point(data=sites3y,aes(Longitude,Latitude), colour="black", fill="pink", pch=21, cex=5)+
  geom_text_repel(data=sites4y, aes(Longitude, Latitude, label=SiteID))+
  geom_text_repel(data=sites3y, aes(Longitude, Latitude, label=SiteID), col="red")+
  # geom_text(data=point.all, aes(Longitude, Latitude, label=SiteID), col="black", fontface = "bold",position=position_jitter(width=3,height=3))+
  # geom_text_repel(data=point.all, aes(Longitude, Latitude, label=SiteID))+
  #scale_colour_manual(values = c("darkgoldenrod1","black","palegoldenrod","midnightblue","lightskyblue2","steelblue3","deepskyblue4"))+
  #scale_colour_viridis(discrete=TRUE,option="magma")+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        #plot.margin =unit(0.5, "cm"),
        plot.title = element_text(size=20, vjust=2, family="sans"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position=c(0.5,0.15),
        legend.direction="none",
        legend.title=element_text(size=16),
        legend.key.size = unit(0.7, "cm"),
        legend.title.align=0.5,
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom",title.position="bottom",ticks = TRUE))
p.map2

ggsave(file="graphics/Sites_map_2.jpeg", p.map2, width=18,height=10, dpi=500)



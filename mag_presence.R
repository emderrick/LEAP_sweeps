library(tidyverse)
library(reshape2)
library(ggplot2)

files<-list.files("Library/Mobile Documents/com~apple~CloudDocs/PhD/Chapter 1/Aim 1A/instrain_output/",recursive = T, pattern=".*genome_info.tsv",full.names = T)

all_mags<-data.frame()

for(i in 1:length(files)){
pond_time_mags<-read.table(files[i],sep="\t",header=T)
timepoint<-gsub(".*profile_output/", "", files[i]) %>% substr(1,4)
pond_time_mags<-cbind(pond_time_mags,timepoint=rep(timepoint,nrow(pond_time_mags)))
all_mags<-rbind(all_mags,pond_time_mags)
}

#reorder columns
all_mags<-all_mags[,c(1,ncol(all_mags),2:(ncol(all_mags)-1))]
all_mags$pond<-all_mags$timepoint%>%substr(1,2)

mag_names<-all_mags$genome%>%unique()
pond_names<-all_mags$timepoint%>%unique()

mag_table<-spread(all_mags[,1:3], key=timepoint, value=coverage)
row.names(mag_table)<-mag_table$genome
mag_table<-mag_table[,-1]
mag_table<-1*(!is.na(mag_table))

times<-colnames(mag_table)%>%substr(4,4)%>%as.numeric()
mag_table<-sweep(mag_table, MARGIN=2, times, `*`)
n_ponds<-apply(mag_table,1,function(x){sum(x!=0)})
mag_order<-order(n_ponds)
mag_table<-mag_table[mag_order,]


plotmag<-melt(mag_table)
plotmag<-plotmag[plotmag$value!=0,]
plotmag$pond<-plotmag$Var2%>%substr(1,2)

plotcolbypond<-ggplot(plotmag, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(pond))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=1, angle=0, vjust=0.3),
                    axis.text.y=element_text(size=1),
                    plot.title=element_text(size=11))

ggsave("plotcolbypond.png", limitsize = FALSE, width=4, height=7, dpi=5000)

plot20<-melt(mag_table[293:313,])
plot20<-plot20[plot20$value!=0,]

plot20mags<-ggplot(plot20, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(value))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plot20mags.png", limitsize = FALSE, width=10, height=5, dpi=500)  

plot50<-melt(mag_table[1:50,])
plot50<-plot50[plot50$value!=0,]

plot50mags<-ggplot(plot50, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(value))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plot50mags.png", limitsize = FALSE, width=5, height=10, dpi=500)                  

plot100<-melt(mag_table[51:100,])
plot100<-plot100[plot100$value!=0,]

plot100mags<-ggplot(plot100, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(value))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plot100mags.png", limitsize = FALSE, width=6, height=10, dpi=500)                  

plot150<-melt(mag_table[101:150,])
plot150<-plot150[plot150$value!=0,]

plot150mags<-ggplot(plot150, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(value))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plot150mags.png", limitsize = FALSE, width=7, height=10, dpi=500)                  

plot200<-melt(mag_table[151:200,])
plot200<-plot200[plot200$value!=0,]

plot200mags<-ggplot(plot200, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(value))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plot200mags.png", limitsize = FALSE, width=8, height=10, dpi=500)                  

plot300<-melt(mag_table[251:300,])
plot300<-plot300[plot300$value!=0,]

plot300mags<-ggplot(plot300, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=as.factor(value))) + 
  labs(x="Pond", y="MAG") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plot300mags.png", limitsize = FALSE, width=8, height=10, dpi=500)                  

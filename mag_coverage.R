library(ggplot2)

ggplot(all_mags, aes(x = genome, y = coverage, colour = pond)) +
  geom_point()+
  labs(x="MAG", y="Coverage") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(axis.text.x=element_text(size=5, angle=0, vjust=0.3),
        axis.text.y=element_text(size=7),
        plot.title=element_text(size=11))
ggsave("plotcoverage.png", limitsize = FALSE, width=20, height=4, dpi=500)  

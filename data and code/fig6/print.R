library("ggplot2")
library("pheatmap")
library("ggbreak") 
data = read.table("/print_abundance.csv", header=TRUE, sep=",")
colors = rev(c('#eeeeee', '#6fa8dc', '#e06666'))
data1 = reshape2::melt(data,id.vars = "group",variable.name = "habitat", value.name = "abundance")
data1$group = factor(data1$group)
data1$habitat = factor(data1$habitat)
p = ggplot(data=data1, aes(x= habitat, y=abundance, fill=group,width=0.4)) + geom_bar(stat="identity",position="stack",color="black") + theme_classic()
p1 = p + theme(axis.text.x = element_text(angle=30, vjust = 0.5))
p2 = p1 + scale_fill_manual(values=colors) + scale_y_continuous(expand = c(0,0))
p3 = p2 + labs(x="Habitat",y="Relative Abundance",
       fill="cluster type",title="具有IAA降解基因簇的微生物在不同生境的相对丰度")
p4 = p3 + theme(
    text=element_text(size=17),
    plot.title = element_text(hjust = 0.5,vjust = 1), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,
                                       vjust = 0.5),
    legend.title=element_text(size=12), 
legend.text=element_text(size=12))
p5 = p4 + scale_y_break(c(5,85),#截断位置及范围
                space = 0.5,#间距大小
                scales = 2.5)#上下显示比例，大于1上面比例大，小于1下面比例大 
ggsave(filename = "genus1.pdf",
       p5,
       width=10,
       heigh=8)
#######################################################################
library(tidyverse)
library(scatterpie)
library(reshape2)
data = read.table("print_point.csv", header=TRUE, sep=",")
library(RColorBrewer)
library(scales)
color_ce= c("#e06666","#6fa8dc")
names(color_ce)=c("iac","iad")
data$Order=ifelse(data$Order=="Puniceispirillales",-1,
                    ifelse(data$Order =="Propionibacteriales",-2,
                           ifelse(data$Order =="Sphingomonadales",-3,
                                  ifelse(data$Order =="Burkholderiales",-4,
									     ifelse(data$Order =="Moraxellales",-5,
                                                ifelse(data$Order =="Pseudomonadales",-6,								  
                                                       ifelse(data$Order =="Enterobacterales",-7,
                                                              ifelse(data$Order =="Caulobacterales",-8,
                                                                     ifelse(data$Order =="Rhizobiales",-9,
																	        ifelse(data$Order =="Actinomycetales",-10,
																			       ifelse(data$Order =="Rhodobacterales",-11,-12)))))))))))

data$environment=ifelse(data$environment=="animal gut",3,
                      ifelse(data$environment=="human gut",6,
                             ifelse(data$environment=="human oral",9,
                                    ifelse(data$environment=="Aquatic",12,
                                           ifelse(data$environment=="human skin",15,
                                                  ifelse(data$environment=="Soil",18,
                                                         ifelse(data$environment=="Wastewater",21,
														        ifelse(data$environment=="Plants",24,
														               ifelse(data$environment=="shoot",27,
														                      ifelse(data$environment=="root",30,31))))))))))
ggplot(xlim = c(0,35), ylim = c(-12,0)) + 
geom_scatterpie_legend(data$r1, x=2, y=-6,n=4,labeller = function(x) round(x*x*400)) +
geom_hline(yintercept = -11,linetype=2,color = "grey") +
geom_hline(yintercept = -10,linetype=2,color = "grey") +
geom_hline(yintercept = -9,linetype=2,color = "grey") +
geom_hline(yintercept = -8,linetype=2,color = "grey") +
geom_hline(yintercept = -7,linetype=2,color = "grey") +
geom_hline(yintercept = -6,linetype=2,color = "grey") +
geom_hline(yintercept = -5,linetype=2,color = "grey") +
geom_hline(yintercept = -4,linetype=2,color = "grey") +
geom_hline(yintercept = -3,linetype=2,color = "grey") +
geom_hline(yintercept = -2,linetype=2,color = "grey") +
geom_hline(yintercept = -1,linetype=2,color = "grey") +
geom_vline(xintercept = 3, linetype=2,color = "grey") +
geom_vline(xintercept = 6, linetype=2,color = "grey") +
geom_vline(xintercept = 9, linetype=2,color = "grey") +
geom_vline(xintercept = 12, linetype=2,color = "grey") +
geom_vline(xintercept = 15, linetype=2,color = "grey") +
geom_vline(xintercept = 18, linetype=2,color = "grey") +
geom_vline(xintercept = 21, linetype=2,color = "grey") +
geom_vline(xintercept = 24, linetype=2,color = "grey") +
geom_vline(xintercept = 27, linetype=2,color = "grey") +
geom_vline(xintercept = 30, linetype=2,color = "grey") +
geom_scatterpie(aes(x=environment, y=Order, r=r1), data=data,cols=colnames(data)[3:4],color=NA) + 
 #geom_hline(yintercept = -4.5,linetype=5)+
 #geom_vline(xintercept = 7.5,linetype=5)+
 #geom_scatterpie_legend(data$r, x=3.5, y=-1.3,n=4,labeller = function(x) round(2^(x*25)))+ 
  coord_equal()+
  scale_fill_manual(values = color_ce)+
  scale_x_continuous("",expand = c(0,0),limits = c(1,35),breaks = c(3,6,9,12,15,18,21,24,27,30),labels = c("animal gut","human gut","human oral","Aquatic","human skin","Soil","Wastewater","Plants","shoot","root"))+
  scale_y_continuous("",expand = c(0,0),limits = c(-12,0),breaks = -11:-1,labels = c("Rhodobacterales","Actinomycetales","Rhizobiales","Caulobacterales","Enterobacterales","Pseudomonadales","Moraxellales","Burkholderiales","Sphingomonadales","Propionibacteriales","Puniceispirillales"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )
ggsave("point.pdf",width = 20,height = 15,units = "cm")

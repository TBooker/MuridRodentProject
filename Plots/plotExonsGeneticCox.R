rm(list=ls())
library(ggplot2)
library(ggpubr)
library(reshape2)
TommyTheme <-theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.title = element_text(size =15),
    legend.text = element_text(size =13, face = 'italic'),
    strip.text.y = element_text(size = 15)
  )
JC<-function(raw){
  -0.75*(log(1-(4*raw)/3))
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


Mmc     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmc.CoxDistance.csv')
Mmc$lab <- 'M. m. castaneus'
Mmc$rat <- Mmc$pi/JC(Mmc$div2)
Mmc$pi_ref <- Mmc$rat/ mean( Mmc[ (abs(Mmc$dist) > 1500) & abs(Mmc$dist) < 2500  , ]$rat ) 


Mmdg<-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmd_ger.CoxDistance.csv')
Mmdg$lab <- 'M. m. domesticus - Germany'
Mmdg$rat <- Mmdg$pi/JC(Mmdg$div2)
Mmdg$pi_ref <- Mmdg$rat/ mean( Mmdg[ (abs(Mmdg$dist) > 1500) & abs(Mmdg$dist) < 2500  , ]$rat ) 

Mmdf    <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmd_fra.CoxDistance.csv')
Mmdf$lab <- 'M. m. domesticus - France'
Mmdf$rat <- Mmdf$pi/JC(Mmdf$div2)
Mmdf$pi_ref <- Mmdf$rat/ mean( Mmdf[ (abs(Mmdf$dist) > 1500) & abs(Mmdf$dist) < 2500  , ]$rat ) 

Mmdi     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmd_ira.CoxDistance.csv')
Mmdi$lab <- 'M. m. domesticus - Iran'
Mmdi$rat <- Mmdi$pi/JC(Mmdi$div2)
Mmdi$pi_ref <- Mmdi$rat/ mean( Mmdi[ (abs(Mmdi$dist) > 1500) & abs(Mmdi$dist) < 2500  , ]$rat ) 

Mmma     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmm_afg.CoxDistance.csv')
Mmma$lab <- 'M. m. musculus - Afganistan'
Mmma$rat <- Mmma$pi/JC(Mmma$div2)
Mmma$pi_ref <- Mmma$rat/ mean( Mmma[ (abs(Mmma$dist) > 1500) & abs(Mmma$dist) < 2500  , ]$rat )

Mmmc     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmm_cze.CoxDistance.csv')
Mmmc$lab <- 'M. m. musculus - Czech Rep.'
Mmmc$rat <- Mmmc$pi/JC(Mmmc$div2)
Mmmc$pi_ref <- Mmmc$rat/ mean( Mmmc[ (abs(Mmmc$dist) > 1500) & abs(Mmmc$dist) < 2500  , ]$rat )

Mmmk     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmm_kaz.CoxDistance.csv')
Mmmk$lab <- 'M. m. musculus - Kazakhstan'
Mmmk$rat <- Mmmk$pi/JC(Mmmk$div2)
Mmmk$pi_ref <- Mmmk$rat/ mean( Mmmk[ (abs(Mmmk$dist) > 1500) & abs(Mmmk$dist) < 2500  , ]$rat )


Ms     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Ms.CoxDistance.csv')
Ms$rat <- Ms$pi/JC(Ms$div2)
Ms$pi_ref <- Ms$rat/ mean( Ms[ (abs(Ms$dist) > 1500) & abs(Ms$dist) < 2500  , ]$rat )
Ms$lab <- 'M. spretus'

Mm<- rbind(Mmc, Mmdi, Mmdf, Mmdg, Mmmc, Mmmk, Mmma, Ms)
str(Mm)
Mm 

pdf('Work/MuridRodentTroughs/Plots/exons_pi_genetic_Cox.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = pi, col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = pi, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(pi))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,2000))+
  TommyTheme
dev.off()

pdf('Work/MuridRodentTroughs/Plots/exons_piRef_genetic_Cox.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>1,], aes(x= dist, y = pi_ref, col = lab))+
 # geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = pi_ref, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(pi/pi[Ref]))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
#  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,3000))+
  TommyTheme
dev.off()

png('Work/MuridRodentTroughs/Plots/exons_piRef_genetic_Cox.png', width =600, height = 400)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = pi_ref, col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = pi_ref, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(pi/pi[Ref]))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,2000))+
  TommyTheme
dev.off()

pdf('Work/MuridRodentTroughs/Plots/exons_divFam_genetic_Cox.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = 100*JC(div1), col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = 100*JC(div1), col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(italic(d[famulus])*' (%)'), limits = c(0,5))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,2000))+
  TommyTheme
dev.off()

pdf('Work/MuridRodentTroughs/Plots/exons_divRat_genetic_Cox.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = 100*JC(div2), col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = 100*JC(div2), col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(italic(d[Rat])*' (%)'), limits = c(10,20))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.00001,2000))+
  TommyTheme
dev.off()


ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = 1-(pi/thW), col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = pi, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(pi))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,50))+
  TommyTheme





png('Work/TBooker.github.io/img/research/exons_divFam_genetic_Cox.png', width =600, height = 400)
ggplot(data = Mm[Mm$dist>1,], aes(x= dist, y = pi_ref, col = lab))+
  geom_line(data = Mm[Mm$dist<-1,], aes(x= abs(dist), y = pi_ref, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(pi/pi[Ref]))+
  scale_colour_manual('',values = cbPalette)+
  theme_bw()
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,2000))+
  TommyTheme
dev.off()



ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = sites, col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = sites, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_y_continuous(expression(pi))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from Exon ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(1,2000))+
  TommyTheme


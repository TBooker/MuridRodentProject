rm(list=ls())

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
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



Mmc     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmc.castaneusDistance.csv')
Mmc$lab <- 'M. m. castaneus'
Mmc$rat <- Mmc$pi/JC(Mmc$div2)
Mmc$pi_ref <- Mmc$rat/ mean( Mmc[ (abs(Mmc$dist) > 150) & abs(Mmc$dist) < 250  , ]$rat ) 
str(Mmc)


Mmdg<-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmd_ger.castaneusDistance.csv')
Mmdg$lab <- 'M. m. domesticus - Germany'
Mmdg$rat <- Mmdg$pi/JC(Mmdg$div2)
Mmdg$pi_ref <- Mmdg$rat/ mean( Mmdg[ (abs(Mmdg$dist) > 150) & abs(Mmdg$dist) < 250  , ]$rat ) 

Mmdf    <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmd_fra.castaneusDistance.csv')
Mmdf$lab <- 'M. m. domesticus - France'
Mmdf$rat <- Mmdf$pi/JC(Mmdf$div2)
Mmdf$pi_ref <- Mmdf$rat/ mean( Mmdf[ (abs(Mmdf$dist) > 150) & abs(Mmdf$dist) < 250  , ]$rat ) 

Mmdi     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmd_ira.castaneusDistance.csv')
Mmdi$lab <- 'M. m. domesticus - Iran'
Mmdi$rat <- Mmdi$pi/JC(Mmdi$div2)
Mmdi$pi_ref <- Mmdi$rat/ mean( Mmdi[ (abs(Mmdi$dist) > 150) & abs(Mmdi$dist) < 250  , ]$rat ) 

Mmma     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmm_afg.castaneusDistance.csv')
Mmma$lab <- 'M. m. musculus - Afganistan'
Mmma$rat <- Mmma$pi/Mmma$div2
Mmma$pi_ref <- Mmma$rat/ mean( Mmma[ (abs(Mmma$dist) > 150) & abs(Mmma$dist) < 250  , ]$rat )

Mmmc     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmm_cze.castaneusDistance.csv')
Mmmc$lab <- 'M. m. musculus - Czech Rep.'
Mmmc$rat <- Mmmc$pi/Mmmc$div2
Mmmc$pi_ref <- Mmmc$rat/ mean( Mmmc[ (abs(Mmmc$dist) > 150) & abs(Mmmc$dist) < 250  , ]$rat )

Mmmk     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmm_kaz.castaneusDistance.csv')
Mmmk$lab <- 'M. m. musculus - Kazakhstan'
Mmmk$rat <- Mmmk$pi/Mmmk$div1
Mmmk$pi_ref <- Mmmk$rat/ mean( Mmmk[ (abs(Mmmk$dist) > 150) & abs(Mmmk$dist) < 250  , ]$rat )


Ms     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Ms.castaneusDistance.csv')
Ms$rat <- Ms$pi/Ms$div2
Ms$pi_ref <- Ms$rat/ mean( Ms[ (abs(Ms$dist) > 150) & abs(Ms$dist) < 250  , ]$rat )
Ms$lab <- 'M. spretus'






Mm<- rbind(Mmc, Mmdi, Mmdf, Mmdg, Mmmc, Mmmk, Mmma, Ms)

pdf('Work/MuridRodentTroughs/Plots/CNEs_pi_genetic_castaneus.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = pi, col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = pi, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  geom_vline(aes(xintercept = 2.5), alpha=0.3, lty = 2, col = 'red')+
  scale_y_continuous(expression(italic(d[Rat])*'%'), limits = c(10,20))+
  scale_y_continuous(expression(italic(pi)))+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from CNE ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.0001,200), breaks = seq(0,200,30))+
  TommyTheme
dev.off()


pdf('Work/MuridRodentTroughs/Plots/CNEs_piRef_genetic_castaneus.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = pi_ref, col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = pi_ref, col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  geom_vline(aes(xintercept = 2.5), alpha=0.3, lty = 2, col = 'red')+
  scale_y_continuous(expression(italic(d[Rat])*'%'), limits = c(10,20))+
  scale_y_continuous(expression(italic(pi/pi[Ref])))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_colour_manual('',values = cbPalette)+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from CNE ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.000,200), breaks = seq(0,200,30))+
  TommyTheme
dev.off()


pdf('Work/MuridRodentTroughs/Plots/CNEs_divFam_genetic_castaneus.pdf', width =8.0, height = 5.5)

ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = 100*JC(div1), col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = 100*JC(div1), col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  geom_vline(aes(xintercept = 2.5), alpha=0.3, lty = 2, col = 'red')+
  scale_y_continuous(expression(italic(d[famulus])*'%'), limits = c(0,5))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_colour_manual('',values = cbPalette)+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from CNE ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.000,200), breaks = seq(0,200,30))+
  TommyTheme
dev.off()

pdf('Work/MuridRodentTroughs/Plots/CNEs_divRat_genetic_castaneus.pdf', width =8.0, height = 5.5)
ggplot(data = Mm[Mm$dist>0,], aes(x= dist, y = 100*JC(div2), col = lab))+
  geom_line(data = Mm[Mm$dist<0,], aes(x= abs(dist), y = 100*JC(div2), col = lab), alpha = 0.75)+
  geom_line(alpha = 0.75)+
  geom_vline(aes(xintercept = 2.5), alpha=0.3, lty = 2, col = 'red')+
  guides(colour = guide_legend(override.aes = list(size=3)))+
    scale_y_continuous(expression(italic(d[Rat])*'%'), limits = c(10,20))+
  scale_colour_manual('',values = cbPalette)+
  theme_bw()+
  scale_x_sqrt(expression('Genetic Distance from CNE ('*rho*' = '*italic('4N'[e]*'r')*')'), limits = c(0.000,200), breaks = seq(0,200,30))+
  TommyTheme
dev.off()

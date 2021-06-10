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
    strip.text.x = element_text(size = 15),
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
Mmmc$rat <- Mmmc$pi/JC(Mmmc$div2)1
Mmmc$pi_ref <- Mmmc$rat/ mean( Mmmc[ (abs(Mmmc$dist) > 1500) & abs(Mmmc$dist) < 2500  , ]$rat )

Mmmk     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Mmm_kaz.CoxDistance.csv')
Mmmk$lab <- 'M. m. musculus - Kazakhstan'
Mmmk$rat <- Mmmk$pi/JC(Mmmk$div2)
Mmmk$pi_ref <- Mmmk$rat/ mean( Mmmk[ (abs(Mmmk$dist) > 1500) & abs(Mmmk$dist) < 2500  , ]$rat )


Ms     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/exons/Ms.CoxDistance.csv')
Ms$rat <- Ms$pi/JC(Ms$div2)
Ms$pi_ref <- Ms$rat/ mean( Ms[ (abs(Ms$dist) > 1500) & abs(Ms$dist) < 2500  , ]$rat )
Ms$lab <- 'M. spretus'

Mm_exons<- rbind(Mmc, Mmdi, Mmdf, Mmdg, Mmmc, Mmmk, Mmma, Ms)


### Now read in the CNEs

Mmc     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmc.coxDistance.csv')
Mmc$lab <- 'M. m. castaneus'
Mmc$rat <- Mmc$pi/JC(Mmc$div2)
Mmc$pi_ref <- Mmc$rat/ mean( Mmc[ (abs(Mmc$dist) > 150) & abs(Mmc$dist) < 250  , ]$rat ) 
str(Mmc)


Mmdg<-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmd_ger.coxDistance.csv')
Mmdg$lab <- 'M. m. domesticus - Germany'
Mmdg$rat <- Mmdg$pi/JC(Mmdg$div2)
Mmdg$pi_ref <- Mmdg$rat/ mean( Mmdg[ (abs(Mmdg$dist) > 150) & abs(Mmdg$dist) < 250  , ]$rat ) 

Mmdf    <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmd_fra.coxDistance.csv')
Mmdf$lab <- 'M. m. domesticus - France'
Mmdf$rat <- Mmdf$pi/JC(Mmdf$div2)
Mmdf$pi_ref <- Mmdf$rat/ mean( Mmdf[ (abs(Mmdf$dist) > 150) & abs(Mmdf$dist) < 250  , ]$rat ) 

Mmdi     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmd_ira.coxDistance.csv')
Mmdi$lab <- 'M. m. domesticus - Iran'
Mmdi$rat <- Mmdi$pi/JC(Mmdi$div2)
Mmdi$pi_ref <- Mmdi$rat/ mean( Mmdi[ (abs(Mmdi$dist) > 150) & abs(Mmdi$dist) < 250  , ]$rat ) 

Mmma     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmm_afg.coxDistance.csv')
Mmma$lab <- 'M. m. musculus - Afganistan'
Mmma$rat <- Mmma$pi/Mmma$div2
Mmma$pi_ref <- Mmma$rat/ mean( Mmma[ (abs(Mmma$dist) > 150) & abs(Mmma$dist) < 250  , ]$rat )

Mmmc     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmm_cze.coxDistance.csv')
Mmmc$lab <- 'M. m. musculus - Czech Rep.'
Mmmc$rat <- Mmmc$pi/Mmmc$div2
Mmmc$pi_ref <- Mmmc$rat/ mean( Mmmc[ (abs(Mmmc$dist) > 150) & abs(Mmmc$dist) < 250  , ]$rat )

Mmmk     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Mmm_kaz.coxDistance.csv')
Mmmk$lab <- 'M. m. musculus - Kazakhstan'
Mmmk$rat <- Mmmk$pi/Mmmk$div1
Mmmk$pi_ref <- Mmmk$rat/ mean( Mmmk[ (abs(Mmmk$dist) > 150) & abs(Mmmk$dist) < 250  , ]$rat )


Ms     <-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/CNEs/Ms.coxDistance.csv')
Ms$rat <- Ms$pi/Ms$div2
Ms$pi_ref <- Ms$rat/ mean( Ms[ (abs(Ms$dist) > 150) & abs(Ms$dist) < 250  , ]$rat )
Ms$lab <- 'M. spretus'


Mm_cnes<- rbind(Mmc, Mmdi, Mmdf, Mmdg, Mmmc, Mmmk, Mmma, Ms)

Mm_exons$class <- 'Protein-Coding Exons'
Mm_cnes$class <- 'CNEs'

Mm<- rbind(Mm_exons, Mm_cnes)
Mm<- Mm[Mm$dist >=1,]
Mm<- Mm[Mm$dist <=2500,]

names(Mm)
mouse_melt <- melt(Mm, id = c("dist",   "thW",    "div1",   "div2",   "sites",  "lab",    "rat", "class"))

mouse_melt$variable <- factor(mouse_melt$variable, levels = levels(mouse_melt$variable), labels = c('italic(pi)', 'pi/pi[italic(Ref)]'))

mouse_melt$class <- as.factor(mouse_melt$class)

mouse_melt$class <- factor(mouse_melt$class, levels = c('Protein-Coding Exons','CNEs'), labels = c(expression("Protein-Coding~Exons"),expression("CNEs")))

# 
# ggplot(data = Mm, aes(x= dist, y = pi_ref, col = lab))+
#   geom_line()+
#   scale_y_continuous(expression(pi/pi[Ref]))+
#   scale_colour_manual('',values = cbPalette)+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   theme_bw()+
#   facet_grid(~class, scales = 'free')+
#   scale_x_continuous(expression('Genetic Distance from Element ('*rho*' = '*italic('4N'[e]*'r')*')'))+
#   TommyTheme

troughs <- ggplot(data = mouse_melt, aes(x= dist, y = value, col = lab))+
  geom_line()+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  facet_grid(variable~class, scales = 'free', labeller = label_parsed)+
  scale_x_continuous(expression('Genetic Distance from Element ('*rho*' = '*italic('4N'[e]*'r')*')'))+
  TommyTheme+
  theme(
    axis.title.y = element_blank()
  )

pdf('~/Work/MuridRodentTroughs/Plots/troughs.Cox.pdf', width = 10, height = 7)
print(troughs)
dev.off()

Mm$sites<-Mm$sites/1e6

mouse_melt_classic <- melt(Mm, id = c("dist",   "thW",   "div2" ,  "lab", "class",'pi_ref'))

mouse_melt_classic$variable <- factor(mouse_melt_classic$variable, levels = c('pi','div1','rat','sites'), labels = c('italic(pi)', 'italic(d[famulus])','italic(pi/d[famulus])', 'Sites~(Mbp)'))

mouse_melt_classic$class <- as.factor(mouse_melt_classic$class)

mouse_melt_classic$class <- factor(mouse_melt_classic$class, levels = c('Protein-Coding Exons','CNEs'), labels = c(expression("Protein-Coding~Exons"),expression("CNEs")))


troughs_extra <- ggplot(data = mouse_melt_classic, aes(x= dist, y = value, col = lab))+
  geom_line()+
  scale_colour_manual('',values = cbPalette)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()+
  facet_grid(variable~class, scales = 'free', labeller = label_parsed)+
  scale_x_continuous(expression('Genetic Distance from Element ('*rho*' = '*italic('4N'[e]*'r')*')'))+
  TommyTheme+
  theme(
    axis.title.y = element_blank()
  )

pdf('~/Work/MuridRodentTroughs/Plots/trough_extras.Cox.pdf', width = 10, height = 15)
print(troughs_extra)
dev.off()


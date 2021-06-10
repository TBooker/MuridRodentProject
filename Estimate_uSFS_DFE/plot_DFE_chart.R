rm(list = ls())
library(ggplot2)
library(reshape2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

TommyTheme <-theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.title = element_text(size =15),
    legend.text = element_text(size =13, face = 'italic'),
    strip.text.y = element_text(size = 15),
    legend.text.align = 0
  )

## Let's plot the DFE inferred for different site classes and taxa

Mmc_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmc.dDFE.csv', header = F)
Mmc_0$lab <- 'M. m. castaneus - India'
Mmd_fra_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmd_fra.dDFE.csv', header = F)
Mmd_fra_0$lab <- 'M. m. domesticus - France'
Mmd_ger_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmd_ger.dDFE.csv', header = F)
Mmd_ger_0$lab <- 'M. m. domesticus - Germany'
Mmd_ira_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmd_ira.dDFE.csv', header = F)
Mmd_ira_0$lab <- 'M. m. domesticus - Iran'
Mmm_afg_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmm_afg.dDFE.csv', header = F)
Mmm_afg_0$lab <- 'M. m. musculus - Afganistan'
Mmm_cze_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmm_cze.dDFE.csv', header = F)
Mmm_cze_0$lab <- 'M. m. musculus - Czech Rep.'
Mmm_kaz_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Mmm_kaz.dDFE.csv', header = F)
Mmm_kaz_0$lab <- 'M. m. musculus - Kazakhstan'
Ms_0 <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/0-fold/Ms.dDFE.csv', header = F)
Ms_0$lab <- 'Mus spretus'


zeroFoldRaw = rbind(Mmc_0,Mmd_fra_0, Mmd_ger_0, Mmd_ira_0, Mmm_afg_0, Mmm_cze_0, Mmm_kaz_0, Ms_0)

names(zeroFoldRaw) <- c('Taxa', 'Type','Range','Point','Lower','Upper', 'lab')

zeroFold <- melt(zeroFoldRaw, id = c('Taxa', 'Type','Point','Lower','Upper', 'lab'))

Mmc_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE//Mmc.DFE.csv', header = F)
Mmc_c$lab <- 'M. m. castaneus - India'
Mmd_fra_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Mmd_fra.DFE.csv', header = F)
Mmd_fra_c$lab <- 'M. m. domesticus - France'
Mmd_ger_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Mmd_ger.DFE.csv', header = F)
Mmd_ger_c$lab <- 'M. m. domesticus - Germany'
Mmd_ira_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Mmd_ira.DFE.csv', header = F)
Mmd_ira_c$lab <- 'M. m. domesticus - Iran'
Mmm_afg_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Mmm_afg.DFE.csv', header = F)
Mmm_afg_c$lab <- 'M. m. musculus - Afganistan'
Mmm_cze_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Mmm_cze.DFE.csv', header = F)
Mmm_cze_c$lab <- 'M. m. musculus - Czech Rep.'
Mmm_kaz_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Mmm_kaz.DFE.csv', header = F)
Mmm_kaz_c$lab <- 'M. m. musculus - Kazakhstan'
Ms_c <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/CNE/Ms.DFE.csv', header = F)
Ms_c$lab <- 'Mus spretus'


CNE_Raw = rbind(Mmc_c,Mmd_fra_c, Mmd_ger_c, Mmd_ira_c, Mmm_afg_c, Mmm_cze_c, Mmm_kaz_c, Ms_c)

names(CNE_Raw) <- c('Taxa', 'Type','Range','Point','Lower','Upper', 'lab')

CNE <- melt(CNE_Raw, id = c('Taxa', 'Type','Point','Lower','Upper', 'lab'))


Mmc_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmc.DFE.csv', header = F)
Mmc_u$lab <- 'M. m. castaneus - India'
Mmd_fra_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmd_fra.DFE.csv', header = F)
Mmd_fra_u$lab <- 'M. m. domesticus - France'
Mmd_ger_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmd_ger.DFE.csv', header = F)
Mmd_ger_u$lab <- 'M. m. domesticus - Germany'
Mmd_ira_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmd_ira.DFE.csv', header = F)
Mmd_ira_u$lab <- 'M. m. domesticus - Iran'
Mmm_afg_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmm_afg.DFE.csv', header = F)
Mmm_afg_u$lab <- 'M. m. musculus - Afganistan'
Mmm_cze_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmm_cze.DFE.csv', header = F)
Mmm_cze_u$lab <- 'M. m. musculus - Czech Rep.'
Mmm_kaz_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Mmm_kaz.DFE.csv', header = F)
Mmm_kaz_u$lab <- 'M. m. musculus - Kazakhstan'
Ms_u <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/plot/UTR/Ms.DFE.csv', header = F)
Ms_u$lab <- 'Mus spretus'


UTR_Raw = rbind(Mmc_u,Mmd_fra_u, Mmd_ger_u, Mmd_ira_u, Mmm_afg_u, Mmm_cze_u, Mmm_kaz_u, Ms_u)

names(UTR_Raw) <- c('Taxa', 'Type','Range','Point','Lower','Upper', 'lab')

UTR <- melt(UTR_Raw, id = c('Taxa', 'Type','Point','Lower','Upper', 'lab'))


allElements = rbind(zeroFold, CNE, UTR)

nameLabels <-  c(expression(italic("M. m. castaneus")), 
                 expression(italic("M. m. domesticus")*" - France") ,  
                 expression(italic("M. m. domesticus")*" - Germany") ,  
                 expression(italic("M. m. domesticus") *" - Iran") , 
                 expression(italic("M. m. musculus")*" - Afghanistan"), 
                 expression(italic("M. m. musculus")*" - Czech Rep."), 
                 expression(italic("M. m. musculus")*" - Kazakhstan"), 
                 expression(italic("M. spretus")))


allElements$Type<- factor(allElements$Type, levels = levels(allElements$Type), labels = c('0-fold Nonsynonymous', 'CNEs','UTRs'))


dfe_plot <- ggplot(data = allElements, aes(x = value, y = Point, fill = lab))+
  geom_bar(stat='identity', position = 'dodge', alpha = 0.65)+
  geom_errorbar(aes(ymin = Lower, ymax = Upper, col = lab),width=.2, position=position_dodge(.9))+ 
  scale_y_continuous('Proportion of Mutations in Range', limits = c(0,1))+
  scale_x_discrete(expression('Strength of Purifying Selection (2'*italic(N[e]*s[d])*')'))+
  scale_fill_manual('',
                    values = cbPalette,
                    labels = nameLabels)+
  scale_colour_manual('',
                      values = cbPalette,
                      labels = nameLabels)+
  facet_grid(Type~.,  labeller = label_context)+
  theme_bw()+
  TommyTheme+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf('~/Work/MuridRodentTroughs/Plots/DFE_plot.pdf', width = 11, height = 9)
print(dfe_plot)
dev.off()

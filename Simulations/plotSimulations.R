rm(list=ls())

library(ggplot2)
## Read in all the simulation data and 

exons_eq_cast<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc.exons.eq.genetic.castMap.csv')
exons_eq_cast$map <- 'LD-based'
exons_eq_cast$demography <- 'Constant Size'

exons_neq_cast<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc.exons.neq.genetic.castMap.csv')
exons_neq_cast$map <- 'LD-based'
exons_neq_cast$demography <- '3-epoch'

exons_eq_cox<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc.exons.eq.genetic.coxMap.csv')
exons_eq_cox$map <- 'Pedigree-based'
exons_eq_cox$demography <- 'Constant Size'

exons_neq_cox<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc.exons.neq.genetic.coxMap.csv')
exons_neq_cox$map <- 'Pedigree-based'
exons_neq_cox$demography <- '3-epoch'

exons <- rbind(exons_eq_cast, exons_neq_cast, exons_eq_cox,exons_neq_cox)
exons$Element <- 'Protein-Coding Exons'

cnes_eq_cast<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc_3epoch_castMap/Mmc.CNE.eq.genetic.castaneus.csv')
cnes_eq_cast$map <- 'LD-based'
cnes_eq_cast$demography <- 'Constant Size'

cnes_neq_cast<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc_3epoch_castMap/Mmc.CNE.neq.genetic.castaneus.csv')
cnes_neq_cast$map <- 'LD-based'
cnes_neq_cast$demography <- '3-epoch'

cnes_eq_cox<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc_3epoch_coxMap_SNPs/Mmc.CNE.eq.genetic.Cox.csv')
cnes_eq_cox$map <- 'Pedigree-based'
cnes_eq_cox$demography <- 'Constant Size'

cnes_neq_cox<-read.csv('~/work/MouseAdap/Simulations/simulationRuns/Mmc_3epoch_coxMap_SNPs/Mmc.CNE.neq.genetic.Cox.csv')
cnes_neq_cox$map <- 'Pedigree-based'
cnes_neq_cox$demography <- '3-epoch'

cnes <- rbind(cnes_eq_cast, cnes_neq_cast)
cnes <- rbind(cnes_eq_cast, cnes_neq_cast, cnes_eq_cox, cnes_neq_cox)
cnes$Element <- 'CNEs'

all_sims <- rbind(exons, cnes)
all_sims$Element<- as.factor(all_sims$Element)
all_sims$demography<- as.factor(all_sims$demography)
levels(all_sims$demography)


all_sims <- na.omit(all_sims[all_sims$dist>=1,])


ggplot(data = all_sims, aes(x = abs(dist), y = pi/0.01, col = map))+
  geom_line()+
  facet_grid(demography~Element, scales = 'free')+
  scale_y_continuous(expression(italic('B = ') * pi/pi[0]))+
  theme_bw()


rm(list = ls())
### Script to bootstrap the estiamtes of fitness change and get the ratio of adaptation in the two classes of sites.

## I'll use a truncated normal distribution for the sampling so that we don't get negative values

library(msm)

## This function gets the fitness change (delta W) for a particular set of selection parameters
delta_W <- function(Nes_1, Nes_2, pa_1, pa_2, sites){
  s1_squared <- (Nes_1/426200/2)^2
  s2_squared <- (Nes_2/426200/2)^2

  dW1 <- s1_squared * pa_1 * sites * 1e6 * 5.4e-9   
  dW2 <- s2_squared * pa_2 * sites * 1e6 * 5.4e-9   
  return( dW1 + dW2 )

}

  ### LD-based castaneus map
exon_2Nes_mean_1 = 19200
exon_2Nes_se_1 = 1560
    
exon_2Nes_mean_2 = 228
exon_2Nes_se_2 = 28.7

exon_pa_mean_1 = 0.000009
exon_pa_se_1 = 0.000001

exon_pa_mean_2 = 0.00098
exon_pa_se_2 = 0.00014


cne_2Nes_mean_1 = 448
cne_2Nes_se_1 = 62.1

cne_2Nes_mean_2 = 24.2
cne_2Nes_se_2 = 6.5

cne_pa_mean_1 = 0.00030
cne_pa_se_1 = 0.000069

cne_pa_mean_2 = 0.0078
cne_pa_se_2 = 0.0023

exon_sites = 24
cne_sites = 54

    
  
  
castaneus_point_estimate = delta_W(exon_2Nes_mean_1, exon_2Nes_mean_2, exon_pa_mean_1, exon_pa_mean_2, exon_sites)/delta_W(cne_2Nes_mean_1, cne_2Nes_mean_2, cne_pa_mean_1, cne_pa_mean_2, cne_sites)


exon_2Nes_1_boots <- rtnorm( 10000, exon_2Nes_mean_1, exon_2Nes_se_1, lower = 0)
exon_2Nes_2_boots <- rtnorm( 10000, exon_2Nes_mean_2, exon_2Nes_se_2, lower = 0)

exon_pa_1_boots <- rtnorm( 10000, exon_pa_mean_1, exon_pa_se_1, lower = 0)
exon_pa_2_boots <- rtnorm( 10000, exon_pa_mean_2, exon_pa_se_2, lower = 0)

cne_2Nes_1_boots <- rtnorm( 10000, cne_2Nes_mean_1, cne_2Nes_se_1, lower = 0)
cne_2Nes_2_boots <- rtnorm( 10000, cne_2Nes_mean_2, cne_2Nes_se_2, lower = 0)

cne_pa_1_boots <- rtnorm( 10000, cne_pa_mean_1, cne_pa_se_1, lower = 0)
cne_pa_2_boots <- rtnorm( 10000, cne_pa_mean_2, cne_pa_se_2, lower = 0)


castaneus_ratio  = delta_W(exon_2Nes_1_boots, exon_2Nes_2_boots, exon_pa_1_boots, exon_pa_2_boots, exon_sites)/delta_W(cne_2Nes_1_boots, cne_2Nes_2_boots, cne_pa_1_boots, cne_pa_2_boots, cne_sites)

hist(castaneus_ratio, breaks = 50, xlim = c(0,100), main = expression("LD-based map based on data from" * italic(" M. m. castaneus")))
abline(v =1 , col = "red", lty = 2)
abline(v =point_estimate , col = "black", lwd = 2)

###### Cox map

exon_2Nes_mean_1 = 6173
exon_2Nes_se_1 = 2645

exon_2Nes_mean_2 = 208
exon_2Nes_se_2 = 105

exon_pa_mean_1 = 0.000008
exon_pa_se_1 = 0.000005

exon_pa_mean_2 = 0.00035
exon_pa_se_2 = 0.00018


cne_2Nes_mean_1 = 1910
cne_2Nes_se_1 = 673

cne_2Nes_mean_2 = 7.0
cne_2Nes_se_2 = 3.5

cne_pa_mean_1 = 0.000013
cne_pa_se_1 = 0.000006

cne_pa_mean_2 = 0.0178
cne_pa_se_2 = 0.0129

exon_sites = 24
cne_sites = 54




cox_point_estimate = delta_W(exon_2Nes_mean_1, exon_2Nes_mean_2, exon_pa_mean_1, exon_pa_mean_2, exon_sites)/delta_W(cne_2Nes_mean_1, cne_2Nes_mean_2, cne_pa_mean_1, cne_pa_mean_2, cne_sites)


exon_2Nes_1_boots <- rtnorm( 1000, exon_2Nes_mean_1, exon_2Nes_se_1, lower = 0)
exon_2Nes_2_boots <- rtnorm( 1000, exon_2Nes_mean_2, exon_2Nes_se_2, lower = 0)

exon_pa_1_boots <- rtnorm( 1000, exon_pa_mean_1, exon_pa_se_1, lower = 0)
exon_pa_2_boots <- rtnorm( 1000, exon_pa_mean_2, exon_pa_se_2, lower = 0)

cne_2Nes_1_boots <- rtnorm( 1000, cne_2Nes_mean_1, cne_2Nes_se_1, lower = 0)
cne_2Nes_2_boots <- rtnorm( 1000, cne_2Nes_mean_2, cne_2Nes_se_2, lower = 0)

cne_pa_1_boots <- rtnorm( 1000, cne_pa_mean_1, cne_pa_se_1, lower = 0)
cne_pa_2_boots <- rtnorm( 1000, cne_pa_mean_2, cne_pa_se_2, lower = 0)

cox_ratio  = delta_W(exon_2Nes_1_boots, exon_2Nes_2_boots, exon_pa_1_boots, exon_pa_2_boots, exon_sites)/delta_W(cne_2Nes_1_boots, cne_2Nes_2_boots, cne_pa_1_boots, cne_pa_2_boots, cne_sites)

hist(cox_ratio  ,main = expression("LD-based map based on data from" * italic(" M. m. castaneus")))
abline(v =1 , col = "red", lty = 2)
abline(v =point_estimate , col = "black", lwd = 2)


both<- data.frame( "Pedigree.based.map" = cox_ratio, "LD.based.map" = castaneus_ratio)

library(ggplot2)
library(reshape2)

options(scipen=999)

cox_point_estimate
quantile(cox_ratio, c(0.025,0.975))

castaneus_point_estimate
quantile(castaneus_ratio, c(0.025,0.975))

lineHolder <- data.frame(variable = c("Pedigree.based.map", "LD.based.map"), pointEstimate = c(cox_point_estimate, castaneus_point_estimate))

ggplot(data = melt(both))+
  geom_density( aes(x = value), fill = "grey")+
  scale_x_log10(expression(Delta*W["Exons"]/Delta*W["CNEs"]))+
  facet_wrap(~variable)+
  geom_vline(xintercept = 1)+
  geom_vline(data = lineHolder, aes(xintercept = pointEstimate), col = "red", lwd = 1.3, lty = 3)+
  theme_bw()



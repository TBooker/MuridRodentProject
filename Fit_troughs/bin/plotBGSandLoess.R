rm(list=ls())

library(ggplot2)

exon_cast<-read.csv('~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/Mmc.exons.eq.genetic.castMap.csv')

exon_cast$B<-exon_cast$pi/0.01
exon_cast<-exon_cast[order(abs(exon_cast$dist)),]

loessMod10.exon_cast <- loess(B ~ abs(dist), data=exon_cast, span=0.3, weights = sites) # 30% smoothing span
smoothed10.exon_cast <- predict(loessMod10.exon_cast, se=T) 

plot(exon_cast$B, x=abs(exon_cast$dist), type="p", xlim=c(0.01,3000))
lines(smoothed10.exon_cast$fit, x=abs(exon_cast$dist), col="red")

exonD_cast <- data.frame(cbind(smoothed10.exon_cast$fit, smoothed10.exon_cast$s, abs(exon_cast$dist), abs(exon_cast$B)))
exonD_cast$lab<- 'Protein-Coding Exons'
exon_cast$Bsmooth <- smoothed10.exon_cast$fit
exon_cast$map <- 'castaneus - LD-based'
exon_cast$lab <- 'Protein-Coding Exons'

write.csv(exon_cast, file = '~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/Exons.GeneticDistance.castaneus.LogBins.Loess.csv')

#########
#########
#########

exon_cox<-na.omit(read.csv('~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/Mmc.exons.eq.genetic.coxMap.csv'))

exon_cast$dist
exon_cox$dist

exon_cox$B<-exon_cox$pi/0.01
#exon_cox<-exon_cox[exon_cox$dist>0.0000001,]
exon_cox<-exon_cox[order(abs(exon_cox$dist)),]

loessMod10.exon_cox <- loess(B ~ dist, data=exon_cox, span=0.3, weights = sites) # 30% smoothing span
smoothed10.exon_cox <- predict(loessMod10.exon_cox, se=T) 

plot(exon_cox$B, x=abs(exon_cox$dist), type="p", xlim=c(0.01,3000))
lines(smoothed10.exon_cox$fit, x=abs(exon_cox$dist), col="red")

exonD_cox <- data.frame(cbind(smoothed10.exon_cox$fit, smoothed10.exon_cox$s, abs(exon_cox$dist), abs(exon_cox$B)))
exonD_cox$lab<- 'Protein-Coding Exons'
exon_cox$Bsmooth <- smoothed10.exon_cox$fit
exon_cox$map <- 'Cox - Pedigree-based'
exon_cox$lab <- 'Protein-Coding Exons'

write.csv(exon_cox, file = '~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/Exons.GeneticDistance.cox.LogBins.Loess.csv')

exons<- rbind( exon_cox, exon_cast)


#########
#########
#########


cne_cast<-na.omit(read.csv('~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/Mmc.CNE.eq.genetic.castaneus.csv'))

cne_cast$B<-cne_cast$pi/0.01
cne_cast<-cne_cast[order(abs(cne_cast$dist)),]

loessMod10.cne_cast <- loess(B ~ abs(dist), data= cne_cast, span=0.3, weights = sites) # 30% smoothing span
smoothed10.cne_cast <- predict(loessMod10.cne_cast, se= T) 

plot(cne_cast$B, x=abs(cne_cast$dist), type="p", xlim=c(0.01,300))
lines(smoothed10.cne_cast$fit, x=abs(cne_cast$dist), col="red")

cneD_cast <- data.frame(cbind(as.numeric(smoothed10.cne_cast$fit), as.numeric(smoothed10.cne_cast$s), as.numeric(abs(cne_cast$dist)), as.numeric(abs(cne_cast$B))))
cneD_cast$lab <- 'Conserved Non-Coding Elements'

cne_cast$Bsmooth <- smoothed10.cne_cast$fit
cne_cast$map <- 'castaneus - LD-based'
cne_cast$lab <- 'CNEs'

write.csv(cne_cast, file = '~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/CNEs.GeneticDistance.castaneus.LogBins.Loess.csv')

#########
#########
#########


cne_cox<-na.omit(read.csv('~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/Mmc.CNE.eq.genetic.Cox.csv'))
cne_cox$B<-cne_cox$pi/0.01
cne_cox<-cne_cox[order(abs(cne_cox$dist)),]

loessMod10.cne_cox <- loess(B ~ abs(dist), data= cne_cox, span=0.3, weights = sites) # 30% smoothing span
smoothed10.cne_cox <- predict(loessMod10.cne_cox, se= T) 

plot(cne_cox$B, x=abs(cne_cox$dist), type="p", xlim=c(0.01,300))
lines(smoothed10.cne_cox$fit, x=abs(cne_cox$dist), col="red")

cneD_cox <- data.frame(cbind(as.numeric(smoothed10.cne_cox$fit), as.numeric(smoothed10.cne_cox$s), as.numeric(abs(cne_cox$dist)), as.numeric(abs(cne_cox$B))))
cneD_cox$lab <- 'Conserved Non-Coding Elements'

cne_cox$Bsmooth <- smoothed10.cne_cox$fit
cne_cox$map <- 'Cox - Pedigree-based'
cne_cox$lab <- 'CNEs'

write.csv(cne_cox, file = '~/Work/MuridRodentTroughs/dataAnalysis/troughFitting/B/CNEs.GeneticDistance.cox.LogBins.Loess.csv')


cnes <- rbind(cne_cast, cne_cox)
# 
# 
# pdf('~/PhD/Coding/Estimate_selection_from_troughs/BGS/BGSplotLoess.pdf', height = 8, width = 10)
# 
# ggplot(data = cneD_cast, aes(x= X3, y = X4))+
#   geom_point(lwd = 0.8)+
#   geom_line(aes(x= X3, y = X1), col = 'red')+
# #  geom_ribbon(aes(x= X3, ymax = X1+(2*X2), ymin = X1-(2*X2)), fill = 'red', alpha = 0.35)+
#   geom_point(data = exonD_cox, aes(x= X3, y = X4),lwd = 0.8)+
#   geom_line(data = exonD_cox, aes(x= X3, y = X1), col = 'red')+
# #  geom_ribbon(data = exonD_cast, aes(x= X3, ymax = X1+(2*X2), ymin = X1-(2*X2)), fill = 'red', alpha = 0.35)+
#   facet_grid( ~ lab, scales = 'free_x')+
#   xlab(expression('Distance from Element (4'*N[e]*'r)'))+
#   ylab('B')+
#   theme_bw()+
#   theme(
#     axis.title.x = element_text(size=14,angle=0),
#     axis.title.y = element_text(size=20,angle=0,vjust=0.5, face = 'italic'),
#     axis.text.x = element_text(size=12,angle=0),
#     axis.text.y = element_text(size=12,angle=0),
#     legend.text = element_text(size =13),
#     strip.text.x = element_text(size = 15)
#   )
# dev.off()

all_sites<- rbind( exons, cnes)
all_sites <- all_sites[all_sites$dist > 1,]

BGS<- ggplot(data = all_sites, aes(x= dist, y = B))+
  geom_point(lwd = 0.8, aes(col = map), alpha = 0.5)+
  geom_line(aes(x= dist, y = Bsmooth, col = map))+
  facet_grid(. ~ lab, scales = 'free_x')+
  scale_color_brewer("Recombination Map",palette = "Dark2")+
  xlab(expression('Distance from Element ('*italic('4'*N[e]*'r')*')'))+
  ylab(expression(italic("B = "* pi/pi[0])))+
  scale_y_continuous(limits = c(0.7,1))+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    legend.title = element_text(size =15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15)
  )

pdf('~/Work/MuridRodentTroughs/Plots/makeTheBGS.pdf', width = 10, height = 4)
print(BGS)
dev.off()

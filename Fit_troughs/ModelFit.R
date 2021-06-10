rm(list=ls())

library(ggplot2)
library(wesanderson)
options(scipen=10000)
TommyTheme <-   
  theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,vjust=0.5, angle = 0),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    strip.text.y = element_text(size = 13, face = 'italic'),
    strip.text.x = element_text(size = 13, face = 'italic'),
    legend.text = element_text(size = 13, face = 'italic')
  )


pal <- c(wes_palette(n = 2, name = "Royal1", type = "continuous"),'black')

  ex<-read.csv('/Users/s0784966/Work/MuridRodentTroughs/dataAnalysis/troughFitting/modelOutput/exons.cox.0.0081.2.csv')
ex$element<-'Protein-Coding Exons'

cn<-read.csv('/Users/s0784966/Work/MuridRodentTroughs/dataAnalysis/troughFitting/modelOutput/CNEs.cox.0.0091.2.csv')
cn$element<-'Conserved Non-Coding Elements'

ex2 <- ex
ex2$distance <- ex2$distance*-1
z<- rbind(ex,cn)

#z$element <- factor(z$element, levels = c('Protein-Coding Exons', 'Conserved Non-Coding Elements'), labels = c('Protein-Coding Exons', 'Conserved Non-Coding Elements'))
#
#cairo_pdf('/Users/s0784966//PhD/Coding/Estimate_selection_from_troughs/mouse_analysis/ModelFit.pdf', width = 8, height = 4)
cairo_pdf('/Users/s0784966/UBC/Talks/EcoEvo/ModelFit_Cox.pdf', width = 8, height = 4)

ggplot(data = z, aes(x = distance, y = pi, col = 'Observed'))+
  geom_line(lwd = 1.1)+
  geom_line(aes(x = distance, y = BGS, col = 'Background Selectio'), lwd = 1.05, lty = 'dashed')+
  geom_line(aes(x = distance, y = fitted, col = 'Fitted'), lwd = 1.2, lty = 'dashed')+
  xlab(expression('Distance to Exon (4'*N[e]*'r)'))+
  facet_grid(~element,scales = 'free_x')+
  scale_color_manual('',values = pal)+
  ylab(expression(pi/pi[0]))+
  theme_bw()+
  TommyTheme

dev.off()


cairo_pdf('/Users/s0784966/Work/MuridRodentTroughs/dataAnalysis/troughFitting/ModelFit_Cox.pdf', width = 10, height = 4)

ggplot(data = z, aes(x = distance, y = pi, col = 'Observed'))+
  geom_line(lwd = 1.1)+
  geom_line(aes(x = distance, y = BGS, col = 'B (From simulations)'), lwd = 1.05, lty = 'dashed')+
  geom_line(aes(x = distance, y = fitted, col = 'Model fit'), lwd = 1.2, lty = 'dashed')+
  xlab(expression('Distance from Element ('*italic(4*N[e]*'r)')))+
  facet_grid(~element,scales = 'free_x')+
  scale_color_manual('',values = pal)+
  ylab(expression(pi/pi[0]))+
  theme_bw()+
  TommyTheme

dev.off()



### Make a figure with just the exons for research statement
cairo_pdf('/Users/s0784966/UBC/Applications/Loftus-Hills/fig/exonsModelFit.pdf', width = 7, height = 4)

ggplot(data = z[z$element == "Protein-Coding Exons",], aes(x = distance, y = pi, col = 'Observed'))+
  geom_line(lwd = 1.1)+
  geom_line(aes(x = distance, y = BGS, col = 'Background selection'), lwd = 1.05, lty = 'dashed')+
  geom_line(aes(x = distance, y = fitted, col = 'Background selection +\n Selective Sweeps'), lwd = 1.2, lty = 'dashed')+
  xlab(expression('Distance to Exon (4'*N[e]*'r)'))+
  scale_color_manual('',values = pal)+
  ylab(expression(pi/pi[0]))+
  theme_bw()+
  TommyTheme
dev.off()

rm(list= ls())

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Let's plot the summary stats for each taxa and funcitonal element.

cnes <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/CNEs.csv')
cnes$elements <- 'CNEs'
cneflanks <- read.csv('Work/MuridRodentTroughs/dataAnalysis/dfe/CNEflanks.csv')
cneflanks$elements <- 'CNE flanks'
utr_3 <- read.csv('Work/MuridRodentTroughs/DFE/3p_UTR.csv')
utr_3$elements <- "3' UTR"
utr_5 <- read.csv('Work/MuridRodentTroughs/DFE/5p_UTR.csv')
utr_5$elements <- "5' UTR"

data<- rbind( cnes, cneflanks, utr_3, utr_5)
data$statistic_clean <- factor(data$statistic, levels = levels(data$statistic), labels = c('Divergence', 'Diversity',"Sites", "Tajima's D", "Upwards Skew"))


str(data)
ggplot(data = data, aes( x = value, y= taxa, col = taxa) )+
  geom_point(cex = 2.5, shape = 'triangle')+
  geom_errorbarh(aes(xmin=lower, xmax=upper), height = 0)+
  scale_color_manual(values = cbPalette,
                      name = element_blank(), guide = guide_legend(reverse = TRUE)) +
  facet_grid(elements ~statistic_clean, scales = 'free_x')+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_blank(),
    axis.ticks.y =element_blank(),
    legend.title = element_text(size =15),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

cneDFE<- read.csv('Work/MuridRodentTroughs/DFE/CNE.DFE.csv')

cneDFE[cneDFE$stat=='S_d',]$est<- log10(abs(cneDFE[cneDFE$stat=='S_d',]$est))
cneDFE[cneDFE$stat=='S_d',]$upper<- log10(abs(cneDFE[cneDFE$stat=='S_d',]$upper))
cneDFE[cneDFE$stat=='S_d',]$lower<- log10(abs(cneDFE[cneDFE$stat=='S_d',]$lower))

ggplot(data = cneDFE, aes( x= est, y = taxa , col = taxa))+
    geom_point()+
  geom_errorbarh(aes(xmin=lower, xmax=upper), height = 0)+
  facet_wrap(~stat,scales = 'free')+
  scale_x_continuous('')+
  scale_color_manual(values = cbPalette,
                     name = element_blank(), guide = guide_legend(reverse = TRUE)) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_blank(),
    axis.ticks.y =element_blank(),
    legend.title = element_text(size =15),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

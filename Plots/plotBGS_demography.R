rm(list = ls())

mmc_eq_cast <- read.csv("~/Work/MuridRodentTroughs/Simulations/AnalysisFiles/Mmc.exons.eq.genetic.castMap.csv")
mmc_eq_cast$map <- "LD-based map"
mmc_eq_cast$population <- "Constant Size"
mmc_eq_cast$pi_ref <- mmc_eq_cast$pi/0.01
mmc_neq_cast <- read.csv("~/Work/MuridRodentTroughs/Simulations/AnalysisFiles/Mmc.exons.neq.genetic.castMap.csv")
mmc_neq_cast$map <- "LD-based map"
mmc_neq_cast$population <- "3-epoch"
mmc_neq_cast$pi_ref <- mmc_neq_cast$pi/0.004195

mmc_eq_cox <- read.csv("~/Work/MuridRodentTroughs/Simulations/AnalysisFiles/Mmc.exons.eq.genetic.coxMap.csv")
mmc_eq_cox$map <- "Cox map"
mmc_eq_cox$population <- "Constant Size"
mmc_eq_cox$pi_ref <- mmc_eq_cox$pi/0.01
mmc_neq_cox <- read.csv("~/Work/MuridRodentTroughs/Simulations/AnalysisFiles/Mmc.exons.neq.genetic.coxMap.csv")
mmc_neq_cox$map <- "Cox map"
mmc_neq_cox$population <- "3-epoch"
mmc_neq_cox$pi_ref <- mmc_neq_cox$pi/0.004195

mmc <- rbind(mmc_eq_cast, mmc_eq_cox, mmc_neq_cast, mmc_neq_cox)

library(ggplot2)

BGS<-ggplot(data = mmc, aes(x= abs(dist), y = pi_ref, col = population))+
  geom_smooth(span = 0.2)+
#  geom_point(data = mmc_neq, aes(x= abs(dist), y = pi/0.004), col = "red")+
  xlab(expression('Distance from Element ('*italic('4'*N[e]*'r')*')'))+
  ylab(expression(italic("B = "* pi/pi[0])))+
  scale_x_continuous(limits = c(1,3000))+
  facet_grid(.~map, scales = "free_y")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15)
  )

pdf('~/Work/MuridRodentTroughs/Plots/BGS_demography.pdf', width = 10, height = 4)
print(BGS)
dev.off()




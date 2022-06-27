
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
library(yaml)
library(stringr)

mcol <-
  wes_palette("Darjeeling1", 5, type = "continuous")[c(1:3, 5)]

lh_csv <- read.csv("results/lh/lhd.csv", row.names = 1, as.is = T)

sigma_config <-
  read_yaml("config/config.yaml")$selfing_rates_backward_in_time
lh_csv$tsigma <-
  sapply(lh_csv$demographic_scenario, function(demography_index) {
    tsigma_string <- sigma_config[[demography_index + 1]][[2]][[2]]
    tsigma_string <- str_replace(tsigma_string, "[^0-9.]", "")
    return(as.numeric(tsigma_string))
  })

lhdf <-
  lh_csv %>%
  # subset(infmodel %in% c("Constant", "OneTransition")) %>%
  pivot_wider(names_from = infmodel, values_from = lh, -infile)

# lhr
lhdf$OneTransition_vs_Constant <- lhdf$OneTransition - lhdf$Constant
lhdf$Free_vs_OneTransition <- lhdf$Free - lhdf$OneTransition
lhdf$Free_vs_Constant <- lhdf$Free - lhdf$Constant

p_lh <- lhdf %>%
  ggplot() +
  geom_point(
    data = lhdf,
    aes(tsigma, OneTransition_vs_Constant, col = "OneTransVsConst"),
    position = position_jitter(width = 0.015),
    alpha=0.4
  ) +
  stat_summary(data = lhdf,
               mapping = aes(tsigma, OneTransition_vs_Constant, col = "OneTransVsConst"),
               fun.data = mean_se,
               geom = "line") +
  stat_summary(data = lhdf,
               mapping = aes(tsigma, OneTransition_vs_Constant, col = "OneTransVsConst"),
               fun.data = mean_se) +
  geom_point(
    data = lhdf,
    aes(tsigma, Free_vs_OneTransition, col = "Free_vs_OneTransition"),
    position = position_jitter(width = 0.015),
    alpha=0.4
  ) +
  stat_summary(data = lhdf,
               mapping = aes(tsigma, Free_vs_OneTransition, col = "Free_vs_OneTransition"),
               fun.data = mean_se,
               geom = "line") +
  stat_summary(data = lhdf,
               mapping = aes(tsigma, Free_vs_OneTransition, col = "Free_vs_OneTransition"),
               fun.data = mean_se) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = mcol) +
  theme(legend.position = "bottom")

p_lh_no_stat_summary <- lhdf %>%
  ggplot() +
  geom_point(
    data = lhdf,
    aes(tsigma, OneTransition_vs_Constant, col = "OneTransVsConst"),
    position = position_jitter(width = 0.015),
    # alpha=0.4
  ) +
  geom_point(
    data = lhdf,
    aes(tsigma, Free_vs_OneTransition, col = "Free_vs_OneTransition"),
    position = position_jitter(width = 0.015),
    # alpha=0.4
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = mcol) +
  theme(legend.position = "bottom")


df <- read.csv("results/esmc/point_estim_for_OneTrans.csv", row.names = 1)

p_tsigma_ <- df %>%
  ggplot(aes(tsigma, tsigma_estim, col=infmodel))+
  geom_abline()+
  geom_point(position = position_jitterdodge(jitter.width = 0.02, dodge.width = 0.0))+
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = mcol) +
  theme(legend.position = "bottom")

plot_grid(p_tsigma_, p_lh, p4, ncol = 1, align = "hv", labels = "AUTO")

p_tsigma_+
  theme(aspect.ratio = 1)
  
p_lh +
  theme(aspect.ratio = 1)

p_lh_no_stat_summary+
  theme(aspect.ratio = 1)


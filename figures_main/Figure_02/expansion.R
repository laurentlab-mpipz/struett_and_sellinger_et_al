setwd("teSMC_exp/")


sigma_config <-
  read_yaml("config/config.yaml")$selfing_rates_backward_in_time


df <-
  read.csv("results/esmc/point_estim_for_OneTrans.csv", row.names = 1)


message("Pop size changes are hard coded for the plot")
df$t_pop <- sapply(df$demography, function(x) {
  if (x <= 21) {
    tpop <- 10010
  } else if (x >= 22) {
    tpop <- 40010
  }
  
  return(tpop)
})

tpop1 <- 10010
df1 <- df %>%
  subset(t_pop == tpop1)
tpop2 <- 40010
df2 <- df %>%
  subset(t_pop == tpop2)


p_tsigma_expansion_tNe10k <- df1 %>%
  ggplot(aes(tsigma, tsigma_estim, col = infmodel)) +
  geom_vline(xintercept = tpop1, col = "gray") +
  geom_abline() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.0),
             alpha = 0.8,
             shape = 16) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = mlimits) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = mlimits) +
  scale_color_manual(values = mcol) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

show(p_tsigma_expansion_tNe10k)


p_tsigma_expansion_tNe40k <- df2 %>%
  ggplot(aes(tsigma, tsigma_estim, col = infmodel)) +
  geom_vline(xintercept = tpop2, col = "gray") +
  geom_abline() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.0),
             alpha = 0.8,
             shape = 16) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = mlimits) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel,
                     limits = mlimits) +
  scale_color_manual(values = mcol) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1
  )

show(p_tsigma_expansion_tNe40k)


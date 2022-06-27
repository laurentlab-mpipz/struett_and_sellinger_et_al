setwd("teSMC_const/")

sigma_config <-
  read_yaml("config/config.yaml")$selfing_rates_backward_in_time


df <-
  read.csv("results/esmc/point_estim_for_OneTrans.csv", row.names = 1)


p_tsigma_constant <- df %>%
  ggplot(aes(tsigma, tsigma_estim, col = infmodel)) +
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


show(p_tsigma_constant)


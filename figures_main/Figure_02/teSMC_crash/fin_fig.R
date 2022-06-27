
# library(tidyr)

source("plot_gof_by_mse.R")
source("plot_lh_and_true_vs_estim.R")

df <- df1  # time pop size change 1
# df <- df2  # time pop size change 2

names(df)[names(df) == "infmodel"] = "model_optim"

logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})

t_pop_size_change <- df$t_pop %>% unique()
stopifnot(length(t_pop_size_change) == 1)

# main
set.seed(1234)
p_main <- df %>%
  ggplot(aes(tsigma, tsigma_estim, col = model_optim)) +
  geom_vline(xintercept = t_pop_size_change, col="black", size=0.1)+
  geom_abline() +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.06,
                                    dodge.width = 0.0),
    shape = 16,
    size = 1.7,
    alpha = 0.5
  ) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_color_manual(values = mcol) +
  labs(
    x = expression("t"[sigma] ~ " [log"[10] ~ "]"),
    y = expression(hat("t"[sigma]) ~ " [log"[10] ~ "]"),
    col = "Optimization"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", size = 0.9),
    axis.line = element_blank(),
    axis.ticks.length = unit(-3, units = "pt")
  )

pdf("main_figure_tsigma_estim_tpop1.pdf")
# pdf("main_figure_tsigma_estim_tpop2.pdf")
show(p_main)
dev.off()

p_lh <- lhdf %>%
  ggplot() +
  geom_point(
    data = lhdf,
    aes(tsigma, OneTransition_vs_Constant, col = "OneTransVsConst"),
    position = position_jitter(width = 0.015),
    shape = 16,
    size = 1.7,
    alpha = 0.9
  ) +
  geom_point(
    data = lhdf,
    aes(tsigma, Free_vs_OneTransition, col = "Free_vs_OneTransition"),
    position = position_jitter(width = 0.015),
    shape = 16,
    size = 1.7,
    alpha = 0.9
  ) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_color_manual(
    values = mcol,
    labels = c("Free/Transition", "Transition/Constant")
    ) +
  labs(
    x = expression("t"[sigma] ~ " [log"[10] ~ "]"),
    y = expression("likelihood-ratio [log"[10] ~ "]"),
    col = "Model-LHR"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", size = 0.9),
    axis.line = element_blank(),
    axis.ticks.length = unit(-3, units = "pt")
  )

pdf("lhr.pdf")
show(p_lh)
dev.off()


p_demography <- rmse_csv %>%
  subset(param == "pop") %>%
  subset(infmodel %in% c("Constant", "OneTransition")) %>%
  ggplot(aes(tsigma, mse, col = infmodel)) +
  geom_point(position = position_jitter(width = 0.1), shape = 16, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_color_manual(values = mcol) +
  labs(
    x = expression("t"[sigma] ~ " [log"[10] ~ "]"),
    y = "RMSE",
    col = "Optimization",
    title = "Demography"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", size = 0.9),
    axis.line = element_blank(),
    axis.ticks.length = unit(-3, units = "pt")
  )

p_sigma <- rmse_csv %>%
  subset(param == "sigma") %>%
  subset(infmodel %in% c("Constant", "OneTransition")) %>%
  ggplot(aes(tsigma, mse, col = infmodel)) +
  geom_point(position = position_jitter(width = 0.1), shape = 16, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se) +
  scale_x_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_color_manual(values = mcol) +
  labs(
    x = expression("t"[sigma] ~ " [log"[10] ~ "]"),
    y = "RMSE",
    col = "Optimization",
    title = "Selfing rate"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", size = 0.9),
    axis.line = element_blank(),
    axis.ticks.length = unit(-3, units = "pt")
  )

p_both <- plot_grid(p_demography, p_sigma, ncol=1, align = T, labels = "auto")


pdf("RMSE_demography_and_sigma.pdf", height = 12)
show(p_both)
dev.off()


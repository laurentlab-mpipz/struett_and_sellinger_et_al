

# setwd("~/Desktop/some_data_analysis/tsabc_tables_txt/teSMC/")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
library(yaml)
library(stringr)

mcol <-
  wes_palette("Darjeeling1", 5, type = "continuous")[c(1:3, 5)]

path_to_csv <- "results/mse/table.csv"

rmse_csv <- read.csv(path_to_csv, row.names = 1)

# translate demography to time of transition from config file
sigma_config <- read_yaml("config/config.yaml")$selfing_rates_backward_in_time
rmse_csv$tsigma <- sapply(rmse_csv$demography, function(demography_index) {
  tsigma_string <- sigma_config[[demography_index+1]][[2]][[2]]
  tsigma_string <- str_replace(tsigma_string,"[^0-9.]","" )
  return(as.numeric(tsigma_string))
})


p1 <- rmse_csv %>%
  # subset(mse < 6) %>%
  subset(infmodel %in% c("Constant", "OneTransition")) %>%
  ggplot(aes(demography, mse, col = infmodel)) +
  geom_point(position = position_jitter(width = 0.1), shape = 16, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se) +
  facet_grid(param ~ ., scales = "free_y") +
  scale_color_manual(values = mcol) +
  labs(title = "Parameter inference by GOF")+
  theme(legend.position = "bottom")

p2 <- rmse_csv %>%
  # subset(mse < 6) %>%
  subset(infmodel %in% c("Constant", "OneTransition")) %>%
  ggplot(aes(tsigma, mse, col = infmodel)) +
  # geom_point(position = position_jitter(width = 0.1), shape = 16, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se) +
  facet_grid(param ~ ., scales = "free_y") +
  scale_x_continuous(trans = "log10")+
  scale_color_manual(values = mcol) +
  labs(title = "Parameter inference by GOF")+
  theme(legend.position = "bottom")

# same as p2 but with data points
p3 <- rmse_csv %>%
  # subset(mse < 6) %>%
  subset(infmodel %in% c("Constant", "OneTransition")) %>%
  ggplot(aes(tsigma, mse, col = infmodel)) +
  geom_point(position = position_jitter(width = 0.02), shape = 16, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se) +
  facet_grid(param ~ ., scales = "free_y") +
  scale_x_continuous(trans = "log10")+
  scale_color_manual(values = mcol) +
  labs(title = "Parameter inference by GOF")+
  theme(legend.position = "bottom")

# same as p3 but panels next to each other
p4 <- rmse_csv %>%
  # subset(mse < 6) %>%
  subset(infmodel %in% c("Constant", "OneTransition")) %>%
  ggplot(aes(tsigma, mse, col = infmodel)) +
  geom_point(position = position_jitter(width = 0.02), shape = 16, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se) +
  facet_wrap(.~ param , scales = "free_y", ncol = 2) +
  scale_x_continuous(trans = "log10")+
  scale_color_manual(values = mcol) +
  labs(title = "Parameter inference by GOF")+
  theme(legend.position = "bottom")

plot_grid(p1, p2, ncol = 2, labels = "AUTO", align = T)+
  # theme(aspect.ratio = 0.707)+
  theme()

p2
p3
p4

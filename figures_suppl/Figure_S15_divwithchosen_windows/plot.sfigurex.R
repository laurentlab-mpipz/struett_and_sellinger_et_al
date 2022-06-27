library(tidyverse)
library(cowplot)
library(wesanderson)
theme_set(theme_cowplot())


region_list <- list()
region_list[[1]] <- c(6500000,	7500000)
region_list[[2]] <- c(16000000,	17000000)
region_list[[3]] <- c(6000000,	7000000)
region_list[[4]] <- c(13000000,	14000000)
region_list[[5]] <- c(1500000,	2500000)


plot_list <- list()
plot_counter <- 0
for (i in 1:5) {
  plot_counter <- plot_counter + 1
  plot_list[[plot_counter]] <-
    read.csv("get_div_for_paper/diversity_sliding_window.csv") %>%
    mutate(POS = X0, DIV = X1) %>%
    select(!starts_with("X")) %>%
    subset(chr == i) %>%
    ggplot(aes(POS, DIV)) +
    geom_rect(
      xmin = region_list[[i]][1],
      xmax = region_list[[i]][2],
      ymin = 0,
      ymax = Inf,
      col = "gray50",
      fill = "gray50"
    ) +
    geom_line() +
    scale_x_continuous(
      breaks = seq(0, 30000000, 5000000),
      labels = c("0", "", "10", "", "20", "", "30")
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.012, 0.001),
      labels = c(
        "0",
        "",
        "",
        "0.003",
        "",
        "",
        "0.006",
        "",
        "",
        "0.009",
        "",
        "",
        "0.012"
      )
    ) +
    # facet_grid(chr~.)+
    theme(aspect.ratio = 0.707/1.7)
}


p <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  align = T,
  labels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5")
)

pdf("SFIGURE.pdf", width = 12, useDingbats = F)
show(p)
dev.off()





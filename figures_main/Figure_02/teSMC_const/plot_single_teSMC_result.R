library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)
library(tidyr)

mcol <-
  wes_palette("Darjeeling1", 5, type = "continuous")[c(1:3, 5)]

logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})

df <- readRDS("results/esmc/table.rds")
df$t_sigma <- factor(df$t_sigma, levels = sort(unique(df$t_sigma)))

{
  rdf <- df %>%
    subset(inference_model %in% c("Constant", "OneTransition"))
  
  pop_const = 40000
  
  tsigma1 = 10000
  tsigma2 = 40000
  
  ylim_dem <- c(1e4, 3.5e5)
  mlimits = c(60, NA)
  
  mtheme <-
    theme(aspect.ratio = 1,
          legend.position = "none")
  
  s1 <-  rdf %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t + 1, sigma, col = inference_model, fill = file)) +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma1, 1e6),
        sigma = c(0.99, 0.1, 0.1),
        inference_model = as.factor("true"),
        file = as.factor("true")
      ),
      col = "black",
      size = 1
    ) +
    geom_step(alpha = 0.6) +
    scale_x_continuous(
      trans = "log",
      limits = mlimits,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    scale_y_continuous(limits = c(0, 1)) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    theme(legend.position = "none") +
    labs(x = "time [log10]",
         y = expression(sigma))
  
  d1 <- rdf %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t + 1, pop_size, col = inference_model, fill = file)) +
    geom_step(alpha = 0.6) +
    geom_vline(xintercept = tsigma1) +
    geom_hline(yintercept = pop_const) +
    scale_x_continuous(
      trans = "log",
      limits = mlimits,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_y_continuous(
      trans = "log",
      limits = ylim_dem,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    labs(x = "time [log10]",
         y = "N")
  
  
  s2 <-  rdf %>%
    subset(t_sigma == tsigma2) %>%
    ggplot(aes(t + 1, sigma, col = inference_model, fill = file)) +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma2, 1e6),
        sigma = c(0.99, 0.1, 0.1),
        inference_model = as.factor("true"),
        file = as.factor("true")
      ),
      col = "black",
      size = 1
    ) +
    geom_step(alpha = 0.6) +
    # geom_vline(aes(xintercept = tsigma2)) +
    scale_x_continuous(
      trans = "log",
      limits = mlimits,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = mcol) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    theme(legend.position = "none") +
    labs(x = "time [log10]",
         y = expression(sigma))
  
  
  d2 <- rdf %>%
    subset(t_sigma == tsigma2) %>%
    ggplot(aes(t + 1, pop_size, col = inference_model, fill = file)) +
    geom_step(alpha = 0.6) +
    geom_vline(xintercept = tsigma2) +
    geom_hline(yintercept = pop_const) +
    scale_x_continuous(
      trans = "log",
      limits = mlimits,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_y_continuous(
      trans = "log",
      limits = ylim_dem,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    labs(x = "time [log10]",
         y = "N")
  
  
  p <- plot_grid(s1,
                 s2,
                 d1,
                 d2,
                 ncol = 2,
                 align = T,
                 labels = "auto")
  
}

pdf("examplatory_inferred_dem_sigma_fun.pdf")
show(p)
dev.off()

if (FALSE) {
  {
  rdf <- df %>%
    subset(inference_model %in% c("Constant", "OneTransition"))
  
  tsigma1 = 3e4
  
  ylim_dem <- c(1e4, 3.5e5)
  mlimits = c(60, NA)
  
  mtheme <-
    theme(aspect.ratio = 1,
          legend.position = "none")
  
  s_ <-  rdf %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t + 1, sigma, col = inference_model, fill = file)) +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma1, 1e6),
        sigma = c(0.99, 0.1, 0.1),
        inference_model = as.factor("true"),
        file = as.factor("true")
      ),
      col = "black",
      size = 1
    ) +
    geom_step(alpha = 0.6) +
    scale_x_continuous(
      trans = "log",
      limits = mlimits,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    scale_y_continuous(limits = c(0, 1)) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    theme(legend.position = "none") +
    labs(x = "time [log10]",
         y = expression(sigma))
  
  d_ <- rdf %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t + 1, pop_size, col = inference_model, fill = file)) +
    geom_step(alpha = 0.6) +
    geom_vline(xintercept = tsigma1) +
    scale_x_continuous(
      trans = "log",
      limits = mlimits,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_y_continuous(
      trans = "log",
      limits = ylim_dem,
      breaks = logbreak,
      labels = loglabel
    ) +
    scale_color_manual(values = mcol) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    labs(x = "time [log10]",
         y = "N")
  
  p_ <- plot_grid(s_, d_,
                 ncol = 2,
                 align = T,
                 labels = "auto")
  
  show(p_)
  }
}

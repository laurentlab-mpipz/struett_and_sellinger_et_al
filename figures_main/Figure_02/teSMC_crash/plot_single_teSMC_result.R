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
df$t_pop <- sapply(df$demographic_scenario, function(x) {
  if (x <= 21) {
    tpop <- 10010
  } else if (x >= 22) {
    tpop <- 40010
  }
  
  return(tpop)
})
message("the demography is devided by pop size change time, the times are hardcoded")

# first t_pop
df_total <- df
df1 <- df_total %>% subset(t_pop == 10010)
df2 <- df_total %>% subset(t_pop == 40010)

rdf_total <- df_total %>%
  subset(inference_model %in% c("Constant", "OneTransition"))

{
  df <- df1
  rdf <- df %>%
    subset(inference_model %in% c("Constant", "OneTransition"))
  
  pop_1 = 200000
  pop_2 = 40000
  tpop = df$t_pop %>% unique()
  
  tsigma1 = 10000
  tsigma2 = 40000
  
  ylim_dem <- c(2e3, 3e6)
  mlimits = c(60, NA)
  
  mtheme <-
    theme(aspect.ratio = 1,
          legend.position = "none")
  
  s1 <-  rdf %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t + 1, sigma, col = inference_model, fill = file)) +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma1, Inf),
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
    geom_vline(xintercept = tsigma1, col="gray50") +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tpop, Inf),
        pop_size = c(pop_1, pop_2, pop_2),
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
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma2, Inf),
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
    geom_vline(xintercept = tsigma2, col="gray50") +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tpop, Inf),
        pop_size = c(pop_1, pop_2, pop_2),
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
  
  p_fin1 <- plot_grid(s1,
                 s2,
                 d1,
                 d2,
                 ncol = 2,
                 align = T,
                 labels = "auto")
}
show(p_fin1)

{
  df <- df2
  rdf <- df %>%
    subset(inference_model %in% c("Constant", "OneTransition"))
  
  pop_1 = 200000
  pop_2 = 40000
  tpop = df$t_pop %>% unique()
  
  tsigma1 = 10000
  tsigma2 = 40000
  
  ylim_dem <- c(2e3, 3e6)
  mlimits = c(60, NA)
  
  mtheme <-
    theme(aspect.ratio = 1,
          legend.position = "none")
  
  s1 <-  rdf %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t + 1, sigma, col = inference_model, fill = file)) +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma1, Inf),
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
    geom_vline(xintercept = tsigma1, col="gray50") +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tpop, Inf),
        pop_size = c(pop_1, pop_2, pop_2),
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
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tsigma2, Inf),
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
    geom_vline(xintercept = tsigma2, col="gray50") +
    geom_step(
      data = data.frame(
        t = c(min(rdf$t[rdf$t >= mlimits[1]]), tpop, Inf),
        pop_size = c(pop_1, pop_2, pop_2),
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
  
  p_fin2 <- plot_grid(s1,
                      s2,
                      d1,
                      d2,
                      ncol = 2,
                      align = T,
                      labels = "auto")
}
show(p_fin2)

pdf("examplatory_inferred_dem_sigma_fun_tpop0.5.pdf")
show(p_fin1)
dev.off()

pdf("examplatory_inferred_dem_sigma_fun_tpop2.0.pdf")
show(p_fin2)
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

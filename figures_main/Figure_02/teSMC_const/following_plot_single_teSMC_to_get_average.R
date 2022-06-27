library(tidyverse)
library(colorspace)

source("plot_single_teSMC_result.R")

complete_df <- rdf %>% subset(t_sigma %in% c(tsigma1, tsigma2))

# get all time points
{
  full_df_list <- list()
  full_df_index <- 0
  uni_times <- unique(complete_df$t) %>% sort()
  for (tsigma in unique(complete_df$t_sigma)) {
    ddf <- complete_df %>%
      subset(t_sigma == tsigma)
    
    for (m in unique(ddf$inference_model)) {
      subdf <- ddf %>%
        subset(inference_model == m)
      
      for (f in unique(ddf$file)) {
        full_df_index <- full_df_index + 1
        subdf <- ddf %>%
          subset(file == f)
        
        this_subdf_list = list()
        list_index = 0
        
        cat(full_df_index, "of",
            length(unique(ddf$file)) * length(unique(ddf$inference_model)) * length(unique(complete_df$t_sigma)),
            m, dim(subdf),
            "\t\t\t\t\r")
        
        for (this_t in uni_times) {
          list_index = list_index + 1
          which_interval <-
            findInterval(this_t, subdf$t, all.inside = T)
          this_time_row <- subdf[which_interval, ]
          this_time_row$t <- this_t
          this_subdf_list[[list_index]] <- this_time_row
        }
        
        this_subdf <- do.call("rbind.data.frame", this_subdf_list)
        full_df_list[[full_df_index]] <- this_subdf
      }
    }
    full_df <- do.call("rbind.data.frame", full_df_list)
  }
}

# get averages by group
stat_df <- full_df %>%
  group_by(inference_model, t_sigma, t, demography) %>%
  summarize(
    mean_sigma = mean(sigma),
    mean_pop_size = mean(pop_size),
    sd_sigma = sd(sigma),
    sd_pop_size = sd(pop_size),
    quant025_pop_size = quantile(pop_size, 0.025),
    quant975_pop_size = quantile(pop_size, 0.975),
    quant025_sigma = quantile(sigma, 0.025),
    quant975_sigma = quantile(sigma, 0.975)
  )

p_new <- stat_df %>%
  subset(t_sigma == tsigma1) %>%
  ggplot(aes(t+1, mean_pop_size, col=inference_model, fill=inference_model,
             xmin=t+1,
             xmax=lead(t)+1,
             ymin=quant025_pop_size,
             ymax=quant975_pop_size)) +
  geom_rect(
    alpha = 0.5,
    col = NA
  ) +
  geom_step(size=1) +
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
  scale_fill_manual(values = lighten(mcol, amount = 0.85)) +
  # facet_grid(t_sigma ~ .) +
  # mtheme +
  theme(aspect.ratio = 0.3,
        legend.position = "none") +
  labs(x = "time [log10]",
       y = "N")

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
  
  d1 <- stat_df %>%
    subset(t_sigma == tsigma1) %>%
    ggplot(aes(t+1, mean_pop_size, col=inference_model, fill=inference_model,
               xmin=t+1,
               xmax=lead(t)+1,
               ymin=quant025_pop_size,
               ymax=quant975_pop_size)) +
    geom_rect(
      alpha = 0.4,
      col = NA
    ) +
    geom_step(size=1) +
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
    scale_fill_manual(values = lighten(mcol, amount = 0.8)) +
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
  
  
  d2 <- stat_df %>%
    subset(t_sigma == tsigma2) %>%
    ggplot(aes(t+1, mean_pop_size, col=inference_model, fill=inference_model,
               xmin=t+1,
               xmax=lead(t)+1,
               ymin=quant025_pop_size,
               ymax=quant975_pop_size)) +
    geom_rect(
      alpha = 0.4,
      col = NA
    ) +
    geom_step(size=1) +
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
    scale_fill_manual(values = lighten(mcol, amount = 0.8)) +
    # facet_grid(t_sigma ~ .) +
    mtheme +
    labs(x = "time [log10]",
         y = "N")
  
  
  p <- plot_grid(d1,
                 d2,
                 s1,
                 s2,
                 ncol = 2,
                 align = T,
                 labels = "auto")
  
  }


pdf("examplatory_inferred_dem_sigma_fun_with_95conf.pdf")
show(p)
dev.off()

row1 <- plot_grid(d1+theme(aspect.ratio = 0.707/2,
                           axis.title.x = element_blank()), labels = c("a"))
row2 <- plot_grid(d2+theme(aspect.ratio = 0.707/2,
                           axis.title.x = element_blank()), labels = c("b"))
row3 <- plot_grid(s1, s2, nrow=1, labels = c("c", "d"))
p_alt <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1, 1, 1.1)
)

pdf("examplatory_inferred_dem_sigma_fun_with_95conf_alt.pdf", width = 7)
show(p_alt)
dev.off()


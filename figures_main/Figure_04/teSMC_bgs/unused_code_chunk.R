bind two lh methods together
df <- df_bw %>%
  pivot_longer(-c(rep, fname),
               names_to = "time_window",
               values_to = "xi") %>%
  separate(
    col = "fname",
    sep = "_+",
    into = c(
      "torm1",
      "torm2",
      "dtype",
      "optim_method",
      "optim_model",
      "torm3",
      "torm4",
      "torm5",
      "torm6",
      "dfe",
      "torm7",
      "tsigma",
      "torm8",
      "torm9",
      "masked",
      "torm10"
    )
  )

# remove no content columns
df[, grep("torm", names(df))] <- NULL

# correct mask information
df$masked <- df$masked == "masked"

# find the correct rescaled time based on the time window
rescaled_time_windows <- df %>%
  subset(dtype == "t")  # subset for time windows
df <- df %>%
  subset(dtype != "t")  # subset for population size and selfing


rescaled_time <- pbapply(df, 1, function(x) {
  rep_ <- x[["rep"]]
  dtype_ <- x[["dtype"]]
  optim_method_ <- x[["optim_method"]]
  optim_model_ <- x[["optim_model"]]
  dfe_ <- x[["dfe"]]
  tsigma_ <- x[["tsigma"]]
  masked_ <- x[["masked"]]
  time_window_ <- x[["time_window"]]
  xi_ <- x[["xi"]]
  
  rescaled_time <- (
    rescaled_time_windows %>%
      subset(rep == as.numeric(rep_)) %>%
      subset(optim_method == optim_method_) %>%
      subset(optim_model == optim_model_) %>%
      subset(dfe == dfe_) %>%
      subset(tsigma == tsigma_) %>%
      subset(masked == masked_) %>%
      subset(time_window == time_window_)
  )[["xi"]]
  stopifnot(length(rescaled_time) == 1)
  
  return(rescaled_time)
})
df$rescaled_time <- rescaled_time

# safety copy
if (!exists("df_safetycopy")) {
  df_safetycopy <- df
} else {
  if (FALSE) {
    # only run manually if at all
    df <- df_safetycopy  # only run manually if at all
  }
}


# reformat data table
df$time_window <- as.numeric(df$time_window)
names(df)[names(df) == "xi"] <- "value"
df$value[df$dtype == "p"] <- 10**(df$value[df$dtype == "p"])  # Thibaut always talks about log10 values
df$optim_method <- factor(df$optim_method, levels = c("LH", "BW"))
df$tsigma <- as.numeric(df$tsigma)

df_for_abc_combin <- df

# plot value setup
logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})
msbreaks <- seq(0, 1, 0.1)
mslabels <- sapply(msbreaks, function(x) {
  if (x %in% c(0, 1))
    return(as.character(x))
  else
    return("")
})
mcol <- wes_palette("FantasticFox1")
xlim <- c(1e4, NA)
ylim <- c(NA, NA)
step_alpha <- 1
vertical.lines <- df$rescaled_time %>% unique() %>% sort()


# plot into 1 plot
{
  
  subdf <- df %>%
    subset(optim_model %in% c("OT"))
  
  pp <- subdf %>%
    subset(dtype == "p") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      y = value,
      fill = as.factor(rep),
      col = dfe,
      linetype = optim_model
    )) +
    geom_step(alpha = step_alpha, size=1) +
    facet_grid(tsigma ~ masked, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = ylim
    ) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y="population size [log10]") +
    theme(
      # legend.position = "none"
    )
  
  sp <- subdf %>%
    subset(dtype == "s") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      y = value,
      fill = as.factor(rep),
      col = dfe,
      linetype = optim_model
    )) +
    geom_step(alpha = step_alpha, size=1) +
    facet_grid(tsigma ~ masked, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = msbreaks, 
      labels = mslabels) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y = "selfing rate") +
    theme(
      # legend.position = "none"
    )
  
  p <- plot_grid(pp, sp, ncol = 1, align = "T") +
    theme(aspect.ratio = 1)
}
show(p)

pdf("ara_results_teSMC.pdf", width = 7, height = 7)
show(p)
dev.off()

png("ara_results_teSMC.png", width = 7, height = 7, units = "in", res = 300)
show(p)
dev.off()


# manually plot abc results
message("abc results are manually provided, Relicts are still missing")
abcpop <- cbind.data.frame(
  value = c(315243, 641509, 641509, 307288, 466664, 466664),
  t = c(0, 394741, Inf, 0, 460915, Inf),
  pop = c("CEU", "CEU", "CEU", "IBnr", "IBnr", "IBnr"),
  dtype = "p"
)

abcsigma <- cbind.data.frame(
  value = c(0.995, 0.102, 0.102, 0.993, 0.078, 0.078),
  t = c(0, 552920, Inf, 0, 607514, Inf),
  pop = c("CEU", "CEU", "CEU", "IBnr", "IBnr", "IBnr"),
  dtype = "s"
)

add_times <- function(df, n=1e4) {
  df_pop_list <- list()
  pop_ix <- 1
  for (my_p in unique(df$pop)) {
    sub <- df %>% subset(pop == my_p)
    
    # for each time point add the value
    x_time=seq(min(sub$t), maxInf(max(sub$t)),
               (maxInf(max(sub$t)) - min(sub$t))/n)
    
    values <- numeric(length(x_time))
    for (my_t_ix in 1:length(x_time)) {
      my_t <- x_time[my_t_ix]
      values[my_t_ix] <- sub$value[max(which(sub$t <= my_t))]
    }
    
    sub <- rbind.data.frame(
      sub,
      cbind.data.frame(
        value = values,
        t = x_time,
        pop = unique(sub$pop),
        dtype = unique(sub$dtype)
      )
    )
    df_pop_list[[pop_ix]] <- sub
    pop_ix <- pop_ix + 1
  }
  
  df <- do.call(rbind.data.frame, df_pop_list)
  return(df)
}

abcpop <- add_times(abcpop)
abcsigma <- add_times(abcsigma)

abc <- rbind.data.frame(abcpop, abcsigma)

xlim <- c(1e4, NA)
ylim <- c(7e4, 1e6)

{
  
  subdf <- abc
  
  pp <- subdf %>%
    subset(dtype == "p") %>%
    ggplot(aes(
      x = t + 1,
      y = value,
      col = pop
    )) +
    geom_step(alpha = step_alpha, size=1) +
    # facet_grid(. ~ optim_model, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = ylim
    ) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y="population size [log10]") +
    theme(
      # legend.position = "none"
    )
  
  sp <- subdf %>%
    subset(dtype == "s") %>%
    ggplot(aes(
      x = t + 1,
      y = value,
      col = pop
    )) +
    geom_step(alpha = step_alpha, size=1) +
    # facet_grid(. ~ optim_model, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = msbreaks, 
      labels = mslabels) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y = "selfing rate") +
    theme(
      # legend.position = "none"
    )
  
  p <- plot_grid(pp, sp, ncol = 1, align = "T") +
    theme(aspect.ratio = 0.707)
}

show(p)

pdf("ara_results_tsABC.pdf", width = 7, height = 7)
show(p)
dev.off()

png("ara_results_tsABC.png", width = 7, height = 7, units = "in", res = 300)
show(p)
dev.off()


# combine two methods into 1 plot
head(df_for_abc_combin); names(df_for_abc_combin)

head(abc)
names(abc) <- c("value", "rescaled_time", "pop", "dtype")
abc$rep = 1
abc$optim_method = "abc"
abc$optim_model = "6par"
abc$time_window = 1

df <- rbind.data.frame(abc, df_for_abc_combin)

xlim <- c(1e4, 1e6)
ylim <- c(NA, NA)


# plot into 1 plot
{
  
  subdf <- df %>%
    subset(optim_model %in% c("OT", "6par"))
  
  subdf$optim_model[subdf$optim_model == "OT"] <- "teSMC"
  subdf$optim_model[subdf$optim_model == "6par"] <- "tsABC"
  
  subdf$optim_model <- factor(subdf$optim_model, levels = c("tsABC", "teSMC"))
  
  pp <- subdf %>%
    subset(dtype == "p") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      y = value,
      fill = as.factor(rep),
      col = pop,
      linetype = optim_model
    )) +
    geom_step(alpha = step_alpha, size=1) +
    # facet_grid(. ~ optim_model, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = ylim
    ) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y="population size [log10]") +
    theme(
      # legend.position = "none"
    )
  
  sp <- subdf %>%
    subset(dtype == "s") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      xintercept = rescaled_time + 1,
      y = value,
      fill = as.factor(rep),
      col = pop,
      linetype = optim_model
    )) +
    geom_vline(xintercept = vertical.lines, col="gray95", size=0.06) +
    geom_step(alpha = step_alpha, size=1) +
    # facet_grid(. ~ optim_model, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = msbreaks, 
      labels = mslabels) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y = "selfing rate") +
    theme(
      # legend.position = "none"
    )
  
  p <- plot_grid(pp, sp, ncol = 1, align = "T") +
    theme(aspect.ratio = 0.707)
}

pdf("ara_results_both_methods.pdf", width = 7, height = 7)
show(p)
dev.off()

png("ara_results_both_methods.pdf.png", width = 7, height = 7, units = "in", res = 300)
show(p)
dev.off()
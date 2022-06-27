library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)
library(pbapply)


# return mio if INF
maxInf <- function(x) {
  if (is.infinite(x)) {
    value <- 1e6
  } else {
    value <- x
  }
  return(value)
}


fnames_bw <-
  list.files(path = "Arabidopsis_results_good_window",
             pattern = "*.csv",
             full.names = T)

csvFiles2df <- function(flist) {
  df <- do.call(rbind.data.frame, lapply(flist, function(x) {
    t <- read.csv(x)
    names(t) <- c("rep", 1:(ncol(t) - 1))
    t$fname <- strsplit(x, "/+", perl = TRUE)[[1]][2]
    return(t)
  }))
  
  return(df)
}

# read files with added filenames
df_bw <- csvFiles2df(fnames_bw)

# bind two lh methods together
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
      "pop",
      "torm6"
    )
  )

# remove no content columns
df[, grep("torm", names(df))] <- NULL

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
  pop_ <- x[["pop"]]
  time_window_ <- x[["time_window"]]
  xi_ <- x[["xi"]]
  
  rescaled_time <- (
    rescaled_time_windows %>%
      subset(rep == as.numeric(rep_)) %>%
      subset(optim_method == optim_method_) %>%
      subset(optim_model == optim_model_) %>%
      subset(pop == pop_) %>%
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
# show(p)

pdf("ara_results_teSMC.pdf", width = 7, height = 7)
show(p)
dev.off()

png("ara_results_teSMC.png", width = 7, height = 7, units = "in", res = 300)
show(p)
dev.off()


# manually plot abc results
message("abc results are manually provided, Relicts are still missing")
abcpop <- cbind.data.frame(
  value = c(340366, 462241.1, 462241.1, 268730.7, 517442.4, 517442.4, 393592.9, 384689.9, 384689.9),
  t = c(0, 573224.4, Inf, 0, 498924.1, Inf, 0, 786032.6, Inf),
  pop = c("CEU", "CEU", "CEU", "IBnr", "IBnr", "IBnr", "Relict", "Relict", "Relict"),
  dtype = "p"
)

abcsigma <- cbind.data.frame(
  value = c(0.9959358, 0.1006572, 0.1006572, 0.9933578, 0.09106993, 0.09106993, 0.9970796, 0.08647879, 0.08647879),
  t = c(0, 707995.2, Inf, 0, 756976, Inf, 0, 592321.1, Inf),
  pop = c("CEU", "CEU", "CEU", "IBnr", "IBnr", "IBnr", "Relict", "Relict", "Relict"),
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


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
library(yaml)
library(stringr)
library(pbapply)

outfile <- "MAN_FIG_2.pdf"

this_working_dir <- paste0(
  "/Users/struett/Dropbox/phd/manuscripts/transitioning_to_selfing/",
  "figure_finals_main/Figure_02/good_window"
)


setwd(this_working_dir)

fnames <-
  list.files(path = "Figure_2_Good_window/",
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


df <- csvFiles2df(fnames)

df <- df %>%
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
      "scenario",
      "torm5",
      "pod",
      "torm7"
    )
  )


# remove no content columns
df[, grep("torm", names(df))] <- NULL


# get tsgima true
sigma_config_const <-
  read_yaml("config_constant.yaml")$selfing_rates_backward_in_time
sigma_config_crash <-
  read_yaml("config_crash.yaml")$selfing_rates_backward_in_time
sigma_config_expansion <-
  read_yaml("config_expansion.yaml")$selfing_rates_backward_in_time


df$tsigma <-
  pbapply(df, 1, function(x) {
    demography_index <- as.numeric(x["pod"])
    scenario <- x["scenario"]
    
    if (scenario == "constant") {
      tsigma_string <- sigma_config_const[[demography_index + 1]][[2]][[2]]
      tsigma_string <- str_replace(tsigma_string, "[^0-9.]", "")
    } else if (scenario == "crash") {
      tsigma_string <- sigma_config_crash[[demography_index + 1]][[2]][[2]]
      tsigma_string <- str_replace(tsigma_string, "[^0-9.]", "")
    } else if (scenario == "expansion") {
      tsigma_string <- sigma_config_expansion[[demography_index + 1]][[2]][[2]]
      tsigma_string <- str_replace(tsigma_string, "[^0-9.]", "")
    } else {
      stop("unknown demography")
    }

    return(as.numeric(tsigma_string))
  })


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
  scenario_ <- x[["scenario"]]
  pod_ <- x[["pod"]]
  time_window_ <- x[["time_window"]]
  xi_ <- x[["xi"]]
  
  rescaled_time <- (
    rescaled_time_windows %>%
      subset(rep == as.numeric(rep_)) %>%
      subset(optim_method == optim_method_) %>%
      subset(optim_model == optim_model_) %>%
      subset(scenario == scenario_) %>%
      subset(pod == pod_) %>%
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
mcol <- c(
  wes_palette("Chevalier1")[1],
  "gray"
)
xlim <- c(1e4, NA)
ylim <- c(NA, NA)
step_alpha <- 1
vertical.lines <- df$rescaled_time %>% unique() %>% sort()

df$pod <- as.numeric(df$pod)

# get tsigma_estim
tsigma_estim <- pbapply(df, 1, function(x) {
  rep_ <- as.numeric(x[["rep"]])
  dtype_ <- x[["dtype"]]
  optim_method_ <- x[["optim_method"]]
  optim_model_ <- x[["optim_model"]]
  scenario_ <- x[["scenario"]]
  pod_ <- as.numeric(x[["pod"]])
  tsigma_ <- as.numeric(x[["tsigma"]])
  
  df_sub <- df %>%
    subset(rep == as.numeric(rep_)) %>%
    subset(dtype == "s") %>%  # neeeds to be s always
    subset(optim_method == optim_method_) %>%
    subset(optim_model == optim_model_) %>%
    subset(scenario == scenario_) %>%
    subset(pod == pod_) %>%
    subset(tsigma == tsigma_)
  stopifnot(nrow(df_sub) == 40)
  
  tsigma_estim_dem <- df_sub$value
  tsigma_estim_time <- df_sub$rescaled_time
  
  tsigma_estim <- tsigma_estim_time[min(which(tsigma_estim_dem <= 0.5))]
  
  return(tsigma_estim)
})

df$tsigma_estim <- tsigma_estim

range(df$tsigma_estim, na.rm = T)
mlimits <- c(6e2, 3e5)

# safety copy 2
if (!exists("df_safetycopy2")) {
  df_safetycopy2 <- df
} else {
  if (FALSE) {
    # only run manually if at all
    df <- df_safetycopy2  # only run manually if at all
  }
}


# reduce data set to tsigma
dfreduced <- df[,c("rep", "dtype", "optim_method", "optim_model", "scenario", "pod", "tsigma", "tsigma_estim")]
dfreduced <- dfreduced[!duplicated(dfreduced),]
df <- dfreduced

{
  p_tsigma_constant <- df %>%
    subset(scenario == "constant") %>%
    ggplot(aes(tsigma, tsigma_estim, col=optim_model)) +
    geom_abline() +
    geom_point(position = position_jitter(width = 0.02,height = 0),
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
  
  message("Pop size changes are hard coded for the plot")
  df$t_pop <- sapply(df$pod, function(x) {
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
  
  p_tsigma_crash_tNe10k <- df1 %>%
    subset(scenario == "crash") %>%
    ggplot(aes(tsigma, tsigma_estim, col=optim_model)) +
    geom_vline(xintercept = tpop1, col = "gray") +
    geom_abline() +
    geom_point(position = position_jitter(width = 0.02,height = 0),
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
  
  show(p_tsigma_crash_tNe10k)
  
  
  p_tsigma_crash_tNe40k <- df2 %>%
    subset(scenario == "crash") %>%
    ggplot(aes(tsigma, tsigma_estim, col=optim_model)) +
    geom_vline(xintercept = tpop2, col = "gray") +
    geom_abline() +
    geom_point(position = position_jitter(width = 0.02,height = 0),
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
  
  show(p_tsigma_crash_tNe40k)
  
  message("Pop size changes are hard coded for the plot")
  df$t_pop <- sapply(df$pod, function(x) {
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
    subset(scenario == "expansion") %>%
    ggplot(aes(tsigma, tsigma_estim, col=optim_model)) +
    geom_vline(xintercept = tpop1, col = "gray") +
    geom_abline() +
    geom_point(position = position_jitter(width = 0.02,height = 0),
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
    subset(scenario == "expansion") %>%
    ggplot(aes(tsigma, tsigma_estim, col=optim_model)) +
    geom_vline(xintercept = tpop2, col = "gray") +
    geom_abline() +
    geom_point(position = position_jitter(width = 0.02,height = 0),
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
  
  
}


plot_list <- list(
  p_tsigma_constant,
  p_tsigma_expansion_tNe10k,
  p_tsigma_expansion_tNe40k,
  p_tsigma_crash_tNe10k,
  p_tsigma_crash_tNe40k
)

linesize <- 0.5


plot_list <- lapply(plot_list, function(x) {
  return(
    x + theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.border = element_rect(color = "black", size = linesize),
      axis.line = element_blank(),
      axis.ticks = element_line(size = linesize * 0.8),
      axis.ticks.length = unit(0.015, "in"),
      plot.margin = margin(8, 8, 8, 8, unit = "pt"),
      legend.position = "none"
    )
  )
})

### multipanel
p <- plot_grid(
  plot_list[[1]],
  NULL,
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  plot_list[[5]],
  ncol = 2,
  # labels = c("A", "", "B", "C", "D", "E"),
  align = TRUE
)

mwidth <- 4.4
mheight <- 6.8

pdf(outfile,
    useDingbats = FALSE,
    height = mheight,
    width = mwidth)
show(p)
dev.off()

png(
  paste0(outfile, ".png"),
  width = mwidth,
  height = mheight,
  res = 600,
  units = "in"
)
show(p)
dev.off()

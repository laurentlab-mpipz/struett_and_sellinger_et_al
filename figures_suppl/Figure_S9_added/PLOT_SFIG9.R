library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
library(yaml)
library(stringr)
library(pbapply)

outfile <- "MAN_FIG_S9.pdf"

this_working_dir <- paste0(
  "~/Dropbox/phd/manuscripts/transitioning_to_selfing/",
  "figure_finals_suppl/Figure_S9_added/"
)

setwd(this_working_dir)



# plot specs
mcol <- c(wes_palette("Chevalier1")[2],
          wes_palette("Chevalier1")[1])
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
mlimits <- c(6e2, 3e5)


fnames <-
  list.files(path = "Map_teSMC/",
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

df <- csvFiles2df(fnames) %>%
  na.omit() %>%
  pivot_longer(-c(rep, fname),
               names_to = "time_window",
               values_to = "xi") %>%
  separate(
    col = "fname",
    sep = "_+|\\.",
    into = c(
      "torm1",
      "torm2",
      "torm3",
      "torm4",
      "dtype",
      "torm5",
      "torm6",
      "torm7",
      "torm8"
    )
  )


mint <- 1e2
maxt <- 5e4
# create table for true
df_true <- rbind.data.frame(
  data.frame(
    rep = NA,
    dtype = "p",
    time_window = NA,
    value = c(5, 5),
    rescaled_time = c(mint, maxt)
  ),
  data.frame(
    rep = NA,
    dtype = "s",
    time_window = NA,
    value = c(0.9, 0.9),
    rescaled_time = c(mint, maxt)
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
  time_window_ <- x[["time_window"]]
  xi_ <- x[["xi"]]
  rescaled_time <- (
    rescaled_time_windows %>%
      subset(rep == as.numeric(rep_)) %>%
      # subset(scenario == scenario_) %>%
      # subset(dtype == dtype_) %>%
      subset(time_window == as.numeric(time_window_))
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
df$value[df$dtype == "p"] <-
  10 ** (df$value[df$dtype == "p"])  # Thibaut always talks about log10 values


# add true values
df$istrue <- F
df_true$istrue <- T

df <- rbind.data.frame(df,
                       df_true)

df$istrue <- factor(df$istrue, levels = c(T, F))


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
mcol <- c("gray",
          wes_palette("Chevalier1")[1])
xlim <- c(5e3, NA)
ylim <- c(NA, NA)
step_alpha <- 1
vertical.lines <-
  (df %>% subset(istrue == F))$rescaled_time  %>% unique() %>% sort()


{
  
  vertical.lines <- df$rescaled_time %>% unique() %>% sort()
  
  pp <- df %>%
    subset(dtype == "p") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      y = value,
      col = as.factor(rep),
      fill = rep,
      size = istrue
    )) +
    geom_step(alpha = step_alpha, size = 1) +
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
    # scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y = "population size [log10]") +
    theme(legend.position = "none")
  
  sp <- df %>%
    subset(dtype == "s") %>%
    ggplot(
      aes(
        x = rescaled_time + 1,
        xintercept = rescaled_time + 1,
        y = value,
        fill = rep,
        col = as.factor(rep)
      )
    ) +
    geom_vline(xintercept = vertical.lines,
               col = "gray95",
               size = 0.06) +
    geom_step(alpha = step_alpha, size = 1) +
    # facet_grid(. ~ optim_model, scales = "free_y") +
    scale_x_continuous(
      trans = "log",
      breaks = logbreak,
      labels = loglabel,
      limits = xlim
    ) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = msbreaks,
                       labels = mslabels) +
    # scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y = "selfing rate") +
    theme(legend.position = "none")
  
  p <- plot_grid(sp, pp, ncol = 1, align = "T") +
    theme(aspect.ratio = 2 * 0.707)
  
  show(p)
  
  }

my_final_figure <- p

pdf(outfile,
    width = 8.8,
    height = 9,
    useDingbats = F)
show(my_final_figure)
dev.off()

png(
  paste0(outfile, ".png"),
  width = 8.8,
  height = 9,
  units = "in",
  res = 600
)
show(my_final_figure)
dev.off()

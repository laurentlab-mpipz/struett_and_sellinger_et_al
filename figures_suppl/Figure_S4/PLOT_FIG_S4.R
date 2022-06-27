library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(wesanderson)
library(yaml)
library(stringr)
library(pbapply)

outfile <- "MAN_FIG_S4.pdf"

this_working_dir <- paste0(
  "/Users/struett/Dropbox/phd/manuscripts/transitioning_to_selfing/",
  "figure_finals_suppl/Figure_S4/"
)

setwd(this_working_dir)

# plot specs
mcol <- c(
  wes_palette("Chevalier1")[2],
  wes_palette("Chevalier1")[1]
)
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
  list.files(path = "Best_case_convergence_teSMC",
             pattern = "*.csv",
             full.names = T)

fnames_rec <- fnames[grep("rec", fnames)]
fnames_sigma <- fnames[grep("selfing", fnames)]

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
df_sigma <- csvFiles2df(fnames_sigma) %>% na.omit()
df_rec <- csvFiles2df(fnames_rec) %>% na.omit()


# bind two methods together
df <- rbind.data.frame(
  df_sigma %>%
    pivot_longer(-c(rep, fname),
                 names_to = "time_window",
                 values_to = "xi") %>%
    separate(
      col = "fname",
      sep = "_+",
      into = c(
        "torm1",
        "scenario",
        "torm2",
        "torm3",
        "torm4",
        "torm5",
        "torm6",
        "dtype",
        "torm7",
        "torm8",
        "torm9",
        "torm10"
      )
    ),
  df_rec %>%
    pivot_longer(-c(rep, fname),
                 names_to = "time_window",
                 values_to = "xi") %>%
    separate(
      col = "fname",
      sep = "_+",
      into = c(
        "torm1",
        "scenario",
        "torm2",
        "torm3",
        "torm4",
        "torm5",
        "torm6",
        "dtype",
        "torm7",
        "torm8",
        "torm9",
        "torm10"
      )
    )
)


mint <- 1e2
maxt <- 5e4
# create table for true
df_true <- rbind.data.frame(
  data.frame(
    rep = NA,
    scenario = "rec",
    dtype = "r",
    time_window = NA,
    value = c(8, 2, 2),
    rescaled_time = c(mint, 5000,  maxt)
  ),
  data.frame(
    rep = NA,
    scenario = "rec",
    dtype = "p",
    time_window = NA,
    value = c(1e4, 1e3, 1e5, 1e5),
    rescaled_time = c(mint, 1000 , 10000, maxt)
  ),
  data.frame(
    rep = NA,
    scenario = "selfing",
    dtype = "s",
    time_window = NA,
    value = c(.8, .2, .2),
    rescaled_time = c(mint, 5000,  maxt)
  ),
  data.frame(
    rep = NA,
    scenario = "selfing",
    dtype = "p",
    time_window = NA,
    value = c(1e4, 1e3, 1e5, 1e5),
    rescaled_time = c(mint, 1000 , 10000, maxt)
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
  scenario_ <- x[["scenario"]]
  dtype_ <- x[["dtype"]]
  time_window_ <- x[["time_window"]]
  xi_ <- x[["xi"]]
  rescaled_time <- (
    rescaled_time_windows %>%
      subset(rep == as.numeric(rep_)) %>%
      subset(scenario == scenario_) %>%
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
df$value[df$dtype == "p"] <- 10**(df$value[df$dtype == "p"])  # Thibaut always talks about log10 values


# add true values
df$istrue <- F
df_true$istrue <- T

df <- rbind.data.frame(
  df,
  df_true
)

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
mcol <- c(
  "gray",
  wes_palette("Chevalier1")[1]
)
xlim <- c(1e2, 1e5)
ylim <- c(1e2, 1e6)
step_alpha <- 1
vertical.lines <- (df %>% subset(istrue == F))$rescaled_time  %>% unique() %>% sort()

{
  subdf <- df %>%
    subset(scenario == "selfing")
  
  vertical.lines <- subdf$rescaled_time %>% unique() %>% sort()
  
  pp <- subdf %>%
    subset(dtype == "p") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      y = value,
      col = as.factor(istrue),
      size = istrue
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
      legend.position = "none"
    )
  
  sp <- subdf %>%
    subset(dtype == "s") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      xintercept = rescaled_time + 1,
      y = value,
      col = as.factor(istrue)
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
      legend.position = "none"
    )
  
  p_sigma <- plot_grid(sp, pp, ncol = 1, align = "T") +
    theme(aspect.ratio = 2* 0.707)
  
  show(p_sigma)
  
}



{
  subdf <- df %>%
    subset(scenario == "rec")
  
  vertical.lines <- subdf$rescaled_time %>% unique() %>% sort()
  
  ppr <- subdf %>%
    subset(dtype == "p") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      y = value,
      col = as.factor(istrue)
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
      legend.position = "none"
    )
  
  spr <- subdf %>%
    subset(dtype == "r") %>%
    ggplot(aes(
      x = rescaled_time + 1,
      xintercept = rescaled_time + 1,
      y = value,
      col = as.factor(istrue)
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
      # trans = "log",
      limits = c(0, 10),
      breaks = seq(0, 10, 2),
      # labels = loglabel
      ) +
    scale_color_manual(values = mcol) +
    labs(x = "generations before present [log10]", y = "recombination rate") +
    theme(
      legend.position = "none"
    )
  
  p_rec <- plot_grid(spr, ppr, ncol = 1, align = "T") +
    theme(aspect.ratio = 2* 0.707)
  
  show(p_rec)
  
}


mwidth <- 4.4
mheight <- 6.8

pdf(paste0("selfing.", outfile),
    useDingbats = FALSE,
    height = mheight,
    width = mwidth)
show(p_sigma)
dev.off()

png(
  paste0("selfing.", outfile, ".png"),
  width = mwidth,
  height = mheight,
  res = 600,
  units = "in"
)
show(p_sigma)
dev.off()



pdf(paste0("rec", outfile),
    useDingbats = FALSE,
    height = mheight,
    width = mwidth)
show(p_rec)
dev.off()

png(
  paste0("rec", outfile, ".png"),
  width = mwidth,
  height = mheight,
  res = 600,
  units = "in"
)
show(p_rec)
dev.off()




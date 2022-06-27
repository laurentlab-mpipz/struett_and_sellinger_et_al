
getwd()
# setwd("~/Desktop/tables_txt/")

library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)
library(scales)  # show colors

mcol = wes_palette("Darjeeling1")[1:2]
mcol = c(mcol[1], "gray90", mcol[2])
show_col(mcol)
colfunc<-colorRampPalette(mcol)
show_col(colfunc(9))
mcol = colfunc(9)[c(4:9)]
mcol = colfunc(9)[c(3, 5:9)]
show_col(mcol)

# setwd("/Users/struett/Desktop/tables_txt/")

df <- read.csv("model_choice_pod_aggregation_5.table.csv", as.is = T) # 5 because rho/theta = 5
df$X = NULL
df$path = NULL
df$ident_0 = NULL
df$param_0 = NULL
df$param_1 = NULL
df$param_2 = NULL
df$param_3 = NULL
df$param_4 = NULL
df$rand_0 = NULL
df$marginal_density_csigma = NULL
df$marginal_density_tsigma = NULL

# if not sufficient amount of pls
df <- df %>% subset(regression_method %in% c("mnlogistic", "rejection"))
df$regression_method <- factor(df$regression_method, levels = c("rejection", "mnlogistic"))

# explify sscomp
df$sscomp <- sapply(df$sscomp, function(x) {
  switch(
    as.character(x),
    "1" = {
      return("sfs/ld")
    },
    "2" = {
      return("tl_true")
    },
    "3" = {
      return("tm_true")
    },
    "4" = {
      return("tm_win")
    },
    "5" = {
      return("sfs/ld/tm_win")
    }
  )
})
df$sscomp <- factor(df$sscomp, levels = c(
  "sfs/ld", "tl_true", "tm_true", "tm_win", "sfs/ld/tm_win"
))

df <- na.omit(df)

discretize_bayes_factor <- function(bf, do.factor = T) {
  if (bf <= 10 ** 0)
    dbf = "Negative"
  else if (bf <= 10 ** 0.5)
    dbf = "Barely worth mentioning"
  else if (bf <= 10 ** 1)
    dbf = "Substantial"
  else if (bf <= 10 ** 1.5)
    dbf = "Strong"
  else if (bf <= 10 ** 2)
    dbf = "Very strong"
  else if (bf > 10 ** 2)
    dbf =  "Decisive"
  else
    stop(paste0(c("unknown bayes: ", as.character(bf)), collapse = ""))
  
  dbf <- factor(
    dbf,
    levels = c(
      "Negative",
      "Barely worth mentioning",
      "Substantial",
      "Strong",
      "Very strong",
      "Decisive"
    )
  )
  
  if (do.factor) {
    return(dbf)
  } else {
    return(as.numeric(dbf)-1)
  }
}

df$discrete_bayes_factor <-
  sapply(df$bayes_factor, discretize_bayes_factor)

# df %>%
#   # subset(regression_method %in% "mnlogistic") %>%
#   ggplot(aes(x = as.factor(param_5),
#              y = as.numeric(discrete_bayes_factor),
#              fill = regression_method)) +
#   geom_point(aes(col = regression_method),
#              position = position_jitter(width = 0.2, height = 0.1))+
#   facet_grid(sscomp~pls)+
#   scale_y_continuous(breaks = 0:6)+
#   scale_color_manual(values = c("darkgreen", "hotpink3"))+
#   theme(
#     axis.text.x = element_text(angle = 60, hjust = 1)
#   )

df %>%
  # subset(regression_method == "mnlogistic") %>%
  # subset(pls <= 30) %>%
  subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win", "tl_true")) %>%
  ggplot(aes(x = as.factor(param_5), fill=discrete_bayes_factor))+
  geom_bar(position = "fill")+
  facet_grid(sscomp~regression_method+pls)+
  scale_fill_manual(values = mcol)+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

# subset percentage calculation
data <- df %>%
  subset(sscomp %in% c("sfs/ld", "tm_win", "sfs/ld/tm_win"))

data$N <- NA
data$n <- NA
for (a in unique(data$regression_method)) {
  for (b in unique(data$param_5)) {
    for (c in unique(data$sscomp)) {
      for (d in unique(data$pls)) {
        N <- data %>%
          subset(regression_method == a) %>%
          subset(param_5 == b) %>%
          subset(sscomp == c) %>%
          subset(pls == d) %>%
          nrow()
        data$N[which(data$regression_method == a &
                       data$param_5 == b &
                       data$sscomp == c &
                       data$pls == d)] <- N
        
        for (e in unique(data$discrete_bayes_factor)) {
          n <- data %>%
            subset(regression_method == a) %>%
            subset(param_5 == b) %>%
            subset(sscomp == c) %>%
            subset(pls == d) %>%
            subset(discrete_bayes_factor == e) %>%
            nrow()
          data$n[which(
            data$regression_method == a &
              data$param_5 == b &
              data$sscomp == c &
              data$pls == d &
              data$discrete_bayes_factor == e
            
          )] <- n
        }
      }
    }
  }
}
data$bayes_factor = NULL

data <- data %>%
  group_by(param_5, regression_method, pod, sscomp, pls, discrete_bayes_factor, N) %>%
  summarize(prop=n()/N) %>%
  ungroup() %>%
  distinct()

logbreak <- sapply(10**(-10:10), function(x) {
  x*(1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10**(-10:10))
    return(as.character(log10(x)))
  else
    return("")
})

# add zeros
zcounter = 0
for (a in unique(data$regression_method)) {
  for (b in unique(data$param_5)) {
    for (c in unique(data$sscomp)) {
      for (d in unique(data$pls)) {
        for (e in unique(data$discrete_bayes_factor)) {
          n <- data %>%
            subset(regression_method == a) %>%
            subset(param_5 == b) %>%
            subset(sscomp == c) %>%
            subset(pls == d) %>%
            subset(discrete_bayes_factor == e) %>%
            nrow()
          if(n == 0) {
            zcounter = zcounter + 1
            data = rbind.data.frame(
              data,
              data.frame(
                param_5=b,
                regression_method=a,
                pod=NA,
                sscomp=c,
                pls=d,
                discrete_bayes_factor=e,
                N=100,
                prop=0
              )
            )
          }
        }
      }
    }
  }
}
print(zcounter)

data$n = NULL

# relabel the pls
data$pls = sapply(data$pls, function(x) {
  if (x == 0) {
    r = "full"
  } else if (x == 5) {
    r = "5-pls"
  } else if (x == 10) {
    r = "10-pls"
  } else if (x == 20) {
    r = "20-pls"
  } else if (x == 30) {
    r = "30-pls"
  } else {
    print(x)
    stop("you are completely useless")
  }
  return(r)
})
data$pls = factor(data$pls, levels = c("full", "5-pls", "10-pls", "20-pls", "30-pls"))

data %>%
  subset(regression_method == "mnlogistic") %>%
  subset(pls != "30-pls") %>%
  # subset(discrete_bayes_factor == "Negative") %>%
  ggplot(aes(x=param_5, y=prop*100, fill=discrete_bayes_factor))+
  geom_area()+
  geom_hline(yintercept = 5, col="white", size=0.1)+
  geom_hline(yintercept = 20, col="white", size=0.1)+
  geom_hline(yintercept = 50, col="white", size=0.1)+
  geom_hline(yintercept = 80, col="white", size=0.1)+
  geom_hline(yintercept = 95, col="white", size=0.1)+
  geom_vline(xintercept = 50000, col="white", size=0.1)+
  # facet_grid(sscomp~regression_method+pls)+
  facet_grid(sscomp~pls)+
  scale_x_continuous(trans = "log", breaks = logbreak, labels = loglabel)+
  scale_y_continuous(breaks = c(5, 95), labels=c("0%", "100%"))+
  scale_fill_manual(values = mcol)+
  theme(
    # axis.text.x = element_text(angle = 60, hjust = 1),
    aspect.ratio = 1,
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    # axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(-0.1, "lines"),
    panel.spacing = unit(-0.9, "lines")
  )+
  labs(x="time after transition [log10]",
       fill="model support",
       main="dada")

data$param_5 %>% unique() %>% sort()

mu_r_five_data <- data %>%
  subset(regression_method == "mnlogistic") %>%
  subset(pls != "30-pls")



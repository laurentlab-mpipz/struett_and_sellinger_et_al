

# Read the config file and the RDS file from Thibaut's teSMC, create plot

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(yaml)
library(wesanderson)
library(dplyr)

# save.image()
# stop(paste0(rep("#", 600), collapse = ""))
# setwd(dir = paste0("/Users/struett/MPIPZ/netscratch/dep_tsiantis/grp_laurent/",
#                    "struett/temp_mspts/generate_ts_data/",
#                    collapse = "")
# )
# load(file = ".RData")

# wes_palettes %>% names
col_true = wes_palette(names(wes_palettes[6]))[1]
col_inf = wes_palette(names(wes_palettes[7]))[5]
mcol <- c(col_true, col_inf)

# read config file as yaml
N = yaml.load_file(snakemake@params[["path_to_config"]])[[
  "population_sizes_backward_in_time"]][[1+as.numeric(
    snakemake@wildcards$demography)]] %>% data.frame %>% t %>% data.frame(
      row.names = NULL); colnames(N) <- c("N", "t")
sigma = yaml.load_file(snakemake@params[["path_to_config"]])[[
  "selfing_rates_backward_in_time"]][[1+as.numeric(
    snakemake@wildcards$demography)]] %>% unlist %>% matrix(ncol=2) %>%
  data.frame; colnames(sigma) <- c("sigma", "t")

# read estimated demographies
infile_list = snakemake@input[["data"]]
result_list = list()
df_list = list()
for (i in 1:length(infile_list)) {
# for (i in c(1,3,10)) {
  result = readRDS(infile_list[i])
  
  # further params
  mu_estim = result$mu
  L = result$L
  real_mu = as.numeric(snakemake@wildcards$mutation_rate)
  
  df_list[[i]] <- cbind.data.frame(
    LH = result$LH, #
    Tc = result$Tc, # determined piecewise constant time frames
    Xi = result$Xi, # estimated Ne
    mu = result$mu, # estimated mutation rate
    beta = result$beta, # estimated (or not) seed bank
    sigma = result$sigma, # estimated sigma
    rho = result$rho, # estimated rho
    rep = as.numeric(strsplit(infile_list[i], split="/|\\.|_")[[1]][which(
      strsplit(infile_list[i], split="/|\\.|_")[[1]] == "rep")+1]),
    real_mu = real_mu,
    L = L,
    mu_estim = mu_estim
  )
}

df <- do.call("rbind.data.frame", df_list)

# add time intervals to true df according to the estimated time frames
accomodate_times_to_estimated_times <- function(df, times, time_scaling) {
  # assuming two columns: 1. value, 2. time
  lines_to_add = list()
  i = 0
  for (t in times) {
    # rescale the time
    t = t * time_scaling
    # find corresponding time in df
    if (! t %in% df$t ) {
      i = i + 1
      val = df[max(which(df$t < t)),1]
      lines_to_add[[i]] = c(val, t)
    }
  }
  a = do.call("rbind.data.frame", lines_to_add)
  names(a) <- names(df)
  
  df <- rbind.data.frame(df, a) %>% arrange(., t)
  
  return(df)
}
N = accomodate_times_to_estimated_times(N, unique(df$Tc), mu_estim/real_mu)
sigma = accomodate_times_to_estimated_times(sigma, unique(df$Tc), mu_estim/real_mu)

# adjust format
df$rep = as.factor(df$rep)

## plot N, r, sigma over time
p1 <- ggplot()+
  geom_step(data = N,
            aes(x=t,y=N),
            col=mcol[1],
            size = 2)+
  geom_step(data = df,
            aes(x=Tc*(mu_estim/real_mu), y=Xi*0.5*mu_estim/real_mu, fill=rep), 
            col=mcol[2])+
  scale_x_continuous(trans = "log10",
                     limits = c(3e2, 18*N$N[1]),
                     breaks = 10**(1:10))+
  scale_y_continuous(trans = "log10")+
  labs(
    # title=expression('N'[e] * ' over time'),
    x = "scaled time",
    y = expression('N'[e]))+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

p2 <- ggplot()+
  geom_step(data = sigma,
            aes(x=t,y=sigma),
            col=mcol[1],
            size=2)+
  geom_step(data = df,
            aes(x=Tc*(mu_estim/real_mu), y=sigma, fill=rep), 
            col=mcol[2])+
  scale_x_continuous(trans = "log10",
                     limits = c(3e2, 18*N$N[1]),
                     breaks = 10**(1:10))+
  scale_y_continuous(limits = c(0, 1))+
  labs(
    # title=expression(sigma * ' over time'),
    x = "scaled time",
    y = expression(sigma))+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

p3 <- ggplot()+
  geom_step(data = df,
            aes(x=Tc*(mu_estim/real_mu),y=rho,fill=rep),
            col=mcol[2])+
  scale_x_continuous(trans = "log10",
                     limits = c(3e2, 18*N$N[1]),
                     breaks = 10**(1:10))+
  # scale_y_continuous(limits = c(0, 1))+
  labs(
    # title=expression(rho * ' over time'),
    x = "scaled time",
    y = expression(rho))+
  theme(legend.position = "bottom")

title <- ggdraw()+
  draw_label("Simulated vs. inferred (teSMC) demography",
             fontface = 'bold', x = 0, hjust = 0, vjust = 1)

rwth = -0.1 # e.g. -0.2, negative values shrinkt the distance between plots
plot_list = list(title, NULL, p1, NULL, p2, NULL, p3)
pt <- plot_grid(
  plotlist = plot_list,
  align = T,
  ncol = 1,
  rel_heights = c(0.3, rwth, rep(c(1, rwth), length(plot_list)))
)

# plot N over time
cat("target file: ", snakemake@output[["pdf"]], "\n")
pdf(snakemake@output[["pdf"]], useDingbats = F)
show(pt)
dev.off()



library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)

mcol = wes_palette("Darjeeling1")[c(1, 2, 3, 5)]

# # save image for debugging
# save.image(); stop("dadada+++++++++++++++++++")
# setwd("../..")
# load(".RData")

infile = snakemake@input$pesti
outfile = snakemake@output$table
outplot = snakemake@output$plot

df <- read.csv(infile)
df$X = NULL

r_list = list()
r_ix = 0
for (a in unique(df$infmodel)) {
  for (b in unique(df$demography)) {
    for (c in unique(df$sample_size)) {
      for (d in unique(df$mutation_rate)) {
        for (e in unique(df$recombination_rate)) {
          for (f in unique(df$mutation_rate)) {
            dfsub <- df %>%
              subset(infmodel == a) %>%
              subset(demography == b) %>%
              subset(sample_size == c) %>%
              subset(mutation_rate == d) %>%
              subset(recombination_rate == e)
            
            stopifnot(length(unique(dfsub$t_sigma_true))==1)
            
            rb = 1/nrow(dfsub) * sum((dfsub$t_sigma_inferred - dfsub$t_sigma_true)/dfsub$t_sigma_true)
            rmse = sqrt(1/nrow(dfsub) * sum(((dfsub$t_sigma_inferred - dfsub$t_sigma_true)/dfsub$t_sigma_true)**2))
            
            r = data.frame(a, b, c, d, e, rb, rmse, unique(dfsub$t_sigma_true))
            r_ix = r_ix + 1
            r_list[[r_ix]] = r
          }
        }
      }
    }
  }
}
r = do.call("rbind.data.frame", r_list) %>% data.frame
colnames(r) = c("infmodel", "demography", "sample_size", "mutation_rate", "recombination_rate", "rb", "rmse", "true_value")

write.csv(x = r, file = outfile, row.names = FALSE)

r$infmodel <- factor(r$infmodel, levels=c("Free", "OneTransition", "GivenTransition", "Constant"))

p_rb <- r %>%
  ggplot(aes(as.factor(true_value), rb*true_value, fill=infmodel))+
  geom_col(position = "dodge")+
  facet_grid(infmodel~., scales = "free_y")+
  scale_fill_manual(values = mcol)+
  labs(x=expression('t'[sigma]), title = "Bias")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_rmse <- r %>%
  ggplot(aes(as.factor(true_value), (rmse)*true_value, fill=infmodel))+
  geom_col(position = "dodge")+
  facet_grid(infmodel~., scales = "free_y")+
  scale_fill_manual(values = mcol)+
  labs(x=expression('t'[sigma]), title = "Precision")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_corr <- r %>%
  ggplot(aes(as.factor(true_value), (rmse-abs(rb))*true_value, fill=infmodel))+
  geom_col(position = "dodge")+
  facet_grid(infmodel~., scales = "free_y")+
  scale_fill_manual(values = mcol)+
  labs(x=expression('t'[sigma]), title = "Unbiased precision")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p <- plot_grid(p_rb, p_rmse, p_corr, nrow = 1, align = T, labels = "AUTO")

pdf(outplot, width = 12)
show(p)
dev.off()

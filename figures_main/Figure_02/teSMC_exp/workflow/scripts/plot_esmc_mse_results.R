
## save image for debugging
# save.image(); stop("dadada+++++++++++++++++++")
# setwd("../..")
# load(".RData")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)

sigma_d = snakemake@config["selfing_rates_backward_in_time"][[1]]
mcol = wes_palette("Darjeeling1")[c(1,2,3,5)]

df <- read.csv(snakemake@input$df, row.names = 1, as.is = T)
df$param[df$param == "pop"] = "demography"

get_transition_time_from_demogr_scenario <- function(d, sigma_d) {
  t = sigma_d[(d*4+1):(d*4+4)][4]
  return(t)
}

df$tr <- sapply(df$demography, function(x) {
  return(get_transition_time_from_demogr_scenario(x, sigma_d))
})

p <- ggplot(df, aes(x=as.factor(tr), y=mse, fill=infmodel,
               group=interaction(rep, infmodel)))+
  geom_col(position = "dodge")+
  facet_grid(param~., scales = "free_y")+
  scale_fill_manual(values = mcol)+
  theme(
    aspect.ratio = 0.707/2,
    axis.text.x = element_text(angle = 90)
  )

pdf(snakemake@output$pdf)
show(p)
dev.off()

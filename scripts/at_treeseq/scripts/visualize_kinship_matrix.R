
setwd("~/Desktop/tsinfer_on_1135genomes_at/")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(reshape2)

pd <- reticulate::import("pandas")

infile = "1135_kinship_matrix.gzip"

df <- pd$read_pickle(infile)
df$sample_row = rownames(df)

a <- df %>%
  gather(key="sample_col", value="kinship", -sample_row)
a$sample_row = as.numeric(a$sample_row)
a$sample_col = as.numeric(a$sample_col)

p <- a %>%
  ggplot(aes(x=sample_col, y=sample_row, fill=kinship))+
  geom_tile()+
  scale_fill_viridis_c()
p









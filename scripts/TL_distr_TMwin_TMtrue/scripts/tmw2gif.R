
# make gif from tl_true

library(dplyr)
library(purrr) 
library(magick)

args <- commandArgs(trailingOnly=TRUE)
outfile = args[1]
infiles = args[2:length(args)]

# read files
list_of_files <- infiles

# sort files
my_order <- sapply(list_of_files, function(x) {
  as.numeric(strsplit(strsplit(x, "time")[[1]][2], "_")[[1]][2])
})
sorted_files = list_of_files[order(my_order)]

message("sorted files:")
message(paste(sorted_files, "\n"))

# create gif
tm_win_gif <- 
  sorted_files %>% 
  map(image_read) %>% 
  image_join() %>% 
  image_animate(delay = c(
    75, # 1
    rep(10, 39),
    50,
    rep(10, 89),
    50)) %>% 
  image_write(outfile, quality = 100)




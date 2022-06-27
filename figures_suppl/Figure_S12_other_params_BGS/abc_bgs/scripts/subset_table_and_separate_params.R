
library(reticulate); pd <- import("pandas")
#library(data.table)

source("scripts/rfunctions.R")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)

# explify args
transformer_file = args[1]
sumstat_file = args[2]
outfile_prefix = args[3]

# define output file
outfile_sumstat = paste0(outfile_prefix, ".sumstat", collapse="")

# read files 
transformer <- read.table(transformer_file, head=F, sep="\t", as.is=T)
sumstats <- pd$read_pickle(sumstat_file)

# save.image("/netscratch/dep_tsiantis/grp_laurent/struett/temp_tsabc/mchoice/tsabc3/TEMP.RData");
# stop("#######################################################################")
# load("TEMP.RData")

# explify names of columns to transform
colnames_to_be_transformed = transformer[,1]

# subset wanted sumstats
sumstats_of_interest <- dplyr::select(sumstats,
	all_of(colnames_to_be_transformed))

# extract params
params_all <- obtain_columns_with_wanted_sumstats(sumstats, "param")

# write subsetted data frames (params, sumstats) into separate files as txt
write.table(cbind(params_all, sumstats_of_interest),
	outfile_sumstat, quote = F, sep = "\t", row.names=F)
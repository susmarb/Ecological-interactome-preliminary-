#!/bin/R

# This script performs centered log-ratio tranformations

# define arguments
args <- commandArgs(TRUE)
data_in <- args[1]
data_out <- args[2]

packs <- c("tidyverse", "compositions")
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}

# input table has MAGs (or equivalent) as one column called "variable", and dates as the rest of the columns
# Done once, from relative abundance to centered log-ratio transformations
read_tsv(data_in) %>%
 gather(date, relab, -variable) %>%
 spread(variable, relab) %>%
 column_to_rownames("date") %>% as.matrix() %>%
 clr() %>%
 write.table(., data_out, sep = "\t", quote = F, col.names = T, row.names = T)

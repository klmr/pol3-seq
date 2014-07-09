#!/usr/bin/env Rscript
filename = commandArgs(trailingOnly = TRUE)[1]
data = read.table(filename, sep = '\t', header = FALSE, row.names = 1)
data = data[nrow(data) - 5, , drop = FALSE]
summary(data[, 1])

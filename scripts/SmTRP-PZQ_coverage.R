#!/usr/bin/env Rscript
# Title: SmTRP-PZQ_coverage.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-02-11
# Modified in: 



#==========#
# Comments #
#==========#

# Plot read depth of the SmTRP-PZQ (Smp_246790) gene



#==========#
# Versions #
#==========#

# v0.0 - 2021-02-11: creation



#===================#
# Packages required #
#===================#

suppressMessages({
    library("magrittr")
})



#===========#
# Variables #
#===========#

# Folders
graph_fd  <- "../graphs/"
result_fd <- "../results/2-Coverage/"

cov_file <- paste0(result_fd, "SmTRP-PZQ.cov")



#=================#
# Data processing #
#=================#

cat("Loading and processing data. This may take some time...\n")

mydata <- read.delim(cov_file, header = FALSE, stringsAsFactors = FALSE)

# Position along the gene
mypos  <- mydata[2, 2:ncol(mydata)] %>% as.numeric()

# Summaryzing read depth
myrd <- mydata[3:nrow(mydata),2:ncol(mydata)]
myrd <- apply(myrd, 1, as.numeric)
myrd <- rowSums(myrd)
myrd <- runmed(myrd,101)



#=========#
# Figures #
#=========#

cat("Drawing graph...\n")

if (! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

png(paste0(graph_fd, "SmTRP-coverage.png"), width=50*72)
# Point graph
plot(myrd ~ mypos, xlab = "Position on the gene (bp)", ylab = "Read depth", ylim = c(10, max(myrd)), xlim = range(mypos), main = "Smp_246790", pch = 20, col = "grey", log = "y")

# # Smoothed line
# lw1 <- loess(myrd ~ mypos, span=0.05)
# lines(mypos, predict(lw1, mypos), col = "red", lwd = 3)
dev.off()

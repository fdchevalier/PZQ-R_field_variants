#!/usr/bin/env Rscript
# Title: variant_read_depth.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-02-01
# Modified in: 



#==========#
# Comments #
#==========#

# v0.0 - 2021-02-01: creation



#======================#
# Packages and options #
#======================#

cat("Loading packages...\n")

suppressMessages({
    library("Gviz")
    library("Biostrings")
    library("magrittr")
})

options(ucscChromosomeNames=FALSE)



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd    <- "../data/"
graph_fd   <- "../graphs/"
result_fd  <- "../results/"

mut_fl    <- paste0(result_fd, "2-mutations\ of\ interest/samples_list.tsv")
genome_fl <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa")



#=================#
# Data processing #
#=================#

cat("Processing data...\n")

if (! dir.exists(graph_fd)) { dir.create(graph_fd) }

# Load genome
sm_genome <- readDNAStringSet(genome_fl)

# Update chromosome names by removing extra information
names(sm_genome) <- sapply(names(sm_genome), function(x) strsplit(x, " ")) %>% sapply(., "[[", 1)

# Number of columns
no_col <- max(count.fields(mut_fl, sep = "\t"))
mymut <- read.delim(mut_fl, header = FALSE, col.names = 1:no_col)

# Genome track
sTrack <- SequenceTrack(sm_genome)

for (i in 1:nrow(mymut)) {
    m <- mymut[i, 1]
    myspl <- mymut[i, 3:no_col] %>% .[ . != "" ]

    for (s in myspl) {
        mybam <- dir(paste0(data_fd, "libraries/", s), pattern = ".*recal.bam$", full.names = TRUE)

        mypos   <- mymut[i, 2]
        myshift <- 100

        afrom <- mypos - myshift
        ato   <- mypos + myshift
        
        alTrack <- AlignmentsTrack(mybam, isPaired = TRUE)
        ht      <- HighlightTrack(alTrack, start = mypos, width = 0.5, chromosome = "SM_V7_3", inBackground = FALSE, alpha = 0.2)
        pdf(paste0(graph_fd, m, " - ", s, ".pdf"), width=10)
        plotTracks(c(ht, sTrack), chromosome = "SM_V7_3", from = afrom, to = ato, main = paste0(m, " - ", s))
        dev.off()
    }
}



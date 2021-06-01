#!/usr/bin/env Rscript
# Title: SmTRP-PZQ_coverage.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-02-11
# Modified in: 2021-05-31



#==========#
# Comments #
#==========#

# Plot read depth of the SmTRP-PZQ (Smp_246790) exons



#==========#
# Versions #
#==========#

# v0.1 - 2021-05-31: limit plotting to exons
# v0.0 - 2021-02-11: creation



#===================#
# Packages required #
#===================#

suppressMessages({
    library("magrittr")
    library("rtracklayer")
    library("Sushi")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("coord_intersect.R")



#===========#
# Variables #
#===========#

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/2-Coverage/"

cov_file <- paste0(result_fd, "SmTRP-PZQ.cov")

# GFF file
mygff.fl <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
mygff <- readGFF(mygff.fl)

# Samples to select
sel_spl <- c("Sm")

# Filtering parameters
mydp  <- 10
mycov <- 0.8 # must be between ]0; 1]

myseed <- 1453225



#=================#
# Data processing #
#=================#

cat("Loading and processing data. This may take some time...\n")

mydata <- read.delim(cov_file, header = FALSE, stringsAsFactors = FALSE)

# Position along the gene
mypos  <- mydata[2, 2:ncol(mydata)] %>% as.numeric()

# Sample selection
myspl  <- mydata[,1] %>% grep(sel_spl, .)
mydata <- mydata[c(1:2, myspl), ]

# Summaryzing read depth
myrd <- mydata[3:nrow(mydata),2:ncol(mydata)]
myrd <- apply(myrd, 1, as.numeric)
# myrd <- rowSums(myrd)
myrd <- rowMeans(myrd)
myrd <- runmed(myrd, 101)
myrd <- data.frame(t(mydata[1:2, -1]), myrd, stringsAsFactors = FALSE)

# Gene coordinates
a <- mygff[, "Parent"] %>% grep("Smp_246790.", .)
mygff_a <- mygff[a,]

b <- mygff_a[, "type"] == "CDS"
mygff_b <- mygff_a[b, ]

mygff_b[, "score"] <- 0

# Generate bed-like coordinate table
mybed <- export(mygff_b, format="BED") %>% unique() %>% sapply(., function(x) strsplit(x, "\\t")) %>% do.call("rbind", .) %>% unname() %>% as.data.frame()
mybed[, 2:3] <- apply(mybed[, 2:3], 2, as.numeric)
mybed <- mybed[order(mybed[,2]),]


mydata_t   <- t(mydata)
mydata_t_b <- coord_intersect(mydata_t, mybed)
cov_exons  <- lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist

for (i in 1:length(cov_exons)) {cat("Coverage of exon ", i,": ", cov_exons[i], "\n", sep="") }

# Exon read depth only
myrd_bed <- coord_intersect(myrd, mybed)


# Other stats
## Total number of exons
((lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist) > 0) %>% sum()
## Average number of sample per exon
lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist %>% mean()
## Standard error
(lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist %>% sd()) / sqrt(length(mydata_t_b))



#=========#
# Figures #
#=========#

cat("Drawing graph...\n")

if (! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

# Gene version
png(paste0(graph_fd, "SmTRP-coverage.png"), width=50*72, height=10*72)
layout(matrix(1:2, ncol=1), height = c(0.8, 0.2))
# Point graph
plot(myrd[,2:3], xlab = "", ylab = "Read depth", ylim = c(0, max(myrd[, 3])), xlim = range(myrd[, 2]), main = "Smp_246790", pch = 20, col = "grey", xaxt = "n", bty = "n") # log = "y", axes = FALSE)


par(mar = c(4, 5, 0, 2)+0.1)
plot(mypos, xlim = c(mypos[1], mypos[length(mypos)]), ylim=c(-.25, .25), type="n", bty="n", axes=FALSE, ann=FALSE)
my.y <- par("usr")[4]*0.2

for (e in 1:nrow(mybed)) {
    rect(mybed[e, 2], -my.y, mybed[e, 3], my.y, col="black", border=NA)

    if (e < nrow(mybed)) {
        lines(c(mybed[e, 3], (mybed[e+1, 2] -  mybed[e, 3]) / 2 +  mybed[e, 3]), c(0, my.y), col = "black")
        lines(c((mybed[e + 1, 2] -  mybed[e, 3]) / 2 + mybed[e, 3], mybed[e + 1, 3]), c(my.y, 0), col = "black")
    }
}


myrg <- range(mybed[,2:3])
mylg <- myrg %>% log10() %>% floor() %>% min() - 1
myat <- ((myrg / 10^mylg) %>% round()) * 10^mylg

axis (1, xaxp = c(myat, 8)) #, cex.axis=mycex.axis)
title(xlab = "Position on the gene (bp)")

dev.off()

# Exon version
png(paste0(graph_fd, "SmTRP-coverage_exons.png"), width=50*72, height=10*72)

layout(matrix(1:2, ncol=1), height = c(0.8, 0.2))

myclr <- rainbow(nrow(mybed))
set.seed(myseed)
#myclr <- sample(myclr, length(myclr))
myclr <- rep("black", length(myclr))

# Point graph
gene_lg <- lapply(myrd_bed, function(x) nrow(x)) %>% unlist() %>% sum()
plot(NULL, xlab = "", ylab = "Read depth", ylim = c(0, max(myrd[, 3])), xlim = c(1, gene_lg), main = "Smp_246790", xaxt = "n", bty = "n") # log = "y", axes = FALSE)

mypos2 <- 1
for (e in 1:length(myrd_bed)) {
    mypos1 <- mypos2
    mypos2 <- nrow(myrd_bed[[e]]) + mypos1 - 1
    mypos_tmp <- mypos1:mypos2

    points(myrd_bed[[e]][, 3] ~ mypos_tmp, pch = 20, col = myclr[e])
    if (e < length(myrd_bed)) { abline(v = mypos2, lty = 3) }
}
   
par(mar = c(4, 5, 0, 2)+0.1)
plot(NULL, xlim = c(1, gene_lg), ylim=c(-.25, .25), type="n", bty="n", axes=FALSE, ann=FALSE)
my.y <- par("usr")[4]*0.2


mypos2 <- 1
for (e in 1:nrow(mybed)) {
    mypos1 <- mypos2
    mypos2 <- nrow(myrd_bed[[e]]) + mypos1 - 1

    rect(mypos1, -my.y, mypos2, my.y, col=myclr[e], border=NA)
    if (e < length(myrd_bed)) { abline(v = mypos2, lty = 3) }
}

dev.off()


# Alternative exon version
png(paste0(graph_fd, "SmTRP-coverage_exons.png"), width=50*72, height=10*72)

myclr <- rainbow(nrow(mybed))
set.seed(myseed)
#myclr <- sample(myclr, length(myclr))
myclr <- rep("black", length(myclr))

# Point graph
gene_lg <- lapply(myrd_bed, function(x) nrow(x)) %>% unlist() %>% sum()
plot(NULL, xlab = "", ylab = "Read depth", ylim = c(-10, max(myrd[, 3])), xlim = c(1, gene_lg), main = "Smp_246790", xaxt = "n", bty = "n") # log = "y", axes = FALSE)
# my.y <- par("usr")[4]*0.2
my.y <- -10 * 0.2

mypos2 <- 1
for (e in 1:length(myrd_bed)) {
    mypos1 <- mypos2
    mypos2 <- nrow(myrd_bed[[e]]) + mypos1 - 1
    mypos_tmp <- mypos1:mypos2

    points(myrd_bed[[e]][, 3] ~ mypos_tmp, pch = 20, col = myclr[e])

    rect(mypos1, -10-my.y, mypos2, -10+my.y, col=myclr[e], border=NA)
    if (e < length(myrd_bed)) { abline(v = mypos2, lty = 3) }
}

dev.off()

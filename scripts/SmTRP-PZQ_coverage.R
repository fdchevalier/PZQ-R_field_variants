#!/usr/bin/env Rscript
# Title: SmTRP-PZQ_coverage.R
# Version: 0.5
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-02-11
# Modified in: 2021-11-23



#==========#
# Comments #
#==========#

# Plot read depth of the SmTRP-PZQ (Smp_246790) exons



#==========#
# Versions #
#==========#

# v0.5 - 2021-11-23: correct file path
# v0.4 - 2021-08-08: add panel B / clean code
# v0.3 - 2021-08-06: improve plot
# v0.2 - 2021-06-09: add high frequency mutations / add new plot layout
# v0.1 - 2021-05-31: limit plotting to exons
# v0.0 - 2021-02-11: creation



#===================#
# Packages required #
#===================#

suppressMessages({
    library("magrittr")
    library("rtracklayer")
    library("Sushi")
    library("rsvg")
    library("grImport2")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("functions/coord_intersect.R")
source("functions/line2user.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/"

cov_file <- paste0(result_fd, "2-coverage/SmTRP-PZQ.cov")

# GFF file
mygff.fl <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
mygff    <- readGFF(mygff.fl)

# Variant file
var_fl <- paste0(result_fd, "1-reports/PZQ-R_field_stringent_Smp_246790.5.flt.norm.tsv")
var    <- read.delim(var_fl)

svg_fl <- paste0(graph_fd, "Sm.TRPM_PZQ structure + mutation.svg")

# High frequency variants
min_rd <- 25
min_af <- 0.05

# Resistant mutation
res_pos <- 737367

# Samples to select
sel_spl <- c("Sm")

# Filtering parameters
mydp  <- 10
mycov <- 0.8 # must be between ]0; 1]



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

cat("\nNumber of samples showing ", mycov * 100, "% of exon coverage with a minimum read depth of ", mydp, ":\n", sep="")
for (i in 1:length(cov_exons)) {cat("\t- Exon ", i,": ", cov_exons[i], "\n", sep="") }

# Exon read depth only
myrd_bed <- coord_intersect(myrd, mybed)


# Other stats
## Total number of exons
tot_cov <- ((lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist) > 0) %>% sum()
## Average number of samples per exon
avg_spl <- lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist %>% mean()
## Standard error
se <- (lapply(mydata_t_b, function(x) ((x[,-(1:2)] >= mydp) %>% colSums() >= nrow(x) * mycov) %>% sum() ) %>% unlist %>% sd()) / sqrt(length(mydata_t_b))

cat("\nStatistics of global coverage:\n", sep="")
cat("\t- Total number of exons covered: ", tot_cov, "\n", sep="")
cat("\t- Average number of samples per exon (mean ± s.e.): ", avg_spl, " ± ", se, "\n", sep="")


# Select high frequency variants + resistant variant
var <- var[var[,4] == "CDS" & var[, 8] != "-", ]
var <- var[rowSums(var[,9:10]) > min_rd, ]
var <- var[var[,10] / rowSums(var[,9:10]) > 0.05 | var[,1] == res_pos, 1:10]

# Select exon variants
var2   <- cbind("SM_V7_3", var)
var2_b <- coord_intersect(var2, mybed)



#=========#
# Figures #
#=========#

cat("Drawing graphs...\n")

if (! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

# Gene version
png(paste0(graph_fd, "SmTRP-coverage.png"), width=50*72, height=10*72)
layout(matrix(1:2, ncol=1), height = c(0.8, 0.2))
# Point graph
plot(myrd[,2:3], xlab = "", ylab = "Read depth", ylim = c(0, max(myrd[, 3])), xlim = range(myrd[, 2]) %>% as.numeric(), main = "Smp_246790", pch = 20, col = "grey", xaxt = "n", bty = "n") # log = "y", axes = FALSE)


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
pdf(paste0(graph_fd, "Supp. Fig. S9 - SmTRP-coverage_exons.pdf"), width = 17, height = 8)

par(mar = c(5, 4, 4, 0) + 0.1)

layout(matrix(1:3, ncol = 1))

myclr <- rep("grey", length(myrd_bed))
mycex <- 0.85

# Point graph
gene_lg1 <- lapply(myrd_bed[1:floor(length(myrd_bed)/2)], function(x) nrow(x)) %>% unlist() %>% sum()
gene_lg2 <- lapply(myrd_bed[1:length(myrd_bed)], function(x) nrow(x)) %>% unlist() %>% sum()

plot(NULL, xlab = "", ylab = "Read depth", ylim = c(-10, max(myrd[, 3])), xlim = c(1, gene_lg1), xaxt = "n", bty = "n") # log = "y", axes = FALSE)
# my.y <- par("usr")[4]*0.2
my.y <- -10 * 0.2
mtext(LETTERS[1], side=3, line=2, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

mypos2 <- 1
for (e in 1:floor(length(myrd_bed)/2)) {
    mypos1 <- mypos2
    mypos2 <- nrow(myrd_bed[[e]]) + mypos1 - 1
    mypos_tmp <- mypos1:mypos2

    # points(myrd_bed[[e]][, 3] ~ mypos_tmp, pch = 20, col = myclr[e])
    lines(myrd_bed[[e]][, 3] ~ mypos_tmp, pch = 20, col = myclr[e])

    rect(mypos1, -10-my.y, mypos2, -10+my.y, col=myclr[e], border=NA)
    if (e < length(myrd_bed)) { abline(v = mypos2, lty = 3) }

    if (dim(var2_b[[e]])[1] > 0) {
        var_tmp <- var2_b[[e]]
        for (l in 1:nrow(var_tmp)) {
            if (var_tmp[l, 2] == res_pos) { mybx_clr <- "red" } else {mybx_clr <- "blue" }
            v <- which(myrd_bed[[e]][, 2] == var_tmp[l, 2]) + mypos1 - 1
            rect(v-10, -10-my.y, v+10, -10+my.y, col = mybx_clr, border = NA)
            text(v, -20-my.y, labels = var_tmp[l, 9], adj = 1, srt = 90, cex = mycex, xpd = TRUE)
        }
    }

    # Add exon number
    text(mean(mypos_tmp), par("usr")[4]*1.15, labels = e, xpd = TRUE)
}

plot(NULL, xlab = "", ylab = "Read depth", ylim = c(-10, max(myrd[, 3])), xlim = c(gene_lg1, gene_lg2), xaxt = "n", bty = "n") # log = "y", axes = FALSE)
mypos2 <- gene_lg1
for (e in ceiling(length(myrd_bed)/2):length(myrd_bed)) {
    mypos1 <- mypos2
    mypos2 <- nrow(myrd_bed[[e]]) + mypos1 - 1
    mypos_tmp <- mypos1:mypos2

    # points(myrd_bed[[e]][, 3] ~ mypos_tmp, pch = 20, col = myclr[e])
    lines(myrd_bed[[e]][, 3] ~ mypos_tmp, pch = 20, col = myclr[e])

    rect(mypos1, -10-my.y, mypos2, -10+my.y, col=myclr[e], border=NA)
    if (e < length(myrd_bed)) { abline(v = mypos2, lty = 3) }

    if (dim(var2_b[[e]])[1] > 0) {
        var_tmp <- var2_b[[e]]
        for (l in 1:nrow(var_tmp)) {
            if (var_tmp[l, 2] == res_pos) { mybx_clr <- "red" } else {mybx_clr <- "blue" }
            v <- which(myrd_bed[[e]][, 2] == var_tmp[l, 2]) + mypos1 - 1
            rect(v-10, -10-my.y, v+10, -10+my.y, col = mybx_clr, border = NA)
            text(v, -20-my.y, labels = var_tmp[l, 9], adj = 1, srt = 90, cex = mycex, xpd = TRUE)
        }
    }
    
    # Add exon number
    text(mean(mypos_tmp), par("usr")[4]*1.15, labels = e, xpd = TRUE)
}


plot.new()
mtext(LETTERS[2], side=3, line=2, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

mysvg <- rsvg_svg(svg_fl) 
mysvg <- readPicture(rawToChar(mysvg))
grid.picture(mysvg, 0.5, 0.15, width=0.5) # Coordinates are relative coordinates at the scale of the device and not the plot

dev.off()

# System
options("repos" = c(CRAN = "https://cloud.r-project.org/"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref="7ea9440", upgrade="never")

# Newer version of dplyr for rlang
install_version("Rcpp",       version = "1.0.0",   upgrade = "never")
install_version("pkgconfig",  version = "2.0.2",   upgrade = "never")
install_version("pillar",     version = "1.3.1",   upgrade = "never")
install_version("tibble",     version = "2.0.0",   upgrade = "never")
install_version("tidyselect", version = "0.2.5",   upgrade = "never")
install_version("dplyr",      version = "0.8.0",   upgrade = "never")

# Gviz
install_version("foreign",      version = "0.8-76", upgrade = "never")
install_version("latticeExtra", version = "0.6-28", upgrade = "never")
install_bioc("3.9/Gviz", upgrade="never")

# Sushi
install_version("zoo",    version = "1.8-9", upgrade = "never")
install_bioc("3.8/Sushi", upgrade="never")

# SVG
install_version("rsvg",      version = "2.1.2", upgrade = "never")
install_version("grImport2", version = "0.2-0", upgrade = "never")

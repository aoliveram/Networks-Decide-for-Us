# Install dependencies for Networks-Decide-for-Us

# Install CRAN packages
install.packages(c("igraph", "doParallel", "ggplot2", "dplyr", "readr", "patchwork", "intergraph", "cluster", "devtools"))

# Install netdiffuseR from the specific branch
devtools::install_github("USCCANA/netdiffuseR", ref = "stochastic-transmission")

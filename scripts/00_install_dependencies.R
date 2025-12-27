# Install dependencies for Networks-Decide-for-Us

# Install CRAN packages
install.packages(c("igraph", "doParallel", "ggplot2", "dplyr", "readr", "patchwork", "intergraph", "cluster", "devtools"), repos = "https://cloud.r-project.org")

# Install netdiffuseR from the specific branch
devtools::install_github("USCCANA/netdiffuseR", ref = "stochastic-transmission", force = TRUE)

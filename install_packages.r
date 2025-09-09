library(renv)
library(usethis)
install("gitcreds")
library(gitcreds)
all_packages <- c(
    'bioc::glmGamPoi',
    'immunogenomics/presto',
    'immunogenomics/harmony',
)

renv::install(packages = all_packages)
renv::snapshot()
here::i_am("scripts/install_packages.r")
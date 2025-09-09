chooseCRANmirror()
1
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE)
.libPaths()

install.packages("languageserver") 
install.packages("IRkernel")
system.file('kernelspec', package = 'IRkernel')
IRkernel::installspec()

# Install mamba:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh

conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n jupyter jupyter
conda activate jupyter  # activate our environment

IRkernel::installspec()
#


# from within R:

# install.packages("devtools")
# devtools::install_github("IRkernel/IRkernel")
# system.file('kernelspec', package = 'IRkernel')
# The last line should give you the location of Jupyter will need to find the kernel. Mine was /home/ubuntu/R/x86_64-pc-linux-gnu-library/4.0/IRkernel/kernelspec

# From the command line:

# inspect the path that you receive when you were in R. There should be a .json file in it.
# jupyter kernelspec list (run this to be sure that jupyter is in your path, you should see information about the current available kernels.
# jupyter kernelspec install /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/IRkernel/kernelspec --name 'R' --user 
#(you will use path that you received while working in R which could be different)
# jupyter kernelspec list (this list should now include R)
# restart jupyter


########################################

# GitHub Installation token if necessary to update
usethis::create_github_token()
usethis::edit_r_environ()

# Installation packages
install.packages("BiocManager")
install.packages("remotes")
install.packages('devtools')

# General R and plotting 
install.packages('here')
install.packages('tidyverse')
install.packages('ggplot2')
install.packages('scales')
remotes::install_github("thomasp85/patchwork")
install.packages('cowplot')
install.packages('gridExtra')
install.packages('ggrepel')
install.packages('stringr')
install.packages('VennDiagram')
install.packages('pheatmap')
install.packages('viridis')
install.packages("svglite")
install.packages('UpSetR')
install.packages('ggpubr')
install.packages('heatmap3')
install.packages("knitr")
install.packages("httpgd")
remotes::install_github("nx10/httpgd")

# Single Cell Analysis Packages
install.packages('Seurat')
install.packages('scRepertoire')
BiocManager::install("scRepertoire")
packageVersion("scRepertoire")
install.packages('circlize')
install.packages('scCustomize')
BiocManager::install('SingleR')
BiocManager::install('celldex')
BiocManager::install('UCell')
# devtools::install_github('scplotter')

# DEG, pathway enrichment and visualization packages
BiocManager::install('DESeq2')
BiocManager::install('clusterProfiler')
install.packages('DOSE')
BiocManager::install('pathview')
BiocManager::install('org.Mm.eg.db')
install.packages('enrichplot')
install.packages('msigdbr')
install.packages('gprofiler2')
devtools::install_github('immunogenomics/presto')


install.packages("shiny")

# Other/unknown packages

# options("download.file.method"="wget")
# options(BioC_mirror = "http://bioconductor.org")
install.packages(c("HGNChelper", "openxlsx", "data.tree"))
BiocManager::install(c("ComplexHeatmap"), force=TRUE)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages(c("fields"))
BiocManager::install("Nebulosa")
install.packages('tidymodels')
devtools::install_github('YuLab-SMU/ggtree')
BiocManager::install("apeglm")
a
install.packages('packrat')
install.packages('arrow')
install.packages('sf')
BiocManager::install('glmGamPoi')
a
remotes::install_github("bnprks/BPCells/r") ##Errors
install.packages('colorspace')
install.packages('ggridges')

# Installing cellhashR
BiocManager::install(c("DropletUtils", "demuxmix","SingleCellExperiment", "S4Vectors"))
a
#Latest version:
devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE)
#devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'seuratVersion', dependencies = TRUE)

#Cluster identification packages
install.packages("RSQLite")

#Annotation Hub
BiocManager::install("AnnotationHub", force=TRUE)
a
BiocManager::install("BiocFileCache", force=TRUE)
BiocManager::install("ensembldb", force=TRUE)


## Seurat Develop Version Installation 
lb = .libPaths()
.libPaths <- c('./package_version_folder/', lb)
.libPaths

new_path <- './package_version_folder/'
dir.create('./package_version_folder/')

old_paths <- .libPaths
old_paths

.libPaths(c('./package_version_folder/',lb))
.libPaths(new = c('./package_version_folder/'))

#.libPaths() <- c('./package_version_folder/', lb)
.libPaths <- 'a'

library(devtools)
devtools::install_github('satijalab/seurat@develop', args = c('--library="./package_version_folder/"'))
devtools::install_github('satijalab/seurat@develop', force = TRUE)
.libPaths(lb)
.libPaths(c(lb,'./package_version_folder/'))
.libPaths()


BiocManager::install("scrapper", force=TRUE)
a
BiocManager::install("immunarch", force=TRUE)
BiocManager::install("immunarch", force=TRUE)
remotes::install_github(c("BorchLab/immApex", "BorchLab/scRepertoire"), force = TRUE)

BiocManager::install("demuxmix", force=TRUE)
a
remotes::install_github('satijalab/azimuth', ref = 'master')
devtools::install_github('satijalab/seurat-data', 'seurat5')
1

# Install Signac and dependences
setRepositories(ind  = 1:3)
install.packages("Signac")
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

# Doublet finder
BiocManager::install("scDblFinder")
a
devtools::install_github('immunogenomics/harmony')
#renv
install.packages("renv")

# GIT Credentials
install.packages("gitcreds")

install.packages("usethis")

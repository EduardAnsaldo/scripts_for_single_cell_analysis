chooseCRANmirror()
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE)
Sys.setenv(GITHUB_PAT='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')
.libPaths()

install.packages("languageserver") 
install.packages("IRkernel")
IRkernel::installspec()
install.packages("shiny")

install.packages("Seurat")
Yes
install.packages("ggplot2")
# remotes::install_github("thomasp85/patchwork")

install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = FALSE))
    install.packages("BiocManager")
BiocManager::install("scRepertoire", force=TRUE)
packageVersion("scRepertoire")
#devtools::install_github("ncborcherding/scRepertoire")

#BiocManager::install(version='devel')
#BiocManager::install("scRepertoire")
.
# options("download.file.method"="wget")
# options(BioC_mirror = "http://bioconductor.org")

install.packages(c("HGNChelper", "openxlsx", "data.tree"))
# ???????????????????? install.packages("scater")

BiocManager::install(c("ComplexHeatmap", "DESeq2"), force=TRUE)
n
# install.packages("ROCR")

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

install.packages(c("fields", "Stat2Data", "tidyverse", "circlize","scCustomize", "remotes"))


#pkgbuild::check_build_tools(debug = TRUE)

# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
#local({options(repos = BiocManager::repositories())})


# Installing cellhashR
BiocManager::install(c("DropletUtils", "demuxmix","SingleCellExperiment", "S4Vectors"))
install.packages('devtools')
#Latest version:
devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE)
#devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'seuratVersion', dependencies = TRUE)

#remove.packages('Seurat')
#package <- 'https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.4.0.tar.gz'
#install.packages(package, repos=NULL, type='source')

#Cluster identification packages
install.packages("RSQLite")
BiocManager::install('msigdb', force=TRUE)
BiocManager::install('clusterProfiler', force=TRUE)

#Topic Modeling Packages
#install.packages("fastTopics")



#Annotation Hub
BiocManager::install("AnnotationHub", force=TRUE)
BiocManager::install("BiocFileCache", force=TRUE)
BiocManager::install("ensembldb", force=TRUE)

install.packages("pheatmap")

install.packages("svglite")

usethis::create_github_token()

usethis::edit_r_environ()

install.packages('devtools')
install.packages('remotes')

devtools::install_github('immunogenomics/presto')

remotes::install_github("thomasp85/patchwork")
# remotes::install_github('LTLA/scuttle')

BiocManager::install("Nebulosa")
Yes
BiocManager::install("SingleR")
BiocManager::install("celldex")

install.packages('tidymodels')

n

update.packages('ps')

BiocManager::install("pathview")
BiocManager::install("org.Mm.eg.db", force = T)

# devtools::install_github("YuLab-SMU/clusterProfiler")

install.packages('msigdbr')

BiocManager::install("UCell", force = T)

install.packages('VennDiagram')
install.packages('UpSetR')

devtools::install_github('YuLab-SMU/ggtree')

#gitcreds::gitcreds_set()

install.packages('ggpubr')
# install.packages('rstatix')

BiocManager::install("apeglm", force = T)

install.packages('gprofiler2')

#devtools::install_github('ncborcherding/scRepertoire')

devtools::install_github('pwwang/scplotter')
3

install.packages('packrat')

devtools::install_github('satijalab/seurat@develop')

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

install.packages('arrow')

install.packages('sf')

#install.packages('BiocManager')
BiocManager::install('glmGamPoi', force = TRUE)
a
Yes

remotes::install_github("bnprks/BPCells/r")

install.packages('colorspace')

install.packages('ggridges')

install.packages('heatmap3')

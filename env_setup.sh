## Create environment by running these commands in a terminal (depends on device and conda/mamba versions installed)
## This approach seems to work best for me
## The environment is called COSMX_RNA_env

## Not all packages are required, so only install what you need

CONDA_SUBDIR=osx-64 mamba create -n COSMX_RNA_env --yes \
  jupyter r-base r-essentials r-tidyverse r-irkernel r-patchwork r-seurat r-magick \
  r-rcolorbrewer r-viridis r-plotly radian r-languageserver r-devtools r-usethis r-gert libgit2 r-enrichr \
  r-tibble r-tidyr r-cluster r-clustree r-pheatmap r-ggplot2 \
  numpy scipy pandas seaborn matplotlib scikit-learn plotly tqdm \
  bioconductor-biocparallel \
  -c conda-forge -c bioconda

conda activate COSMX_RNA_env

mamba install -c conda-forge r-ggrastr

R
install.packages("harmony")
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
remotes::install_github("prabhakarlab/Banksy")
remotes::install_github("satijalab/seurat-wrappers")
q()


R 
# Make sure Matrix and irlba are both up to date (otherwise versioning issues cause prcomp_irlba to error out)
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
BiocManager::install("sparseMatrixStats")

# Install Insitutype:
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
devtools::install_github("https://github.com/Nanostring-Biostats/InSituCor")

install.packages('BiocManager')
BiocManager::install('glmGamPoi')

install.packages("remotes")
options(download.file.method = "curl")
options(timeout = 600)
remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
                         subdir = "_code/scPearsonPCA", ref = "Main")

remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", 
                        subdir = "_code/HieraType", ref = "Main")

BiocManager::install("ComplexHeatmap")
BiocManager::install("DESeq2")

install.packages("Cairo")
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

q()

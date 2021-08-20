<a href=""><img src="" alt="DOI"></a>
# Dynamic Landscape of Oral Carcinogenesis Reveals Multistage Malignant Transcriptional Programs
Understanding composition changes among tumor initiation and progression is critical to characterize the tumor microenvironment of oral cancer.
## Environment 
Ubuntu 9.3.0
R version 4.0.5	
Python version 3.8.10	

## Install software
### Install R package MAESTRO V1.4.1
    conda config --add channels defaults
    conda config --add channels liulab-dfci
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install mamba -c conda-forge
    mamba create -n MAESTRO maestro=1.4.1 -c liulab-dfci
### Install R package Seurat v2.3.4 	
    source("https://z.umn.edu/archived-seurat")
### Install R package DoubletFinder v2.0
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
### Install R package Monocle v2.8 	
    source("http://bioconductor.org/biocLite.R") 
    biocLite("monocle")	


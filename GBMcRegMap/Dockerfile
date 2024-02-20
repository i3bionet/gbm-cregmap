FROM rocker/r-ver:4.2.1

# Update env
RUN apt-get update && apt-get install -y  cmake gdal-bin git-core libcurl4-openssl-dev libgdal-dev libgeos-dev libgeos++-dev libgit2-dev libicu-dev libpng-dev libproj-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc python zlib1g-dev libxt-dev && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/

# Copy spatstat.core binaries and install
#RUN mkdir /spatstat
#COPY spatstat.core_2.4-4.tar.gz /spatstat/spatstat.core_2.4-4.tar.gz

# Copy BLCARegMap to Docker container
RUN mkdir /GBM_RegMap_pkg
COPY GBMRegMap_0.0.1.tar.gz /GBM_RegMap_pkg/GBMRegMap_0.0.1.tar.gz
COPY inst/ /inst

RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
#RUN Rscript -e 'remotes::install_bioc("depmap",upgrade="never", version = "1.14.0")'
RUN Rscript -e 'remotes::install_version("RcppTOML",upgrade="never", version = "0.1.7")'
RUN Rscript -e 'remotes::install_version("RColorBrewer",upgrade="never", version = "1.1-3")'
RUN Rscript -e 'remotes::install_version("stringr",upgrade="never", version = "1.4.1")'
RUN Rscript -e 'remotes::install_version("scales",upgrade="never", version = "1.2.1")'
RUN Rscript -e 'remotes::install_version("rlang",upgrade="never", version = "1.0.5")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.3.6")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.2")'
RUN Rscript -e 'remotes::install_version("tidyr",upgrade="never", version = "1.2.1")'
RUN Rscript -e 'remotes::install_version("purrr",upgrade="never", version = "0.3.4")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.0.10")'
RUN Rscript -e 'remotes::install_version("ggrepel",upgrade="never", version = "0.9.1")'
RUN Rscript -e 'remotes::install_version("gtools",upgrade="never", version = "3.9.3")'
RUN Rscript -e 'remotes::install_version("gplots",upgrade="never", version = "3.1.3")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.2")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.0")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("waiter",upgrade="never", version = "0.2.5")'
RUN Rscript -e 'remotes::install_version("visNetwork",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("survminer",upgrade="never", version = "0.4.9")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.7.3")'
RUN Rscript -e 'remotes::install_version("shinythemes",upgrade="never", version = "1.2.0")'
RUN Rscript -e 'remotes::install_version("shinybusy",upgrade="never", version = "0.3.1")'
## Install spatstat
#RUN Rscript -e 'utils::install.packages("/spatstat/spatstat.core_2.4-4.tar.gz", repos = NULL, force = TRUE)'
#RUN Rscript -e 'remotes::install_version("Seurat",upgrade="never", version = "4.1.1")'
#Specify Seurat version 4.3.0
RUN Rscript -e 'remotes::install_github("satijalab/seurat@v4.3.0")' 
RUN Rscript -e 'remotes::install_version("rlist",upgrade="never", version = "0.4.6.2")'
RUN Rscript -e 'remotes::install_bioc("limma",upgrade="never", version = "3.52.3")'
RUN Rscript -e 'remotes::install_version("leaflet",upgrade="never", version = "2.1.1")'
RUN Rscript -e 'remotes::install_bioc("iterClust",upgrade="never", version = "1.18.0")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.3")'
RUN Rscript -e 'remotes::install_version("e1071",upgrade="never", version = "1.7-11")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.25")'
#RUN Rscript -e 'remotes::install_bioc("depmap",upgrade="never", version = "1.10.0")'
RUN Rscript -e 'remotes::install_bioc("CoRegNet",upgrade="never", version = "1.34.0")'
RUN Rscript -e 'remotes::install_bioc("ComplexHeatmap",upgrade="never", version = "2.12.1")'
RUN Rscript -e 'remotes::install_version("formattable",upgrade="never", version = "0.2.1")'
RUN Rscript -e 'remotes::install_version("whereami",upgrade="never", version = "0.1.9")'
RUN Rscript -e 'remotes::install_bioc("TCGAbiolinks",upgrade="never", version = "2.24.3")'
RUN Rscript -e 'remotes::install_version("reticulate",upgrade="never", version = "1.26")'
RUN Rscript -e 'remotes::install_version("umap",upgrade="never", version = "0.2.9.0")'
RUN Rscript -e 'remotes::install_bioc("GEOquery",upgrade="never", version = "2.66.0")'
RUN Rscript -e 'remotes::install_version("reactable",upgrade="never", version = "0.4.4")'
## Install BLCARegmap
RUN Rscript -e 'utils::install.packages("/GBM_RegMap_pkg/GBMRegMap_0.0.1.tar.gz", repos = NULL, force = TRUE)'
#RUN mkdir /build_zone
#ADD . /build_zone
#WORKDIR /build_zone
#RUN R -e 'remotes::install_local(upgrade="never")'
#RUN rm -rf /build_zone
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');GBMRegMap::run_app()"
EXPOSE 80

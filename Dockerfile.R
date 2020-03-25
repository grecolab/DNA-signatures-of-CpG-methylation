
FROM openanalytics/r-base

MAINTAINER Tobias Verbeke "tobias.verbeke@openanalytics.eu"

# system libraries of general use
# RUN apt-get update && apt-get install -y \
#     sudo \
#     pandoc \
#     pandoc-citeproc \
#     libcurl4-gnutls-dev \
#     libcairo2-dev \
#     libxt-dev \
#     libssl-dev \
#     libssh2-1-dev \
#     libssl1.0.0

# system library dependency for the euler app
# RUN apt-get update && apt-get install -y \
#    libmpfr-dev

# basic shiny functionality
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'readr', 'DT', 'tibble', 'gplots',\
                              'dendextend', 'foreach', 'doParallel', 'XML', 'methods', 'BiocManager'), repos='https://cloud.r-project.org/')"

# install dependencies of the euler app
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19','IlluminaHumanMethylation450kanno.ilmn12.hg19'), version = '3.8')"

# copy the app to the image
RUN mkdir /root/euler
COPY euler /root/euler

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/euler')"]
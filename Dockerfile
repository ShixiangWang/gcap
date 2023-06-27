# FROM rocker/r-ver:4.3
# https://github.com/rocker-org/rocker-versioned2
# https://hub.docker.com/r/rocker/r-ver/tags

# https://rocker-project.org/use/extending.html#conda-forge
# https://hub.docker.com/r/continuumio/miniconda3
FROM continuumio/miniconda3:latest

LABEL \
    maintainer="Shixiang Wang" \
    email="wangsx1@sysucc.org.cn" \
    description="Docker Image for GCAP (Gene-level Circular Amplicon Prediction)" \
    org.label-schema.license="Non-Commercial Academic License (c) Shixiang Wang" \
    org.label-schema.vcs-url="https://github.com/ShixiangWang/gcap"

RUN apt update && apt install -y build-essential zip cmake libcairo2-dev &&\
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/* &&\
    conda install mamba -n base -c conda-forge &&\
    mamba create -n cancerit -c bioconda -c conda-forge cancerit-allelecount &&\
    mamba clean -yaf

# Install GCAP & deploy it
# XGBOOST should be <1.6
# The default path for conda in the container is /opt/conda
RUN mamba install -y -c conda-forge -c bioconda r-base=4.2 r-remotes r-biocmanager r-tidyverse r-sigminer sequenza-utils samtools tabix  &&\
    mamba clean -yaf &&\
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz", repos = NULL)' &&\
    R -e 'BiocManager::install("ShixiangWang/copynumber", update = FALSE, force = TRUE)' &&\
    R -e 'BiocManager::install("ShixiangWang/facets", update = FALSE, force = TRUE)' &&\
    R -e 'BiocManager::install("ShixiangWang/ascat@v3-for-gcap-v1", subdir = "ASCAT", dependencies = TRUE, update = FALSE)' &&\
    R -e 'BiocManager::install("ShixiangWang/gcap", dependencies = TRUE, update = FALSE)' &&\
    R -e 'gcap::deploy()' &&\
    cd /opt/conda/lib/R/library/facets/extcode/ &&\
    g++ -std=c++11 -I/opt/conda/include snp-pileup.cpp -L/opt/conda/lib -lhts -Wl,-rpath=/opt/conda/lib -o snp-pileup &&\
    rm -rf /tmp/downloaded_packages

# Deploy
WORKDIR $HOME
ENTRYPOINT [ "gcap" ]
CMD [ "--help" ]

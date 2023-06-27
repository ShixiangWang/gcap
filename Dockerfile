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

RUN apt update && apt install -y build-essential zip cmake &&\
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*
    conda install mamba -n base -c conda-forge &&\
    mamba create -n cancerit -c bioconda -c conda-forge cancerit-allelecount &&\
    mamba clean -yaf

# Install GCAP & deploy it
# XGBOOST should be <1.6
# The default path for conda in the container is /opt/conda
RUN mamba install -y -c conda-forge -c bioconda r-base=4.3 r-remotes r-biocmanager sequenza-utils samtools tabix  &&\
    mamba clean -yaf &&\
    R -e 'BiocManager::install("ShixiangWang/ascat@v3-for-gcap-v1", subdir = "ASCAT", dependencies = TRUE)' &&\
    R -e 'remotes::install_github("ShixiangWang/gcap", dependencies = TRUE)' &&\
    R -e 'gcap::deploy()' &&\
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz", repos = NULL)' &&\
    cd /data3/wsx/R/x86_64-pc-linux-gnu-library/4.2/facets/extcode/ &&\
    g++ -std=c++11 -I/opt/conda/include snp-pileup.cpp -L/opt/conda/lib -lhts -Wl,-rpath=/opt/conda/lib -o snp-pileup &&\
    rm -rf /tmp/downloaded_packages

# Deploy
WORKDIR $HOME
ENTRYPOINT [ "gcap" ]
CMD [ "--help" ]

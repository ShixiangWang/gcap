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

RUN conda install mamba -n base -c conda-forge &&\
    conda create -n cancerit -c bioconda -c conda-forge cancerit-allelecount &&\
    conda clean -yaf

# Install GCAP & deploy it
RUN mamba install -y -c conda-forge r-base=4.3 r-remotes &&\
    mamba clean -yaf &&\
    R -e 'remotes::install_github("ShixiangWang/gcap@v1.1.2", dependencies = TRUE)' &&\
    R -e 'gcap::deploy()' &&\
    rm -rf /tmp/downloaded_packages

# Install 
WORKDIR $HOME
ENTRYPOINT [ "gcap" ]
CMD [ "--help" ]
  

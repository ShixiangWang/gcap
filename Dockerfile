FROM continuumio/miniconda3:latest

# Avoid Docker overlay2 hardlink issues with conda package extraction.
# Without this, mamba produces "Cannot find a valid extracted directory cache"
# for random packages during install. See: conda/conda#12156, conda/conda#11314
ENV CONDA_ALWAYS_COPY=true

LABEL \
    maintainer="Shixiang Wang" \
    email="wangsx1@sysucc.org.cn" \
    description="Docker Image for GCAP (Gene-level Circular Amplicon Prediction)" \
    org.label-schema.license="Non-Commercial Academic License (c) Shixiang Wang, Qi Zhao; All commercial use is strictly prohibited and requires a commercial use licence." \
    org.label-schema.vcs-url="https://github.com/ShixiangWang/gcap" \
    org.opencontainers.image.source="https://github.com/ShixiangWang/gcap" \
    org.opencontainers.image.description="Docker Image for GCAP (Gene-level Circular Amplicon Prediction)"

RUN apt update && apt install -y build-essential zip cmake libcairo2-dev &&\
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/* &&\
    conda install mamba -n base -c conda-forge -y &&\
    mamba clean -yaf

# Install GCAP & deploy it
# XGBOOST should be <1.6
# The default path for conda in the container is /opt/conda
RUN mamba install -y -c conda-forge -c bioconda \
        r-base=4.3 python=3.10 \
        r-remotes r-biocmanager r-tidyverse r-sigminer \
        sequenza-utils samtools tabix cancerit-allelecount &&\
    mamba clean -yaf &&\
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz", repos = NULL)' &&\
    R -e 'BiocManager::install("jokergoo/GetoptLong", update = FALSE, force = TRUE)' &&\
    R -e 'BiocManager::install("ShixiangWang/copynumber", update = FALSE, force = TRUE)' &&\
    R -e 'BiocManager::install("ShixiangWang/facets", update = FALSE, force = TRUE)' &&\
    R -e 'BiocManager::install("VanLoo-lab/ascat", subdir = "ASCAT", dependencies = TRUE, update = FALSE)' &&\
    R -e 'BiocManager::install("ShixiangWang/gcap", dependencies = TRUE, update = FALSE)' &&\
    cd /opt/conda/lib/R/library/facets/extcode/ &&\
    g++ -std=c++11 -I/opt/conda/include snp-pileup.cpp -L/opt/conda/lib -lhts -Wl,-rpath=/opt/conda/lib -o snp-pileup &&\
    ls -alh &&\
    R -e 'gcap::deploy()' &&\
    rm -rf /tmp/downloaded_packages &&\
    ls -alh ~ &&\
    echo "Deploy success!"

# Entrypoint wrapper: activate conda base environment before running gcap
# The exec form ENTRYPOINT does not go through bash, so conda is not auto-activated.
# This wrapper sources conda first, then runs gcap with all arguments forwarded.
RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh &&\
    echo '. /opt/conda/etc/profile.d/conda.sh && conda activate base' >> /usr/local/bin/entrypoint.sh &&\
    echo 'exec gcap "$@"' >> /usr/local/bin/entrypoint.sh &&\
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]

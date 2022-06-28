FROM ubuntu:22.04

RUN apt-get update && apt-get install -y wget \
    build-essential \
    clang \
    git \ 
    gdb \
    gfortran \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get autoremove -y

ENV PATH="/usr/local/miniconda3/bin:${PATH}"
ARG PATH="/usr/local/miniconda3/bin:${PATH}"
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /usr/local/miniconda3 -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

RUN conda install -c conda-forge -c intel cmake mkl-static mkl-devel intel-openmp

ENV LD_LIBRARY_PATH /usr/local/miniconda3/lib:$LD_LIBRARY_PATH
ARG LD_LIBRARY_PATH /usr/local/miniconda3/lib:$LD_LIBRARY_PATH

ENTRYPOINT bash
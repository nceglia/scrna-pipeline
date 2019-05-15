FROM nceglia/cellranger-base:latest AS cellranger
FROM ubuntu:18.04

MAINTAINER Nicholas Ceglia <nickceglia@gmail.com>

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev libffi-dev build-essential software-properties-common  \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip


RUN python3 -m pip install tensorflow
RUN python3 -m pip install PyYAML
RUN python3 -m pip install matplotlib
RUN python3 -m pip install pandas
RUN python3 -m pip install rpy2
RUN python3 -m pip install networkx
RUN python3 -m pip install dill
RUN python3 -m pip install azure
RUN python3 -m pip install scanorama
RUN python3 -m pip install scanpy
RUN mkdir /results
RUN mkdir /runs
RUN apt-get install -y git
RUN git clone https://github.com/shahcompbio/scvis.git /codebase/scvis/
RUN cd /codebase/scvis && python3 setup.py install
RUN git clone https://github.com/nceglia/scrna-pipeline.git /codebase/SCRNApipeline/
RUN git clone https://github.com/shahcompbio/pypeliner.git /codebase/pypeliner/
RUN cd /codebase/pypeliner && python3 setup.py install
ENV GITHUB_PAT=b4e0f61f39f253efd3a8d6358528db8bd775b6a2
RUN apt-get clean && apt-get update && apt-get install -y locales
RUN locale-gen en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y r-base
RUN apt-get install -y libcurl4-openssl-dev
RUN R -e "install.packages('httr')"
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('httr')"
RUN R -e "BiocManager::install('scran')"
RUN R -e "BiocManager::install('scater')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "BiocManager::install('EnsDb.Hsapiens.v86')"
RUN R -e "devtools::install_github('Irrationone/cellassign')"
RUN R -e "devtools::install_github('kieranrcampbell/clonealign')"


COPY --from=cellranger /codebase/cellranger-3.0.2 /codebase/cellranger-3.0.2
COPY lsf.template /codebase/cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/lsf.template
CMD ["/bin/bash"]

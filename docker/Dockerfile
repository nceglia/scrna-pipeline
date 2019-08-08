FROM nceglia/cellranger-base:latest AS cellranger
FROM ubuntu:18.04

MAINTAINER Nicholas Ceglia <nickceglia@gmail.com>

RUN apt-get clean && apt-get update && apt-get install -y locales
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev libffi-dev build-essential software-properties-common libcurl4-openssl-dev libssl-dev apt-transport-https r-base  \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

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
RUN python3 -m pip install tzlocal
RUN python3 -m pip install --user scipy==1.2.1

RUN apt-get install -y git
RUN git clone https://github.com/shahcompbio/scvis.git /codebase/scvis/
RUN cd /codebase/scvis && python3 setup.py install
RUN git clone https://github.com/nceglia/scrna-pipeline.git /codebase/SCRNApipeline/
RUN git clone https://github.com/shahcompbio/pypeliner.git /codebase/pypeliner/
RUN cd /codebase/pypeliner && python3 setup.py install

RUN apt-get install -y r-base

RUN R -e "install.packages('httr')"
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('Seurat')"

RUN R -e "BiocManager::install('scran')"
RUN R -e "BiocManager::install('scater')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "BiocManager::install('EnsDb.Hsapiens.v86')"

RUN R -e "install.packages('XML')"
RUN R -e "install.packages('tensorflow')"
RUN git clone https://github.com/Irrationone/cellassign.git /codebase/cellassign/
RUN cd /codebase; R CMD build --no-build-vignettes cellassign; R CMD INSTALL cellassign_0.99.2.tar.gz
RUN git clone https://github.com/kieranrcampbell/clonealign.git /codebase/clonealign
RUN cd /codebase; R CMD build --no-build-vignettes clonealign; R CMD INSTALL clonealign_0.99.0.tar.gz

RUN mkdir /results
RUN mkdir /runs

COPY --from=cellranger /codebase/cellranger-3.0.2 /codebase/cellranger-3.0.2
COPY lsf.template /codebase/cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/lsf.template
RUN cd /codebase/SCRNApipeline/; git pull;
RUN python3 -m pip install --user scipy==1.2.1
RUN cd /codebase/SCRNApipeline/; git pull;
CMD ["/bin/bash"]

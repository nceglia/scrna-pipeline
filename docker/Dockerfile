FROM rocker/verse:latest

MAINTAINER Nicholas Ceglia <nickceglia@gmail.com>

RUN apt-get clean && apt-get update && apt-get install -y locales && apt-get install -y libhdf5-dev
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y install python3-dev
RUN mkdir /codebase

RUN apt-get install -y git
RUN apt-get install -y python3-setuptools
RUN git clone https://github.com/nceglia/scrna-pipeline.git /codebase/scrna-pipeline/
RUN git clone https://github.com/shahcompbio/pypeliner.git /codebase/pypeliner/
RUN cd /codebase/pypeliner && python3 setup.py install

RUN apt-get install -y curl
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python3 get-pip.py

RUN python3 -m pip install PyYAML
RUN python3 -m pip install matplotlib
RUN python3 -m pip install pandas
RUN python3 -m pip install rpy2
RUN python3 -m pip install networkx
RUN python3 -m pip install dill
RUN python3 -m pip install scanpy
RUN python3 -m pip install tzlocal
RUN python3 -m pip install scipy

ENV PYTHONPATH=/codebase/scrna-pipeline/

COPY hrd.yaml /codebase/hrd.yaml
COPY exhaustion.yaml /codebase/exhaustion.yaml
COPY markers.yaml /codebase/markers.yaml
COPY markers_positive.yaml /codebase/markers_positive.yaml
COPY markers_negative.yaml /codebase/markers_negative.yaml
COPY run /bin/run
COPY run_integration /bin/run_integration

RUN chmod 777 /bin/run
RUN chmod 777 /bin/run_integration

RUN python3 -m pip install numpy

RUN R -e "install.packages('tensorflow')"

RUN apt-get install -y python3-venv
RUN /root/.virtualenvs/r-reticulate/bin/python /get-pip.py
RUN R -e "BiocManager::install('multtest',suppressUpdates=TRUE)"
RUN R -e "install.packages('mutoss')"
RUN R -e "install.packages('Seurat')"
RUN R -e "BiocManager::install('scater')"
RUN R -e "BiocManager::install('scran')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "tensorflow::install_tensorflow(extra_packages='tensorflow-probability')"
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('Irrationone/cellassign',dep=FALSE)"
RUN R -e "BiocManager::install('DropletUtils')"

COPY hrd.yaml /codebase/hrd.yaml
RUN rm -rf /codebase/scrna-pipeline/
RUN git clone https://github.com/nceglia/scrna-pipeline.git /codebase/scrna-pipeline/
ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib/"
CMD ["/bin/bash"]

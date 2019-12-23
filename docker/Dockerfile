FROM nceglia/base-scrna-r:latest

MAINTAINER Nicholas Ceglia <nickceglia@gmail.com>

RUN apt-get clean && apt-get update && apt-get install -y locales
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN mkdir /codebase

RUN apt-get install -y git
RUN git clone https://github.com/nceglia/scrna-pipeline.git /codebase/scrna-pipeline/
RUN git clone https://github.com/shahcompbio/pypeliner.git /codebase/pypeliner/
RUN cd /codebase/pypeliner && python3 setup.py install

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

CMD ["/bin/bash"]

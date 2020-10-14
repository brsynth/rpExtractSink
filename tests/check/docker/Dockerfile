ARG IMAGE
FROM ${IMAGE}

COPY tests/check/environment.yml environment-test.yml
COPY extras/environment-run.yml .

RUN conda update -n base -c defaults conda
RUN conda install -y -c conda-forge conda-merge
RUN conda-merge environment-run.yml environment-test.yml > environment.yml
RUN conda env create -f environment.yml

ARG HOME
WORKDIR ${HOME}/tests

ADD tests .

ARG PKG
ADD ${PKG} ${PKG}

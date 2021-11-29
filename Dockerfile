FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

WORKDIR "/analysis"

ENV TMPDIR="/scratch"
RUN mkdir -p $TMPDIR

COPY environment.yml .

RUN conda config --add channels bioconda && conda config --add channels conda-forge && \
    conda install -c conda-forge -n base mamba && \
    mamba env update -f environment.yml -n base && \
    conda clean -a

RUN wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz && \
    gunzip usearch11.0.667_i86linux32.gz && chmod +x usearch11.0.667_i86linux32 && \
    mv usearch11.0.667_i86linux32 $CONDA_PREFIX/bin/usearch

COPY src src
COPY setup.py .
COPY README.md .

RUN python -m pip install .

ENTRYPOINT ["amplicon_benchmark"]
CMD ["-h"]
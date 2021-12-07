# Amplicon data benchmarking

A python package to benchmark taxonomic assignment tools for amplicon databases

## Description
This tool was designed to split a fasta file into `test` and `train` datasets 
in order to evaluate the performance of various taxonomic assigners. 

## Installing

### From GitHub
1. Clone the repo: `git clone https://github.com/johnne/amplicon_benchmark`
2. Create the environment: `mamba env create -f environment.yml`
3. Activate the environment: `conda activate amplicon_benchmark`
4. Install the package: `python -m pip install --use-feature=in-tree-build .`

### From DockerHub
1. `docker pull johnne/amplicon_benchmark:latest`

#### With Singularity
1. `singularity build amplicon_benchmark.sif docker://johnne/amplicon_benchmark:latest`

## Searching for PCR products

To use the `search_pcr` function you must have `usearch` installed in your `$PATH`.

**Note**: If you're using the docker image then `usearch` is already installed and you can 
skip the steps below.

Download your version from [this page](https://drive5.com/usearch/download.html).

Example on MAC:
1. `curl https://drive5.com/downloads/usearch11.0.667_i86osx32.gz -o usearch.gz`
2. `gunzip -c usearch.gz > $CONDA_PREFIX/bin/usearch`
3. `chmod +x $CONDA_PREFIX/bin/usearch`
4. `rm usearch.gz`
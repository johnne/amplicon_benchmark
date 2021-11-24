# Amplicon data benchmarking

A python package to benchmark taxonomic assignment tools for amplicon databases

## Installing

1. Clone the repo: `git clone https://github.com/johnne/amplicon_benchmark`
2. Create the environment: `mamba env create -f environment.yml`
3. Activate the environment: `conda activate amplicon_benchmark`
4. Install the package: `python -m pip install .`

## Searching for PCR products

To use the `search_pcr` function you must have `usearch` installed in your `$PATH`.
Download your version from [this page](https://drive5.com/usearch/download.html).

Example on MAC:
1. `curl https://drive5.com/downloads/usearch11.0.667_i86osx32.gz -o usearch.gz`
2. `gunzip -c usearch.gz > $CONDA_PREFIX/bin/usearch`
3. `chmod +x $CONDA_PREFIX/bin/usearch`
4. `rm usearch.gz`
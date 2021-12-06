from importlib_resources import files, as_file
import yaml
from argparse import ArgumentParser
import logging
import sys
import os
import amplicon_benchmark

from snakemake import snakemake
from snakemake.utils import available_cpu_count
from snakemake.utils import validate


def populate_config(schema):
    """
    Initializes an empty config dictionary with default values taken
    from config schema.
    
    :param schema: JSON schema with default settings 
    :return: dict with config settings
    """
    config = {}
    validate(config, schema)
    return config


def update_config(config, args):
    """
    Updates the config dictionary with settings given at the CLI
    
    :param config: config dictionary 
    :param args: arguments from command line execution
    :return: config dict with values updated
    """
    if args.primers:
        config["search_pcr"]["primers"] = args.primers
        config["search_pcr"]["run_search"] = True
    else:
        config["search_pcr"]["run_search"] = False
    if args.fastafile:
        config["fastafile"] = args.fastafile
    if args.infofile:
        config["infofile"] = args.infofile
    if args.minamp:
        config["search_pcr"]["minamp"] = args.minamp
    if args.maxamp:
        config["search"]["maxamp"] = args.maxamp
    if args.threads:
        config["threads"] = args.threads
    return config


class SnakemakeError(Exception):
    pass


def run(args):
    """
    Runs the workflow with arguments
    
    :param args: Arguments from command line execution 
    :return: 
    """
    # Read default settings from schema
    schema_source = files(amplicon_benchmark).joinpath("config.schema.yaml")
    with as_file(schema_source) as schema_path:
        config = populate_config(schema_path)
    # Update config
    config = update_config(config, args)
    # Write to tempfile
    configfile = f".config.{os.getpid()}.yml"
    with open(configfile, 'w') as fh:
        yaml.dump(config, fh)
    
    snakefile_source = files(amplicon_benchmark).joinpath("Snakefile")
    with as_file(snakefile_source) as snakefile_path:
        forcerun = []
        if args.force:
            forcerun = args.targets
        success = snakemake(
            snakefile_path, targets=args.targets, dryrun=args.dryrun,
            cores=args.cores, configfiles=[configfile],
            cluster_config=args.cluster_config, workdir=args.workdir,
            printshellcmds=args.printshellcmds, unlock=args.unlock,
            forcerun=forcerun,
        )
        os.remove(configfile)
        return success


def main():
    parser = ArgumentParser()
    prog_opts = parser.add_argument_group("program-options")
    prog_opts.add_argument("-f", "--fastafile", type=str, required=True,
                        help="Fasta file with reference sequences to generate"
                             "training and test data from")
    prog_opts.add_argument("-i", "--infofile", type=str, required=True,
                        help="Tab-delimited file with taxonomic information for "
                             "records in fastafile")
    prog_opts.add_argument("--primers", type=str,
                        help="If specified, the usearch -search_pcr function will "
                             "run to identify amplicon sequences in the database")
    prog_opts.add_argument("--minamp", type=int,
                           help="Minimum amplicon length when running search_pcr"
                                " (including primer)")
    prog_opts.add_argument("--maxamp", type=int,
                           help="Maximum amplicon length when running search_pcr"
                                " (including primer")
    prog_opts.add_argument("--threads", type=int,
                           help="Max threads to run with")
    snakemake_opts = parser.add_argument_group("snakemake-options")
    snakemake_opts.add_argument("targets", nargs='*', default=[],
                        help="File(s) to create or steps to run. If omitted, "
                             "the full pipeline is run.")
    snakemake_opts.add_argument("-n", "--dryrun", action="store_true",
                        help="Only print what to do, don't do anything [False]")
    snakemake_opts.add_argument("-j", "--cores", type=int, default=4,
                        help="Number of cores to run with [4]")
    snakemake_opts.add_argument("-F", "--force", action="store_true",
                        help="Force workflow run")
    snakemake_opts.add_argument("-u", "--unlock", action="store_true",
                        help="Unlock working directory")
    snakemake_opts.add_argument("--cluster-config", type=str,
                        help="Path to cluster config (for running on SLURM)")
    snakemake_opts.add_argument("--workdir", type=str,
                        help="Working directory. Defaults to current dir")
    snakemake_opts.add_argument("-p", "--printshellcmds", action="store_true",
                        help="Print shell commands")
    args = parser.parse_args()
    
    # Run the workflow
    success = run(args)
    if not success:
        raise SnakemakeError()
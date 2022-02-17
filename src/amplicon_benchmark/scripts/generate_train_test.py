#!/usr/bin/env python

from Bio.SeqIO import parse, write as write_fasta
import gzip as gz
import tqdm
import random
import pandas as pd
from argparse import ArgumentParser
import sys
import re


def restrict_to_ids(info, ids):
    # Restrict to ids if given
    if len(ids) > 0:
        i = info.shape[0]
        intersect_ids = set(list(info.index)).intersection(set(ids))
        info = info.loc[intersect_ids]
        sys.stderr.write(
            f"Restricted sampling to {info.shape[0]}/{i} records\n")
    return info


def sample_and_remove_from_taxrank(info, taxrank="bold_id", k=250, ids=None):
    """
    CASE 1
    - Select k genera with at least 2 taxa
    - for one of these taxa select 1 sequence and remove any other
    sequences of the taxrank from the database

    Expectation: Test sequences should be classified at genus level –
    harder since no sequences for same taxrank in db

    :param info: Dataframe of taxonomic information for sequences
    :param taxrank: Taxonomic rank at which to pick sequences
    :param k: Number of genera to randomly select
    :param ids: List of read ids to restrict sampling from
    :return: tuple of lists with picked and excluded sequences
    """
    if ids is None:
        ids = []
    excluded = []
    picked = []
    _info = info.copy()
    # Restrict to ids if given
    info = restrict_to_ids(info, ids)
    # Get genera with at least 2 taxa
    taxa_per_genera = info.groupby(["genus"]).nunique().loc[:, taxrank]
    genera = list(taxa_per_genera.loc[taxa_per_genera >= 2].index)
    sys.stderr.write(f"Found {len(genera)} genera with >= 2 taxa\n")
    # Sample k genera randomly
    sampled_g = random.sample(list(genera), k=min(k, len(list(genera))))
    # For each sampled genus
    for g in tqdm.tqdm(sampled_g, unit=" genera",
                       desc="sampling from genera (CASE 1/3)"):
        # Get taxa for genus
        g_taxa = info.loc[info.genus == g, taxrank].unique()
        # Pick one taxa randomly
        picked_taxa = random.sample(list(g_taxa), k=1)
        # Pick one sequence randomly
        seqs = list(info.loc[info[taxrank] == picked_taxa[0]].index)
        picked_seq = random.sample(seqs, k=1)
        # Exclude the rest
        exclude_seq = seqs
        # Also exclude from original, unrestricted info file
        exclude_seq += list(_info.loc[_info[taxrank] == picked_taxa[0]].index)
        picked += picked_seq
        excluded += exclude_seq
    sys.stderr.write(f"Sampled: {len(picked)}, Excluded: {len(excluded)}\n")
    return picked, excluded


def sample_and_keep_from_taxrank(info, taxrank="bold_id", k=250, ids=None):
    """
    CASE 2
    - Select k genera with at least 1 taxrank with at least 2 sequences,
    - select 1 of these sequences and remove it from the database

    Expectation:
    Test sequences should be classified at genus level – easy since
    other sequences for same taxrank in database

    :param info: Dataframe of taxonomic information for sequences
    :param taxrank: Taxonomic rank at which to pick sequences
    :param k: Number of genera to randomly select
    :param ids: List of read ids to restrict sampling from
    :return tuple of lists with picked and excluded sequences
    """
    # Restrict to ids if given
    if ids is None:
        ids = []
    info = restrict_to_ids(info, ids)
    # Count number of sequences in each genus/taxrank
    seqs_per_taxrank = info.groupby(taxrank).count().iloc[:, 1]
    # Filter to taxa with at least 2 sequences
    filtered_taxa = list(seqs_per_taxrank.loc[seqs_per_taxrank >= 2].index)
    # Filter to genera matching filtered_taxa
    filtered_genera = info.loc[info[taxrank].isin(filtered_taxa), "genus"]
    # Remove NAs and make unique
    genera = list(
        set(filtered_genera.loc[filtered_genera == filtered_genera].values))
    sys.stderr.write(f"Found {len(genera)} genera with >=2 seqs at rank:{taxrank}\n")
    # Sample k random genera from filtered
    sampled_g = random.sample(list(genera), k=min(k, len(list(genera))))
    picked = []
    excluded = []
    for g in tqdm.tqdm(sampled_g, desc="sampling from genera (CASE 2/3)",
                       unit=" genera"):
        seqs = list(info.loc[(info.genus == g) & (
            info[taxrank].isin(filtered_taxa)), taxrank].index)
        picked_seq = random.sample(seqs, k=1)
        picked += picked_seq
        excluded += picked_seq
    sys.stderr.write(f"Sampled: {len(picked)}, Excluded: {len(excluded)}\n")
    return picked, excluded


def sample_genera_exclusively(info, k=500, ids=None):
    """
    CASE 3
    - For k genera, select one random sequence each
    - remove all other sequences for the genera from the database

    Expectation: Test sequences should be classified as unknown at genus
    level

    :param info: Dataframe of taxonomic information for sequences
    :param k: Number of genera to randomly select
    :param ids: List of read ids to restrict sampling from
    :return tuple of lists with picked and excluded sequences
    """
    if ids is None:
        ids = []
    _info = info.copy()
    # Restrict to ids if given
    info = restrict_to_ids(info, ids)
    genera = info.loc[info.genus == info.genus].genus.unique()
    sys.stderr.write(f"Found {len(genera)} remaining genera\n")
    sampled_g = random.sample(list(genera), k=min(k, len(list(genera))))
    picked = []
    excluded = []
    for g in tqdm.tqdm(sampled_g, desc="sampling from genera (CASE 3/3)",
                       unit=" genera"):
        # get sequences in genus
        seqs = info.loc[info.genus == g]
        #TODO: Check whether families are still included
        #families = info.loc[info.genus == g, "family"].unique()
        # pick one sequence randomly
        picked_seq = random.sample(list(seqs.index), k=1)
        # add all sequences for this genera to the excluded list
        excluded += list(set(_info.loc[_info.genus==g].index))
        picked += picked_seq
    sys.stderr.write(f"Sampled: {len(picked)}, Excluded: {len(excluded)}\n")
    return picked, excluded


def ranks_unclassified(desc):
    """
    Returns number of unclassified ranks in description
    :param desc: 
    :return:
    """
    regex = re.compile(".+(_[X]+)$")
    return sum(
        [1 if x else 0 for x in [regex.match(s) for s in desc.split(";")]])


def read_records(f, no_unclassified, strip_string="centroid="):
    """
    Reads fasta file with sequences

    :param f: Fasta file
    :param strip_string: string to be removed from each record id
    :return: Dictionary of records
    """
    records = {}
    openfn = open
    skipped = 0
    if f.endswith(".gz"):
        openfn = gz.open
    with openfn(f, 'rt') as fhin:
        for record in tqdm.tqdm(parse(fhin, "fasta"),
                                desc=f"Reading fasta file {f}", unit=" records"):
            i = (record.id).replace(strip_string, "")
            if no_unclassified:
                if ranks_unclassified(record.description) > 1:
                    skipped+=1
                    continue
            records[i] = record
    if skipped > 0:
        sys.stderr.write(f"Excluded {skipped} records with unclassified ranks\n")
    return records


def read_info(f, read_ids):
    """
    Read info file with taxonomic information

    :param f: Input file
    :param read_ids: Read ids to store
    :return: Dataframe with stored info
    """
    sys.stderr.write("Reading info file\n")
    df = pd.read_csv(f, sep="\t", header=0, index_col=0)
    return df.loc[read_ids]


def write_records(records, info, fhrecs, fhinfo, tag="", header=True):
    """
    Writes fasta and info records to files

    :param records: Dictionary of records to write
    :param info: Dataframe with info for records
    :param fhrecs: File handle for writing records
    :param fhinfo: File handle for writing info
    :param tag: Prefix
    :return:
    """
    write_fasta(records, fhrecs, "fasta")
    if tag != "":
        info = info.assign(
            tag=pd.Series([tag] * info.shape[0], index=info.index))
    info.to_csv(fhinfo, sep="\t", header=header)


def write_output(fastafile, infofile, records, info, case1, case2, case3):
    with open(fastafile, 'w') as fhrecs, open(infofile, 'w') as fhinfo:
        pass
    with open(fastafile, 'a') as fhrecs, open(infofile, 'a') as fhinfo:
        write_records([records[x] for x in case1], info.loc[case1], fhrecs,
                      fhinfo, "CASE1", header=True)
        write_records([records[x] for x in case2], info.loc[case2], fhrecs,
                      fhinfo, "CASE2", header=False)
        write_records([records[x] for x in case3], info.loc[case3], fhrecs,
                      fhinfo, "CASE3", header=False)


def read_ids(f, strip_string="centroid="):
    l = []
    with open(f, 'r') as fhin:
        for line in fhin:
            l.append(line.rstrip().replace(strip_string, ""))
    return l


def main():
    parser = ArgumentParser()
    parser.add_argument("fasta", type=str, help="Sequence file")
    parser.add_argument("info", type=str, help="Info file")
    parser.add_argument("--test_seqs",
                        help="Supply a separate file fasta file with "
                             "sequences that should be written for the"
                             " test data")
    parser.add_argument("--testfasta", default="test.fasta", type=str,
                        help="Test fasta file")
    parser.add_argument("--testinfo", default="test_info.tsv", type=str,
                        help="Test info file")
    parser.add_argument("--trainfasta", default="train.fasta", type=str,
                        help="Train fasta file")
    parser.add_argument("--traininfo", default="train_info.tsv", type=str,
                        help="Train info file")
    parser.add_argument("--taxrank", type=str, default="bold_id",
                        help="Taxonomic rank below genus to sample based on")
    parser.add_argument("--num_genera1", type=int, default=250,
                        help="Number of genera to sample from in case1")
    parser.add_argument("--num_genera2", type=int, default=250,
                        help="Number of genera to sample from in case2")
    parser.add_argument("--num_genera3", type=int, default=500,
                        help="Number of genera to sample from in case3")
    parser.add_argument("--no-unclassified", action="store_true",
                        help="Only use sequences classified at all ranks, "
                             "i.e. no ranks should end in '_X*'")
    args = parser.parse_args()

    records = read_records(args.fasta, args.no_unclassified)
    # Read test records from separate file if specified and different
    # from the fastafile
    if args.test_seqs != args.fasta:
        test_records = read_records(args.test_seqs, args.no_unclassified)
    else:
        test_records = records
        
    # Read info for sequences, limit to records given in fastafile 
    info = read_info(args.info, list(records.keys()))
    restrict_ids = list(test_records.keys())
    
    # First pick according to CASE1
    case1_picked, case1_excluded = sample_and_remove_from_taxrank(
        info, taxrank=args.taxrank, k=args.num_genera1, ids=restrict_ids)
    sys.stderr.write("\n")
    
    # Then pick according to CASE2
    case2_picked, case2_excluded = sample_and_keep_from_taxrank(
        info.drop(case1_excluded), taxrank=args.taxrank, k=args.num_genera2,
        ids=restrict_ids)
    sys.stderr.write("\n")
    
    # Then pick according to CASE3
    case3_picked, case3_excluded = sample_genera_exclusively(
        info.drop(case1_excluded + case2_excluded), k=args.num_genera3,
        ids=restrict_ids)
    sys.stderr.write("\n")
    
    # Write records for each case
    sys.stderr.write(f"Writing test sequences to {args.testfasta}\n")
    write_output(args.testfasta, args.testinfo, test_records, info, case1_picked,
                 case2_picked, case3_picked)
    
    # Write remaining records to train data
    exclude = case1_picked+case1_excluded+case2_picked+case2_excluded+case3_picked+case3_excluded
    sys.stderr.write(f"{len(exclude)} records excluded from train set\n")
    include = list(set(records.keys()).difference(set(exclude)))
    sys.stderr.write(f"Writing {len(include)} records to train fasta\n")
    with open(args.trainfasta, 'w') as fhrecs, open(args.traininfo, 'w') as fhinfo:
        write_records([records[x] for x in include],
                      info.loc[include], fhrecs, fhinfo, "TRAIN")


if __name__ == "__main__":
    main()

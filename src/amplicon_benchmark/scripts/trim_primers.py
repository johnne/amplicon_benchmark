#!/usr/bin/env python

from Bio.SeqIO import parse
import tqdm
import pandas as pd
from multiprocessing import Pool
from argparse import ArgumentParser


def primer_trim(row):
    query = row[0]
    seq = row[1]["amplicon_seq"]
    start = row[1]["primer1_aln_len"]
    stop = row[1]["primer2_aln_len"]
    trimmed_seq = seq[start:-stop]
    record = f">{query}\n{trimmed_seq}\n"
    return record


def main(args):
    # Read primers
    primers = []
    for record in parse(args.primers, 'fasta'):
        primers.append((record.id, len(record.seq)))            
    df = pd.read_csv(args.hits, sep="\t", header=None, index_col=0, 
        usecols=[0, 4, 5, 6, 7, 8, 9, 10, 11],
        names=["query", "primer1", "primer1_strand", "primer1_aln", 
               "primer2", "primer2_strand", "primer2_aln", "amplicon_len",
               "amplicon_seq"])
    # Filter hits to those that have primer1 and primer2 matches
    df = df.loc[(df.primer1 == primers[0][0]) & (df.primer2 == primers[1][0])]
    df.drop(["primer1","primer2", "primer1_strand", "primer2_strand", "amplicon_len"], inplace=True, axis=1)
    # Get lengths
    df["primer1_aln_len"] = df["primer1_aln"].str.len()
    df["primer2_aln_len"] = df["primer2_aln"].str.len()
    # Remove rows with aligned primer lengths that differ from primers
    df = df.loc[(df["primer1_aln_len"]==primers[0][1]) & (df["primer2_aln_len"]==primers[1][1])]
    df.drop(["primer1_aln","primer2_aln"], axis=1, inplace=True)
    with Pool(processes=args.threads) as pool:
        trimmed_records = list(tqdm.tqdm(pool.imap(primer_trim, df.iterrows()),
            unit=" seqs", desc=f"trimming primers with {args.threads} threads",
            total=df.shape[0], ncols=100))
    with open(args.outfile, 'w') as fhout:
        for r in trimmed_records:
            fhout.write(r)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("hits", help="Hits output from search_pcr function")
    parser.add_argument("primers", help="File with primers used in search_pcr")
    parser.add_argument("outfile", help="File with trimmed sequences")
    parser.add_argument("--threads", type=int, help="Threads to run with", 
                        default=1)
    args = parser.parse_args()
    main(args)
import argparse
from Bio import SeqIO
import pyrodigal
import gzip
import os
import collections
import logging

import pyhmmer.easel

HMM_FILE = "combined_marker_hmms_latest_version.HMM"

def get_options():

    parser = argparse.ArgumentParser(description='Get marker genes from HMMER profiles',
                                     prog='hmm_overlap')

    # input options
    parser.add_argument('genome', help='Fasta file')
    return parser.parse_args()

def read_fasta(input_file):
    """Reads a FASTA file, gzipped or not, and returns the sequences."""
    if input_file.endswith(".gz"):
        with gzip.open(input_file, "rt") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    else:
        with open(input_file, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    return sequences

def get_basename(file_path):
    """Gets base name of file"""
    base_name = os.path.basename(file_path)
    file_name_without_extension, _ = os.path.splitext(base_name)
    return file_name_without_extension

def iter_aa(filename):
    """Trains and runs Pyrodigal on input sequences."""

    logging.info("Predicting genes")
    sequences = read_fasta(filename)

    # train the orf_finder
    orf_finder = pyrodigal.GeneFinder(meta=True, closed=False)
    #orf_finder = pyrodigal.GeneFinder()
    #orf_finder.train("TTAATTAATTAA".join([str(record.seq) for record in sequences]))

    # iterate over records and generate gene sequences, writing sequences to file
    for record in sequences:
        # find genes in contig
        orfs = orf_finder.find_genes(str(record.seq))
        for gene in orfs:
            yield gene.translate()


if __name__ == "__main__":
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S',
        force=True)

    # Check input ok
    args = get_options()

    logging.info("Reading HMMs")
    with pyhmmer.plan7.HMMFile(HMM_FILE) as hmm_file:
        hmms = list(hmm_file)

    queries = []
    for idx, query_gene in enumerate(iter_aa(args.genome)):
        queries.append(pyhmmer.easel.TextSequence(name = str(idx).encode(), sequence = query_gene).digitize(pyhmmer.easel.Alphabet.amino()))

    logging.info("Running hmmsearch")
    Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])
    results = []
    for hits in pyhmmer.hmmsearch(hmms, queries, bit_cutoffs="trusted"):
        cog = hits.query.name.decode()
        for hit in hits:
            if hit.included:
                results.append(Result(hit.name.decode(), cog, hit.score))

    logging.info("Filtering results")
    best_results = {}
    keep_cog = set()
    for result in results:
        if result.cog in best_results:
            previous_bitscore = best_results[result.cog].bitscore
            if result.bitscore > previous_bitscore:
                logging.debug(f"{result}\n{previous_bitscore}")
                best_results[result.cog] = result
                keep_cog.add(result.cog)
            elif result.bitscore == previous_bitscore:
                if best_results[result.cog].query != hit.query:
                    keep_cog.remove(result.query)
        else:
            best_results[result.cog] = result
            keep_cog.add(result.cog)

    filtered_results = [best_results[k] for k in sorted(best_results) if k in keep_cog]
    for result in filtered_results:
        print(result.query, "{:.1f}".format(result.bitscore), result.cog, sep="\t")




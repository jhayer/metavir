#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    import argparse

    parser = argparse.ArgumentParser(description="To filter sequences by their\
     length")
    parser.add_argument("-i", "--input_file",
                        help="contigs fasta file from assembly", required=True)
    parser.add_argument("-o", "--output_file",
                        help="fasta file with selected contigs of minimal\
                        length", required=True)
    parser.add_argument("-l", "--length",
                        help="Minimal length to filter the contigs",
                        required=True)

    args = parser.parse_args()

    in_file = args.input_file
    out_file = args.output_file
    length = args.length

    filter_seq_length(in_file, out_file, length)


def filter_seq_length(in_file, out_file, length):
    seq_handle = open(in_file, 'r')
    save_handle = open(out_file, 'w')

    for record in SeqIO.parse(seq_handle, "fasta"):
        new_rec = record
        # seq_length = len(new_rec)
        if len(new_rec) >= int(length):
            SeqIO.write(new_rec, save_handle, "fasta")

    seq_handle.close()
    save_handle.close()


if __name__ == '__main__':
    main()

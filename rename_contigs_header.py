#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    import argparse

    parser = argparse.ArgumentParser(description="To rename contigs headers \
        using a prefix provided by the user")

    parser.add_argument("-i", "--input_file",
                        help="contigs fasta file from assembly", required=True)
    parser.add_argument("-o", "--output_file",
                        help="fasta file with renamed headers", required=True)
    parser.add_argument("-a", "--assembler", default="megahit",
                        help="Assembler used. Options: megahit")
    parser.add_argument("-p", "--prefix",
                        help="Prefix to be added to the contigs headers", required=True)

    args = parser.parse_args()

    in_file = args.input_file
    out_file = args.output_file
    assembler = args.assembler
    prefix = args.prefix

    if args.assembler == "megahit" or args.assembler == "spades":
        rename_headers(in_file, out_file, assembler, prefix)
    else:
        raise Exception("Unknown assembler: '%s'" % args.assembler)


def rename_headers(in_file, out_file, ass, prefix):

    seq_handle = open(in_file, 'r')
    save_handle = open(out_file, 'w')

    for record in SeqIO.parse(seq_handle, "fasta"):
        new_rec = record
        if ass == "megahit":
            new_rec = megahit_rename(new_rec, prefix)
        else:
            if ass == "spades":
                new_rec = spades_rename(new_rec, prefix)
            else:
                print("Doing nothing, function not implemented yet for this assembler")

        SeqIO.write(new_rec, save_handle, "fasta")

    seq_handle.close()
    save_handle.close()


def megahit_rename(new_rec, prefix):
    # new_rec = record
    header = new_rec.description
    if " " in header:
        head_list = header.split(" ")
        # header = header.replace(" ", "_")
    # new_header = prefix+'_'+header
    new_header = prefix+'_'+head_list[0]+'_'+head_list[3]
    new_rec.id = new_header
    new_rec.description = ""
    return(new_rec)

def spades_rename(new_rec, prefix):
    header = new_rec.description

    new_header = prefix+'_'+header
    new_rec.id = new_header
    new_rec.description = ""
    return(new_rec)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO, SeqRecord


def parse_kaiju(kaiju_file, tax_dic):

    with open(kaiju_file, 'r') as f:
        for line in f:
            line_tab = line.split('\t')
            if line_tab[2].rstrip() in tax_dic:
                tax_dic[line_tab[2].rstrip()].append(line_tab[1])


def retrieve_seq(fastq_file, tax_dic, out_file):
    # fasta_save = open(out_file, 'w')
    fasta_dic = {taxid: [] for taxid in tax_dic.keys()}

    fastq_handle = open(fastq_file, 'r')
    in_file_format = "fastq"
    # if the user wants to input a fasta file
    if fastq_file.endswith(".fasta"):
        in_file_format = "fasta"

    for record in SeqIO.parse(fastq_handle, in_file_format):

        for tax, reads in tax_dic.items():
            if record.id in reads:
                new_rec = record
                header = record.id+'_'+tax
            #    new_rec = SeqRecord.SeqRecord(record.seq, id=header, description="")
                new_rec.id = header
                fasta_dic[tax].append(new_rec)

    if  out_file.endswith(".fasta"):
        # writing the fasta file (without erasing previous seqs)
        with open(out_file, 'a') as out:
            for taxid in fasta_dic.keys():
                SeqIO.write(fasta_dic[taxid], out, "fasta")
    else:
        # writing the output in fastq
        with open(out_file, 'a') as out:
            for taxid in fasta_dic.keys():
                SeqIO.write(fasta_dic[taxid], out, "fastq")


def taxa_id_listing(taxid_list):
    if taxid_list.endswith(".ids"):
        tax_dic={}
        file_r = open(taxid_list, "r")
        for line in file_r:
            if line != "":
                for taxid in line.split(','):
                    if taxid != "":
                        tax_dic.update( {taxid: []})
    else:tax_dic = {taxid: [] for taxid in taxid_list.split(',')}
    return (tax_dic)

def outfile_test(out_file_name,count):
    if out_file_name.endswith('.fasta'):
        out_file= out_file_name.strip(".fasta")+'_File%d.fasta'%count
        return (out_file)
    elif out_file_name.endswith('.fastq'):
        out_file= out_file_name.strip(".fastq")+'_File%d.fastq'%count
        return (out_file)
    else:
        print ("Error wrong outfile name, Only .fastq and .fasta extension")
        exit ()

def main():

    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-k", "--kaiju_file", help="Kaiju output file (.out)")
    parser.add_argument("-t", "--tax_id_list", help="list of tax_id to retrieve coma-delimited or a .ids file of taxids coma-delimited")
    parser.add_argument("-f", "--fastq_file", help="one or more file to retrieve sequences from")
    parser.add_argument("-o", "--out_file", help="Output fasta file to create")

    args = parser.parse_args()

    kaiju_file = args.kaiju_file
    taxid_list = args.tax_id_list
    fastq = args.fastq_file
    out_file_name = args.out_file
    
    count = 0 
    for fastq_file in fastq.split(','):
        count += 1
        out_file = outfile_test(out_file_name,count)
        tax_dic = taxa_id_listing(taxid_list)
        parse_kaiju(kaiju_file, tax_dic)
        retrieve_seq(fastq_file, tax_dic, out_file)

if __name__ == '__main__':
    main()

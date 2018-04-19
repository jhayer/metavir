#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from Bio import Entrez, SeqIO, SeqRecord


def checker(branch_name, krona_file, fastq, out_file_name, kaiju_file):

        if not branch_name:
                print ("Give the -b argument")
                exit()

        if len(branch_name) >100:
                print ("Error Invalid branch name")
                exit()

        if not krona_file:
                print ("Give the -f argument")
                exit()

        if not krona_file.endswith(".krona"):
                print ("Error No Krona File submit")
                exit ()

        # for each submitted file, test if it's the expected format
        for fastq_file in fastq.split(','):
                if not fastq_file.endswith('.fastq'):
                        print ("error invalid fastq file  name, must end with fastq")
                        exit ()

        # check if the outfile format is the one expected
        if  out_file_name.endswith('.fasta') or out_file_name.endswith('.fastq'):
                useless=42
        else:
                print ("error invalid out_file name, must end with fastq or fasta")
                exit()

        if not kaiju_file:
                print ('Give the -k argument with a kaiju file (.out)')
                exit()

        if not kaiju_file.endswith('.out'):
                print('Error no Kaiju file submit')
                exit()
        return()



def extract_names(branch_name, krona_file, name_file):

        file_r = open (krona_file, "r")
        good_lines = []
        taxa_names = []

        # search if the given name is present then extract everything after it and split it in a list of taxa names
        for line in file_r:
                if branch_name in line:
                        good_part = (line.split(branch_name,1)[-1]).rstrip('\r\n')
                        taxa = re.split("\t",good_part)
                        for tax in taxa:
                                if tax != "":
                                        taxa_names.append(tax)

        taxa_names = sorted(set(taxa_names))
        file_r.close()

        # output the list of names extracted
        if name_file:
                file_w=open(name_file,"a")
                for tax_n in taxa_names:
                        file_w.write("%s\n"%tax_n)
                file_w.close()

        return (taxa_names)



def extract_id (taxa_names, ids_file, taxid_list):

        raw_taxa = []
        taxa = []
        Entrez.email = 'EMAIL_ADRESS'

        # search for the name in the taxonomy bank and extract the hit
        for name in taxa_names:
                id_handle = Entrez.esearch(db="Taxonomy",term=name,idtype="taxid")
                id_rec = Entrez.read(id_handle)
                raw_taxa.append(id_rec['IdList'])

        # clean the extract id
        taxa = [s.strip("[']") for stri in raw_taxa for s in stri]

        # output the list of extracted ids
        if ids_file:
                file_w=open(ids_file,"a")
                for tax_id in taxa:
                        file_w.write("%s,"%tax_id)
                file_w.close()

        taxid_list.extend(taxa)

        return (taxid_list)



def taxa_id_listing(taxid_list):

        tax_dic = {taxid: [] for taxid in taxid_list}
        return (tax_dic)



def outfile_test(out_file_name,count):

        if out_file_name.endswith('.fasta'):
                out_file= out_file_name.strip(".fasta")+'_File%d.fasta'%count
                return (out_file)
        elif out_file_name.endswith('.fastq'):
                out_file= out_file_name.strip(".fastq")+'_File%d.fastq'%count
                return (out_file)



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

        if  out_file.endswith(".fasta") or out_file.endswith(".fa") :
        # writing the fasta file (without erasing previous seqs)
                with open(out_file, 'a') as out:
                        for taxid in fasta_dic.keys():
                                SeqIO.write(fasta_dic[taxid], out, "fasta")
        else:
        # writing the output in fastq
                with open(out_file, 'a') as out:
                        for taxid in fasta_dic.keys():
                                SeqIO.write(fasta_dic[taxid], out, "fastq")



def main():

        parser = argparse.ArgumentParser(description="")

        parser.add_argument("-k", "--kaiju_file", help="Kaiju output file (.out)")
        parser.add_argument("-f", "--fastq_file", help="one or more file to retrieve sequences from")
        parser.add_argument("-o", "--out_file", help="Output fasta file to create")
        parser.add_argument("-b", "--branch_name", help="One or more names (coma-delemited) of Branches in the .krona file ")
        parser.add_argument("-t", "--krona_file", help="Krona file to retrieve Taxid from")
        parser.add_argument("-i", "--ids_file", help="Output name for the list of ids; Default: None")
        parser.add_argument("-n", "--name_file", help="Output name for the list of names; Default: None")

        args = parser.parse_args()

        branches = args.branch_name
        krona_file = args.krona_file
        kaiju_file = args.kaiju_file
        fastq = args.fastq_file
        out_file_name = args.out_file
        ids_file = args.ids_file
        name_file = args.name_file

        taxid_list = []

        for branch_name in branches.split(','):
                checker(branch_name, krona_file, fastq, out_file_name, kaiju_file)
                taxa_names = extract_names(branch_name, krona_file, name_file)
        taxid_list = extract_id(taxa_names, ids_file, taxid_list)

        tax_dic = taxa_id_listing(taxid_list)
        count = 0
        for fastq_file in fastq.split(','):
                count += 1
                out_file = outfile_test(out_file_name,count)
                parse_kaiju(kaiju_file, tax_dic)
                retrieve_seq(fastq_file, tax_dic, out_file)

if __name__ == '__main__':
        main()

        tax_dic = taxa_id_listing(taxid_list)
        count = 0
        for fastq_file in fastq.split(','):
                count += 1
                out_file = outfile_test(out_file_name,count)
                parse_kaiju(kaiju_file, tax_dic)
                retrieve_seq(fastq_file, tax_dic, out_file)

if __name__ == '__main__':
        main()

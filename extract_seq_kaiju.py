#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from Bio import Entrez, SeqIO, SeqRecord


def checker(branches, krona_file, fastq, out_file_name, kaiju_file):
        if not branches:
                print('please give at least one taxaname to start the search')
                exit()
        for branch_name in branches.split(','):
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
        for seq_file in fastq.split(','):
                input_ext=".fastq"
                if not seq_file.endswith('.fastq'):
                    if seq_file.endswith('.fq'):
                        input_ext=".fq"
                    elif seq_file.endswith('.fa'):
                        input_ext=".fa"
                    elif seq_file.endswith('.fasta'):
                        input_ext=".fasta"
                    else:
                        print ("error invalid fastq file  name, must end with fastq, fq, fasta or fa")
                        exit ()

        # check if the outfile format is the one expected
        if out_file_name:
                if  out_file_name.endswith('.fasta') or out_file_name.endswith('.fastq') or out_file_name.endswith('.fa') or out_file_name.endswith('.fq'):
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
        print('checking done')
        return(input_ext)



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

        taxa_names.append(branch_name)
        taxa_names = sorted(set(taxa_names))
        print('extracting names done')
        file_r.close()

        return (taxa_names)



def extract_id (whole_taxa_names, ids_file, taxid_list):
        raw_taxa = []
        taxa = []
        Entrez.email = 'EMAIL_ADRESS'

        # search for the name in the taxonomy bank and extract the hit
        for name in whole_taxa_names:
                id_handle = Entrez.esearch(db="Taxonomy",term=name,idtype="taxid")
                id_rec = Entrez.read(id_handle)
                raw_taxa.append(id_rec['IdList'])

        # clean the extract id
        taxa = [s.strip("[']") for stri in raw_taxa for s in stri]


        taxid_list.extend(taxa)
        print('extracting id done')
        return (taxid_list)



def output_name_ids(whole_taxa_name, taxid_list, name_file, ids_file):
        # output the list of names extracted
        if name_file:
                file_w=open(name_file,"w")
                for tax_n in whole_taxa_name:
                        file_w.write("%s\n"%tax_n)
                file_w.close()

        # output the list of extracted ids
        if ids_file:
                file_w=open(ids_file,"w")
                for tax_id in taxid_list:
                        file_w.write("%s,"%tax_id)
                file_w.close()
        print('outputing id and name file done')


def taxa_id_listing(taxid_list):

        tax_dic = {taxid: [] for taxid in taxid_list}
        return (tax_dic)



def outfile_test(out_file_name,count,seq_file):
        if not out_file_name:
                # use the name of the fastqfile used
                out_file= seq_file.split(input_ext)[0]+'_extracted.fastq'
                return(out_file)

        if out_file_name.endswith('.fasta'):
                #keep the generic name of the fasta file and add a value for each submitted file
                out_file= out_file_name.split(".fasta")[0]+'_File%d.fasta'%count
                return (out_file)

        if out_file_name.endswith('.fa'):
                #keep the generic name of the fa file and add a value for each submitted file
                out_file= out_file_name.split(".fa")[0]+'_File%d.fa'%count
                return (out_file)

        elif out_file_name.endswith('.fastq'):
                out_file= out_file_name.split(".fastq")[0]+'_File%d.fastq'%count
                return (out_file)



def parse_kaiju(kaiju_file, tax_dic):

        with open(kaiju_file, 'r') as f:
                for line in f:
                        line_tab = line.split('\t')
                        if line_tab[2].rstrip() in tax_dic:
                                tax_dic[line_tab[2].rstrip()].append(line_tab[1])
        print('kaiju parsed')

def retrieve_seq(seq_file, tax_dic, out_file, input_ext):
        # fasta_save = open(out_file, 'w')
        fasta_dic = {taxid: [] for taxid in tax_dic.keys()}

        # check if input is fastq or fasta
        fastq_handle = open(seq_file, 'r')
        if input_ext == ".fastq" or input_ext == ".fq":
            in_file_format = "fastq"
        else:
            in_file_format = "fasta"
        print('seq format retrieved')

        # parse the fastqfile and for each tax search if a read is present in the fastq
        # then it add all the corresponding in a dic to be output
        for record in SeqIO.parse(fastq_handle, in_file_format):
                for tax, reads in tax_dic.items():
                        if record.id in reads:
                                new_rec = record
                                header = record.id+'_'+tax
                        #    new_rec = SeqRecord.SeqRecord(record.seq, id=header, description="")
                                new_rec.id = header
                                fasta_dic[tax].append(new_rec)
        fastq_handle.close()

        if  out_file.endswith(".fasta") or out_file.endswith(".fa"):
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

        desc="Extract sequences from fastq or fasta files that match one or more specified taxa names in Kaiju"
        parser = argparse.ArgumentParser(description=desc)

        parser.add_argument("-k", "--kaiju_file", help="Kaiju output file (.out)", required=True)
        parser.add_argument("-f", "--seq_file", help="One or more fastq or fasta file to retrieve sequences from", required=True)
        parser.add_argument("-o", "--out_file", help="Output fasta/fastq file to create, default: name of fastq submit with '_extracted'")
        parser.add_argument("-b", "--branch_name", help="One or more names (coma-delemited) of the highest taxanames to retrieve in the .krona file ", required=True)
        parser.add_argument("-t", "--krona_file", help="Krona file to retrieve Taxid from (.krona)", required=True)
        parser.add_argument("-i", "--ids_file", help="Output name for the list of Taxa_ids; Default: None")
        parser.add_argument("-n", "--name_file", help="Output name for the list of Taxa_names; Default: None")

        args = parser.parse_args()

        branches = args.branch_name
        krona_file = args.krona_file
        kaiju_file = args.kaiju_file
        fastq = args.seq_file
        out_file_name = args.out_file
        ids_file = args.ids_file
        name_file = args.name_file

        taxid_list = []
        whole_taxa_names = []


        input_ext=checker(branches, krona_file, fastq, out_file_name, kaiju_file)
        # search the list of names for each specified branch
        for branch_name in branches.split(','):
                branch_name = branch_name.replace("+"," ")
                taxa_names = extract_names(branch_name, krona_file, name_file)
                whole_taxa_names.extend(taxa_names)

        whole_taxa_names = sorted(set(whole_taxa_names))
        taxid_list = extract_id(whole_taxa_names, ids_file, taxid_list)
        output_name_ids(whole_taxa_names, taxid_list, name_file, ids_file)

        tax_dic = taxa_id_listing(taxid_list)
        count = 0
        # same kaiju file used for all fastq files
        parse_kaiju(kaiju_file, tax_dic)
        # extract the reads in each specified file
        for fastq_file in fastq.split(','):
                count += 1
                out_file = outfile_test(out_file_name,count,fastq_file)
                retrieve_seq(fastq_file, tax_dic, out_file,input_ext)

if __name__ == '__main__':
        main()

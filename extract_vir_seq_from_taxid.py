#!/usr/bin/python2.7

from Bio.Seq import Seq
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = "")

parser.add_argument("-i", "--virus_seq_file", help="virus sequences file retrieved from kraken")
parser.add_argument("-t", "--tax_id_list", help="list of tax_id to retrieve coma-delimited")
parser.add_argument("-p", "--prefix_output", help="prefix to give for the output files")

args = parser.parse_args()

input_seq_file = args.virus_seq_file
taxid_list = args.tax_id_list
out_prefix = args.prefix_output

vir_seq_handle = open(input_seq_file, 'r')

lst_taxid = taxid_list.split(',')

line=vir_seq_handle.readline()
while line:
	if line.startswith("#") and line.split(" -")[0][1:] in lst_taxid:

		new_file_name = out_prefix+"_"+line.split(" -")[0][1:]+".fa"
		save_seq_file = open(new_file_name, 'w')
		line=vir_seq_handle.readline()
		while line.startswith("#")==False:
			save_seq_file.write(line)
			line=vir_seq_handle.readline()
		save_seq_file.close()
	else:
		line=vir_seq_handle.readline()

vir_seq_handle.close()




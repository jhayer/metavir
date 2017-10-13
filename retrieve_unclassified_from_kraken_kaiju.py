#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="")

parser.add_argument("-k", "--result_file",
                    help="kraken_res file or kaiju out file")
parser.add_argument("-s", "--seq_file",
                    help="sequences file. Can be .fastq or .fasta")
parser.add_argument("-u", "--unclassif_seq_file",
                    help="File to store sequences unclassified by kraken. \
                          Can be .fastq or .fasta")

args = parser.parse_args()

res_file = args.result_file
seq_file = args.seq_file
unclass_file = args.unclassif_seq_file

res_k_handle = open(res_file, 'r')
seq_handle = open(seq_file, 'r')
save_unclass_seq = open(unclass_file, 'w')

# kraken res file parsing: to retrieve the seqid for each viral taxid
unclassif_list = []
for res_line in res_k_handle:
    # for Kaiju because the last column is taxid column
    if res_line.endswith('\n'):
        res_line = res_line.replace("\n", "")
    lst_res = res_line.split('\t')
    # if unclassified
    if lst_res[0] == 'U':
        unclassif_list.append(lst_res[1])
# end kraken res file

# retieval of the sequences
if seq_file.endswith(".fasta"):
    # put the sequences of the fasta input file in a dic
    seq_dic = SeqIO.index(seq_file, "fasta")
else:
    seq_dic = SeqIO.index(seq_file, "fastq")

# write the unclassified sequences fasta/fastq file
for unc_seq in unclassif_list:
    seq_record = seq_dic[unc_seq]
    if unclass_file.endswith(".fasta"):
        formatted_seq = seq_record.format('fasta')
    else:
        formatted_seq = seq_record.format('fastq')
    save_unclass_seq.write(formatted_seq)

res_k_handle.close()
seq_handle.close()
save_unclass_seq.close()

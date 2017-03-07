#!/usr/bin/python2.7

from Bio.Seq import Seq
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = "")

parser.add_argument("-k", "--kraken_res", help="kraken_res file")
parser.add_argument("-r","--kraken_report", help="kraken_report file")
parser.add_argument("-s","--seq_file", help="sequences file")
parser.add_argument("-o","--output_seq_file", help="Output file")
parser.add_argument("-u","--unclassif_seq_file", help="Fasta file to store sequences unclassified by kraken")

args = parser.parse_args()

res_file = args.kraken_res
report = args.kraken_report
seq_file = args.seq_file
out_file = args.output_seq_file
unclass_file = args.unclassif_seq_file

res_k_handle = open(res_file, 'r')
report_handle = open(report, 'r')
seq_handle = open(seq_file, 'r')
save_seq_file = open(out_file, 'w')
save_unclass_seq = open(unclass_file, 'w')

vir_dic = {}
b_vir_section = False

##- kraken report parsing: to retrieve all the viral taxid
for line in report_handle:
	lst_line = line.split('\t')
	tax_name = lst_line[5].replace(" ","")
	tax_name = tax_name.replace("\n","")

	tax_id = lst_line[4]
	if tax_id=='10239':
		b_vir_section = True

	if b_vir_section==True and tax_id=='131567':
		break

	if b_vir_section==True:
		if int(lst_line[2])>0:
			vir_dic[lst_line[4]] = [lst_line[2],tax_name]
##- end kraken report parsing

##-- kraken res file parsing: to retrieve the seqid for each viral taxid
unclassif_list = []
for res_line in res_k_handle:
	lst_res = res_line.split('\t')
	#if classified
	if lst_res[0]=='C':
		if lst_res[2] in vir_dic:
			vir_dic[lst_res[2]].append(lst_res[1])
	else:
		## treatment for unclassified sequences 'U'
		unclassif_list.append(lst_res[1])
		pass
##-- end kraken res file

#####retieval of the sequences
#put the sequences of the fasta input file in a dic
seq_dic = SeqIO.index(seq_file, "fasta")

## write the unclassified sequences fasta file
for unc_seq in unclassif_list:
	seq_record = seq_dic[unc_seq]
	seq_desc = seq_record.description
	sequence = seq_record.seq
	#write the retrieved sequence with description in fasta output
	fasta_seq = str('>'+seq_desc+'\n'+sequence+'\n')
	save_unclass_seq.write(fasta_seq)
##

#### write the fasta file with sequences classified as viruses by kraken
for k, v in vir_dic.iteritems():
	#write a header for the virus: #taxid - virus name - num_seq:X
	vir_header="#"+k+" - "+v[1]+" - num_seq:"+v[0]+"\n"
	save_seq_file.write(vir_header)

	for i, val in enumerate(v): 
		## retrieve all the seqs from the list for each key of the dic
		if i>1:
	
			seq_record = seq_dic[val]
			seq_desc = seq_record.description
			sequence = seq_record.seq
			#write the retrieved sequence with description in fasta output
			fasta_seq = str('>'+seq_desc+'\n'+sequence+'\n')
			save_seq_file.write(fasta_seq)
####

res_k_handle.close()
report_handle.close()
seq_handle.close()
save_seq_file.close()
save_unclass_seq.close()


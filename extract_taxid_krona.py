#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import re
from Bio import Entrez


def checker(branch_name, krona_file,out_file):
	if not branch_name and not krona_file:
		print ("Give at least the -b and -f argument (default output will be ids_extracted.out)")
		exit()
	if not branch_name:
		print ("Give the -b argument")
		exit()
	if len(branch_name)>100:
		print ("Error Invalid branch name")
		exit()
	if not krona_file:
		print ("Give the -f argument")
		exit()	
	if not krona_file.endswith(".krona"):
		print ("Error No Krona File submit")
		exit ()
	if not out_file:
		out_file = "ids_extracted.out"
	return (out_file)

def extract_names(branch_name, krona_file, name_file):
	file_r = open (krona_file, "r")
	good_lines = []
	taxa_names = []
	for line in file_r:
		if branch_name in line:
			good_part = (line.split(branch_name,1)[-1]).rstrip('\r\n')
			taxa = re.split("\t",good_part)
			for tax in taxa:
				if tax != "":
					taxa_names.append(tax)
	taxa_names = sorted(set(taxa_names))
	file_r.close()	
	if name_file:
		file_w=open(name_file,"a")
		for tax_n in taxa_names:			
			file_w.write("%s\n"%tax_n)
		file_w.close()
	return (taxa_names)


def extract_id (taxa_names, out_file):
	raw_taxa = []
	taxa = []
	Entrez.email = 'EMAIL_ADRESS'
	for name in taxa_names:
		id_handle = Entrez.esearch(db="Taxonomy",term=name,idtype="taxid")
		id_rec = Entrez.read(id_handle)
		raw_taxa.append(id_rec['IdList'])
	taxa = [s.strip("[']") for stri in raw_taxa for s in stri]
	file_w=open(out_file,"a")
	for tax_id in taxa:
		file_w.write("%s,"%tax_id)
	file_w.close()
	
	
def main():

    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-b", "--branch_name", help="One or more names (coma-delemited) of Branches in the .krona file ")
    parser.add_argument("-f", "--krona_file", help="Krona file to retrieve Taxid from")
    parser.add_argument("-o", "--out_file", help="Output name for the list of ids; Default: ids_extracted.out")
    parser.add_argument("-n", "--name_file", help="Output name for the list of names; Default: None")

    args = parser.parse_args()

    branches = args.branch_name
    krona_file = args.krona_file
    out_file = args.out_file
    name_file = args.name_file
    
    for branch_name in branches.split(','):
	    out_file = checker(branch_name, krona_file, out_file)
	    taxa_names = extract_names(branch_name, krona_file, name_file)
	    extract_id(taxa_names, out_file)

if __name__ == '__main__':
    main()

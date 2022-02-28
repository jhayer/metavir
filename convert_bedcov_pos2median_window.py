#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import statistics
import sys

def main():
    import argparse
    parser = argparse.ArgumentParser(description="bla bla")
    parser.add_argument("-i", "--input_file", help="bed file with 3 columns: seq_id, position, coverage", required=True)
    parser.add_argument("-o", "--output_file", help="tab bed file with median coverage per window of X positions", required=True)
    parser.add_argument("-w", "--window_size", help="Window size for computing median coverage",required=True)

    args = parser.parse_args()

    in_file = args.input_file
    out_file = args.output_file
    window = int(args.window_size)

    compute_median_per_window(in_file, out_file, window)

def compute_median_per_window(in_file, out_file, window):
    bed_handle = open(in_file, 'r')
    save_handle = open(out_file, 'w')
# NODE_9_length_68152_cov_876.648373	1	96
# NODE_9_length_68152_cov_876.648373	2	122
# NODE_9_length_68152_cov_876.648373	3	124
# NODE_9_length_68152_cov_876.648373	4	150
# NODE_9_length_68152_cov_876.648373	5	151

    per_contig_dic={}
    i=0
    st_pos=1
    contig_name=""
    
    for line in bed_handle:
        l_line=line.split('\t')
        if l_line[0] == contig_name or contig_name=="":
            contig_name=l_line[0]
            #add coverage to list fot start position
            if i<window:
                if st_pos not in per_contig_dic:
                    per_contig_dic[st_pos]=[]
                per_contig_dic[st_pos].append(int(l_line[2].rstrip()))
                i=i+1
            else:
                median=statistics.median(per_contig_dic[st_pos])
                #save the median cov for current contig, next x (window size) position
                save_handle.write(contig_name+"\t"+str(st_pos)+"\t"+str(st_pos+i)+"\t"+str(median)+"\n")
                i=0
                st_pos=int(l_line[1])
                per_contig_dic[st_pos]=[]


        else:
            median=statistics.median(per_contig_dic[st_pos])
            # before changing contig, we save the last record of median cov for the last positions of the previous contig
            save_handle.write(contig_name+"\t"+str(st_pos)+"\t"+str(st_pos+i)+"\t"+str(median)+"\n")
            i=0
            st_pos=int(l_line[1])
            per_contig_dic={}
            per_contig_dic[st_pos]=[]
            # dealing now with new contig
            contig_name=l_line[0]

# Last case to deal with
    median=statistics.median(per_contig_dic[st_pos])
    # save last median cov for last contig last segment
    save_handle.write(contig_name+"\t"+str(st_pos)+"\t"+str(st_pos+i)+"\t"+str(median)+"\n")
    i=0
    st_pos=int(l_line[1])
    per_contig_dic[st_pos]=[]
    # dealing now with new contig
    contig_name=l_line[0]

    bed_handle.close()
    save_handle.close()

if __name__ == '__main__':
    main()

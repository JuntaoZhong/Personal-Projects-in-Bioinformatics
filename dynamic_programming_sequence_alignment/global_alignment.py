#!/usr/bin/env python3
import local_alignment as la
import math
import textwrap
import sys

'''Written by Jimmy Zhong (zhongj2@carleton.edu) at Carleton College 
perform local aligment for 2 sequences in python3, last run Nov 23rd, 2021
Relies on local_alignment.py
'''

def global_backtrack(result, len_seq1):
    path = []
    cur = result[-1][-1]
    print(f"Alignment Score: {cur.score}\n")
    while cur.source != -1:
        path.append(cur.ID)
        cur = result[la.get_row(cur.source, len_seq1)][la.get_col(cur.source, len_seq1)]
    path.reverse()
    return path

def global_algin_sequence(scorer, seq1, seq2):
    '''initialized the first row of the dp array, pass it la.dp_2D_transver() to fill it out,
    then get the back track path using backtrack()'''
    first_row = [la.cell(0, 0)]
    for each in range(len(seq1)): first_row.append(la.cell(each + 1, each + 1))
    dp_array = la.dp_2D_transverse(scorer, [first_row], seq1, seq2, "global")
    backpath = global_backtrack(dp_array, len(seq1))
    alignment = la.print_out_alignment(backpath, seq1, seq2, 60)
    return alignment

if __name__ == "__main__":
    scorer, seq1, seq2 = la.input_check_and_init()
    alignment_result = global_algin_sequence(scorer, seq1.seq.upper(), seq2.seq.upper())
    print('\n'.join(textwrap.wrap(f"top sequence: {seq1.description}", 60)))
    print('\n'.join(textwrap.wrap(f"bottom sequence: {seq2.description}", 60)))
    print("\nOptimal glocal Alignment:")
    for line in alignment_result: print(line)

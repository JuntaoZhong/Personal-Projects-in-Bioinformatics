#!/usr/bin/env python3
# ./local_alignment.py TP53-Cat.fasta TP53-Rat.fasta scoring1.txt 
from Bio import SeqIO
import os
import textwrap
import sys
'''Written by Jimmy Zhong (zhongj2@carleton.edu) at Carleton College 
perform local aligment for 2 sequences in python3, last run Oct 4th, 2021
Homework for Professor Oesper's CS 362 Computational Biology
'''
class cell:
    '''store path and score for each grid of the 2d dp array'''
    def __init__(self, score, ID):
        self.score, self.ID, self.source = score, ID, -1

class scorer:
    '''return the score for match vs mismatch'''
    def __init__(self, match, mismatch, gap):
        #gap panelty and mismatch should be non-positive
        if min(int(mismatch), int(gap))>0:
            sys.exit("mismatch and gap must have a non-positive value (panelty)")
        self.match, self.mismatch, self.gap = int(match), int(mismatch), int(gap)
    def score(self, char1, char2):
        return self.match if char1.upper() == char2.upper() else self.mismatch

def fasta_handler(fasta_file):
    '''return a SeqIO object of a single sequence. Only allows 1 sequence in a fasta file'''
    seq_list = [seq for seq in SeqIO.parse(fasta_file, "fasta")]
    if len(seq_list) != 1: 
        sys.exit(f"We allow exactly 1 sequence in: {fasta_file}, but you have: {len(seq_list)}")
    return seq_list[0]

def init(seq1, seq2, scoring_matrix_path):
    '''read in both sequences files, and the scroring matrix to initialize'''
    with open(scoring_matrix_path, "r") as f: lines = f.readlines()
    for line in lines:
        if not line.startswith("#"): 
            s = line.split(",")
            if len(s) != 3: sys.exit(f"error in the scoring matrix file: {scoring_matrix_path}")
            return scorer(s[0], s[1], s[2]), fasta_handler(seq1), fasta_handler(seq2)

def get_row(ID, seq1_len): return ID//(seq1_len+1) # get row number using cell ID
def get_col(ID, seq1_len): return ID%(seq1_len+1) # get col number using cell ID

def print_out_alignment(path_list, seq1, seq2, per_line):
    ''' given the backtrack path, print out a nicely formated alignment'''
    prev = path_list[0]
    seq1_start = get_col(prev, len(seq1))-1
    seq2_start = get_row(prev, len(seq1))-1
    seq1_char, seq2_char = seq1[seq1_start], seq2[seq2_start]
    print(f"top alignment start pos: {seq1_start+1} \nbottom alignment start pos: {seq2_start+1}\n")
    top, middle, bottom = seq1_char, "|" if seq1_char == seq2_char else " ", seq2_char
    for step in path_list[1:]:
        seq1_char = seq1[get_col(step, len(seq1))-1]
        seq2_char = seq2[get_row(step, len(seq1))-1]
        difference = step - prev
        if difference == 1:
            bottom += "-"
            middle += " "
            top += seq1_char
        elif difference%(len(seq1)+1) == 0:
            bottom += seq2_char
            middle += " "
            top += "-"
        elif difference%(len(seq1)+1) == 1:
            bottom += seq2_char
            middle += "|" if seq1_char == seq2_char else " "
            top += seq1_char
        else: sys.exit("unknown error occurs at the traceback step")
        prev = step
    pretty_align = [f'''{"".join(top[i:i+per_line])}\n{"".join(middle[i:i+per_line])}\
    \n{"".join(bottom[i:i+per_line])}\n''' for i in range(0, len(top), per_line)]
    return pretty_align

def backtrack(result, len_seq1, max_col, max_row):
    '''find the path to the start of the match given the cell with maximum score'''
    path = []
    cur = result[max_row][max_col]
    while cur.source != -1:
        path.append(cur.ID)
        cur = result[get_row(cur.source, len_seq1)][get_col(cur.source, len_seq1)]
    path.reverse()
    return path

def dp_2D_transverse(scorer, dp_array, seq1, seq2, global_alignment = False):
    '''fill out the dp_array using Smith–Waterman local sequence alignment algorithm'''
    if min(len(seq1), len(seq2)) < 3: 
        sys.exit("sequence must be longer than 3bp (1 amino acid codon)!")
    for row in range(len(seq2)):
        last_row = dp_array[-1]
        new_row = [cell(0, last_row[-1].ID + 1)] # first col of the row
        for col in range(len(seq1)):
            choices = [last_row[col].score + scorer.score(seq1[col], seq2[row]), 
                last_row[col+1].score + scorer.gap, new_row[-1].score + scorer.gap, 0]
            if global_alignment: 
                choices = choices[:-1]
            max_score = max(choices)
            max_index = choices.index(max_score) # return the first index of max_score
            cur_cell = cell(max_score, new_row[-1].ID + 1) # take previous ID + 1 as new ID
            if max_index == 0: cur_cell.source = last_row[col].ID 
            elif max_index == 1: cur_cell.source = last_row[col+1].ID
            elif max_index == 2: cur_cell.source = new_row[-1].ID
            elif max_index == 3: pass #local alignment only, source is -1 by default
            else: sys.exit("error during the filling of dp_array")
            new_row.append(cur_cell)
        dp_array.append(new_row)
    return dp_array

def get_max_index(dp_array):
    '''find the max element in the filled out dp_array'''
    max_row, max_col, cur_max = -1, -1, 0
    for r in range(len(dp_array)):
        for c in range(len(dp_array[0])):
            if cur_max < dp_array[r][c].score: 
                max_row, max_col, cur_max = r, c, dp_array[r][c].score
    if cur_max == 0: sys.exit("no alignment: alignment scores are > 0")
    return max_row, max_col, cur_max

def algin_sequence(scorer, seq1, seq2):
    '''initialized the first row of the dp array, pass it dp_2D_transver() to fill it out,
    then get the back track path using backtrack()'''
    first_row = [cell(0, 0)]
    for each in range(len(seq1)): first_row.append(cell(0, each + 1))
    dp_array = dp_2D_transverse(scorer, [first_row], seq1, seq2)
    max_row, max_col, max_score = get_max_index(dp_array)
    print(f"Alignment Score: {max_score}\n")
    backpath = backtrack(dp_array, len(seq1), max_col, max_row)
    alignment = print_out_alignment(backpath, seq1, seq2, 60)
    return alignment

def input_check_and_init():
    if (len(sys.argv) != 4): 
        sys.exit("I need 3 command line inputs: sequence1.fasta, sequence2.fasta, scoring.csv")
    scorer, seq1, seq2 = init(sys.argv[1], sys.argv[2], sys.argv[3])
    # print(seq1, seq2)
    print("-"*20, "local alignment", "-"*20)
    # seq1 and seq2 are SeqIO object from biopython
    if seq1.seq == seq2.seq:
        sys.exit("Two input sequences are identical. Will not perform alignment.")
    return scorer, seq1, seq2

if __name__ == "__main__":
    scorer, seq1, seq2 = input_check_and_init()
    alignment_result = algin_sequence(scorer, seq1.seq.upper(), seq2.seq.upper())
    print('\n'.join(textwrap.wrap(f"top sequence: {seq1.description}", 60)))
    print('\n'.join(textwrap.wrap(f"bottom sequence: {seq2.description}", 60)))
    print("\nOptimal Local Alignment:")
    for line in alignment_result: print(line)

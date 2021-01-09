# Written by Jimmy Zhong (zhong.juntao.2019@gmail.com) at Carleton College 
# for python3, last run Jan 9th, 2020
# Does a simple sequence alignment with custom alignment scores. 
# Return all possible alignment(s) with the highest alignment score   
# Wrote it to double check the correctness of my Biology hw

import os
import sys
class a_block:
    def __init__(self, score, ID):
        self.score = score
        self.ID = ID
        self.source = []

def print_score_array(array):
    print("\nResult:")
    print("the alignment scores is: " + str(array[-1][-1].score))
    for row in array:
        print_this_row = []
        for each in row:
            if type(each) is str:
                while len(each) < 4:
                    each = each + " "
                print_this_row.append(each)
            elif type(each) is a_block:
                string_score = str(each.score)
                while len(string_score) < 4:
                    string_score = string_score + " "
                print_this_row.append(string_score)
        print(*print_this_row, sep =', ')

def print_path_array(array):
    print("\nHere is the track back path. The first number is the unique ID of this block;")
    print("number(s) after colon shows where this block got its score from (can have multiple source)")
    for row in array:
        print_this_row = []
        for each_block in row:
            if type(each_block) is str:
                while len(each_block) < 10:
                    each_block = each_block + " "
                print_this_row.append(each_block)
            elif type(each_block) is a_block:
                this_cell = str(each_block.ID) + ":"
                source_list = each_block.source
                for each_source in source_list:
                    this_cell = this_cell + str(each_source.ID) + " "
                while len(this_cell) < 10:
                    this_cell = this_cell + " "
                print_this_row.append(this_cell)
                #print(len(source_list))
        print(*print_this_row, sep=", ")

def print_out_alignment (path_list, seq_1, seq_2):
    line_above = ""
    line_below = ""
    seq_1_index = len(seq_1)-1
    seq_2_index = len(seq_2)-1
    for i in range(1, len(path_list)):
        if (path_list[i-1] - path_list[i]) == 1:
            line_below = line_below + "-"
            line_above = line_above + seq_1[seq_1_index]
            seq_1_index = seq_1_index - 1
        elif (path_list[i-1] - path_list[i]) == (len(seq_1) + 1):
            line_above = line_above + "-"
            line_below = line_below + seq_2[seq_2_index]
            seq_2_index = seq_2_index - 1
        elif (path_list[i-1] - path_list[i]) == (len(seq_1) + 2):
            line_above = line_above + seq_1[seq_1_index]
            line_below = line_below + seq_2[seq_2_index]
            seq_1_index, seq_2_index = seq_1_index - 1, seq_2_index - 1
        else:
            print("error!!!")
    print(line_above[::-1] + "\n"+ line_below[::-1] + "\n")

def recursive_path_helper(path, cur, des, seq_1, seq_2):
    this_path = path.copy()
    this_path.append(cur.ID)
    if cur == des:
        print("\nthe path on the matrix is:" + str(this_path[::-1]) + ", and the alignment is: ")
        print_out_alignment (this_path, seq_1, seq_2)
    elif len(cur.source) == 0:
        print("I should never come here, error!!!!!!")
        return
    else:
        for source in cur.source:
            recursive_path_helper(this_path, source, des, seq_1, seq_2)

def print_backtrack(array, seq_1, seq_2):
    start = array[-1][-1]
    destination = array[1][1]
    recursive_path_helper([], start, destination, seq_1, seq_2)

def dp_2D_transverse(array, match, mismatch, gap):
    if len(array) < 3:
        print("yooo, check your sequence length, they are too short!")
        return
    for row in range(2, len(array)):
        for col in range(2, len(array[0])):
            #case 2, mismatch or gap
            src_score = [array[row-1][col-1].score + mismatch, array[row-1][col].score + gap, array[row][col-1].score + gap]
            #case 1, match! 
            if array[0][col] == array[row][0]:
                src_score.append(array[row-1][col-1].score + match)
            max_score = max(src_score)

            array[row][col] = a_block(max_score, array[row][col-1].ID + 1)
            #find its source
            src_list =[]
            for i in range(len(src_score)):
                if max_score == src_score[i]:
                    src_list.append(i)
            for a_source in src_list:
                if a_source == 0 or a_source == 3:
                    array[row][col].source.append(array[row-1][col-1])
                elif a_source == 1:
                    array[row][col].source.append(array[row-1][col])
                elif a_source == 2:
                    array[row][col].source.append(array[row][col-1])
                

#gap panelty and mismatch should be non-positive
def algin_sequence(match, mismatch, gap, seq_1, seq_2):
    #initialize the 2-D Array
    array_2D = []
    first_row = [letter for letter in seq_1]
    first_row = [' ', ' '] + first_row
    array_2D.append(first_row)
    second_row = [' ']
    for i in range(len(seq_1) + 1):
        second_row.append(a_block(0+i*gap, i))
        if i > 0:
            second_row[i+1].source.append(second_row[i])
    array_2D.append(second_row)

    block_above = array_2D[1][1]
    for i in range(len(seq_2)):
        this_row = ['u_d']*(len(seq_1)+2)
        this_row[0], this_row[1] = seq_2[i], a_block(0+(i+1)*gap, (i+1)*(len(seq_1) + 1))
        this_row[1].source.append(block_above)
        block_above = this_row[1]
        array_2D.append(this_row)
    dp_2D_transverse(array_2D, match, mismatch, gap)
    print_score_array(array_2D)
    print_path_array(array_2D)
    print_backtrack(array_2D, seq_1, seq_2)

def main():
    score = [2, -1, -2]
    print("Parameter set up: Align_score, Misalign_punishment, Gap_punishment (default: 2, -1, -2)")
    input_score = input("Press Enter to use default, else, give 3 integer score separate by comma:")
    if len(input_score) != 0:
        input_score = input_score.replace(" ", "")
        score = input_score.split(",")
    print("alignment based on Align: " + str(score[0]) + ", Misalign: " + str(score[1]) + ", Gap: " + str(score[2]))

    input_seq1 = input("what's your first sequence? DNA sequences only have ATCG or AUCG: ")
    input_seq2 = input("Second sequence:")

    print("aligning: " + input_seq1.upper() + " and " + input_seq2.upper())
    algin_sequence(int(score[0]), int(score[1]), int(score[2]), input_seq1.upper(), input_seq2.upper())


if __name__ == "__main__":
    main()
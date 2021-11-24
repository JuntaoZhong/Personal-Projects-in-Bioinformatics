import sequence as s
import math
import BuildPairList as build
from Bio.SubsMat import MatrixInfo as matlist

aminos = list("ARNDBCEQZGHILKMFPSTWYV")
blosum62 = matlist.blosum62
seqs = [s.chicken(), s.human(), s.macaque(), s.cat(), s.crab_eating_macque(), s.macaw()]

seq_len = len(seqs[0]) # sequences have the same lengths, aligned.

def get_score(char1, char2):
    return -blosum62[(char1, char2)] if (char1, char2) in blosum62 else -blosum62[(char2, char1)]

def permute(seed):
    '''Given a list of characters, return all of its permutation
    for a list of length 6, there will be 6! or 720 permutaions
    '''
    if len(seed) < 2: 
        return seed
    result = [] 
    for i in range(len(seed)):
       this = seed[i]
       rest = seed[:i] + seed[i+1:]
       for p in permute(rest):
           result.append(this+p)
    return result

def _amino_sankoff(btree, root, traceback = False):
    '''Running Sankoff Algorithm (Dynamic Programming tree) 
    on distance matrix of amino acids
    '''
    res = []
    while btree:
        l, r = btree.pop(0)
        res.append([l,r])
        for amino in aminos:
            amino_l = amino_r = math.inf
            for a in aminos:
                score = get_score(a, amino)
                amino_l = min(l.distance[a]+score, amino_l)
                amino_r = min(r.distance[a]+score, amino_r)
            l.parent.distance[amino] = int(amino_l + amino_r)
        # if traceback_mode:
    if traceback:
        return res,root
    else:
        return root

def get_single_char_diff(root, node_list):
    ''' tree is build based on a specific 
        permutation and a specific tree shape
    '''
    done_root = _amino_sankoff(node_list, root)
    min_dist = math.inf
    for key in root.distance.keys():
        min_dist = min(min_dist, done_root.distance[key])
    return min_dist

def build_single_char_tree(permutation, tree_shape, cur_index):
    ''' Give a specific permutation (list), a tree shape, 
        and the current index on the sequence (from 0 to end),
        build a tree to pass into get_single_char_diff()
    '''
    hextuple = [seqs[int(order)][cur_index] for order in permutation] # a 6-character list, from the i th possition of taxa 1-6, order by the permuation input
    parsimony_score = build.tree(hextuple, tree_shape)
    return parsimony_score

def sankoff_analysis():
    ''' ASSUME THERE ARE 6 TAXA!
    Loop from all 6 tree shapes possible for 6 taxa, 
    each tree shape has be tested with 720 permutations.
    return the most parsimony tree + permute combo.
    The record of all parimony trees is in result.txt
    '''
    permutations = permute(list("012345"))

    with open("result.txt", "w") as f: 
        f.write("parsimony score, tree type, permutation\n")
    
    min_parsimony_score = min_tree_type = math.inf
    min_permutation = permutations[0]
    counter = 0
    for tree_type in range(1, 7):
        for per in permutations:
            counter += 1
            sum_score = 0
            for i in range(seq_len): # sc = single character
                sc_tree_root, sc_tree_list = build_single_char_tree(per, tree_type, cur_index=i)
                sum_score += get_single_char_diff(sc_tree_root, sc_tree_list)
            
            if counter % 10 == 1:
                print("trail #: ", counter, " current best: ", min_parsimony_score, "this: ", sum_score)
            
            if sum_score <= min_parsimony_score:
                min_parsimony_score = sum_score
                min_tree_type = tree_type
                min_permutation = per
                
                with open("result.txt", "a") as f:
                    f.write(f"{min_parsimony_score}, {min_tree_type}, {per}\n")
        
        print("done tree type: ", tree_type)
    return min_parsimony_score, min_tree_type, min_permutation

if __name__ == "__main__":
    min_parsimony_score, min_tree_type, min_permutation = sankoff_analysis()

    print(min_parsimony_score, min_tree_type, min_permutation)

    # result： 78 1 053124
    # or: 78 4 053124
    name = ["chiken", "human", "monkey", "cat", "crab-eating monkey", "macaw"]
    for index in "053124":
        print(name[int(index)])

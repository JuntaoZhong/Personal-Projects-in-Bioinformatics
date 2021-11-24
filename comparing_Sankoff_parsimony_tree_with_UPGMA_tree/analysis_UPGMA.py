import sequence as s
import UPGMA

chicken, human, monkey, cat, monkey2, macaw = s.chicken(), s.human(), s.macaque(), s.cat(), s.crab_eating_macque(), s.macaw()

seq = [chicken, human, monkey, cat, monkey2, macaw]
dist_matrix = []
for i in range(len(seq)):
    row = []
    for j in range(len(seq)):
        hamming_dist = s.compare_seq(seq[i], seq[j])
        row.append(hamming_dist)
    dist_matrix.append(row)

print(dist_matrix)
mat_obj = UPGMA.UPGMA(dist_matrix, names = ["chiken", "human", "monkey", "cat", "crab-eating monkey", "macaw"])
print('*' * 60)
print('The Phylogenetic tree is:\n')
mat_obj.cluster().print_tree()


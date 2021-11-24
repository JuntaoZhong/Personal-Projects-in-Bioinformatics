import sequence as s
import BuildPairList as build
import analysis_Sankoff as san

node_seq_file = "sankofftree.txt"
seqs = [s.chicken(), s.human(), s.macaque(), s.cat(), s.crab_eating_macque(), s.macaw()]

print(seqs)
def _Traversal(root):
    '''
    Traverse back the given tree and fill each nodes with the assinged character from root to leaves by recursive methods. 
    '''  
    if root.char == "root":
        root.char = min(root.distance,key=root.distance.get)
    if root.left and root.right:
        char = root.char
        for i in root.left.distance.keys():
            for j in root.right.distance.keys():
                score = root.left.distance[i] + root.right.distance[j]+san.get_score(i, char)+san.get_score(j, char)
                if score == root.distance[char]:
                  root.left.char = i
                  root.right.char = j       
        _Traversal(root.left) 
        _Traversal(root.right)

def sankoff_analysis(leaves, treetype):
    '''
    After known which tree type is the optimal phlogenetic tree type and optimal sequence of leaves notation, traverse the anaylsis the whole tree and create the output for sankoff's algorithm in sankofftree.txt.
    '''  
    res_string=""
    root_list = []
    tree_list = []

    for i in range(len(seqs[0])): 
      re_tree_root, re_tree_list = san.build_single_char_tree(leaves, treetype, cur_index=i)
      tree,root = san._amino_sankoff(re_tree_list, re_tree_root, True)
      root_list.append(root)
      tree_list.append(tree) 

    for root in root_list:
      _Traversal(root)
    
    for i in range(len(tree_list[0])):
      for j in range(0,2):
        res_string+="node"+str(i)+"."+str(j)+":\n"  
        for k in range(len(tree_list)):
          res_string += tree_list[k][i][j].char
        res_string += "\n"
    res_string+="root:\n"    
    for i in range(len(root_list)):
      res_string+= root_list[i].char   
    with open(node_seq_file, "w") as f: 
        f.write(res_string)

def SankoffDistCal(cur, seq_dict):
    '''
    Calculate the edge lengths of the phylogenetic tree Sankoff's algorithm output.
    '''
    if cur.left and cur.right:
        cur_seq = seq_dict[cur.char]
        l_dist = s.compare_seq(cur_seq, seq_dict[cur.left.char])
        r_dist = s.compare_seq(cur_seq, seq_dict[cur.right.char])
        cur.distance = (l_dist, r_dist)
        if cur.left.left:
            SankoffDistCal(cur.left, seq_dict)
        if cur.right.left:
            SankoffDistCal(cur.right, seq_dict)

if __name__ == "__main__":
    sankoff_analysis("421350",4)

    with open(node_seq_file, "r") as f:
        lines = f.readlines()
    seq_dict = {}
    for i in range(len(lines)//2):
        node_name = lines[2*i].strip().strip(":").strip("node")
        seq = lines[2*i + 1].strip()
        seq_dict[node_name] = seq

    node_list = ["0.0", "0.1", "2.1", "3.1", "1.0", "1.1", "2.0", "3.0", "4.0", "4.1", "root"] # abcdefghij
    root = build.tree(node_list, 4)[0] # don't need tree list
    SankoffDistCal(root, seq_dict)
    root.print_tree()

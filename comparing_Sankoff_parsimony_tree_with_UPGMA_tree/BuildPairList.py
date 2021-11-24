import BinaryTree as BT
import math
import sys

def _build_amino_dist(leaf):
    aminos = ["A","R","N","D","B","C","E","Q","Z","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    leaf.distance = {a:math.inf for a in aminos}
    if leaf.char in aminos:
        leaf.distance[leaf.char] = 0

def list_add_amino_dist(leaf_list):
    [_build_amino_dist(leaf) for leaf in leaf_list]

def create_all_nodes(hextuple, only_leaf = False):
    '''
    Helper function that initialize the leaves by assinging characters to them according to the given sequence. 
    '''
    if len(hextuple) == 11: # already have all the nodes
        a, b, c, d, e, f = [BT.BinaryTree(i+1, char=hextuple[i]) for i in range(6)]
        g,h,i,j = [BT.BinaryTree(-1, char=hextuple[i]) for i in range(6,10)]
        k = BT.BinaryTree(-1, char=hextuple[10])
        return a,b,c,d,e,f,g,h,i,j,k

    a, b, c, d, e, f = [BT.BinaryTree(i+1, char=hextuple[i]) for i in range(6)]
    g,h,i,j = [BT.BinaryTree(-1, char="inner") for i in range(4)]
    k = BT.BinaryTree(-1, char="root")
    list_add_amino_dist([a,b,c,d,e,f,g,h,i,j,k])
    if only_leaf: 
        return [a,b,c,d,e,f]
    else:
        return a,b,c,d,e,f,g,h,i,j,k

def tree(hextuple, tree_type):
    '''
    Create all possible trees types that have 6 leaves and assign their leaves by the given hextuple.
    '''
    if not (len(hextuple) == 6 or len(hextuple) == 11): 
        sys.exit("must have 6 taxa, or all 11 nodes")
    if tree_type == 1: return _tree1(hextuple)
    elif tree_type == 2: return _tree2(hextuple)
    elif tree_type == 3: return _tree3(hextuple)
    elif tree_type == 4: return _tree4(hextuple)
    elif tree_type == 5: return _tree5(hextuple)
    elif tree_type == 6: return _tree6(hextuple)
    else: sys.exit("there are only 6 shapes for trees with 6 taxa")

def _tree1(hextuple):
    leaf_list = create_all_nodes(hextuple, "only_leaf")
    left = leaf_list[0]
    tuple_list = []
    for i in range(1,6):
        right = leaf_list[i]
        tuple_list.append((left, right))
        p = BT.BinaryTree(-1, left = left, right = right, char = "inner")
        _build_amino_dist(p)
        left.parent = p
        right.parent = p
        left = p
    left.char == "root"
    return left, tuple_list
    
def _tree5(hextuple):
    leaf_list = create_all_nodes(hextuple, "only_leaf")
    tuple_list = []
    while len(leaf_list) > 1:
        left = leaf_list.pop(0)
        right = leaf_list.pop(0)
        p = BT.BinaryTree(-1, left = left, right = right, char = "inner")
        _build_amino_dist(p)
        leaf_list.append(p)
        left.parent = p
        right.parent = p
        tuple_list.append((left,right))
    return leaf_list[0], tuple_list

def _tree2(hextuple):
    a,b,c,d,e,f,g,h,i,j,k = create_all_nodes(hextuple)
    a.parent = b.parent = g
    c.parent = d.parent = j
    g.parent = j.parent = h
    h.parent = e.parent = i
    i.parent = f.parent = k
    g.left = a
    g.right = b
    j.left = c
    j.right = d
    h.left = g
    h.right = j 
    i.left = h
    i.right = e
    k.left = i
    k.right = f
    list2 = [(a,b),(c,d),(g,j),(h,e),(i,f)]

    return k,list2

def _tree3(hextuple):
    a,b,c,d,e,f,g,h,i,j,k = create_all_nodes(hextuple)
    a.parent = b.parent = g
    d.parent = e.parent = j
    g.parent = c.parent = h
    h.parent = j.parent = i
    i.parent = f.parent = k
    g.left = a
    g.right = b
    j.left = e
    j.right = d
    h.left = g
    h.right = c 
    i.left = h
    i.right = j
    k.left = i
    k.right = f
    list3 = [(a,b),(d,e),(g,c),(h,j),(i,f)]

    return k,list3

def _tree4(hextuple):
    a,b,c,d,e,f,g,h,i,j,k = create_all_nodes(hextuple)
    a.parent =b.parent = g
    e.parent = f.parent = j
    g.parent = c.parent = h
    h.parent = d.parent = i
    i.parent = j.parent = k
    g.left = a
    g.right = b
    j.left = e
    j.right = f
    h.left = g
    h.right = c 
    i.left = h
    i.right = d
    k.left = i
    k.right = j
    list4 = [(a,b),(e,f),(g,c),(h,d),(i,j)]

    return k,list4

def _tree6(hextuple):
    a,b,c,d,e,f,g,h,i,j,k = create_all_nodes(hextuple)
    a.parent = b.parent = g
    d.parent = e.parent = j
    g.parent = c.parent = h
    j.parent = f.parent = i
    i.parent = h.parent = k
    g.left = a
    g.right = b
    j.left = d
    j.right = e
    h.left = g
    h.left = c
    i.left = j
    i.right = f
    k.left = h
    k.right = i
    list6 = [(a,b),(d,e),(g,c),(j,f),(i,h)]

    return k,list6



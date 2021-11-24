# Comparing and Contrasting Phylogenetic Tree Building Algorithms 

A final project of Professor Layla Oesper's CS362, Computational Biology, at Carleton College, MN, USA. This project implements both UPGMA and Maximum Parsimony Algorithm with Sankoff's Algorithm. 

## Authors

- Zhen Ren
- Jimmy Zhong

## Instructor
- Professor Layla Oesper

## File Introduction
- main.sh
    - the shell script that runs the demo of both UPGMA and Max Parsimony with Sankoff's Algorithm.
- anaylsis_Sankoff.py
    - Implement and run parsimony and Sankoff's Algorithm to find phylogenetic tree.
- anaylsis_UPGMA.py
    - Run UPGMA by calling UPGMA.py to find phylogenetic tree.
- backtrace_Sankoff.py
    - Trace back the tree given by anaylsis_Sankoff.py and find the edge-length of the tree.
- UPGMA.py
    - Code for UPGMA
- Binarytree.py
    - Binary tree construction
- BuildPairList.py
    - Create all possible tree topology given a sequence. 
- sequence.py
    - Dataset of hemoglobin subunit alpha
    - Example sequence: Cat > MVLSAADKSNV KACWGKIGSH AGEYGAEALE RTFCSFPTTK TYFPHFDLSH GSAQVKAHGQ KVADALTQAV AHMDDLPTAM SALSDLHAYK LRVDPVNFKF LSHCLLVTLA CHHPAEFTPA VHASLDKFFS AVSTVLTSKY R
- sankoff_parsimony.txt
    - Output by anaylsis_Sankoff.py. 
    - The phylogenetic tree's parsimony score, leaves order and tree type that Sankoff's algorithm output has. 
- sankofftree.txt
    - Output by backtrace_Sankoff.py.
    - The inner node's sequences of the optimal phylogenetic tree that Sankoff's algorithm outputs. 
- sankoff_shape.txt
    - A hand typed visualization to describe the parsimony tree found by the Sankoff's algorithm in sankoff_parsimony.txt.
    - the numeric distance is later inputted into a newick format to visualize tree quantitatively.

## A New Data Structure for Binary Tree 
- Create a new Binary tree structure using list that facilitate the Sankoff's bottom up dynamic programing. 
- Representing a Binary tree topology by storing pairs of nodes that have the same parent together into a list. With the list, we do not need to transverse the binary tree structure to get inner nodes and leaves. 

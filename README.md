# Personal Projects in Bioinformatics
Some Bioinformatics/Biology related little programs for fun/laziness

# dynamic_programming_sequence_alignment
A program that uses Needleman–Wunsch algorithm to do sequence global alignment and local alignment

# comparing_Sankoff_parsimony_tree_with_UPGMA_tree
Extending on the sequence alignment project. Build phylogenetic tree from amino acid sequences of hemogloblin subunit alpha-a.
Build evolutionary tree of cat, human, chicken, macaw, 2 macaque speices using blosum62 matrix.
Compared and contrast 2 trees build by UPGMA algorithm and by brute force parsimony tree finding (using Sankoff Dynamic Programming to calculate parsimony score). 

# merge_name.py 
(work with Bowties' output, xxx_coverage.txt)

For my Bioinformatics class at Carleton College with amazing prof Rika Anderson, the prokka annotation in our server is always named by a machined code, such as "KGOCCMOE_00001" (from projects with Tara Dataset). These machine codes assure a variety of programs, such as Bowtie2, run smoothly, but if there is no way you will that know "KGOCCMOE_00001" is actually "Argininosuccinate synthase" (much more interesting)?
If you want to do that, put xxxx_coverage.txt and xxxx_ORFs.faa in the same folder as this script, direct to that folder in your terminal and run: python merge_name.py [ORF_coverage.txt] [ORFs.faa]
And there is a file with human readable names, as well as the relative abundance of that corresponding gene! 


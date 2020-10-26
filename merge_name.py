#Turn ORF codes in ORF_coverage.txt into human readable gene names
#input: python merge_names.py [ORF_coverage.txt] [ORFs.faa]
#output: a csv file, [ORF_coverage].csv

#parse out the ORF code(letter+number) and human readable names from .faa（i.e. ERR599142_ORFs.faa)
#then, merge the human readable names from .faa to ERR599142_ORF_coverage.txt

#Insipired by Professor Rika Anderson @Carleton College
#Modified by Jimmy(Juntao) Zhong at her Bioinformatic class. Anyone can use this!!!
import sys
import re
#use package "pandas" to perform merge, if not install: do pip...
#pip install pandas
import pandas as pd
txtfile = sys.argv[1]
fastafile = sys.argv[2]

txtMatrix = pd.read_csv(txtfile, sep="\t", names=['contig_name', 'Start Coor', 'Stop Coordinate', 'ORF_code', 'Sum of per_base coverage'])

readfasta = open(fastafile).read()
#split with "/n>" because there are "-->" in the name input
fastas = readfasta.split('\n>') #splits into fastas
fastas[0] = fastas[0][1:] #first one does not have '/n', only '>'
fastas = [fasta for fasta in fastas if fasta] #delete empty strings
#print(fastas[0:2])

#2D array OnlyNames.
#OnlyNames[i][0] = [ORF_code], OnlyNames[i][1] = [human_readable_names]. i is each ORF
OnlyNames = []
for fasta in fastas:
    code_and_name = [0]*2
    split_code = fasta.split(" ", 1)
    code_and_name[0] = split_code[0]
    split_ORF_name = split_code[1].splitlines()
    code_and_name[1] = split_ORF_name[0]
    OnlyNames.append(code_and_name)
    #print(fasta)

df_OnlyNames = pd.DataFrame(data = OnlyNames, columns=['ORF_code','human_readable_names'])
#print(df_OnlyNames)
#print(txtMatrix)

#merge step
OutMatrix = pd.merge(txtMatrix, df_OnlyNames, how='outer', on=['ORF_code'])

print("Nothing went wrong (so it's successful?); here is the first and last 3 code-to-human_readable_names in the .faa file:")
print(OnlyNames[-4:-1])
print(OnlyNames[0:3])

#write out file
txtfilename = txtfile.rstrip('.txt')
outfile = txtfilename + '.csv'
OutMatrix.to_csv(outfile, header=True)
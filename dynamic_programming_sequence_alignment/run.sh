#!/bin/bash
./local_alignment.py TP53-Cat.fasta TP53-Human.fasta scoring.txt
./global_alignment.py TP53-Cat.fasta TP53-Human.fasta scoring.txt

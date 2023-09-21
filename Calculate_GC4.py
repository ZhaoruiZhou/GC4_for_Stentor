#!/usr/bin/env python3
# -*- coding= UTF-8 -*-
# By Jerry Zhou 2023
# This script is used to calculate the GC4 of CDS
# usage: python Calculate_gc4.py input.fasta > GC4.txt
import sys
cds=sys.argv[1]
gencode = {
    'ATA':'', 'ATC':'', 'ATT':'', 'ATG':'',
    'ACA':'@', 'ACC':'*', 'ACG':'*', 'ACT':'@',
    'AAC':'', 'AAT':'', 'AAA':'', 'AAG':'',
    'AGC':'', 'AGT':'', 'AGA':'', 'AGG':'',
    'CTA':'@', 'CTC':'*', 'CTG':'*', 'CTT':'@',
    'CCA':'@', 'CCC':'*', 'CCG':'*', 'CCT':'@',
    'CAC':'', 'CAT':'', 'CAA':'', 'CAG':'',
    'CGA':'@', 'CGC':'*', 'CGG':'*', 'CGT':'@',
    'GTA':'@', 'GTC':'*', 'GTG':'*', 'GTT':'@',
    'GCA':'@', 'GCC':'*', 'GCG':'*', 'GCT':'@',
    'GAC':'', 'GAT':'', 'GAA':'', 'GAG':'',
    'GGA':'@', 'GGC':'*', 'GGG':'*', 'GGT':'@',
    'TCA':'@', 'TCC':'*', 'TCG':'*', 'TCT':'@',
    'TTC':'', 'TTT':'', 'TTA':'', 'TTG':'',
    'TAC':'', 'TAT':'', 'TAA':'', 'TAG':'',
    'TGC':'', 'TGT':'', 'TGA':'', 'TGG':''} #the 4-fold degenerate site == A/T ==@  && G/C ==*
def translate(dna):   #In this way, the GC content at the 4-fold degenerate site can be calculated
    amino_acid_sequence = ""
    for start in range(0,len(dna) - 2, 3):
        stop = start + 3
        codon = dna[start:stop]
        aa = gencode.get(codon.upper(),'X')
        amino_acid_sequence = amino_acid_sequence + aa
    return(amino_acid_sequence)

##############################################################
fa=cds

def readFa(fa):
    with open(fa,'r') as FA:
        seqName,seq='',''
        while 1:
            line=FA.readline()
            line=line.strip('\n')
            if (line.startswith('>') or not line) and seqName:
                yield((seqName,seq))
            if line.startswith('>'):
                seqName = line.split(" ")[0]
                seq=''
            else:
                seq+=line
            if not line:break

for seq in readFa(fa):
    cdsname=list(seq)[0]
    dna=list(seq)[1]
    fasta=translate(dna)
    at=fasta.count("@")
    gc=fasta.count("*")
    gc4=(gc/(gc+at))
    print(str(gc4))

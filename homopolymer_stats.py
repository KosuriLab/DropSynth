import csv
import Bio
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import re
import math

##################################
#INPUTS:
chip_file = "chip_out/DHFR_chip_3.csv"

bases = ['A','T','G','C']

minlength = 4
maxlength = 13

check_final_RC = True
##################################
polyA = dict()
polyT = dict()
polyG = dict()
polyC = dict()

csv_file_handle = open(chip_file)
csv_input_file = csv.DictReader(csv_file_handle)
for crec in csv_input_file:
    seqdesigned = crec["Sequence"]
    if str(seqdesigned).count('A') > str(seqdesigned).count('T') and check_final_RC:
        newoligo = str(Seq(seqdesigned,generic_dna).reverse_complement())
    else:
        newoligo = seqdesigned
    for n in range(minlength,maxlength):
        for ba in bases:
            search_str = ba * n
            #print(m)
            matchindx = newoligo.find(search_str)
            
            if matchindx != -1 and newoligo[matchindx+n] != ba:
                #print('match '+search_str +' in ' + crec["Name"])
                #print(newoligo)
                
                #found homopolymer of length n
                if ba == 'A':
                    polyA[n]= polyA.get(n, 0) + 1
                    
                    #print(y2)
                elif ba == 'T':
                    polyT[n]= polyT.get(n, 0) + 1
                    
                elif ba == 'G':
                    polyG[n]= polyG.get(n, 0) + 1
                    #print(y1)
                    
                elif ba == 'C':
                    polyC[n]= polyC.get(n, 0) + 1
                    
csv_file_handle.close()

print("polyA")
for aa in polyA:
    print(''+str(aa)+' '+str(polyA[aa]))
print("polyT")
for aa in polyT:
    print(''+str(aa)+' '+str(polyT[aa]))
print("polyG")
for aa in polyG:
    print(''+str(aa)+' '+str(polyG[aa]))
print("polyC")
for aa in polyC:
    print(''+str(aa)+' '+str(polyC[aa]))
    
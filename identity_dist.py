import Bio
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import time
import json
from Bio.Alphabet import IUPAC
import pylab
import csv
from Bio import Entrez
import os
import itertools
import math
from Bio import pairwise2

def getFastaSeqs(filename):
    fastqseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastqseqs.append(seqrec)
    handle.close()
    return fastqseqs

if __name__ == '__main__':
    start_time = time.time()
    
    ##################
    #seq in file
    seq_input_file = "db/DHFR_IP_targetlist_parsed_4_5_res_oligo_toFASTA.csv"

    # output file
    seq_indent_out_file = "db_ident/DHFR_IP_identity_dist.csv"

    ##################
    csv_file_handle = open(seq_input_file)
    csv_input_file = csv.DictReader(csv_file_handle)
    
    uniqseqdict = dict()
    accList = []
    seqList = []
    #fieldnames = ['Lib','AccID','Length','Source','TaxID','Taxa1','Taxa2','Taxa3','Taxa4','Definition','Sequence']
    
    current_seq_counter = 0
    for crec in csv_input_file:
        currentseq = crec["Sequence"]
        if not (currentseq in uniqseqdict):
            
            #get AccID
            accID = crec["AccID"]
            uniqseqdict[currentseq] = accID
            seqList.append(currentseq)
            accList.append(accID)
            
            
    csv_file_handle.close()
    
    print(str(len(seqList))+" unique seqs in CSV file.")
    
    pairsList = list(itertools.combinations(seqList, 2))
    
    print(str(len(pairsList)) +" combinations in pairsList.")
    
    csvfile = open(seq_indent_out_file, 'w')
    fieldnames = ['AccID1','AccID2','Identity']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    aligncount = 1
    for cpair in pairsList:
        alignment = pairwise2.align.globalxx(cpair[0], cpair[1], one_alignment_only=1, score_only=1)
        len1 = lengthDict[cpair[0]]
        len2 = lengthDict[cpair[1]]
        longer_prot_length = max(len1,len2)
        scaledid = alignment/longer_prot_length
        writer.writerow({'AccID1':str(uniqseqdict[str(cpair[0])]), 'AccID2':str(uniqseqdict[str(cpair[1])]), 'Identity':str(scaledid)})
        
        if aligncount % 5000 == 0:
            print(str(aligncount))
            print(str(alignment))
            
        aligncount += 1
        
    csvfile.close()
    print("--- %s seconds ---" % (time.time() - start_time))
    
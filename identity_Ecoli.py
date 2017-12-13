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
    
    ##########################
    #seq input
    seq_input_file = "db/DHFR_IP_targetlist_parsed_4_5_res_oligo_toFASTA.csv"
    #indent out file
    seq_indent_output_file = "db_ident/DHFR_IP_identity_Ec_Hu.csv"
    
    ##########################
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
    
    
    csvfile = open(seq_indent_output_file, 'w')
    fieldnames = ['AccID1','AccID2','IdentityEcoli','IdentityHuman']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    aligncount = 1
    for cpair in seqList:
        alignmentEc = pairwise2.align.globalxx('misliaalavdrvigmenampwnlpadlawfkrntlnkpvimgrhtwesigrplpgrkniilssqpgtddrvtwvksvdeaiaacgdvpeimvigggrvyeqflpkaqklylthidaevegdthfpdyepddwesvfsefhdadaqnshsycfeilerr', cpair, one_alignment_only=1, score_only=1)
        alignmentHu = pairwise2.align.globalxx('mvgslncivavsqnmgigkngdlpwpplrnefryfqrmtttssvegkqnlvimgkktwfsipeknrplkgrinlvlsrelkeppqgahflsrslddalklteqpelankvdmvwivggssvykeamnhpghlklfvtrimqdfesdtffpeidlekykllpeypgvlsdvqeekgikykfevyeknd', cpair, one_alignment_only=1, score_only=1)
        writer.writerow({'AccID1':str(uniqseqdict[str(cpair)]), 'IdentityEcoli':str(alignmentEc), 'IdentityHuman':str(alignmentHu)})
        
        if aligncount % 5000 == 0:
            print(str(aligncount))
            
        aligncount += 1
        
    csvfile.close()
    print("--- %s seconds ---" % (time.time() - start_time))
    
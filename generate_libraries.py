import Bio
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import time
from Bio.Alphabet import IUPAC
import pylab
import csv
import math
import random
import hashlib
from random import randrange
random.seed(15643243242342759)

##################################
#FUNCTIONS:
def getFastaSeqs(filename):
    fastqseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastqseqs.append(seqrec)
    handle.close()
    return fastqseqs

if __name__ == '__main__':
    start_time = time.time()
    
    #####################################
    ########## OPTIONS ##################
    #####################################

    min_nt = 440 #min length (4 oligo assemblies)
    max_nt_4 = 528 #for 230mer 4 oligo max is 531nt
    max_nt = 669 #max length (5 oligo assemblies)
    
    #database of sequences for E. coli DHFR Blast results
    ecoli_BLASTP_file = "seq_input/BlastP_Ecoli_DHFR.fasta"

    #where to save the output:
    output_file4 = "db/IP_protein_after_padding_to_384_4olgio.fasta"
    output_file5 = "db/IP_protein_after_padding_to_384_5olgio.fasta"
    output_file4ex = "db/IP_protein_after_padding_to_384_4olgio_extra.fasta"
    output_file5ex = "db/IP_protein_after_padding_to_384_5olgio_extra.fasta"

    #this is our database of sequences
    csv_file_handle = open("db/DHFR_IP_targetlist_parsed_4_5_res_oligo_toFASTA.csv")
    csv_input_file = csv.DictReader(csv_file_handle)
    
    #this contains a list of IDs for which the oligo split fails:
    IDs_to_remove = "db/IDs_to_remove.csv"

    #####################################
    ######### / OPTIONS #################
    #####################################

    #list of all unique sequences (seq is key, ID is value):
    uniqseqdict = dict()
    #list of all unique sequences for NCBI Blast seqs (seq is key, ID is value):
    uniqseqdict2 = dict()
    #list of all IDs:
    accList = []
    #list of all seqs:
    seqList = []
    #4 oligo assembly records:
    lib4 = []
    #4 oligo assembly records (resistant):
    lib4res = []
    #5 oligo assembly records:
    lib5 = []
    #5 oligo assembly records (resistant):
    lib5res = []
    #dictionary with seq as key and length as value:
    lengthDict = dict()
    #extra 4 oligo assembly records from NCBI Blast:
    added4seq = []
    #extra 5 oligo assembly records from NCBI Blast:
    added5seq = []
    #4 oligo assembly counter:
    o4count = 0
    #4 oligo assembly counter (resistant):
    o4rescount = 0
    #5 oligo assembly counter:
    o5count = 0
    #5 oligo assembly counter (resistant):
    o5rescount = 0
    
    with open(IDs_to_remove, 'r') as f:
            reader = csv.reader(f)
            lib_ident_list = list(reader)
    if len(lib_ident_list) > 0:
        badIDsarray = lib_ident_list[0]
    else:
        badIDsarray = []
    
    print("Bad IDs (not added to libs):")
    print(badIDsarray)
    
    current_seq_counter = 0
    #sort each sequence into different lists
    for crec in csv_input_file:
        #get the current sequence
        currentseq = crec["Sequence"].upper()
        accID = crec["AccID"]
        seqlength = len(currentseq)*3 #in bp
        
        #make sure it's unique
        #not a bad sequence
        #correct length
        #screen for ambigous a.a. "X"
        if (not (currentseq in uniqseqdict) and not (accID in badIDsarray) and 
            (seqlength >min_nt) and (seqlength <=max_nt) and not ('X' in currentseq)):
            
            recseq = Seq(currentseq.upper(),IUPAC.protein)
            record = SeqRecord(recseq,id=accID)
            uniqseqdict[currentseq] = accID
            seqList.append(currentseq)
            accList.append(accID)
            lengthDict[currentseq] = len(currentseq)

            #update counters and record lists
            if (crec["Lib"]== '4'):
                o4count += 1
                lib4.append(record)
            elif crec["Lib"]== '5':
                o5count += 1
                lib5.append(record)
            elif crec["Lib"]== '4res':
                o4rescount += 1
                lib4res.append(record)
            elif crec["Lib"]== '5res':
                o5rescount += 1
                lib5res.append(record)
        else:
            if (currentseq in uniqseqdict):
                print("Duplicate sequence: "+accID)
            if (accID in badIDsarray):
                print("Seq in Bad IDs list: "+accID)
            if ('X' in currentseq):
                print("Ambigous X a.a.: "+accID)

            
    csv_file_handle.close()
    
    #print stats
    print(str(len(seqList))+" unique seqs in CSV file.")
    print(str(o4count)+" - 4 oligo seqs in CSV file.")
    print(str(o4rescount)+" - 4 oligo (resistant) seqs in CSV file.")
    print(str(o5count)+" - 5 oligo seqs in CSV file.")
    print(str(o5rescount)+" - 5 oligo (resistant) seqs in CSV file.")
    
    #how many extra seuqnces to add from NCBI Blast list? for 4 oligo libs
    toadd_4 = 500
    #toadd_4 = math.ceil((o4count + o4rescount)/384)*384-(o4count+o4rescount)
    print(str(toadd_4)+" - 4 oligo seqs to add.")

    #how many extra seuqnces to add from NCBI Blast list? for 5 oligo libs
    toadd_5 = 500
    #toadd_5 = math.ceil((o5count + o5rescount)/384)*384-(o5count+o5rescount)
    print(str(toadd_5)+" - 5 oligo seqs to add.")
    
    #load e.coli DHFR Blast fasta:
    readseqs = getFastaSeqs(ecoli_BLASTP_file)
    count = 1
    
    #unused seqs (4 oligo):
    added4seqextra = []
    #unused seqs (5 oligo):
    added5seqextra = []

    for currentseq in readseqs:
        #check unique seq
        #dont use PDB sequences
        #dont use seq with ambigous a.a.
        #start with M
        
        if (not (currentseq.seq in uniqseqdict) and not ('X' in currentseq.seq) and
            not (currentseq.seq in uniqseqdict2) and not ('pdb' in currentseq.id) and
            not (currentseq.id in badIDsarray) and (currentseq.seq[0] == 'M')):
            currentseq_size = len(currentseq.seq)*3
            uniqseqdict2[currentseq.seq] = currentseq.id
            #check length:
            if (currentseq_size <= max_nt_4) and (currentseq_size > min_nt):# and (count4added < toadd_4):
                added4seq.append(currentseq)
            elif (currentseq_size <= max_nt) and (currentseq_size > max_nt_4):# and (count5added < toadd_5):
                added5seq.append(currentseq)
    
    #add resistant sequences (these will end up in Lib 13 for 4 oligo and Lib 15 for 5 oligo):
    for currentseq in lib4res:
        lib4.append(currentseq)
    for currentseq in lib5res:
        lib5.append(currentseq)
        
    count4added = 0
    count5added = 0

    #add new sequences, randomly, for 4 oligo lib:
    while count4added < toadd_4:
        random_index = randrange(0,len(added4seq))
        temp_rec = added4seq[random_index]
        if not (temp_rec.seq in uniqseqdict):
            lib4.append(temp_rec)
            count4added += 1
            uniqseqdict[temp_rec.seq] = temp_rec.id
    
    #add new sequences, randomly, for 5 oligo lib:
    while count5added < toadd_5:
        random_index = randrange(0,len(added5seq))
        temp_rec = added5seq[random_index]
        if not (temp_rec.seq in uniqseqdict):
            lib5.append(temp_rec)
            count5added += 1
            uniqseqdict[temp_rec.seq] = temp_rec.id
    
    #put extra unused sequences into a file:
    for currentseq in added4seq:
        if not currentseq.seq in uniqseqdict:
            added4seqextra.append(currentseq)
    for currentseq in added5seq:
        if not currentseq.seq in uniqseqdict:
            added5seqextra.append(currentseq)
    
    #save 4 and 5 oligo constructs:
    
    output_handle4 = open(output_file4, "w")
    SeqIO.write(lib4, output_handle4, "fasta")
    output_handle4.close()
    print("Wrote "+str(len(lib4))+" records for 4 oligo.")
    
    output_handle5 = open(output_file5, "w")
    SeqIO.write(lib5, output_handle5, "fasta")
    output_handle5.close()
    print("Wrote "+str(len(lib5))+" records for 5 oligo.")
    
    #save unused seqs:
    output_handle4ex = open(output_file4ex, "w")
    SeqIO.write(added4seqextra, output_handle4ex, "fasta")
    output_handle4ex.close()
    
    output_handle5ex = open(output_file5ex, "w")
    SeqIO.write(added5seqextra, output_handle5ex, "fasta")
    output_handle5ex.close()
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
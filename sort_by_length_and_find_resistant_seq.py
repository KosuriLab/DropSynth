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
    
    ##################################
    #VARIABLES:
    target_file = "seq_input/M20160921ACFE4208EAFA842A78A1B3BA7138A93DB4D25CK.list"
    seq_input_file = "db/DHFR_IP_targetlist_parsed_4_5_res_oligo.csv"
    entrez_dir = "folA_entrez/"
    min_nt = 300
    max_nt_4 = 531 #for 230mer 4 oligo max is 531nt
    max_nt = 670
    ##################################
    
    fastaparsed = []
    
    f = open(target_file, 'r')
    
    length_list_nt = []
    length_list_parsed_nt = []
    parsed_fasta = []
    parsed_accID = []
    uniqseqdict = dict()
    #save parsed seqs here:
    csvfile = open(seq_input_file, 'w')
    fieldnames = ['Lib','AccID','Length','Source','TaxID','Taxa1','Taxa2','Taxa3','Taxa4','Definition','Sequence']
    
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    counter_no_source = 0
    tax_count = 0
    
    current_seq_counter = 0
    for line in f:
        #get AccID
        accID = line.rstrip('\n')

        filename = entrez_dir+accID+".xml"
        if os.path.isfile(filename):#does XML file exist
            handle = open(filename, 'r')
            record = Entrez.read(handle)
            handle.close()
            
            currentlengthnt = int(record[0]["GBSeq_length"])*3
            currentseq = record[0]["GBSeq_sequence"]
            
            #make sure seq is in correct range and that it is unique
            if (currentlengthnt < max_nt and currentlengthnt > min_nt) and not (currentseq in uniqseqdict) and not ('X' in currentseq.upper()):
                
                #save this for later lookup to make sure all are unique
                uniqseqdict[currentseq] = accID
                GBSeq_feature_table = record[0]["GBSeq_feature-table"]
                
                recnum=0
                source_flag = 0
                taxid = []
                for items_feat_table in GBSeq_feature_table:
                    if items_feat_table["GBFeature_key"] == "source":
                        source_flag = 1
                        for items_feat_quals in items_feat_table["GBFeature_quals"]:
                            if items_feat_quals["GBQualifier_name"] == "db_xref":
                                #this is the taxonomic ID
                                taxid = items_feat_quals["GBQualifier_value"].split(":")[1]
                                tax_count += 1
                if source_flag == 0:
                    counter_no_source += 1
                
                #parsed_fasta.append(currentseq)
                parsed_accID.append(accID)
                length_list_parsed_nt.append(currentlengthnt)
                desc_split = record[0]["GBSeq_definition"].split('[')
                
                #will this be a 4 or 5 oligo assembly
                if currentlengthnt < max_nt_4:
                    LibStr = '4'
                else:
                    LibStr = '5'
                
                #is this seq resistant?
                if ('trimethoprim' in record[0]["GBSeq_definition"]) or ('resista' in record[0]["GBSeq_definition"]):
                    LibStr += 'res'
                
                #find taxa info
                taxa_split = record[0]["GBSeq_taxonomy"].split(';',4)
                taxa_length = len(taxa_split)
                if taxa_length == 1:
                    taxa_split_1 = record[0]["GBSeq_taxonomy"]
                    taxa_split_2 = ""
                    taxa_split_3 = ""
                    taxa_split_4 = ""
                elif taxa_length == 2:
                    taxa_split_1 = taxa_split[0].lstrip(' ')
                    taxa_split_2 = taxa_split[1].lstrip(' ')
                    taxa_split_3 = ""
                    taxa_split_4 = ""
                elif taxa_length == 3:
                    taxa_split_1 = taxa_split[0].lstrip(' ')
                    taxa_split_2 = taxa_split[1].lstrip(' ')
                    taxa_split_3 = taxa_split[2].lstrip(' ')
                    taxa_split_4 = ""
                elif taxa_length > 3:
                    taxa_split_1 = taxa_split[0].lstrip(' ')
                    taxa_split_2 = taxa_split[1].lstrip(' ')
                    taxa_split_3 = taxa_split[2].lstrip(' ')
                    taxa_split_4 = taxa_split[3].lstrip(' ')
                
                #save all info to csv file
                writer.writerow({'Lib':str(LibStr),'AccID':str(accID),'Length':str(currentlengthnt), 'Source': record[0]["GBSeq_source"],'TaxID':str(taxid), 'Taxa1': taxa_split_1, 'Taxa2': taxa_split_2, 'Taxa3': taxa_split_3, 'Taxa4': taxa_split_4, 'Definition': desc_split[0],'Sequence':str(currentseq)})
                recseq = Seq(currentseq.upper(),IUPAC.protein)
                record = SeqRecord(recseq,id=accID)
                parsed_fasta.append(record)
        current_seq_counter += 1
        if current_seq_counter % 200 == 0: #update every 200 seqs.
            print("Seq counter: " + str(current_seq_counter), end='\r')
    
    csvfile.close()
    print(str(len(parsed_fasta))+" seqs after parsing.")
    print(str(counter_no_source)+" seqs with no source found.")
    print(str(tax_count)+" seqs with taxID")
    print(str(min(length_list_parsed_nt))+" min length.")
    print(str(max(length_list_parsed_nt))+" max length.")
    
    
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
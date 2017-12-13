#python 3
#updated 10.28.2016 Calin Plesa
from Bio import SeqIO
import random
random.seed(15643243242342759)
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
from Bio.Restriction import *
import csv
from Bio import Entrez
import os
import time

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print seqrecord.id
        #print len(seqrecord)
    #print len(records)
    return records

######################################################################
Entrez.email = "plesa@ucla.edu"     # Always tell NCBI who you are

#these are the input files, one per 384 lib:
input_files = ["db_oligo/DHFR_Lib01_4oligo.proteins","db_oligo/DHFR_Lib02_4oligo.proteins",
                "db_oligo/DHFR_Lib03_4oligo.proteins","db_oligo/DHFR_Lib04_4oligo.proteins",
                "db_oligo/DHFR_Lib05_4oligo.proteins","db_oligo/DHFR_Lib06_4oligo.proteins",
                "db_oligo/DHFR_Lib07_4oligo.proteins","db_oligo/DHFR_Lib08_4oligo.proteins",
                "db_oligo/DHFR_Lib09_4oligo.proteins","db_oligo/DHFR_Lib10_4oligo.proteins",
                "db_oligo/DHFR_Lib11_4oligo.proteins","db_oligo/DHFR_Lib12_4oligo.proteins",
                "db_oligo/DHFR_Lib13_4oligo.proteins","db_oligo/DHFR_Lib14_5oligo.proteins",
                "db_oligo/DHFR_Lib15_5oligo.proteins"]

output_folder = 'db_oligo'

#output file
csv_info_file_out = "db/DHFR_Lib_Info_All.csv"

# number of oligos to use to split genes in each lib:
num_oligos = [4,4,4,4,4,4,4,4,4,4,4,4,4,5,5]

#generate two different codon optimizations for each ptoetin?
gen_two_codon = True

entrez_dir = "folA_entrez/"

#files for BLAT analysis
protein_all_unique_aa_fileout_name = "db/DHFR_aa_proteins_LibAll_unique.fasta"
protein_all_unique_nt_fileout_name = "db/DHFR_nt_LibAll_unique.fasta"

######################################################################
csvfile = open(csv_info_file_out, 'w')
fieldnames = ['AccID','LibNumCodon1','LibNumCodon2','NumOligos','Lengthnt','LengthntREBufPrim','Source','TaxID','Taxa1','Taxa2','Taxa3','Taxa4','Definition','Sequence','c12ident','ntCodon1','ntCodon2','ntCodon1REBufPrim','ntCodon2REBufPrim']
#fieldnames = ['AccID','Length','Sequence']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()
total_lib_per_codon = len(input_files)
lib_counter = 1
file_extension = "proteins"
seqlist_for_unique = []
seqdictuniq = dict()
seqdictuniqlib = dict()
masterunique_seq_rec_nt = []
masterunique_seq_rec_aa = []
for libfile in input_files:
    print("Lib "+str(lib_counter)+":")
    #grab all record in the library:
    aarecords = obtainFastaSequences(input_files[lib_counter-1])
    
    gene_nt_fileout = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'full_wRE_wPrim.genes')
    gene_nt_no_primer_or_RE_fileout = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'full_nRE_nPrim.genes')
    gene_nt_no_primer_fileout = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'full_wRE_noPrim.genes')
    protein_nt_fileout = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'proteins')
    
    ntrecord1 = obtainFastaSequences(gene_nt_no_primer_or_RE_fileout)
    ntrecord1REBufPrim = obtainFastaSequences(gene_nt_fileout)
    
    if gen_two_codon:
        gene_nt_fileout2 = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'codon2.full_wRE_wPrim.genes')
        gene_nt_no_primer_or_RE_fileout2 = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'codon2.full_nRE_nPrim.genes')
        gene_nt_no_primer_fileout2 = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'codon2.full_wRE_noPrim.genes')
        protein_nt_fileout2 = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'codon2.proteins')
        
        ntrecord2 = obtainFastaSequences(gene_nt_no_primer_or_RE_fileout2)
        ntrecord2REBufPrim = obtainFastaSequences(gene_nt_fileout2)
        proteinnt_ident_fileout = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'lib_ident.csv')
        with open(proteinnt_ident_fileout, 'r') as f:
            reader = csv.reader(f)
            lib_ident_list = list(reader)
    
    final_records_gene_level = []
    final_records_gene_level_no_primer_or_RE = []
    final_records_gene_level_no_primer = []
    final_records_protein_level = []
    
    pct_ident_array = []
    
    if gen_two_codon:
        final_records_gene_level2 = []
        final_records_gene_level_no_primer_or_RE2 = []
        final_records_gene_level_no_primer2 = []
        final_records_protein_level2 = []
    
    #split oligo files here:
    splitoligo1_file = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'oligos')
    if gen_two_codon:
        splitoligo2_file = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace(file_extension,'codon2.oligos')
        
    rec_count = 1
    total_rec_count = len(aarecords)
    
    #fieldnames = ['AccID','LibNumCodon1','LibNumCodon2','NumOligos','Length','Source','TaxID','Taxa1','Taxa2','Taxa3','Taxa4','Definition','Sequence','ntCodon1','ntCodon2']
    
    tax_count = 0
    #for each protein record:
    for ii in range(len(aarecords)):
        caarec = aarecords[ii]
        accID=caarec.id
        if caarec.seq in seqlist_for_unique:
            print("Duplicate sequence found, ID:"+accID+". Seq: "+caarec.seq)
            print("Matches with "+seqdictuniq[caarec.seq]+" in Lib "+str(seqdictuniqlib[caarec.seq]))
        else:
            seqlist_for_unique.append(caarec.seq)
            seqdictuniq[caarec.seq] = accID
            seqdictuniqlib[caarec.seq] = lib_counter
            masterunique_seq_rec_aa.append(caarec)
            
            filename = entrez_dir+accID+".xml"
            if os.path.isfile(filename):#does XML file exist
                handle = open(filename, 'r')
                record = Entrez.read(handle)
                handle.close()
                #print("ii: "+str(ii)+". len(ntrecord1):"+str(len(ntrecord1)))
                nt1rec = ntrecord1[ii]
                nt2rec = ntrecord2[ii]
                
                nt1REBufPrimrec = ntrecord1REBufPrim[ii]
                nt2REBufPrimrec = ntrecord2REBufPrim[ii]
                
                nt1rec.id = nt1rec.id + ";C1;"+"Lib"+str(lib_counter)
                nt2rec.id = nt2rec.id + ";C2;"+"Lib"+str(lib_counter+15)
                
                masterunique_seq_rec_nt.append(nt1rec)
                masterunique_seq_rec_nt.append(nt2rec)
                
                reclengthnt = len(caarec)*3
                
                #########################################################################
                GBSeq_feature_table = record[0]["GBSeq_feature-table"]
                
                recnum=0
                source_flag = 0
                taxid = []
                for items_feat_table in GBSeq_feature_table:
                    if items_feat_table["GBFeature_key"] == "source":
                        source_flag = 1
                        for items_feat_quals in items_feat_table["GBFeature_quals"]:
                            if items_feat_quals["GBQualifier_name"] == "db_xref":
                                taxid = items_feat_quals["GBQualifier_value"].split(":")[1]
                                tax_count += 1
                if source_flag == 0:
                    counter_no_source += 1
            
                desc_split = record[0]["GBSeq_definition"].split('[')
            
            
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
                #print("libcount: "+str(lib_counter)+". rec_count: "+str(rec_count))
                #print(lib_ident_list[0][rec_count-1])
                #print("lib_counter:"+str(lib_counter))
                #print("rec_count:"+str(rec_count))
                #print(lib_ident_list[0][rec_count-1])
                #print(desc_split[0])
                #print(record[0]["GBSeq_source"])
                #ntrecord1REBufPrim
                #print(str(nt1REBufPrimrec.seq))
                writer.writerow({'AccID':str(caarec.id),'LibNumCodon1':str(lib_counter),'LibNumCodon2':str(lib_counter+total_lib_per_codon),'NumOligos':str(num_oligos[lib_counter-1]),'Lengthnt':str(reclengthnt),'LengthntREBufPrim':str(len(nt1REBufPrimrec.seq)), 'Source': record[0]["GBSeq_source"],'TaxID':str(taxid), 'Taxa1': taxa_split_1, 'Taxa2': taxa_split_2, 'Taxa3': taxa_split_3, 'Taxa4': taxa_split_4, 'Definition': desc_split[0],'Sequence':str(caarec.seq),'c12ident':lib_ident_list[0][rec_count-1],'ntCodon1':str(nt1rec.seq),'ntCodon2':str(nt2rec.seq),'ntCodon1REBufPrim':str(nt1REBufPrimrec.seq),'ntCodon2REBufPrim':str(nt2REBufPrimrec.seq)})
                rec_count += 1
            else:
                print("Error: could not open XML file for: "+accID+" Trying to download from Entrez now.")
                if not 'pdb' in accID:
                    # Downloading...
                    net_handle = Entrez.efetch(db="protein",id=accID,retmode="xml")
                    out_handle = open(filename, "w")
                    out_handle.write(net_handle.read())
                    out_handle.close()
                    net_handle.close()
                    print("Saved "+accID)
                    time.sleep(0.5)
    
    print("Processed "+str(rec_count-1)+" records.")
    
    lib_counter += 1

csvfile.close()

#save files with all sequences for use with BLAT:
protein_all_unique_aa_fileout = open(protein_all_unique_aa_fileout_name,'w')
SeqIO.write(masterunique_seq_rec_aa, protein_all_unique_aa_fileout, "fasta")
protein_all_unique_nt_fileout = open(protein_all_unique_nt_fileout_name,'w')
SeqIO.write(masterunique_seq_rec_nt, protein_all_unique_nt_fileout, "fasta")

import Bio
import csv
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Restriction import *
import time
import os
from Bio import pairwise2

start_time = time.time() #time the run

##################################
#FUNCTIONS:
def getFastaSeqs(filename):
    fastaseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastaseqs.append(seqrec)
    handle.close()
    return fastaseqs

##################################
#INPUTS:
inputfiles = ['db_oligo/DHFR_Lib01_4oligo.proteins','db_oligo/DHFR_Lib02_4oligo.proteins','db_oligo/DHFR_Lib03_4oligo.proteins',
                'db_oligo/DHFR_Lib04_4oligo.proteins','db_oligo/DHFR_Lib05_4oligo.proteins','db_oligo/DHFR_Lib06_4oligo.proteins',
                'db_oligo/DHFR_Lib07_4oligo.proteins','db_oligo/DHFR_Lib08_4oligo.proteins','db_oligo/DHFR_Lib09_4oligo.proteins',
                'db_oligo/DHFR_Lib10_4oligo.proteins','db_oligo/DHFR_Lib11_4oligo.proteins','db_oligo/DHFR_Lib12_4oligo.proteins',
                'db_oligo/DHFR_Lib13_4oligo.proteins','db_oligo/DHFR_Lib14_5oligo.proteins','db_oligo/DHFR_Lib15_5oligo.proteins']

#inputfiles = ['db_oligo/DHFR_Lib15_5oligo.proteins']

barcodes = getFastaSeqs('barcodes/filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction.fasta')

##################################
#OUTPUTS:
csv_output_name = 'chip_out/DHFR_chip_3_proteins.csv'

all_lib_fasta_out = "chip_out/DHFR_all_lib_clean_header.fasta"

#this contains the ID as key and percent ident to ecoli as value:
identity_to_ecoli = dict()

##################################
#VARIABLES:
counter_no_source = 0
tax_count = 0

#list of all unique perfect IDs:
IDlist = []

#all protein records:
buildseqs = []

Ecoliseq = Seq("MISLIAALAVDRVIGMENAMPWNLPADLAWFKRNTLNKPVIMGRHTWESIGRPLPGRKNIILSSQPGTDDRVTWVKSVDEAIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFEILERR",generic_protein)
Ecolirecord = SeqRecord(Ecoliseq,id="NP_414590",description="")

entrez_dir = "folA_entrez/"

##################################
#CODE:

csvfile = open(csv_output_name, 'w')
fieldnames = ['Lib#codon1','Lib#codon2','Construct#','Barcode','PctIdentEcoli','Accession','Source','Definition','TaxID','Taxa1','Taxa2','Taxa3','Taxa4','Sequence']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for fileindex in range(len(inputfiles)):
    constructs = getFastaSeqs(inputfiles[fileindex])#.replace('.oligos','-finaloligos.fasta')

    #loop over all proteins in fasta file
    for index in range(len(constructs)):

        proteinrec = constructs[index]
        buildseqs.append(SeqRecord(proteinrec.seq,id=proteinrec.id,description=""))
        barcode_current = barcodes[index]
        seqid = proteinrec.id
        #print(proteinrec)
        #print(seqid)
        IDlist.append(seqid)

        xmlfilename = entrez_dir+seqid+".xml"

        alignment = pairwise2.align.globalxx(Ecolirecord.seq, proteinrec.seq, one_alignment_only=1, score_only=1)
        if len(Ecolirecord.seq) > len(proteinrec.seq):
            pct_ident = alignment/len(Ecolirecord.seq)
        else:
            pct_ident = alignment/len(proteinrec.seq)

        #does XML file exist:
        if os.path.isfile(xmlfilename):
            handle = open(xmlfilename, 'r')
            record = Entrez.read(handle)
            handle.close()
            #currentseq = record[0]["GBSeq_sequence"]
            GBSeq_feature_table = record[0]["GBSeq_feature-table"]
            
            source_flag = 0
            taxid = []
            for items_feat_table in GBSeq_feature_table:
                #does the record list a source organism:
                if items_feat_table["GBFeature_key"] == "source":
                    source_flag = 1
                    for items_feat_quals in items_feat_table["GBFeature_quals"]:
                        #does it have a taxon ID:
                        if items_feat_quals["GBQualifier_name"] == "db_xref":
                            taxid = items_feat_quals["GBQualifier_value"].split(":")[1]
                            tax_count += 1
            #count the number of records with no source organism:
            if source_flag == 0:
                counter_no_source += 1
            
            #taxonomy split for first four levels:
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
        else:
            print(ID+" XML file not found!")
        #write out info
        #'Lib#codon1','Lib#codon2','Construct#','Barcode','PctIdentEcoli','Accession','Source','Definition','TaxID','Taxa1','Taxa2','Taxa3','Taxa4','Sequence'
        writer.writerow({'Lib#codon1':str(fileindex+1),'Lib#codon2':str(fileindex+16),'Construct#':str(index+1),
                        'Barcode':barcode_current.id,'PctIdentEcoli':str(pct_ident),'Accession': seqid, 
                        'Source': record[0]["GBSeq_source"], 'Definition': record[0]["GBSeq_definition"],
                        'TaxID':str(taxid),'Taxa1': taxa_split_1,'Taxa2': taxa_split_2,'Taxa3': taxa_split_3,
                        'Taxa4': taxa_split_4, 'Sequence': proteinrec.seq})
        print("Construct "+str(index+1)+ " Lib " + str(fileindex+1), end='\r')

csvfile.close()
print("Found TaxID info for "+str(tax_count) +" of "+str(len(IDlist))+" records.")
#######################################################

print("Writing records to "+all_lib_fasta_out)
output_handle5ex = open(all_lib_fasta_out, "w")
SeqIO.write(buildseqs, output_handle5ex, "fasta")
output_handle5ex.close()

#######################################################

print("--- %s seconds ---" % (time.time() - start_time))

#######################################################


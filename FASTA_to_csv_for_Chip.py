from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Restriction import *
import csv

##################################
#FUNCTIONS:
def getOligos(filename):
    constructs = []
    currentconstruct = 'foo'
    for seqrec in SeqIO.parse(filename, "fasta"):
        if not seqrec.id.startswith(currentconstruct):
            currentconstruct = seqrec.id[0:seqrec.id.rfind(';')]
            constructs.append([currentconstruct])
            constructs[-1].append(seqrec)
        else:
            constructs[-1].append(seqrec)
    return constructs

##################################
#INPUTS:
inputfiles = ['db_oligo/DHFR_Lib01_4oligo.oligos','db_oligo/DHFR_Lib02_4oligo.oligos','db_oligo/DHFR_Lib03_4oligo.oligos',
                'db_oligo/DHFR_Lib04_4oligo.oligos','db_oligo/DHFR_Lib05_4oligo.oligos','db_oligo/DHFR_Lib06_4oligo.oligos',
                'db_oligo/DHFR_Lib07_4oligo.oligos','db_oligo/DHFR_Lib08_4oligo.oligos','db_oligo/DHFR_Lib09_4oligo.oligos',
                'db_oligo/DHFR_Lib10_4oligo.oligos','db_oligo/DHFR_Lib11_4oligo.oligos','db_oligo/DHFR_Lib12_4oligo.oligos',
                'db_oligo/DHFR_Lib13_4oligo.oligos','db_oligo/DHFR_Lib14_5oligo.oligos','db_oligo/DHFR_Lib15_5oligo.oligos',
                'db_oligo/DHFR_Lib16_4oligo.oligos','db_oligo/DHFR_Lib17_4oligo.oligos','db_oligo/DHFR_Lib18_4oligo.oligos',
                'db_oligo/DHFR_Lib19_4oligo.oligos','db_oligo/DHFR_Lib20_4oligo.oligos','db_oligo/DHFR_Lib21_4oligo.oligos',
                'db_oligo/DHFR_Lib22_4oligo.oligos','db_oligo/DHFR_Lib23_4oligo.oligos','db_oligo/DHFR_Lib24_4oligo.oligos',
                'db_oligo/DHFR_Lib25_4oligo.oligos','db_oligo/DHFR_Lib26_4oligo.oligos','db_oligo/DHFR_Lib27_4oligo.oligos',
                'db_oligo/DHFR_Lib28_4oligo.oligos','db_oligo/DHFR_Lib29_5oligo.oligos','db_oligo/DHFR_Lib30_5oligo.oligos']


file_out_1 = "chip_out/DHFR_chip_3.csv"
file_out_2 = "chip_out/DHFR_chip_3_raw.csv"

entrez_dir = "folA_entrez/"

##################################

csvfile = open(file_out_1, 'w')
fieldnames = ['Lib#','Construct#','Oligo#','Barcode','Accession','Source','Definition','Sequence']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for fileindex in range(len(inputfiles)):
    constructs = getOligos(inputfiles[fileindex].replace('.oligos','-finaloligos.fasta'))
    first_oligo_counter = 0
    middle_oligo_counter = 0
    last_oligo_counter = 0
    for index in range(len(constructs)):
        oligos = constructs[index][1:]
        oligo_counter=1
        #temp_str = constructs[index][0]
        temp_str = constructs[index][0].split(';')
        seqid = temp_str[0]
        barcode_id = temp_str[1]
        xmlfilename = entrez_dir+seqid+".xml"
        handle = open(xmlfilename, 'r')
        record = Entrez.read(handle)
        handle.close()
        for oligoseqrec in oligos:
            if oligo_counter == 1:
                writer.writerow({'Lib#':str(fileindex+1),'Construct#':str(index+1),'Oligo#':str(oligo_counter),'Barcode':barcode_id,'Accession': seqid, 'Source': record[0]["GBSeq_source"], 'Definition': record[0]["GBSeq_definition"], 'Sequence': oligoseqrec.seq})
            else:
                writer.writerow({'Lib#':'', 'Construct#':'','Oligo#':str(oligo_counter),'Barcode': '','Accession':'', 'Source': '', 'Definition': '', 'Sequence': oligoseqrec.seq})
            oligo_counter += 1

csvfile.close()

csvfile = open(file_out_2, 'w')
fieldnames = ['Name','Sequence']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#writer.writeheader()

for fileindex in range(len(inputfiles)):
    constructs = getOligos(inputfiles[fileindex].replace('.oligos','-finaloligos.fasta'))
    first_oligo_counter = 0
    middle_oligo_counter = 0
    last_oligo_counter = 0
    for index in range(len(constructs)):
        oligos = constructs[index][1:]
        oligo_counter=1
        for oligoseqrec in oligos:
            temp_str = constructs[index][0]+';'+str(oligo_counter)
            writer.writerow({'Name':temp_str, 'Sequence': oligoseqrec.seq})
            oligo_counter += 1

csvfile.close()
#this script does the following:

#for each library
#get oligos
#add typeIIs cutters (BtsI)
#check that no new restriction sites added
#pad the oligos to same length
#add barcode with nicking sites (Nt.BspQI)
#check that no new restriction sites added
#add amplification primers
#check that no new restriction sites added
#output

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Restriction import *
import pylab

##################################
#FUNCTIONS:
def getBuildOligos(filename,num_oligos):
    handle = open(filename,'r')
    buildoligos = []
    countOligos = 0
    localCount = num_oligos
    for line in handle:
        if line.startswith('>'):
            #make sure the previous sequence had the correct number of oligos
            if localCount != num_oligos:
                print("Wrong num oligos. Before: "+str(buildoligos[-1]))
                print("Total oligos: " +str(localCount)+". Expected: "+str(num_oligos))
            buildoligos.append([line.strip()])
            localCount = 0
        elif line != "":
            buildoligos[-1].append(Seq(line.strip(),generic_dna))
            countOligos += 1
            localCount += 1
        
    handle.close()
    return buildoligos

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records

def addBtsI(constructs):
    #add BtsI sites and check that only two exist
    newconstructs = []
    for construct in constructs:
        newconstruct = []
        newconstruct.append(construct[0])
        for i in range(1,len(construct)):
            newoligo = Seq("GCAGTG",generic_dna) + construct[i] + Seq("CACTGC",generic_dna)
            rb = RestrictionBatch([BtsI, BspQI, NdeI, KpnI,EcoRI])
            seqsearch = rb.search(newoligo)
            if i==1:
                #first oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 1 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\naddBtsI: Bad number of restriction sites in first oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
                    
            elif i==(len(construct)-1):
                #last oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=1:# or len(seqsearch[EcoRI])!=0:
                    print("\naddBtsI: Bad number of restriction sites in last oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            else:
                #middle oligos
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\naddBtsI: Bad number of restriction sites in middle oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs

def plotOligo_length(constructs,filename):
    #plot histogram of lengths
    construct_lengths = []
    for construct in constructs:
        for i in range(1,len(construct)):
            construct_lengths.append(len(construct[i]))
    generate_histogram = True
    if generate_histogram:
        data = construct_lengths
        print("min Lev dist: " +str(min(data)) + " max Lev dist: " +str(max(data)))
        pylab.hist(data, bins=(max(data)-min(data)))
        pylab.title("%i oligos length distribution\nfrom %i to %i" \
                    % (len(data),min(data),max(data)))
        pylab.xlabel("Oligo length (nt)")
        pylab.ylabel("Count")
        pylab.savefig(filename+'.png', bbox_inches='tight')
        pylab.savefig(filename+'.pdf', bbox_inches='tight')
        pylab.show()
    
def padOligoCP(constructs,totallength):
    #add bases until oligo is correct length
    newconstructs = []
    for construct in constructs:
        newconstruct = []
        newconstruct.append(construct[0])
        for i in range(1,len(construct)):
            #pad oligo ATGC
            oligo_len = len(construct[i])
            length_to_add = totallength - oligo_len
            full_seq = length_to_add//4
            full_seq_mod = length_to_add % 4
            
            if full_seq_mod == 3:
                last_padding = "ATG"
            elif full_seq_mod == 2:
                last_padding = "AT"
            elif full_seq_mod == 1:
                last_padding = "A"
            else:
                last_padding = ""
            
            full_padding_seq = "ATGC" * full_seq + last_padding
            newoligo =  Seq(full_padding_seq,generic_dna) + construct[i] # add padding between BspQI and BtsI
            rb = RestrictionBatch([BtsI, BspQI, NdeI, KpnI, EcoRI])
            seqsearch = rb.search(newoligo)
            if i==1:
                #first oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 1 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\npadOligo: Bad number of restriction sites in first oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            elif i==(len(construct)-1):
                #last oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=1:# or len(seqsearch[EcoRI])!=0:
                    print("\npadOligo: Bad number of restriction sites in last oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            else:
                #middle oligos
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\npadOligo: Bad number of restriction sites in middle oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs

def addBarcodes(constructs, barcodes):
    newconstructs = []
    for index in range(len(constructs)):
        newconstruct = []
        newconstruct.append(constructs[index][0]+';'+barcodes[index].id)
        for i in range(1,len(constructs[index])):
            newoligo = Seq("GCTCTTCG",generic_dna) + barcodes[index].seq + Seq("CGAAGAGC",generic_dna) + constructs[index][i]
            rb = RestrictionBatch([BtsI, BspQI, NdeI, KpnI,EcoRI])
            seqsearch = rb.search(newoligo)
            if i==1:
                #first oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=2 or len(seqsearch[NdeI])!= 1 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\naddBarcodes: Bad number of restriction sites in first oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
                    
            elif i==(len(constructs[index])-1):
                #last oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=2 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=1:# or len(seqsearch[EcoRI])!=0:
                    print("\naddBarcodes: Bad number of restriction sites in last oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            else:
                #middle oligos
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=2 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\naddBarcodes: Bad number of restriction sites in middle oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs

def addAmpPrimers(constructs,fwdprim,revprim):
    newconstructs = []
    for index in range(len(constructs)):
        newconstruct = []
        newconstruct.append(constructs[index][0])
        for i in range(1,len(constructs[index])):
            newoligo = fwdprim + constructs[index][i] + revprim.reverse_complement()
            rb = RestrictionBatch([BtsI, BspQI, NdeI, KpnI, EcoRI])
            seqsearch = rb.search(newoligo)
            if i==1:
                #first oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=2 or len(seqsearch[NdeI])!= 1 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\naddAmpPrimers: Bad number of restriction sites in first oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
                    
            elif i==(len(constructs[index])-1):
                #last oligo
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=2 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=1:# or len(seqsearch[EcoRI])!=0:
                    print("\naddAmpPrimers: Bad number of restriction sites in last oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            else:
                #middle oligos
                if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=2 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0:# or len(seqsearch[EcoRI])!=0:
                    print("\naddAmpPrimers: Bad number of restriction sites in middle oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                    print(seqsearch)
            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs


def printOligoLibrary(constructs):
    for construct in constructs:
        print(construct[0])
        for oligo in construct[1:]:
            print('\t' + oligo + '\t' + str(len(oligo)))
            

def outputOligoLibrary(constructs,filename,append_or_write,index):
    handle = open(filename,append_or_write)
    if index < 15:
        codon_str = "Codon1"
    else:
        codon_str = "Codon2"
    for construct in constructs:
        basename = construct[0] + ';' + codon_str + ';'
        for index in range(1,len(construct)):
            handle.write(basename+str(index)+'\n')
            handle.write(str(construct[index])+'\n')
    handle.close()

#####################################
########## OPTIONS ##################
#####################################

filenames = ['db_oligo/DHFR_Lib01_4oligo.oligos','db_oligo/DHFR_Lib02_4oligo.oligos','db_oligo/DHFR_Lib03_4oligo.oligos',
                'db_oligo/DHFR_Lib04_4oligo.oligos','db_oligo/DHFR_Lib05_4oligo.oligos','db_oligo/DHFR_Lib06_4oligo.oligos',
                'db_oligo/DHFR_Lib07_4oligo.oligos','db_oligo/DHFR_Lib08_4oligo.oligos','db_oligo/DHFR_Lib09_4oligo.oligos',
                'db_oligo/DHFR_Lib10_4oligo.oligos','db_oligo/DHFR_Lib11_4oligo.oligos','db_oligo/DHFR_Lib12_4oligo.oligos',
                'db_oligo/DHFR_Lib13_4oligo.oligos','db_oligo/DHFR_Lib14_5oligo.oligos','db_oligo/DHFR_Lib15_5oligo.oligos',
                'db_oligo/DHFR_Lib01_4oligo.codon2.oligos','db_oligo/DHFR_Lib02_4oligo.codon2.oligos','db_oligo/DHFR_Lib03_4oligo.codon2.oligos',
                'db_oligo/DHFR_Lib04_4oligo.codon2.oligos','db_oligo/DHFR_Lib05_4oligo.codon2.oligos','db_oligo/DHFR_Lib06_4oligo.codon2.oligos',
                'db_oligo/DHFR_Lib07_4oligo.codon2.oligos','db_oligo/DHFR_Lib08_4oligo.codon2.oligos','db_oligo/DHFR_Lib09_4oligo.codon2.oligos',
                'db_oligo/DHFR_Lib10_4oligo.codon2.oligos','db_oligo/DHFR_Lib11_4oligo.codon2.oligos','db_oligo/DHFR_Lib12_4oligo.codon2.oligos',
                'db_oligo/DHFR_Lib13_4oligo.codon2.oligos','db_oligo/DHFR_Lib14_5oligo.codon2.oligos','db_oligo/DHFR_Lib15_5oligo.codon2.oligos']

filenames_out = ['db_oligo/DHFR_Lib01_4oligo.oligos','db_oligo/DHFR_Lib02_4oligo.oligos','db_oligo/DHFR_Lib03_4oligo.oligos',
                'db_oligo/DHFR_Lib04_4oligo.oligos','db_oligo/DHFR_Lib05_4oligo.oligos','db_oligo/DHFR_Lib06_4oligo.oligos',
                'db_oligo/DHFR_Lib07_4oligo.oligos','db_oligo/DHFR_Lib08_4oligo.oligos','db_oligo/DHFR_Lib09_4oligo.oligos',
                'db_oligo/DHFR_Lib10_4oligo.oligos','db_oligo/DHFR_Lib11_4oligo.oligos','db_oligo/DHFR_Lib12_4oligo.oligos',
                'db_oligo/DHFR_Lib13_4oligo.oligos','db_oligo/DHFR_Lib14_5oligo.oligos','db_oligo/DHFR_Lib15_5oligo.oligos',
                'db_oligo/DHFR_Lib16_4oligo.oligos','db_oligo/DHFR_Lib17_4oligo.oligos','db_oligo/DHFR_Lib18_4oligo.oligos',
                'db_oligo/DHFR_Lib19_4oligo.oligos','db_oligo/DHFR_Lib20_4oligo.oligos','db_oligo/DHFR_Lib21_4oligo.oligos',
                'db_oligo/DHFR_Lib22_4oligo.oligos','db_oligo/DHFR_Lib23_4oligo.oligos','db_oligo/DHFR_Lib24_4oligo.oligos',
                'db_oligo/DHFR_Lib25_4oligo.oligos','db_oligo/DHFR_Lib26_4oligo.oligos','db_oligo/DHFR_Lib27_4oligo.oligos',
                'db_oligo/DHFR_Lib28_4oligo.oligos','db_oligo/DHFR_Lib29_5oligo.oligos','db_oligo/DHFR_Lib30_5oligo.oligos']

#ampF from 00_primer_screen.py output: skpp15 F&R for amplification primers
#these have been further edited manually
ampprimf_file = 'ampprimers-skpp15/skpp15-forward_select_mod.faa'
ampprimersf = obtainFastaSequences(ampprimf_file)

ampylist = []
for ampy in ampprimersf:
    ampylist.append(str(ampy.seq))
print(ampylist)

#ampR (not RC) from 00_primer_screen.py output:
ampprimr_file = 'ampprimers-skpp15/skpp15-reverse_select_mod.faa'
ampprimersr = obtainFastaSequences(ampprimr_file)

ampylist = []
for ampy in ampprimersr:
    ampylist.append(str(ampy.seq))
print(ampylist)

barcodes = obtainFastaSequences('barcodes/filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction.fasta')

# number of oligos to use to split genes in each lib:
num_oligos = [4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5]

constructs_per_lib = [384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384]

#length of payload + BtsaI sites + buffer such that all oligos are full length
length_padded_payload = 172 # 142nt for 200mer with 12nt BC, 172nt for 230mers

file_for_BLAT = "db_oligo/LibAll_finaloligos_noAmp.fasta"

#####################################
######### / OPTIONS #################
#####################################

for index in range(len(filenames)):
    print("Processing Lib"+str(index+1))
    buildoligos = getBuildOligos(filenames[index],num_oligos[index])
    #plotOligo_length(buildoligos,"plots/Oligo_Lib"+str(index)+"_length_in_hist")
    btsIoligos = addBtsI(buildoligos)
    #plotOligo_length(btsIoligos,"plots/Oligo_Lib"+str(index)+"_length_after_btsI_hist")
    paddedOligos = padOligoCP(btsIoligos,length_padded_payload)
    #plotOligo_length(paddedOligos,"plots/Oligo_"+str(index)+"_length_after_padding_hist")
    barcodedoligos = addBarcodes(paddedOligos,barcodes)
    if index == 0:
        outputOligoLibrary(barcodedoligos,file_for_BLAT,'w',index)
    else:
        outputOligoLibrary(barcodedoligos,file_for_BLAT,'a',index)
    finaloligos = addAmpPrimers(barcodedoligos,ampprimersf[index].seq,ampprimersr[index].seq)
    #printOligoLibrary(finaloligos)
    print(filenames_out[index].replace('.oligos','-finaloligos.fasta'))
    outputOligoLibrary(finaloligos,filenames_out[index].replace('.oligos','-finaloligos.fasta'),'w',index)



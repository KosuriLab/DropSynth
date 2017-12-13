from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Restriction import *

##################################
#FUNCTIONS:
def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records

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

def checkConstruct(construct,lengthmax,filenum,ampprimf,ampprimr,barcode,assemblyprimf,assemblyprimr,padding_between_btsaI_ampR,fwdre=[],revre=[]):
    
    oligos = construct[1:]
    count_first = 0
    count_mid = 0
    count_last = 0
    BspQI_length = 8
    BtsaI_length = 6
    amppri_length = len(ampprimf)
    asmpri_length = len(assemblyprimf)
    bc_length = len(barcode)
    btsai_from_end_pos = amppri_length + BtsaI_length - 1
    bspqi_first = amppri_length + BspQI_length + 1
    bspqi_second = amppri_length + BspQI_length + bc_length -2
    asmpriFsearch_min_pos = amppri_length + 2*BspQI_length + bc_length +BtsaI_length-1
    asmpriRsearch_pos = lengthmax - amppri_length - BtsaI_length - asmpri_length
    
    KpnI_min_from_end = lengthmax - amppri_length - BtsaI_length - asmpri_length
    KpnI_min_from_start = len(ampprimf)+ 2*BspQI_length + bc_length + BtsaI_length
    
    #check that lengths all pass
    for oligoseqrec in oligos:
        if len(oligoseqrec.seq)>lengthmax:
            print('Oligo length out of range:')
            print(oligoseqrec.id)
            print(str(oligoseqrec.seq) + '\t' + str(len(oligoseqrec.seq)))
    
    #check that amplification primers are correct
    for oligoseqrec in oligos:
        if ampprimf != str(oligoseqrec.seq[0:len(ampprimf)]) or ampprimr != str(oligoseqrec.seq[-1*len(ampprimr):].reverse_complement()):
            print('Amp primers not found:')
            print(oligoseqrec.id)
            print(str(oligoseqrec.seq))
    
    #check that barcode and BspQI sites are correct
    for oligoseqrec in oligos:
        #first site is 15(amp)+8(BspQI)+1=24
        #second site is 15(amp)+8(BspQI)+12(barcode)-2
        
        if BspQI.search(oligoseqrec.seq) != [bspqi_first, bspqi_second]:
            print("BspQI SITES ARE WRONG")
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(BspQI.search(oligoseqrec.seq))
        #barcode site is 15amp + 8(BspQI)
        barcode_pos = amppri_length + BspQI_length
        if str(oligoseqrec.seq).find(barcode) != barcode_pos:
            print("BARCODE NOT FOUND")
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(barcode)
    
    #check that BtsI sites are correct
    for oligoseqrec in oligos:
    
        btssearch = BtsI.search(oligoseqrec.seq)
        kpnsearch = KpnI.search(oligoseqrec.seq)
        ndesearch = NdeI.search(oligoseqrec.seq)
        asmpriFsearch = str(oligoseqrec.seq).find(assemblyprimf)
        assemblyprimr_seq = Seq(assemblyprimr, generic_dna)
        asmpriRsearch = str(oligoseqrec.seq).find(str(assemblyprimr_seq.reverse_complement()))
        #pos end btsaI from end = amp length + BtsaI_length -1
        if len(btssearch) !=2 or len(oligoseqrec)-btssearch[1] != btsai_from_end_pos: # btssearch[0] != 43:#
            print("End BtsI site bad")
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(barcode)
            print(str(btssearch))
        
        first_oligo = 0
        if asmpriFsearch != -1:
            #this should be the first oligo in assembly
            first_oligo = 1
            count_first += 1
            #min asm search = amplength + 2*bspqilength + bclength+btsai-1
            if asmpriFsearch < asmpriFsearch_min_pos or asmpriFsearch > 100:
                print("Frist oligo Assembly FWD primer wrong")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if len(ndesearch) != 1 or ndesearch[0] != asmpriFsearch + asmpri_length + 3:
                print("First oligo NdeI site bad")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
                print(ndesearch)
            
            if len(kpnsearch) != 0:
                print("First oligo has a KpnI site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
                print(kpnsearch)
            
            if btssearch[0] != asmpriFsearch + 3:
                print("First oligo BtsI site bad")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
        
        last_oligo = 0
        
        if asmpriRsearch != -1: #this should be the last oligo in assembly
            
            last_oligo = 1
            count_last += 1
            
            if first_oligo == 1:
                print("This oligo contains both FWD and REV assembly primers")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if asmpriRsearch != asmpriRsearch_pos:#200mer: 159 #230mer:189 #calc: total length - 15ampPriLength - 6BtsaI -20asmPriLength
                print("Last oligo assembly primer REV wrong")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if len(ndesearch) != 0:
                print("Last oligo has an NdeI site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if len(kpnsearch) != 1 or kpnsearch[0] > KpnI_min_from_end or kpnsearch[0] < KpnI_min_from_start:
                print("Last oligo KpnI site bad")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
        elif first_oligo == 0: # this is middle oligo
            
            count_mid += 1
            if len(ndesearch) != 0:
                print("Middle oligo has an NdeI site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if len(kpnsearch) != 0:
                print("Middle oligo has a KpnI site. Lib_num:"+str(filenum))
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
                
            if asmpriFsearch != -1:
                print("Middle oligo has an Assembly primer FWD site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            	
            if asmpriRsearch != -1:
                print("Middle oligo has an Assembly primer REV site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
                
    return count_first, count_mid, count_last
    
#####################################
########## OPTIONS ##################
#####################################

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


#asmF skpp20 from 00_primer_screen.py output #skpp504F
#Primers for alternate codon versions are offset by len(num_oligos) = 15
assemblyprimf = ['ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG',
                'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG', 'ATCGGGGATGGTAACTAACG']
#enter the reverse primers #skpp504R-rc
#Primers for alternate codon versions are offset by len(num_oligos) = 15
assemblyprimr = ['ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT',
                'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT']

#ampF from 00_primer_screen.py output: skpp15 5##F for amplification primers
ampprimersf =  ['CGCAGGGTCCAGAGT', 'GGGTTCGAGCGGGAG', 'ACTCGACGGCCTCTG',
                'GCGGCACCACAAACT', 'TCCACCGTCGGCAAG', 'GGCGCGCTCTAACAC',
                'AACGCCCAGCCTGTC', 'AGGCACGCTCAACCT', 'CATTGCCGTGCGTGA',
                'CGCCGAGCCGTATGA', 'AGCCCACTTGCCCTC', 'GAGGGCTCCGTTCGT',
                'CCCTCCCACGGACTT', 'CGTCCGCACAAACCC', 'GAGTCTGAGCGGCGT',
                'GCCGGTCCCAACTCT', 'AGTCCAGCGGCTCAC', 'TCTGAGACGGCGAGG',
                'CGGGCGCCTCTTGTT', 'AGGCGCTCATGTGGA', 'CGTGCAATGTGGCGT',
                'GAGAGCCGGCCTGTG', 'GGGCACGCGGTAAGT', 'GCTCGGCCGTAGTGT',
                'ACCTCATGTGGCCGA', 'ACTGATGCGCGGTCT', 'TCCGCGTTCTTGGCT',
                'CAGCACATCCCGCCC', 'GGCACCGTCCTGTCT', 'CCTAACTGCGGGCGT',
                'ACTAGCCCGCGTTCC', 'GGCCTGCGCGTATCT', 'GCGACCCTCCACTGA', 'CGCAGGTACGGGTCT']
                
#ampR (not RC) from 00_primer_screen.py output:
ampprimersr =  ['GTTCGCGCGAAGGAA', 'TAGCGCGCAGAGAGG', 'ACACGCGCGTTGAAG',
                'CGTGGCCTCTGTCCT', 'GGCCGCACCCAGTAG', 'CTCCCTCTCGCAGCA',
                'CCGCGTTGCTGAGTG', 'CCTAGGTCGCACGCA', 'GAGGGTTCCCGCTGA',
                'GCGCATTGGAGGCTG', 'CCAAGCCGGGTTCCA', 'CGGCCAGGTCAGGTC',
                'GGGTCCCTCGTCTCC', 'CCGCATCGTTGACCC', 'GCCTAGCTCGCCTGA',
                'CAGCCATGTCTCGCC', 'CCGCCTTCTAGCCCA', 'AGGACGCCCGTAGTG',
                'AGCGCGATTCAGCCA', 'ACTCAGCAGCGGGAC', 'CGCTGGACTCGTGGT',
                'CACGCAGCCAAACCC', 'TGTGCCGCCAAGACC', 'CGAGTTGTGGCACGG',
                'CCAGTGACGCAGGGA', 'ATGAAGGCGGCAGGT', 'GGGACGTTCGGACCA',
                'CCCTGGTCGCGTCTG', 'CCATGCCCTCCGACT', 'ACGGCGGCCCTAATG',
                'GCCGACAATTCCGCC', 'AGCGTCGCCAAACCC', 'CGTGATCCCGCCAAG', 'CGCGTGGACTTGCTC']

barcodes = obtainFastaSequences('barcodes/filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction.fasta')

oligos_per_construct = [4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5]
constructs_per_lib = [384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384]

padding_between_btsaI_ampR = True

oligo_length = 230

#####################################
######### / OPTIONS #################
#####################################

for fileindex in range(len(inputfiles)):
    constructs = getOligos(inputfiles[fileindex].replace('.oligos','-finaloligos.fasta'))
    first_oligo_counter = 0
    middle_oligo_counter = 0
    last_oligo_counter = 0
    for index in range(len(constructs)):
        add_first, add_mid, add_last = checkConstruct(constructs[index],oligo_length,fileindex+1,ampprimersf[fileindex],ampprimersr[fileindex],str(barcodes[index].seq),assemblyprimf[fileindex],assemblyprimr[fileindex],padding_between_btsaI_ampR)
        first_oligo_counter += add_first
        middle_oligo_counter += add_mid
        last_oligo_counter += add_last
    print("Lib: " + inputfiles[fileindex])
    if first_oligo_counter != constructs_per_lib[fileindex] or middle_oligo_counter != ((oligos_per_construct[fileindex]-2)*constructs_per_lib[fileindex]) or last_oligo_counter != constructs_per_lib[fileindex]:
        print(str(first_oligo_counter) + " start oligos. Expect to have " + str(constructs_per_lib[fileindex]))
        print(str(middle_oligo_counter) + " middle oligos or " + str(middle_oligo_counter/constructs_per_lib[fileindex]) + " per construct. Expect to have "+ str((oligos_per_construct[fileindex]-2)*constructs_per_lib[fileindex])+ " total middle oligos.")
        print(str(last_oligo_counter) + " end oligos. Expect to have "+ str(constructs_per_lib[fileindex]))
    else:
        print("Everything good.")

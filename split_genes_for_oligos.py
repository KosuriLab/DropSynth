#Converts protein sequences into split oligos
#allows multiple attempts for oligo splitting
#python 3.5
#updated 10.20.2016 Calin Plesa

#to profile
#python -m cProfile -o profile_08.txt split_genes_for_oligos.py

from Bio import SeqIO
import random
random.seed(15643243242342759)
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
from Bio.Restriction import *
import Primerselectiontools_py3
import csv
from Bio import pairwise2
import time

##################################
#FUNCTIONS:
def screenRestrictionSites(seqrec):
    recuts_btsI = Restriction.BtsI.search(seqrec.seq)
    recuts_bspQI = Restriction.BspQI.search(seqrec.seq)
    recuts_EcoRI = Restriction.EcoRI.search(seqrec.seq)
    if len(recuts_btsI)>0 or len(recuts_bspQI)>0 or len(recuts_EcoRI)>0: #calin changed
        print('%s\t%s\t%s' % (seqrec.id, recuts_btsI, recuts_bspQI,recuts_EcoRI))
        
#generate random DNA and limit max homopolymer length
#http://stackoverflow.com/questions/21205836/generating-random-sequences-of-dna
def randomDNA(length):
    gsflag = 0
    while gsflag == 0:
        temp_str = ''.join(random.choice('CGTA') for _ in range(length))
        #make sure no homopolymers greater than 5 in length
        if ((temp_str.find('AAAAA') == -1) and (temp_str.find('TTTTT') == -1) and
            (temp_str.find('GGGGG') == -1) and (temp_str.find('CCCCC') == -1)):
            gsflag = 1
    return temp_str

def readCodonUsage(filename="ecoli-codon-usage.txt"):
    handle = open(filename)
    strings = []
    for line in handle:
        splitline = line.split() # breaks by whitespace
        for i in [0,5,10,15]:
            strings.append(splitline[i:i+5])
    
    codons = {}
    for string in strings:
        if string[1] in codons: #string[1] is triplet
            dnacodon = string[0].replace('U','T')
            codons[string[1]].append([dnacodon, float(string[3])]) #string[3] is probability
        else:
            dnacodon = string[0].replace('U','T')
            codons[string[1]] = [[dnacodon, float(string[3])]]
    
    #make a fractional accounting
    for aa in codons:
        runningaacount = 0
        for i in codons[aa]:
            runningaacount += i[1] #add all probabilities
        for i in codons[aa]:
            fraction = i[1]/runningaacount #normalize
            i.append(fraction)
            #if fraction<0.15:
            #    print i
    
    #drop codons less than 15%
    for aa in codons:
        newlist = []
        for i in codons[aa]:
            if i[2]>0.15:
                newlist.append(i)
        codons[aa] = newlist
    
    #do running total accounting
    for aa in codons:
        runningaacount = 0
        runningtotal = 0
        for i in codons[aa]:
            runningaacount += i[1]
        for i in codons[aa]:
            fraction = i[1]/runningaacount
            i[2] = fraction
            runningtotal += fraction
            i.append(runningtotal)

    #take last codon from each and make 1.0
    for aa in codons:
        codons[aa][-1][-1] = 1.0
    
    #make reverse codon dictionary
    revcodons = {}
    for aa in codons:
        for i in codons[aa]:
            revcodons[i[0]] = aa
            
    
    #for aa in codons:
    #    print aa + '\t\t' + str(codons[aa])
    return codons, revcodons

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print seqrecord.id
        #print len(seqrecord)
    #print len(records)
    return records

#do a basic optimization; not checking restriction sites or other sites to avoid
def codonOptimize_two(codons, record, gentwo,furthest_codon):
    sequence = ""
    if gentwo:
        sequence2 = ""
    #print(codons)
    #example of a codons record:
    #'D': [['GAT', 32.2, 0.6276803118908382, 0.6276803118908382], ['GAC', 19.1, 0.3723196881091618, 1.0]]
    #for each amino acid:
    for letter in record.seq:
        #print(str(letter) + '\t' + str(codons[letter][0][0]))
        #pick random number
        rn = random.random()
        for i in codons[letter]:
            if rn < i[-1]:
                #print str(letter) + '\t' + i[0]
                sequence = sequence + i[0]
                if gentwo:
                    sequence2 = sequence2 + furthest_codon[i[0]]
                break
    newseq = Seq(sequence,generic_dna)
    newrecord = SeqRecord(newseq, id=record.id)
    if gentwo:
        newseq2 = Seq(sequence2,generic_dna)
        newrecord2 = SeqRecord(newseq2, id=record.id)
        
        alignment = pairwise2.align.globalxx(newrecord.seq, newrecord2.seq, one_alignment_only=1, score_only=1)
        pct_ident = alignment/len(newrecord.seq)
        
        return newrecord, newrecord2, pct_ident
    else:
        return newrecord

def removeSite(seqrecord, site, codons, revcodons,firstinstance=True):
    #find site
    seqtoremove = site
    if firstinstance:
        index = seqrecord.seq.find(site)
    else:
        index = seqrecord.seq.rfind(site)
    if index == -1:
        if firstinstance:
            index = seqrecord.seq.find(site.reverse_complement())
            seqtoremove = site.reverse_complement()
        else:
            index = seqrecord.seq.rfind(site.reverse_complement())
            seqtoremove = site.reverse_complement()
            
    
    recognitionlength = len(seqtoremove)
    if index == -1: 
        return False
    if index % 3 == 0:
        #first base is first base in codon
        #calculate how many codons can be changed
        
        codonstochange = (((recognitionlength-1) // 3)+1) #py3
        
        for i in range(codonstochange):
            codonstart = i*3 + index
            codontochange = seqrecord.seq[codonstart:codonstart+3]
            for possiblecodon in codons[revcodons[str(codontochange)]]:
                if possiblecodon[0] != str(codontochange):
                    tempseq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart+3:]
                    seqrecord.seq = tempseq
                    if str(seqtoremove) != seqrecord.seq[index:index+recognitionlength]:
                        return True
        print("Could not change restriction site")
        return False
            #find alternatives:
    elif index % 3 == 1:
        #first base is second base in codon
        #calculate how many codons can be changed
        codonstochange = (((recognitionlength) // 3)+1) #py3
        for i in range(codonstochange):
            codonstart = i*3 + index - 1
            codontochange = seqrecord.seq[codonstart:codonstart+3]
            for possiblecodon in codons[revcodons[str(codontochange)]]:
                if possiblecodon[0] != str(codontochange):
                    tempseq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart+3:]
                    seqrecord.seq = tempseq
                    if str(seqtoremove) != seqrecord.seq[index:index+recognitionlength]:
                        return True
        print("Could not change restriction site")
        return False
    elif index % 3 == 2:
        # first base is third base in codon
        #calculate how many codons can be changed
        codonstochange = (((recognitionlength+1) // 3)+1) #py3
        for i in range(codonstochange):
            codonstart = i*3 + index - 2
            codontochange = seqrecord.seq[codonstart:codonstart+3]
            for possiblecodon in codons[revcodons[str(codontochange)]]:
                if possiblecodon[0] != str(codontochange):
                    tempseq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart+3:]
                    seqrecord.seq = tempseq
                    if str(seqtoremove) != seqrecord.seq[index:index+recognitionlength]:
                        return True
        print("Could not change restriction site")
        return False
        

def changeRestrictionSites(seqrecord, codons, revcodons):
    tempseq = Seq(str(seqrecord.seq),generic_dna)
    newseqrecord = SeqRecord(tempseq,id=seqrecord.id,description=seqrecord.description)
    rb = RestrictionBatch([NdeI, EcoRI, KpnI, BtsI, BspQI])
    #print seqrecord.id
    #print seqrecord.seq
    reanalysis = rb.search(newseqrecord.seq)
    #print reanalysis
    for key in reanalysis:
        for i in reanalysis[key]:
            seqkey = Seq(key.site, generic_dna)
            removeSite(newseqrecord, seqkey, codons, revcodons)
    return newseqrecord

def changeRestrictionSites2(seqrecordin, codons, revcodons):
    tempseq = Seq(str(seqrecordin.seq),generic_dna)
    seqrecord = SeqRecord(tempseq,id=seqrecordin.id,description=seqrecordin.description)
    rb = RestrictionBatch([BtsI, BspQI,EcoRI])
    reanalysis = rb.search(seqrecord.seq)
    #print reanalysis
    #if there is BtsI or BspQI or EcoRI, juggle codons to remove site:
    for key in reanalysis:
        for i in reanalysis[key]:
            seqkey = Seq(key.site, generic_dna)
            removeSite(seqrecord, seqkey, codons, revcodons)
    
    rb = RestrictionBatch([NdeI, KpnI])
    reanalysis = rb.search(seqrecord.seq)
    #print reanalysis
    #if there are more than one NdeI or KpnI sites, juggle codons to remove the second site:
    for key in reanalysis:
        if len(reanalysis[key])>1:
            if key==NdeI:
                seqkey = Seq(key.site, generic_dna)
                removeSite(seqrecord, seqkey, codons, revcodons,False)
            else:
                seqkey = Seq(key.site, generic_dna)
                removeSite(seqrecord, seqkey, codons, revcodons)
    return seqrecord

def checkRestrictionSites(seqrecords):
    for seqrecord in seqrecords:
        rb = RestrictionBatch([NdeI, EcoRI, KpnI, BtsI, BspQI])
        print(seqrecord.id)
        print(seqrecord.seq)
        print(rb.search(seqrecord.seq))
        print(len(seqrecord.seq))
        

def checkCodingRegions(aarecord, seqrecord):
    numberOfMatchedProteins = 0
    newproteinseq = seqrecord.seq.translate()
    oldproteinseq = aarecord.seq
    if str(newproteinseq) == str(oldproteinseq):
        numberOfMatchedProteins += 1
    else:
        print(seqrecord.id+" bad codon assignment, translation mismatch")
    #print str(numberOfMatchedProteins) + " proteins with correct aa sequence redesigned"
    return numberOfMatchedProteins

def addBufferREsites(seqrecordin,add_stop_codon):
    forward_buffer = "CAT" # before start codon ATG = NdeI
    if add_stop_codon:
        reverse_buffer = "TAAGGTACC" # KpnI
    else:
        reverse_buffer = "GGTACC" # KpnI
    #print(str(seqrecordin.seq[0:3]).upper())
    if not ("ATG" == str(seqrecordin.seq[0:3]).upper()):
        print("Error! Sequence does not start with ATG. "+str(seqrecordin.id))
    newseq = Seq(forward_buffer + str(seqrecordin.seq) + reverse_buffer,generic_dna)
    record = SeqRecord(newseq,id=seqrecordin.id,description=seqrecordin.description)
    
    return record

def addAssemblyPrimers(seqrecordin, forward_buffer, reverse_buffer, padding_var, padding_length):
    seqlength = len(seqrecordin.seq)
    pad_added_length = 0
    if (padding_var == False) or ((padding_var == True) and (seqlength >= padding_length)):
        tempseq = forward_buffer.seq + seqrecordin.seq + reverse_buffer.seq.reverse_complement()
        #tempseq = Seq(str(forward_buffer + seqrecordin.seq + reverse_buffer),generic_dna)
    else:
        pad_random_seq = randomDNA(padding_length-len(seqrecordin.seq))
        tempseq = forward_buffer.seq + seqrecordin.seq + Seq(pad_random_seq,generic_dna) + reverse_buffer.seq.reverse_complement()
        #tempseq = Seq(str(forward_buffer + seqrecordin.seq + pad_random_seq + reverse_buffer),generic_dna)
        pad_added_length = len(pad_random_seq)
    seqrecord = SeqRecord(tempseq,id=seqrecordin.id,description=seqrecordin.description)
    return seqrecord, pad_added_length


def checkRestrictionSites2(seqrecord):
    
    rb = RestrictionBatch([NdeI, KpnI, EcoRI, BtsI, BspQI])
    seqsearch = rb.search(seqrecord.seq)
    fail_flag = 0
    if len(seqsearch[NdeI])!=1 or len(seqsearch[EcoRI])!=0 or len(seqsearch[KpnI])!=1 or len(seqsearch[BtsI])!=0 or len(seqsearch[BspQI])!=0:
        #print(seqrecord.id+" fail in checkRestrictionSites2")
        #print(seqrecord.seq)
        #print(rb.search(seqrecord.seq))
        #print(len(seqrecord.seq))
        fail_flag = 1
    
    return fail_flag

def checkOligos(oligosin):
    maxolgiolength = 0
    length_in_case_of_error = 1000
    #find length of longest oligo after split:
    for oligo in oligosin:
        if len(oligo) > maxolgiolength:
            maxolgiolength = len(oligo)
            
    checkOligos_msg = ""
    #check that adding BtsaI sites doesn't introduce new restriction sites
    for i in range(1,len(oligosin)):
        newoligo = Seq("GCAGTG",generic_dna) + oligosin[i] + Seq("CACTGC",generic_dna)
        rb = RestrictionBatch([BtsI, BspQI, NdeI, KpnI,EcoRI])
        seqsearch = rb.search(newoligo)
        if i==0:
            #first oligo
            if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 1 or len(seqsearch[KpnI])!=0 or len(seqsearch[EcoRI])!=0:
                checkOligos_msg = "Bad number of restriction sites in first oligo"
                #print(checkOligos_msg)
                #print(newoligo)
                #print(seqsearch)
                maxolgiolength = length_in_case_of_error
                
        elif i==(len(oligosin)-1):
            #last oligo
            if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=1 or len(seqsearch[EcoRI])!=0:
                if checkOligos_msg == "":
                    checkOligos_msg = "Bad number of restriction sites in last oligo"
                else:
                    checkOligos_msg = checkOligos_msg + ", last oligo"
                #print(checkOligos_msg)
                #print(newoligo)
                #print(seqsearch)
                maxolgiolength = length_in_case_of_error
        else:
            #middle oligos
            if len(seqsearch[BtsI])!=2 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0 or len(seqsearch[EcoRI])!=0:
                if checkOligos_msg == "":
                    checkOligos_msg = "Bad number of restriction sites in middle oligo"
                else:
                    checkOligos_msg = checkOligos_msg + ", middle oligo"
                #print(checkOligos_msg)
                #print(newoligo)
                #print(seqsearch)
                maxolgiolength = length_in_case_of_error
    return maxolgiolength, checkOligos_msg

######################################################################
if __name__ == '__main__':
    start_time = time.time()
    codons, revcodons = readCodonUsage()
    
    #####################################
    ########## OPTIONS ##################
    #####################################

    #table for each codon which other codon encodes same aa but is most distant
    furthest_codon = dict([('GGT', 'GGC'), ('GGC', 'GGT'), ('TCG', 'AGC'), ('AGT', 'TCG'), ('AGC', 'AGT'), 
        ('ATG', 'ATG'), ('GCT', 'GCC'), ('GCC', 'GCT'), ('GCA', 'GCG'), ('GCG', 'GCA'), ('CGT', 'CGC'), ('CGC', 'CGT'),
        ('TGT', 'TGC'), ('TGC', 'TGT'), ('CTG', 'CTG'), ('CAA', 'CAG'), ('CAG', 'CAA'), ('CCT', 'CCG'), ('CCA', 'CCT'),
        ('CCG', 'CCA'), ('TAT', 'TAC'), ('TAC', 'TAT'), ('AAT', 'AAC'), ('AAC', 'AAT'), ('TTT', 'TTC'), ('TTC', 'TTT'),
        ('GTT', 'GTC'), ('GTC', 'GTT'), ('GTA', 'GTG'), ('GTG', 'GTA'), ('CAT', 'CAC'), ('CAC', 'CAT'), ('GAA', 'GAG'),
        ('GAG', 'GAA'), ('ACT', 'ACG'), ('ACC', 'ACT'), ('ACG', 'ACC'), ('ATT', 'ATC'), ('ATC', 'ATT'), ('TGG', 'TGG'),
        ('TAA', 'TGA'), ('TGA', 'TAA'), ('AAA', 'AAG'), ('AAG', 'AAA'), ('GAT', 'GAC'), ('GAC', 'GAT')])
    
    #these are the output file base names for each library, one per 384 lib:
    input_files = ["db/DHFR_Lib01_4oligo.fasta","db/DHFR_Lib02_4oligo.fasta",
                    "db/DHFR_Lib03_4oligo.fasta","db/DHFR_Lib04_4oligo.fasta",
                    "db/DHFR_Lib05_4oligo.fasta","db/DHFR_Lib06_4oligo.fasta",
                    "db/DHFR_Lib07_4oligo.fasta","db/DHFR_Lib08_4oligo.fasta",
                    "db/DHFR_Lib09_4oligo.fasta","db/DHFR_Lib10_4oligo.fasta",
                    "db/DHFR_Lib11_4oligo.fasta","db/DHFR_Lib12_4oligo.fasta",
                    "db/DHFR_Lib13_4oligo.fasta","db/DHFR_Lib14_5oligo.fasta",
                    "db/DHFR_Lib15_5oligo.fasta"]
    
    #the source files with all sequences split by number of oligos required
    input_files_by_numoligo = {4:"db/IP_protein_after_padding_to_384_4olgio.fasta",5:"db/IP_protein_after_padding_to_384_5olgio.fasta"}
    
    output_folder = 'db_oligo'
    
    #types of num oligos
    type_num_oligos = [4,5]
    
    #number of libraries for each numOligo type:
    numlibs_for_each_numOligos = {4: 13, 5: 2}
    
    # number of oligos to use to split genes in each lib:
    num_oligos = [4,4,4,4,4,4,4,4,4,4,4,4,4,5,5]
    
    # max space on oligo between BtsaI cut sites (160 for 230mer):
    max_payload_length = 160
    
    #maximum tries to split each gene (default 14):
    max_num_attempts = 14
    
    #generate two different codon optimizations for each protein?
    gen_two_codon = True
    
    #add a stop codon before KpnI?
    add_stop_codon = True
    
    #use padding (add random sequence between final BtsI site and the reverse assembly primer)
    #this makes gel isolation easier since the distribution of lengths is reduced
    padding_var = True
    
    #padding lengths for given number of oligos
    #sequences with lengths below this will get padding
    padding_length = {4: 510, 5: 640}
    
    #block any homopolymer repeats of this length or longer:
    max_homopolymer_repeat_length = 8
    
    #this is a list of all seq IDs which don't work
    IDs_to_remove = "db/IDs_to_remove.csv"
    
    #asmF primer file (only primer 504 used)
    #DO NOT USE skpp20-511 this is used in the anchor oligo!
    #Primers for alternate codon versions are offset by len(num_oligos) = 15
    assemblyprimf_file = 'asmprimers-skpp20/asmF_selected_mod504.fasta'
    assemblyprimf = obtainFastaSequences(assemblyprimf_file)

    #enter the reverse primers (only primer 504 used)
    #Primers for alternate codon versions are offset by len(num_oligos) = 15
    assemblyprimr_file = 'asmprimers-skpp20/asmR_selected_mod504.fasta'
    assemblyprimr = obtainFastaSequences(assemblyprimr_file)
    
    seeds_for_libs = {1: 15643243242342759, 2: 888808080808, 3: 181818181818,
                        4:84939686, 5:4958686,6:94857540,
                        7:984850947, 8:944798879,9:6784394785643782,
                        10:7832956375392,11:92038586762355,12:333554665565,
                        13:83838383838,14:1223723723723,15:2723723732782387}
    
    constructs_per_lib = [384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384]
    
    #####################################
    ######### / OPTIONS #################
    #####################################

    badIDsarray = []
    with open(IDs_to_remove, 'r') as f:
            reader = csv.reader(f)
            lib_ident_list = list(reader)
    
    if len(lib_ident_list) > 0:
        badIDsarray = lib_ident_list[0]
    else:
        badIDsarray = []

    lib_counter = 1
    
    for current_type_num_oligo in type_num_oligos:
        aarec_counter = 0
        #obtain the protein records for this numOligos
        aarecords = obtainFastaSequences(input_files_by_numoligo[current_type_num_oligo])
        print("Starting with "+str(current_type_num_oligo)+" oligo libraries ("+str(len(aarecords))+" records read):")
        for libfile in range(numlibs_for_each_numOligos[current_type_num_oligo]):
            random.seed(seeds_for_libs[lib_counter])
            
            print("###########################################")
            print("Lib "+str(lib_counter)+" ("+str(constructs_per_lib[lib_counter-1])+" records):")
            
            final_records_gene_level = []
            final_records_gene_level_no_primer_or_RE = []
            final_records_gene_level_no_primer = []
            final_records_protein_level = []
            updatedseqrecord = []
            seqrecord = []
            
            pct_ident_array = []
            padding_length_array = []
            
            if gen_two_codon:
                final_records_gene_level2 = []
                final_records_gene_level_no_primer_or_RE2 = []
                final_records_gene_level_no_primer2 = []
                final_records_protein_level2 = []
            
            #write the split oligos here:
            fileout = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','oligos'),'w')
            if gen_two_codon:
                fileout2 = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','codon2.oligos'),'w')
            
            rec_count = 1
            total_rec_count = len(aarecords)
            
            #for each protein record:
            #print(str(constructs_per_lib[lib_counter-1]))
            #print(str(lib_counter-1))
            for ii in range(constructs_per_lib[lib_counter-1]):
                #print(str(ii))
                protein_seq_good_flag = 0
                while protein_seq_good_flag == 0:
                    
                    #get the next protein sequence for this numOligos
                    caarec = aarecords[aarec_counter]
                    num_tries = 1
                    while (num_tries < max_num_attempts):
                        #print("Starting "+caarec.id+" (construct # "+str(rec_count)+" of "+str(total_rec_count)+", Lib "+str(lib_counter)+") try # "+str(num_tries))
                        #generate random weighted codon usage nt record:
                        if gen_two_codon:
                            seqrecord,seqrecord2,pct_ident = codonOptimize_two(codons, caarec, gen_two_codon,furthest_codon)
                        else:
                            seqrecord = codonOptimize_two(codons, caarec, gen_two_codon,furthest_codon)
                        
                        #remove NdeI, EcoRI, KpnI, BtsI, BspQI:
                        updatedseqrecord = changeRestrictionSites(seqrecord, codons, revcodons)
                        if gen_two_codon:
                            updatedseqrecord2 = changeRestrictionSites(seqrecord2, codons, revcodons)
                        
                        #make sure the translation is still the same:
                        numMatch = checkCodingRegions(caarec, updatedseqrecord)
                        if gen_two_codon and numMatch == 1:
                            numMatch = checkCodingRegions(caarec, updatedseqrecord2)
                        
                        #make sure no long homopolymers, 
                        temp_str = str(updatedseqrecord.seq)
                        if gen_two_codon:
                            #concatenate both records together
                            temp_str = str(updatedseqrecord.seq)+str(updatedseqrecord2.seq)
                        gsflag = 0
                        
                        if ((temp_str.find('A' * max_homopolymer_repeat_length) == -1) and (temp_str.find('T' * max_homopolymer_repeat_length) == -1) and
                            (temp_str.find('G' * max_homopolymer_repeat_length) == -1) and (temp_str.find('C' * max_homopolymer_repeat_length) == -1)):
                            gsflag = 1

                        #if everything is good:
                        if numMatch == 1 and gsflag == 1:
                            
                            #add NdeI at start and Stop codon + KpnI at end:
                            bufferedupdatedseqrecord = addBufferREsites(updatedseqrecord,add_stop_codon)
                            if gen_two_codon:
                                bufferedupdatedseqrecord2 = addBufferREsites(updatedseqrecord2,add_stop_codon)
                            RE_fail_flag1 = 0
                            
                            #change codons if BtsaI, BspQI, or EcoRI are present. Change codons if >1 NdeI or KpnI site present.
                            finalupdatedseqrecord = changeRestrictionSites2(bufferedupdatedseqrecord,codons,revcodons)
                            if gen_two_codon:
                                finalupdatedseqrecord2 = changeRestrictionSites2(bufferedupdatedseqrecord2,codons,revcodons)
                            
                            if RE_fail_flag1 == 0:
                                pad_length_added = 0
                                #add the assembly primers and also add padding with random sequence if necessary
                                finalrecwithAssPrimers, pad_length_added = addAssemblyPrimers(finalupdatedseqrecord,assemblyprimf[lib_counter-1],assemblyprimr[lib_counter-1],padding_var,padding_length[num_oligos[lib_counter-1]])
                                if gen_two_codon:
                                    pad_length_added2 = 0
                                    finalrecwithAssPrimers2, pad_length_added2 = addAssemblyPrimers(finalupdatedseqrecord2,assemblyprimf[lib_counter-1+len(num_oligos)],assemblyprimr[lib_counter-1+len(num_oligos)],padding_var,padding_length[num_oligos[lib_counter-1]])
                                
                                RE_fail_flag2 = 0

                                # make sure 1 NdeI site, 1 KpnI site, and none of EcoRI, BtsaI, BspQI:
                                RE_fail_flag2 = checkRestrictionSites2(finalrecwithAssPrimers)
                                if gen_two_codon and RE_fail_flag2 == 0:
                                    RE_fail_flag2 = checkRestrictionSites2(finalrecwithAssPrimers2)
                                
                                #if no extra restriction sites found:
                                if RE_fail_flag2 == 0:
                                    screenRestrictionSites(finalrecwithAssPrimers)
                                    if gen_two_codon:
                                        screenRestrictionSites(finalrecwithAssPrimers2)
                                    
                                    #this is where genes are split into oligos:
                                    #input parameters are:
                                    #seq, oligosizemax, lengthleeway, positionleeway, avgoverlapsize, overlaptemps, deltaGThreshold, selfDimersThreshold, num_of_oligos
                                    oligos = Primerselectiontools_py3.optimizedSplit(finalrecwithAssPrimers.seq,max_payload_length,15,15,20,[58,62],-4,4,num_oligos[lib_counter-1])
                                    if gen_two_codon:
                                        oligos2 = Primerselectiontools_py3.optimizedSplit(finalrecwithAssPrimers2.seq,max_payload_length,15,15,20,[58,62],-4,4,num_oligos[lib_counter-1])
                                    
                                    checkOligos_msg = ""
                                    maxolgiolength,checkOligos_msg = checkOligos(oligos)
                                    
                                    #make sure number of oligos and max oligo size is good:
                                    if (len(oligos) != num_oligos[lib_counter-1]) or (maxolgiolength > max_payload_length):
                                        print(caarec.id+" failed oligo splitting on try # "+str(num_tries)+". "+checkOligos_msg)
                                    else:
                                        if gen_two_codon:
                                            checkOligos_msg = ""
                                            maxolgiolength,checkOligos_msg = checkOligos(oligos2)
                                            
                                            if (len(oligos2) != num_oligos[lib_counter-1]) or (maxolgiolength > max_payload_length):
                                                print(caarec.id+" failed oligo codon2 splitting on try # "+str(num_tries)+". "+checkOligos_msg)
                                            else:
                                                #everything worked
                                                protein_seq_good_flag = 1
                                                aarec_counter += 1
                                                
                                                #write split oligos for codon version 2:
                                                fileout2.write('>' + finalrecwithAssPrimers2.id + '\n')
                                                for oligo in oligos2:
                                                    fileout2.write(str(oligo) + '\n')
                                                
                                                #save all correct records:
                                                final_records_gene_level2.append(finalrecwithAssPrimers2)
                                                final_records_protein_level2.append(caarec)
                                                final_records_gene_level_no_primer_or_RE2.append(updatedseqrecord2)
                                                final_records_gene_level_no_primer2.append(finalupdatedseqrecord2)
                                                
                                                pct_ident_array.append(pct_ident)
                                                padding_length_array.append(pad_length_added)

                                                #write split oligos for codon version 1:
                                                fileout.write('>' + finalrecwithAssPrimers.id + '\n')
                                                for oligo in oligos:
                                                    fileout.write(str(oligo) + '\n')

                                                #save all correct records:
                                                final_records_gene_level.append(finalrecwithAssPrimers)
                                                final_records_protein_level.append(caarec)
                                                final_records_gene_level_no_primer_or_RE.append(updatedseqrecord)
                                                final_records_gene_level_no_primer.append(finalupdatedseqrecord)
                                                
                                                
                                                break
                                        
                                        else: #only one codon
                                            #everything worked
                                            protein_seq_good_flag = 1
                                            aarec_counter += 1

                                            #write split oligos:
                                            fileout.write('>' + finalrecwithAssPrimers.id + '\n')
                                            for oligo in oligos:
                                                fileout.write(str(oligo) + '\n')

                                            #save all correct records:
                                            final_records_gene_level.append(finalrecwithAssPrimers)
                                            final_records_protein_level.append(caarec)
                                            final_records_gene_level_no_primer_or_RE.append(updatedseqrecord)
                                            final_records_gene_level_no_primer.append(finalupdatedseqrecord)
                                            
                                            break
                                else:
                                    print(caarec.id+ " failed checkRestrictionSites2 after adding RE, AsmPri, and buffer seq on try # "+str(num_tries))
                        else:
                            print(caarec.id+ " failed codon generation, translation mismatch or long stretch of homopolymer")
                        num_tries += 1
                    if (num_tries > max_num_attempts-1):
                        protein_seq_good_flag = 0
                        aarec_counter += 1
                        print("*** "+caarec.id+ " reached max number of attempts. Adding seq to badIDsarray.")
                        #this seq cant be split, add it to list of files to remove
                        if not (caarec.id in badIDsarray):
                            badIDsarray.append(caarec.id)
                
                rec_count += 1
            
            print("Processed "+str(rec_count-1)+" records.")
            fileout.close()
            fileout2.close()
            
            #save file for each lib (codon 1) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and assembly primers):
            gene_nt_fileout = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','full_wRE_wPrim.genes'),'w')
            SeqIO.write(final_records_gene_level, gene_nt_fileout, "fasta")
            gene_nt_fileout.close()
            
            #save file for each lib (codon 1) with full gene seq before splitting (with no restriction sites (NdeI & KpnI) and no assembly primers):
            #these seqs have start codon but no stop codon
            gene_nt_no_primer_or_RE_fileout = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','full_nRE_nPrim.genes'),'w')
            SeqIO.write(final_records_gene_level_no_primer_or_RE, gene_nt_no_primer_or_RE_fileout, "fasta")
            gene_nt_no_primer_or_RE_fileout.close()
            
            #save file for each lib (codon 1) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and no assembly primers):
            #these seqs have start codon (part of NdeI) and stop codon (TAA) before KpnI
            gene_nt_no_primer_fileout = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','full_wRE_noPrim.genes'),'w')
            SeqIO.write(final_records_gene_level_no_primer, gene_nt_no_primer_fileout, "fasta")
            gene_nt_no_primer_fileout.close()
            
            #save file for each lib (codon 1 & 2) with full protein seq (a.a.) (these have starting Met)
            protein_nt_fileout = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','proteins'),'w')
            SeqIO.write(final_records_protein_level, protein_nt_fileout, "fasta")
            protein_nt_fileout.close()
            
            #save for codon 2 also:
            if gen_two_codon:
                #save file for each lib (codon 2) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and assembly primers):
                gene_nt_fileout2 = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','codon2.full_wRE_wPrim.genes'),'w')
                SeqIO.write(final_records_gene_level2, gene_nt_fileout2, "fasta")
                gene_nt_fileout2.close()
                
                #save file for each lib (codon 2) with full gene seq before splitting (with no restriction sites (NdeI & KpnI) and no assembly primers):
                #these seqs have start codon but no stop codon
                gene_nt_no_primer_or_RE_fileout2 = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','codon2.full_nRE_nPrim.genes'),'w')
                SeqIO.write(final_records_gene_level_no_primer_or_RE2, gene_nt_no_primer_or_RE_fileout2, "fasta")
                gene_nt_no_primer_or_RE_fileout2.close()
                
                #save file for each lib (codon 2) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and no assembly primers):
                #these seqs have start codon (part of NdeI) and stop codon (TAA) before KpnI
                gene_nt_no_primer_fileout2 = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','codon2.full_wRE_noPrim.genes'),'w')
                SeqIO.write(final_records_gene_level_no_primer2, gene_nt_no_primer_fileout2, "fasta")
                gene_nt_no_primer_fileout2.close()
                
                #save file for each lib (codon 1 & 2) with full protein seq (a.a.) (these have starting Met)
                #redundant
                protein_nt_fileout2 = open(output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','codon2.proteins'),'w')
                SeqIO.write(final_records_protein_level2, protein_nt_fileout2, "fasta")
                protein_nt_fileout2.close()
                
                #save file for each lib with nt identity between two codon versions
                proteinnt_ident_fileout = output_folder+'/'+input_files[lib_counter-1].split('/')[-1].replace('fasta','lib_ident.csv')
                resultFile = open(proteinnt_ident_fileout,'w')
                wr = csv.writer(resultFile, dialect='excel')
                wr.writerow(pct_ident_array)
            
            lib_counter += 1
    
    #save all of the IDs for sequences that failed completely:
    resultFile2 = open(IDs_to_remove,'w')
    wr = csv.writer(resultFile2, dialect='excel')
    wr.writerow(badIDsarray)
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
from Bio import SeqIO
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#screen primer pairs for restriction sites
#asmF: BtsaI + AsmF + NdeI
#asmR: KpnI + AsmR + BtsaI
#ampF: AmpF+ nt.BspQI
#ampR: btsaI + AmpR

#read in fasta files
def getFastaSeqs(filename):
    fastaseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastaseqs.append(seqrec)
    handle.close()
    return fastaseqs

#for the asm priemrs start at 500:
skip_first_on_skpp20 = 500
num_primers_to_screen = 200
num_primers_to_output = 200

#primer input files:
asmF_file = 'asmprimers-skpp20/forward_finalprimers.fasta'
asmR_file = 'asmprimers-skpp20/reverse_finalprimers.fasta'
ampF_file = 'ampprimers-skpp15/skpp15-forward.faa'
ampR_file = 'ampprimers-skpp15/skpp15-reverse.faa'

BtsaI_seq = Seq('GCAGTG',generic_dna)
BspQI_seq = Seq('GCTCTTCNNN',generic_dna)
KpnI_seq = Seq('GGTACC',generic_dna)
NdeI_seq = Seq('CATATG',generic_dna)

asmFseqs = getFastaSeqs(asmF_file)
asmRseqs = getFastaSeqs(asmR_file)
ampFseqs = getFastaSeqs(ampF_file)
ampRseqs = getFastaSeqs(ampR_file)

asmFseqs_pass = []
asmRseqs_pass = []
ampFseqs_pass = []
ampRseqs_pass = []

asmFseqs_pass_seq = []
asmRseqs_pass_seq = []
asmRseqs_pass_seqRC = []
ampFseqs_pass_seq = []
ampRseqs_pass_seq = []
ampRseqs_pass_seqRC = []

ampAllseqs_pass_seq_rec = []
asmFseqs_pass_seq_rec = []
asmRseqs_pass_seq_rec = []

asmFseqs = asmFseqs[skip_first_on_skpp20:skip_first_on_skpp20+num_primers_to_screen]
asmRseqs = asmRseqs[skip_first_on_skpp20:skip_first_on_skpp20+num_primers_to_screen]
ampFseqs = ampFseqs[0:num_primers_to_screen]
ampRseqs = ampRseqs[0:num_primers_to_screen]

#which enzymes to screen for:
rb = RestrictionBatch([BtsI, BspQI, NdeI, KpnI, EcoRI])

# print(BtsI.site)
# print(BspQI.site)
# print(NdeI.site)
# print(KpnI.site)
# print(EcoRI.site)

#asmF: BtsaI + AsmF + NdeI
counter = 0
#loop over AsmF
for primer in asmFseqs:
    #build the sequence
    current_seq = BtsaI_seq + primer.seq + NdeI_seq
    #find illegal sites
    seqsearch = rb.search(current_seq)
    if len(seqsearch[BtsI])!=1 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 1 or len(seqsearch[KpnI])!=0 or len(seqsearch[EcoRI])!=0:
        print(primer.id+ ". asmFseqs: Bad number of restriction sites:")
        #print(newoligo)
        print(seqsearch)
    else:
        #no sites found
        asmFseqs_pass.append(counter)
    counter += 1

#asmR: Stop + KpnI + AsmR + BtsaI
counter = 0
#loop over AsmR
for primer in asmRseqs:
    #build the sequence
    current_seq = Seq('TAA',generic_dna) + KpnI_seq + primer.seq.reverse_complement() + BtsaI_seq.reverse_complement()
    #find illegal sites
    seqsearch = rb.search(current_seq)
    if len(seqsearch[BtsI])!=1 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=1 or len(seqsearch[EcoRI])!=0:
        print(primer.id+ ". asmR: Bad number of restriction sites:")
        #print(newoligo)
        print(seqsearch)
    else:
        #no sites found
        asmRseqs_pass.append(counter)
    counter += 1

#ampF: AmpF+ nt.BspQI
counter = 0
#loop over AmpF
for primer in ampFseqs:
    #build the sequence
    current_seq =  primer.seq + BspQI_seq
    #find illegal sites
    seqsearch = rb.search(current_seq)
    if len(seqsearch[BtsI])!=0 or len(seqsearch[BspQI])!=1 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0 or len(seqsearch[EcoRI])!=0:
        print(primer.id+ ". ampF: Bad number of restriction sites:")
        #print(newoligo)
        print(seqsearch)
    else:
        #no sites found
        ampFseqs_pass.append(counter)
    counter += 1
    
#ampR: btsaI + AmpR
counter = 0
#loop over AmpR
for primer in ampRseqs:
    #build the sequence
    current_seq = BtsaI_seq + primer.seq.reverse_complement()
    #find illegal sites
    seqsearch = rb.search(current_seq)
    if len(seqsearch[BtsI])!=1 or len(seqsearch[BspQI])!=0 or len(seqsearch[NdeI])!= 0 or len(seqsearch[KpnI])!=0 or len(seqsearch[EcoRI])!=0:
        print(primer.id+ ". ampR: Bad number of restriction sites:")
        #print(newoligo)
        print(seqsearch)
    else:
        #no sites found
        ampRseqs_pass.append(counter)
    counter += 1

pass_count = 0
#loop over each Asm primers pair
for i in range(0,num_primers_to_screen):
    #if both primers in pair passed
    if (i in asmFseqs_pass) and (i in asmRseqs_pass) and (pass_count < num_primers_to_output):
        asmFseqs_pass_seq.append(str(asmFseqs[i].seq))
        asmRseqs_pass_seq.append(str(asmRseqs[i].seq))
        asmRseqs_pass_seqRC.append(str(asmRseqs[i].seq.reverse_complement()))
        
        #write pairs to same file:
        asmFR_pair = []
        asmFR_pair.append(asmFseqs[i])
        asmFR_pair.append(asmRseqs[i])
        output_filename = "asmprimers-skpp20/"+ str(asmFseqs[i].id).replace('-F','')+".fasta"
        output_handle = open(output_filename, "w")
        SeqIO.write(asmFR_pair, output_handle, "fasta")
        
        asmFseqs_pass_seq_rec.append(asmFseqs[i])
        asmRseqs_pass_seq_rec.append(asmRseqs[i])
        output_handle.close()
        pass_count += 1

pass_count = 0
#loop over each Amp primers pair
for i in range(0,num_primers_to_screen):
    #if both primers in pair passed
    if (i in ampFseqs_pass) and (i in ampRseqs_pass) and (pass_count < num_primers_to_output):
        ampFseqs_pass_seq.append(str(ampFseqs[i].seq))
        ampRseqs_pass_seq.append(str(ampRseqs[i].seq))
        ampRseqs_pass_seqRC.append(str(ampRseqs[i].seq.reverse_complement()))
        ampAllseqs_pass_seq_rec.append(ampFseqs[i])
        ampAllseqs_pass_seq_rec.append(ampRseqs[i])
        pass_count += 1

#write passed primers to STDOUT:

print("asmFseqs")
print(asmFseqs_pass_seq)

print("asmRseqs (not RC)")
print(asmRseqs_pass_seq)

print("asmRseqs (RC)")
print(asmRseqs_pass_seqRC)

print("ampFseqs")
print(ampFseqs_pass_seq)

print("ampRseqs (not RC)")
print(ampRseqs_pass_seq)

print("ampRseqs (RC)")
print(ampRseqs_pass_seqRC)

#make a file with all amp primer pairs. This has both F and R primers, one after another.
output_handle = open("ampprimers-skpp15/ampAll_200_selected.fasta", "w")
SeqIO.write(ampAllseqs_pass_seq_rec, output_handle, "fasta")
output_handle.close()

#make a file with all AsmF selected.
output_handle = open("asmprimers-skpp20/asmF_selected.fasta", "w")
SeqIO.write(asmFseqs_pass_seq_rec, output_handle, "fasta")
output_handle.close()

#make a file with all AsmR selected.
output_handle = open("asmprimers-skpp20/asmR_selected.fasta", "w")
SeqIO.write(asmRseqs_pass_seq_rec, output_handle, "fasta")
output_handle.close()
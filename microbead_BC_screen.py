from Bio import SeqIO
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def getFastaSeqs(filename):
    fastaseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastaseqs.append(seqrec)
    handle.close()
    return fastaseqs

BC_file = 'barcodes/filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction.fasta'

BCseqs = getFastaSeqs(BC_file)

counter = 0
BCuniq = dict()
for primer in BCseqs:
    if str(primer.seq) in BCuniq:
        print(primer.id + " is a duplicate of " + BCuniq[str(primer.seq)])
    else:
        BCuniq[str(primer.seq)] = primer.id
    counter += 1

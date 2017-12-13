import os, math, subprocess, copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import string
import functools

def cmp(a, b):
    return (a > b) - (a < b) 

def get_class(kls):
    parts = kls.split('.')
    module = ".".join(parts[:-1])
    m = __import__( module )
    for comp in parts[1:]:
        m = getattr(m, comp)            
    return m

def secondaryStructure(primers, cutoff):
    newprimers=[]
    for primer in primers:
        score = calcSecondaryStructure(primer[0])
        if score > cutoff:
            newprimers.append(primer)
    return newprimers

def selfDimers(primers, cutoff):
    newprimers = []
    for primer in primers:
        score = primerdimers(primer[0],primer[0])
        if not score > cutoff:
            newprimers.append(primer)

    return newprimers

def meltingTemperature(primers, templow, temphigh):
    newprimers = []
    for primer in primers:
        Tm = oligoTm(primer[0])
        #print(str(Tm) + '\t' + str(len(primer[0])))
        if Tm > templow and Tm < temphigh:
            avgmeltingtemp = (float(templow) + float(temphigh)) / 2
            deltatemp = Tm - avgmeltingtemp
            primer.append(deltatemp)
            newprimers.append(primer)
            
    return newprimers

def avoidedSeqsInOverlapRegions(seq, primers, seqsToAvoid, coordsOfAvoidedRegionsInSeq):
    # prevents specified sequences from touching the overlap regions
    newprimers = []
    for primer in primers:
        addPrimer = True
        coords = primer[1]
        for avoidedCoords in coordsOfAvoidedRegionsInSeq:
            if (coords[0] >= avoidedCoords[0] and coords[0] <= avoidedCoords[1]) or (coords[1] >= avoidedCoords[0] and coords[1] <= avoidedCoords[1]):
                # if the start point of the primer is inside an avoided region, don't add this primer
                # if the end point of the primer is inside and avoided region, don't add this primer
                addPrimer = False
        if addPrimer:
            newprimers.append(primer)
    return newprimers
            
def evaluateOverlaps(seq, overlapseqs, overlaptemps, deltaGThreshold, selfDimersThreshold):
    overlapseqs = secondaryStructure(overlapseqs, deltaGThreshold)
    overlapseqs = selfDimers(overlapseqs, selfDimersThreshold)
    overlapseqs = meltingTemperature(overlapseqs, overlaptemps[0], overlaptemps[1])
    if len(overlapseqs) > 0:
        #py2: overlapseqs.sort(lambda x, y: cmp(abs(x[2]), abs(y[2])))
        overlapseqs.sort(key=functools.cmp_to_key(lambda x, y: cmp(abs(x[2]), abs(y[2]))))
        finaloverlap = overlapseqs[0]
        return finaloverlap
    else:
        return []

def optimizeOverlap(overlap, seq, lengthleeway, positionleeway, overlaptemps, deltaGThreshold, selfDimersThreshold):
    #First Check Characteristics of all oligos at the same positions with different lengths
    
    newoverlap = overlap[:]
    poslees = list(range(positionleeway*-1,(positionleeway)+1))
    #py2: poslees.sort(lambda x, y: cmp(abs(x), abs(y)))
    poslees.sort(key=functools.cmp_to_key(lambda x, y: cmp(abs(x), abs(y))))
    #print(poslees)
    for pl in poslees:
        for ll in range(lengthleeway+1):
            #print("ll=" + str(ll))
            newoverlap = [overlap[0]+pl, overlap[1]+pl]
            ##print(overlap[0])
            ##print(overlap[1])
            ##print(pl)
            testseqcoords = []
            for i in range(ll+1):
                #print("i="+str(i))
                testseqcoords.append([i,ll-i])
                if not ll==0:
                    testseqcoords.append([i*-1,(ll-i)*-1])
            overlapseqs = []
            #print(testseqcoords)
            for coord in testseqcoords:
                newcoords = [newoverlap[0]+coord[0], newoverlap[1]+coord[1]]
                ##print(newoverlap[0])
                ##print(newoverlap[1])
                #print(str(newcoords) + '\t' + str(overlap[2]))
                if newcoords[1]<=overlap[2]:
                    ##print(newcoords[0])
                    ##print(newcoords[1])
                    overlapseq = seq[newcoords[0]:newcoords[1]]
                    overlapseqs.append([overlapseq, newcoords])
                #print(overlapseqs)
            finaloverlap = evaluateOverlaps(seq, overlapseqs, overlaptemps, deltaGThreshold, selfDimersThreshold)
            if finaloverlap:
                return finaloverlap[1]
        #print("position change")

def optimizedSplit(seq, oligosizemax, lengthleeway, positionleeway, avgoverlapsize, overlaptemps, deltaGThreshold, selfDimersThreshold, num_of_oligos):
    seqsize = len(seq)
    avgoligolength = ((seqsize - ((num_of_oligos+1)*avgoverlapsize))//num_of_oligos) + 1 #py3 changed to floor div
    ##print("seqsize: "+str(seqsize))
    ##print("num_of_oligos: "+str(num_of_oligos))
    #print(avgoligolength)
    #overlaps = unoptimizedOverlaps(seq, avgoverlapsize, oligosizemax-positionleeway)
    optoverlaps = []
    optov = [avgoligolength+avgoverlapsize,avgoligolength+2*avgoverlapsize, oligosizemax]
    ##print("avgoligolength: "+str(avgoligolength))
    ##print("avgoverlapsize: "+str(avgoverlapsize))
    ##print("oligosizemax: "+str(oligosizemax))
    #print('initial seed=' + str(optov))
    #for overlap in overlaps:
    while True:
        ##print("optov[0]: "+str(optov[0]))
        ##print("optov[1]: "+str(optov[1]))
        newoptov = optimizeOverlap(optov, seq, lengthleeway, positionleeway, overlaptemps, deltaGThreshold, selfDimersThreshold)
        optoverlaps.append(newoptov)
        #print(newoptov)
        if not newoptov:
            #no overlaps found suitable
            return optimizedSplit(seq, oligosizemax, lengthleeway, positionleeway, avgoverlapsize, overlaptemps, deltaGThreshold, selfDimersThreshold, num_of_oligos+1)
        if (newoptov[0]+oligosizemax) >= len(seq):
            break #end of sequence reached
        else:
            lengthoffset = 0
            if(len(optoverlaps)>1):
                lengthoffset = ((newoptov[0]-optov[0])*-1) + ((newoptov[1]-newoptov[0])-avgoverlapsize)*-1
            else:
                lengthoffset = ((newoptov[0]-avgoverlapsize)-avgoligolength)*-1 + ((newoptov[1]-newoptov[0])-avgoverlapsize)*-1
    
            optov = [newoptov[1]+avgoligolength+lengthoffset,newoptov[1]+avgoligolength+lengthoffset+avgoverlapsize, newoptov[0]+oligosizemax]
            #print(str(newoptov) + '\t' + str(optov) + '\t' + str(lengthoffset))

    oligos = []
    previndex = 0
    for i in optoverlaps:
        #print('overlap='+str(i) +'\toverlaplen='+ str(i[1]-i[0])+ '\tprimercoords='+str([previndex,i[1]]))
        oligos.append(seq[previndex:i[1]])
        previndex = i[0]
    oligos.append(seq[previndex:])
    
    if len(oligos)>num_of_oligos:
        return optimizedSplit(seq, oligosizemax, lengthleeway, positionleeway, avgoverlapsize, overlaptemps, deltaGThreshold, selfDimersThreshold, num_of_oligos+1)
    
    return oligos
    
def bufferBuildSequences(startingseqrec, forwardprimer, reverseprimer):
    updatedseq = forwardprimer.seq + startingseqrec.seq + reverseprimer.seq.reverse_complement()
    return updatedseq

def checkRestrictionEnzyme(oligos, restrictionsite, restrictionsiteRC):
    for oligo in oligos:
        if restrictionsite in oligo:
            return True
    for oligo in oligos:
        if restrictionsiteRC in oligo:
            return True
    return False

def padOligo(buildseq, oligo, padlength):
    padString = 'T' * padlength
    endpad = Seq(padString, generic_dna)
    newoligo = oligo
    #find oligo
    index = buildseq.find(oligo)
    #add 5' pad
      #if first oligo, add TT
    if index == -1:
        print("An error occurred while padding the oligo %s" % oligo)
        return "ERROR"
    elif index == 0:
        newoligo = endpad + oligo
    else:
        padString = buildseq[index-padlength:index]
        newoligo = padString + newoligo
    #add 3' pad
    #if last oligo, add TT
    index = index + len(oligo)
    if index == len(buildseq):
        newoligo = newoligo + endpad
    else:
        padString = buildseq[index:index+padlength]
        newoligo = newoligo + padString
    
    return newoligo

def bufferOligoSequences(buildseq, oligos, restrictionenzyme, restrictionenzymebuffer, forwardprimer, reverseprimer, plateforwardprimer, platereverseprimer):
    newoligos = []
    for oligo in oligos:
        newoligo = padOligo(buildseq, oligo, restrictionenzymebuffer)
        newoligo = Seq(restrictionenzyme.site, generic_dna)  + newoligo  + Seq(restrictionenzyme.site, generic_dna).reverse_complement()
        newoligo = forwardprimer.seq + newoligo + reverseprimer.seq.reverse_complement()
        newoligo = plateforwardprimer.seq + newoligo + platereverseprimer.seq.reverse_complement()
        newoligos.append(newoligo)
    return newoligos

def oligoTm(seqobj):
    """
    Takes either a SeqRecord object, a Seq object, or a string
    and computes the melting temp based on the NN model (yes?).
    This is Kun's code
    
    CHECK THE NN PARAMETERS
    From Uri Laserson
    """
    
    if isinstance(seqobj,SeqRecord):
        seq = str(seqobj.seq).upper()
    elif isinstance(seqobj,Seq):
        seq = str(seqobj).upper()
    elif isinstance(seqobj,str):
        seq = seqobj.upper()
    
    # set the default tm parameters
    C_primer = 200 # 200nM in standard pcr
    C_Mg = 1.5 #Taq Buffer i
    C_MonovalentIon = 70 #20mM Tris-Cl + 50mM KCL in Taq Buffer  
    C_dNTP = 0.25 #mM
    percentage_DMSO = 0
    percentage_annealed = 50
    
    percentage_annealed = percentage_annealed/100.0
    percentage_DMSO = percentage_DMSO/100.0
    
    # Some constants
    R = 1.987
    deltaH = dict()
    deltaS = dict()
    deltaH =  {"AA": -7.6,  "TT": -7.6, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,"CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,"CG": -10.6,"GC": -9.8, "GG": -8.0, "CC": -8.0, "A": 2.2, "T": 2.2, "G": 0.0, "C": 0.0}
    deltaS = {"AA": -21.3, "TT": -21.3, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,"CG": -27.2, "GC": -24.4, "GG": -19.9, "CC":-19.9, "A": 6.9, "T": 6.9, "G": 0.0, "C": 0.0}
    
    C_SodiumEquivalent = C_MonovalentIon + 120 * math.sqrt(C_Mg-C_dNTP)
    seqLength = len(seq)
    dH = 0.2 + deltaH[str(seq[0])] + deltaH[str(seq[len(seq)-1])]
    dS = -5.7 + deltaS[seq[0]] + deltaS[seq[len(seq)-1]]
    for i in range(0, seqLength - 1):
        dH += deltaH[str(seq[i:i+2])]
        dS +=  deltaS[seq[i:i+2]]
    dS = dS + 0.368 * seqLength * math.log(C_SodiumEquivalent/1000.0)
    # val = math.log(C_primer*(1-percentage_annealed)/percentage_annealed)
    Tm = (dH * 1000) / (dS + R * (math.log(C_primer*(1-percentage_annealed)/percentage_annealed)-21.4164)) - 273.15 - 0.75*percentage_DMSO
    return Tm

def basecompare(base1,base2):
    if base1=="G" and base2=="C":
        return 1
    elif base1=="C" and base2=="G":
        return 1
    elif base1=="A" and base2=="T":
        return 1
    elif base1=="T" and base2=="A":
        return 1
    else:
        return -1

def complementcompare(seq1, seq2):
    if len(seq1)!=len(seq2):
        print("Error: primer sequences in complementcompare are not equal")
    score = 0
    for i in range(len(seq1)):
        if i<10:
            score += basecompare(seq1[i],seq2[((i*-1)-1)])*2
        else:
            score += basecompare(seq1[i],seq2[((i*-1)-1)])

    return score

def primerdimers(primer1, primer2):
    score = []
    for i in range(len(primer1)):
        score.append(complementcompare(primer1[len(primer1)-i-1:len(primer1)], primer2[len(primer1)-i-1:len(primer1)]))
    #print(score)
    return max(score)

def oligolistTm(seqobjlist):
    Tm = []
    for i in seqobjlist:
        Tm.append(oligoTm(i))
    return Tm
    
def retrievePrimerHRregions(therecord, thefeature, updatedstart, updatedend, hrlength, numofprimers):
    #construct list of sequences to use
    updatedstart = updatedstart-hrlength
    updatedend = updatedend+hrlength
    myseq = therecord.seq[updatedstart:updatedend]
    strand = strandChooser(thefeature)
    #reverse complement if necessary
    if strand==1:
        #reverse complement
        myseq = myseq.reverse_complement()
    #construct possible forward subsequences:
    forwardhr = []
    reversehr = []
    for i in range(numofprimers):
        start = i*3
        end = i*3+hrlength
        if strand==1:
            reversehr.append(myseq[start:end])
        else:
            forwardhr.append(myseq[start:end])
    #construct possible reverse subsequences:
    seqlength = len(myseq)
    for i in range(numofprimers):
        start = (seqlength-hrlength-i*3)
        end = (seqlength-(i*3))
        if strand==1:
            forwardhr.append(myseq[start:end])
        else:
            reversehr.append(myseq[start:end])
    #return the two lists
    
    return forwardhr, reversehr

def calcSecondaryStructure(seq):
    cmd = "hybrid-ss-min --NA=DNA -q " + str(seq)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    p.wait()
    dg = p.stdout.read()
    ans = float(dg)
    return ans

def oligolistsecstructure(seqobjlist):
    dg = []
    for i in seqobjlist:
        dg.append(calcSecondaryStructure(i.seq))
    return dg

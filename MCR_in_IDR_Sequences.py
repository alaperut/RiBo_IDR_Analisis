from __future__ import division
import sys
from Bio import SeqIO
import csv,itertools,re, numpy
from localcider.sequenceParameters import SequenceParameters

FastaFile = sys.argv[1]
kappa = float(sys.argv[2]) # Lower limit of charge mixing. 1 is completely segregated and 0 is completely mixed
charge_bal = float(sys.argv[3]) # Fraction of charged residues (FCR) - 0.7
window_size = float(sys.argv[4]) #Length of region -60 in the paper
cut_off = float(sys.argv[5]) # The minimum score to be considered significant -0.6
Result = sys.argv[6] # want results to contain the kappa value, protein ID, start, % of each aa, and end of the motif, file name to include window size, cutoff, kappa charge balance 

def score_polyamp(Sequence, window, cutoff):
    #Sequence: amino acid seq, window: window size, cutoff: the minimum score to be considered as significant
    #returns a string combining all the mixed charge regions and a string of the start and end loci of all mcrs
    polyamph_regions = []
    polyamph_Start = []
    polyamph_End = []
    polyamph_start = 0
    polyamph_end = 0
    last_start = 0
    last_end = 0
    i = 0
    window = int(window)
    for j in range(window//2,1+len(Sequence)-window//2): # double slash rounds to the nearest whole number, scan through the sequence.
        region = ''
        start = j-(window//2)
        end = j+(window//2)
        region = Sequence[start:end]
        pos = region.count('R')+region.count('K') # can look for anything really though, can use this to look for RS regions as well
        neg = region.count('E')+region.count('D')
        score = pos+neg
        if score>=cutoff: # if found a MCR
            if polyamph_start == 0: # if first region found in this sequence
                extended_reg = Sequence[polyamph_start:end]
                if extended_reg.count('R')+extended_reg.count('K')+extended_reg.count('E')+extended_reg.count('D')>=cutoff: # checks to see if extended region is still valid
                    polyamph_start = start # new start to recorded mcr
                    polyamph_end = end # new end to recorded mcr
                    last_start = polyamph_start 
                    last_end = polyamph_end
            elif polyamph_end+5 > start: #if next start is within 5 from last end, join the 2
                if extended_reg.count('R')+extended_reg.count('K')+extended_reg.count('E')+extended_reg.count('D')>=cutoff:
                    polyamph_end = end
                    last_start = polyamph_start
                    last_end = polyamph_end
            else:
                if extended_reg.count('R')+extended_reg.count('K')+extended_reg.count('E')+extended_reg.count('D')>=cutoff: # append another mcr if it's far away in sequence
                    polyamph_regions.append(Sequence[polyamph_start:polyamph_end+1])
                    polyamph_Start.append(polyamph_start)
                    polyamph_End.append(polyamph_end+1)
                    polyamph_start = 0
                    polyamph_end = 0

    if last_end!=0: 
        temp1 = Sequence[last_start:last_end+1]
        if temp1 not in polyamph_regions:
            polyamph_regions.append(temp1)
            polyamph_Start.append(last_start)
            polyamph_End.append(last_end+1)
    return polyamph_regions, polyamph_Start, polyamph_End

def get_polyamp_regions(FastaFile, Result):
    # look for mixed charge regions in each protein. Each entry is preceeded by its accession
    seqs,out,added = [],[],{}
    with open(FastaFile,'r') as FastaFile, open(Result,'w') as Result:
        InputFile = SeqIO.parse(FastaFile,'fasta') # Read FastaFile
        Result.write("accession" + '\t' + "highly charged region" + '\t' + "start" + '\t' + "end" + '\t' + "r/(r+k)" + '\t' + "d/(d+e)" + '\t' +"charge_density" '\t' + "kappa" + '\n')
        for Sequence in InputFile:
            AccessionNumber = Sequence.id
            print(Sequence.seq)
            t2=score_polyamp(str(Sequence.seq),window_size,cut_off*window_size) #s.seq must be part of the SeqIO package, pass into function above
            if len(t2[0])>0:
                for i in range(len(t2[0])):
                    r,k = t2[0][i].count('R'), t2[0][i].count('K')
                    d,e = t2[0][i].count('D'), t2[0][i].count('E')
                    neg = d+e 
                    pos = r+k
                    if neg == 0 or pos == 0:
                        PosNeg = -1
                        NegPos = -1
                    if neg > 0 and pos > 0:
                        PosNeg = pos/neg
                        NegPos = neg/pos
                    if PosNeg>=charge_bal and NegPos>=charge_bal:
                        l = len(t2[0][i])
                        r_ratio = r/pos
                        d_ratio = d/neg
                        region_noX = ''
                        for char in t2[0][i]:
                            if char=='X':
                                region_noX = region_noX+'S'
                            elif char=='S' or char == 'U':
                                region_noX = region_noX+'T'
                            else:
                                region_noX = region_noX+char
                        SeqOb = SequenceParameters(region_noX)
                        kap = SeqOb.get_kappa_X(grp1 = ['E','D'], grp2 = ['K','R']) # get_kappa_X is from localcider SequenceParameters https://github.com/Pappulab/localCIDER/blob/master/localcider/sequenceParameters.py
                        if kap>=kappa: # if the calculated kappa is over the threshold
                            j = i
                            Result.write(str(AccessionNumber) + '\t' + str(t2[0][i]) + '\t' + str(t2[1][i]) + '\t' + str(t2[2][i]) + '\t' + str(r_ratio) + '\t' + str(d_ratio) + '\t' + str((pos+neg)/l) + '\t' + str(kap) + '\n')



get_polyamp_regions(FastaFile, Result)


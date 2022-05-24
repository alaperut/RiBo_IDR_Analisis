import sys
ProteinCodes = open(sys.argv[1], "r")
MetaResults = open(sys.argv[2], "r")
IUPredResults = open(sys.argv[3], "r")
FinalIDRs = open(sys.argv[4], "w")
ProteinCodesCache = ProteinCodes.readlines()
MetaResultsCache = MetaResults.readlines()
IUPredResultsCache = IUPredResults.readlines()
ProteinCodes.close()
MetaResults.close()
IUPredResults.close()
UniprotID = ""
Start = 0
End = 0
StartMet = 0
EndMet = 0
StartPred = 0
EndPred = 0
Test = True
Keep = True

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

for i in ProteinCodesCache:
    UniprotID = i[75:82]
    if i[82]=="-":
        UniprotID = i[75:84]
    for x in MetaResultsCache:
        wrong = UniprotID + "-"
        if UniprotID in x and wrong not in x:
            elem = x.split("\t")
            StartMet = int(elem[1])
            EndMet = int(elem[2])
            Start = StartMet
            End = EndMet
            for y in IUPredResultsCache:
                if UniprotID in y:
                    line = y.split("\t")
                    StartPred = int(line[2])
                    EndPred = int(line[3])
                    Keep = True
                    if StartPred>StartMet:
                        Start = StartPred
                    if EndPred<EndMet:
                        End = EndPred
                    if Start>End:
                        Start = StartMet
                        End = EndMet
                        Keep = False
                    if Keep == True:
                        FinalIDRs.write(UniprotID + '\t' + str(Start) + '\t'+ str(End) + '\n')
                        Start = 0
                        End = 0
FinalIDRs.close()
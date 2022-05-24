import sys
file = sys.argv[1]
new = sys.argv[2] #Just number of IDRs per protein
IDRStartEnd = sys.argv[3]
stringency = int(sys.argv[4])
length = int(sys.argv[5])
IDRlength = 0
FoundIDR = 0
PotentialIDR = 0
IDRcount = 0
TotalIDR = 0
ProteinCount = 0
IDRStartArray = []
IDREndArray = []
ProteinNameArray = []
results = open(file, "r")
newfile = open(new, "w")
IDRStartEnd = open(IDRStartEnd, "w")
AACount = 0
value = 0
z = 0
i = 0
y = 0
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
        

for line in results:
    x = is_int(line[1])
    if len(line)> 6 and x == False: # start a new protein
        if FoundIDR == 0 and IDRcount > 0:  #if there was an IDR in the protein, but not at the end
            newfile.write(str(IDRcount) + '\n')
        if FoundIDR == 0 and IDRcount == 0:
            if z == 1:
                newfile.write(str(IDRcount) + '\n')
            z = 1
        if FoundIDR == 1: # if there's an IDR all the way to the end
            IDRcount = IDRcount + 1
            FoundIDR = 0
            newfile.write(str(IDRcount) + '\n')
            IDREndArray.append(AACount)
        ProteinName = line[0:len(line)-1]
        newfile.write(ProteinName + '\t')
        IDRcount= 0
        AACount = 0
        PotentialIDR = 0
        value = 0
        IDRlength = 0

    if x:
        if len(line) > 40:
            line = line.replace(",", "")
        AACount = AACount + 1
        PreviousValue = value
        value = 0
        if line[1] == "1":
            value= 10
        if len(line) > 2: 
            value = int(line[3])
        if value >= stringency: # Must track the AANumber Start and End sites, current problem- that start site at 1, length doesn't apply
            if FoundIDR == 0 and PreviousValue < stringency:
                IDRStart = AACount
            if FoundIDR == 0 and AACount == 1:
                IDRStart = 1
            IDRlength = IDRlength + 1
            PotentialIDR = 1
            if IDRlength >= length:
                if FoundIDR == 0:
                    IDRStartArray.append(IDRStart)
                    ProteinNameArray.append(ProteinName)
                FoundIDR = 1
        if value < stringency: 
            if PotentialIDR == 1 and IDRlength < length:
                IDRlength = 0
                FoundIDR = 0
            if FoundIDR == 1:
                FoundIDR = 0
                IDRlength = 0
                TotalIDR = TotalIDR + 1
                IDRcount = IDRcount + 1   
                IDREnd = AACount - 1
                IDREndArray.append(IDREnd)
if FoundIDR == 1:
    TotalIDR = TotalIDR + 1
    IDRcount = IDRcount + 1   
    newfile.write(IDRCount + '\n')
if FoundIDR == 0:  #if there was an IDR in the protein, but not at the end
    newfile.write(str(IDRcount) + '\n')
print(IDREndArray)
while i < len(IDRStartArray):
    IDRStartEnd.write(ProteinNameArray[i] + '\t' + str(IDRStartArray[i]) + '\t' + str(IDREndArray[i]) + '\n')
    i = i + 1
 #           ProteinArray.extend("\n" "~ this IDR is" + str(IDRlength) + "amino acids long" "\n")

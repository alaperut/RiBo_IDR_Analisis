import sys
results = sys.argv[1]
ProteinCodes = sys.argv[2]
new = sys.argv[3]
new2 = sys.argv[4]
stringency = int(sys.argv[5])*100
length = int(sys.argv[6])
IDRlength = 0
FoundIDR = 0
PotentialIDR = 0
IDRcount = 0
TotalIDR = 0
ProteinCount = 0
x = 0
results = open(results, "r")
newfile = open(new, "w")
IDRperProtein = open(new2, "w")
ProteinArray = []
IDRArray = []
PotentialIDRArray = []
ProteinIDR = 0
w = 0

file2 = open(ProteinCodes, "r")
f2cache=file2.readlines()
# close the files
file2.close()


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
        

for line in results:
    if line[0] == ">":
        if FoundIDR == 1:
            IDRcount = IDRcount + 1
            TotalIDR = TotalIDR + 1
            FoundIDR = 0
            ProteinArray.extend(IDRArray)
        if FoundIDR == 0 and IDRcount > 0:
            newfile.write("\n".join(str(elem) for elem in ProteinArray))
            ProteinArray = []
            IDRperProtein.write(str(IDRcount)+'\n')
 #           newfile.write ("\n" "- there are" + str(IDRcount) + "IDRs in this protein" "\n")
        if FoundIDR == 0 and IDRcount == 0 and x == 1:
            IDRperProtein.write(str(IDRcount)+'\n')
        ProteinCount = ProteinCount + 1
        IDRlength = 0
        IDRcount= 0
        PotentialIDR = 0
        IDRArray = []
        x = 0
        PotentialIDRArray = []
        ProteinArray.append(line)
        for h in f2cache:
            w = w + 1
            if h[75:82] in line or h[95:101] in line:
                ProteinArray.append(f2cache[w-1][131:134])
        w = 0
        ProteinArray.append("# POS	AMINO ACID	IUPRED SCORE	ANCHOR SCORE")
    if is_int(line[10:14]):
        x = 1
        value = int(line[10:14])
        if value >= stringency:
            PotentialIDRArray.append(line)
            IDRlength = IDRlength + 1
            PotentialIDR = 1
            if IDRlength >= length:
                FoundIDR = 1
                IDRArray.extend(PotentialIDRArray)
                PotentialIDRArray = []
        if value < stringency: 
            if PotentialIDR == 1 and IDRlength < length:
                PotentialIDRArray = []
                IDRlength = 0
            if FoundIDR == 1:
                FoundIDR = 0
                IDRlength = 0
                TotalIDR = TotalIDR + 1
                IDRcount = IDRcount + 1   
                ProteinArray.extend(IDRArray)
                IDRArray = []
                PotentialIDRArray = []
if FoundIDR == 1:
    TotalIDR = TotalIDR + 1
    ProteinIDR = ProteinIDR + 1
    IDRcount = IDRcount + 1   
    ProteinArray.extend(IDRArray)
newfile.write("\n".join(str(elem) for elem in ProteinArray))
IDRperProtein.write(str(IDRcount)+'\n')
newfile.close()
IDRperProtein.close()
 #           ProteinArray.extend("\n" "~ this IDR is" + str(IDRlength) + "amino acids long" "\n")

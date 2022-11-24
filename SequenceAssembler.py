import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

inputFileName = sys.argv[1]
scoreForMatch = int(sys.argv[2])
penaltyForReplace = int(sys.argv[3])
penaltyForIndel = int(sys.argv[4])
outputFileName = sys.argv[5]

def scoreij(f1, f2, i, j):
    if f1[i-1] == f2[j-1]:
        return scoreForMatch
    else:
        return penaltyForReplace

def backTrack(path, f1, f2, maxi, maxj):
    combined = ""
    i = maxi
    j = maxj
    if maxi == len(f1):
        combined = f2[j:]
    else:
        combined = f1[i:]

    while path[i][j] != 'n':
        if path[i][j] == 'm':
            combined = f1[i-1] + combined
            i = i - 1
            j = j - 1
        elif path[i][j] == 'i':
            combined = f1[i-1] + combined
            j = j - 1
        else:
            combined = f2[j-1] + combined
            i = i - 1

    if maxi == len(f1):
        combined = f1[:i] + combined
    else:
        combined = f2[:j] + combined
    return combined

def getScore(f1, f2):
    dp = [[0 for i in range(len(f2)+1)] for j in range(len(f1) + 1)]
    path = [['n' for i in range(len(f2)+1)] for j in range(len(f1) + 1)]
    max_score = -999999999
    maxi = 0
    maxj = 0

    for i in range(1, len(f1)+1):
        for j in range(1, len(f2)+1):
            del_score = dp[i-1][j] + penaltyForIndel
            ins_score = dp[i][j-1] + penaltyForIndel
            mut_score = dp[i-1][j-1] + scoreij(f1, f2, i, j)
            m = max(ins_score, del_score, mut_score)
            dp[i][j] = m
            if m == mut_score:
                path[i][j] = 'm'
            elif m == ins_score:
                path[i][j] = 'i'
            elif m == del_score:
                path[i][j] = 'd'

            if i == len(f1) or j == len(f2):
                if m > max_score:
                    max_score = m
                    maxi = i
                    maxj = j
    combined = ""
    if max_score > -999999999:
        combined = backTrack(path, f1, f2, maxi, maxj)
    
    return max_score, combined

fragments = []
for seq_record in SeqIO.parse(inputFileName, "fasta"):
    fragments.append(str(seq_record.seq))

while True:
    max_score = -999999999
    index1 = -99
    index2 = -99
    comFrag = ""
    
    for j in range(len(fragments)):
        for k in range(j+1, len(fragments)):
            score, frag = getScore(fragments[j], fragments[k])
            if(score > max_score):
                index1 = j
                index2 = k
                max_score = score
                comFrag = frag
    if max_score > 0:
        fragments.pop(index2)
        fragments.pop(index1)
        fragments.append(comFrag)

    if max_score < 0 or len(fragments) == 1:
        break

longestFrag = max(fragments, key=len)
sequences = [SeqRecord(Seq(longestFrag))] # TODO: Add id and description of original input
SeqIO.write(sequences, outputFileName, "fasta")


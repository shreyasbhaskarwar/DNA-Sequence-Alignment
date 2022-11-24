import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq
import random

inputFileName = sys.argv[1]
noOfSequences = int(sys.argv[2])
mutationProbability = float(sys.argv[3])
outputFileName = sys.argv[4]

replaceProbability = 0.8 * mutationProbability
deleteProbability = 0.2 * mutationProbability

# 0 -> Replace, 1 -> Delete, 2 -> No Mutation
def getMutationType():
    weights = [replaceProbability, deleteProbability, 1 - mutationProbability]
    return random.choices(range(3), weights)[0]

def replaceNucleotide(n):
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(n)
    return random.choices(nucleotides)[0]

seqRecord = SeqIO.read(inputFileName, "fasta")
sequences = [seqRecord]

for i in range(noOfSequences-1):
    mutableSeq = MutableSeq(seqRecord.seq)
    for j in range(len(seqRecord.seq)):
        mutationType = getMutationType()
        if mutationType == 0:
            mutableSeq[j] = replaceNucleotide(mutableSeq[j])
        elif mutationType == 1:
            mutableSeq[j] = '-'
        else:
            continue
    # Remove deleted nucleotides
    mutableSeq = mutableSeq.replace("-","")
    #Remove id from description to maintain fasta structure of input
    sequences.append(SeqRecord(Seq(mutableSeq), seqRecord.id+str(i), description = seqRecord.description.replace(seqRecord.id, "").strip()))

SeqIO.write(sequences, outputFileName, "fasta")
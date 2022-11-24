import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
import random
from Bio.Seq import Seq

inputFileName = sys.argv[1]
xMinFragLength = int(sys.argv[2])
yMaxFragLength = int(sys.argv[3])
zAcceptableFragLength = int(sys.argv[4])
outputFileName = sys.argv[5]

lstFragment = []

records = list(SeqIO.parse(inputFileName, "fasta"))

for record in records:
    my_seq = Seq(record.seq)
    recordLength = len(my_seq)
    i = 0
    while(recordLength >= xMinFragLength):

        fragmentSize = random.randint(xMinFragLength,yMaxFragLength)
        fragment = my_seq[i:i+fragmentSize]
        i = i + fragmentSize
        recordLength = recordLength - fragmentSize
        if fragmentSize <= zAcceptableFragLength:
            randomid = str(random.randint(0,100000))
            fragmentId = "SeqID"+randomid
            fragmentDescription = "Description"+randomid
            lstFragment.append(SeqRecord(fragment,fragmentId,"Name",fragmentDescription))

SeqIO.write(lstFragment, outputFileName, "fasta")

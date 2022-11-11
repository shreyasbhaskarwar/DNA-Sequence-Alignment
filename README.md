# DNA Sequence Alignment

The aim of this project is to gain knowledge on sequence alignment by implementing the following three simulations:
1. Sequence generator
2. Sequence partitioning
3. Sequence assembler

## Implementation Details

### Sequence generator

This program generates a given number of sequences from a parent sequence using probability-based mutation.

It takes the following input parameters:

1. Input file name which contains DNA sequence in FASTA format
2. Number of sequences (_k_)
3. Mutation probability in [0:1] interval (_p_)
4. Output file name which will contain mutated DNA sequence in FASTA format

The program will output the nucleotide sequences in FASTA format. The first sequence will be the input sequence. Next _k-1_ sequences will be a copy of the first sequence, but each nucleotide has a _4p/5_ probability to get replaced with a randomly selected nucleotide and a _p/5_ probability to get deleted.

### Sequence partition

This program is used to generate fragments from a given input of mutated sequences.

It takes the following input parameters:

1. Input file name (output of sequence generator)
2. Minimum fragment length
3. Maximum fragment length
4. Maximum acceptable fragment length 
5. Output file name

The program will partition each sequence in the input file into smaller fragments and store them in the output file in FASTA format. The fragment length is chosen randomly and lies between the minimum and maximum fragment length. If the length is more than the acceptable fragment length, the fragment is discarded. While chopping, if there is a small piece left with a length less than the minimum fragment length, it will be discarded.

### Sequence assembler

This program will assemble the fragments of sequences in a file into one long sequence.

It takes the following input parameters:

1. Input file name (output of sequence partition) 
2. Score for match (positive integer)
3. Penalty for replace (negative integer)
4. Penalty for delete/insert (negative integer)
5. Output file name

The program will combine fragments using a greedy strategy as follows:

1. Align all fragments _f<sub>i</sub>_ and _f<sub>j</sub>_ with each other using dovetail alignment and compute the alignment score _v<sub>i,j</sub>_.
2. **Repeat**
   - Merge the two fragments _f<sub>i</sub>_ and _f<sub>j</sub>_ with the largest alignment score _v<sub>i,j</sub>_ into a
new fragment _f'_, if _v<sub>i,j</sub>_ > 0
   - Replace _f<sub>i</sub>_ and _f<sub>j</sub>_ with the new fragment _f′_ and compute align _f′_ with all the
other fragments in the current set using dovetail alignment.
3. **Until** (one fragment is left in the set) OR (the largest alignment score is negative)
4. Write the longest fragment you assembled into output file in FASTA format.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install Biopython library.

```bash
pip install biopython
```

Run sequence generator:
```bash
python SequenceGenerator.py input.fasta 10 0.005 output1.fasta
```

Run sequence partition:
```bash
python SequencePartition.py output1.fasta 100 220 200 output2.fasta
```

Run sequence assembler:
```bash
python SequenceAssembler.py output2.fasta 1 -1 -3 output3.fasta
```

## Authors

- Shreyas Bhaskarwar
- Bhavya Kavdia
- Venkata Shiva Reddy Manchala

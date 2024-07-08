# KMP Parallel Search
## DESCRIPTION

In string computation, the exact pattern matching problem is the problem of
finding all the occurrences of a pattern (string) P, in a text (string) S, where usually P is much
shorter than S. For example the pattern could be the world “stella” and the text the whole Divina
Commedia, or P can be the CCATTGTG motif and the text the human genome.
One strategy to speed up the computation is to create an index on the pattern P and use this index to
scan the text S in a more efficient way.
The Knuth-Morris-Pratt algorithm uses this approach. It first of all builds an index on P and then
uses it to scan S, applying simple rules to the index to decide how to shift the pattern.
Expected outcome: design a Python script that:
1. Takes as input:
a. A FASTA file containing a genome (length ~10kbp), e.g., the genome of SARSCOV-2 (https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)
b. A FASTA file containing a set of short sequences (order of about 100 sequences)
2. Uses KMP to identify in parallel (with multiprocessing) all matches of the short sequences
within the longer genome sequence, and stores the matches, one per line, in a single file in
an appropriate format (information about the matching sequence and location of the match
in the genome).
3. Uses a Lock to manage the access to the output file.

## Requirements

- Python3
- Biopython
- NumPy
- Pandas

```bash
pip install biopython numpy pandas
```

##USAAGE:
To use this script run this comand:
```bash
python kmp_parallel.py <genome_file> <sequences_file> <output_file>
```
where:
<genome_file> is the paths to the genome FASTA file;
<sequences_file> is the path to the sequences FASTA file;
<output_file> is the path of the output file with all the results.

## FUNCTIONS DESCRIPTION:
create_random_sequences(filename, num_sequences, min_len, max_len): Generates a FASTA file with randomic sequences.
create_kmp_table(pattern): Creates the KPM table for the given pattern.
kmp_search(args): Performs the KMP search for the pattern.
parse_arguments(): Parses command line arguments.
load_genome(filename): Loads the genome sequence from a FASTA file.
load_sequences(filename): Loads sequences for searching them from a FASTA file.
main(): Main function to execute the KMP search in parallel and saving all the results.

## EXAMPLE
example for running this code (in the folder containing this project i inserted this 2 files to show a functional example of the program):
```bash
python kmp_parallel.py genome.fasta sequences.fasta output.csv
```

## TESTING:
With the use of 'unittest' is possible to test the creation and table KPM.
Use this comand to run the testing:
```bash
python -m unittest test_kmp.py
```

## GitHub REPOSITORY:
To clone the repository use this comand:
```bash
git clone https://github.com/maribottini/kmp_parallel_search.git
cd kmp_parallel_search

```

## VERSION:
Code created and tested on:
Python 3.10.12

##LICENSE
MIT License

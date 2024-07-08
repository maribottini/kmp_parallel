import sys
import os
import random
import string
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, Lock, Manager

def create_random_sequences(filename, num_sequences=100, min_len=5, max_len=20):
    """Create a random fasta file (with random sequences).

    Args:
        filename (str): The name of the file where to write the sequences into.
        num_sequences (int): The number of random sequences we want to generate.
        min_len (int): The minimum length of the sequences.
        max_len (int): The maximum length of the sequences.
    """
    try:
        nucleotides = ['A', 'T', 'C', 'G']
        with open(filename, 'w') as file:
            for i in range(num_sequences):
                length = random.randint(min_len, max_len)
                sequence = ''.join(random.choices(nucleotides, k=length))
                file.write(f">seq{i+1}\n{sequence}\n")
        print(f"Generated {num_sequences} random sequences in {filename}")
    except IOError as e:
        print(f"Error during the writing into the file {filename}: {e}")

def create_kmp_table(pattern):
    """Creates KMP table for a given pattern.

    Args:
        pattern (str): The pattern string for which to create the KMP table.

    Returns:
        np.ndarray: The KMP table for the pattern.
    """
    length = 0
    kmp_table = np.zeros(len(pattern), dtype=int)
    i = 1
    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            kmp_table[i] = length
            i += 1
        else:
            if length != 0:
                length = kmp_table[length - 1]
            else:
                kmp_table[i] = 0
                i += 1
    return kmp_table

def kmp_search(args):
    """Perform KMP search for the pattern in the text of the genome.

    Args:
        args (tuple): A tuple containing the genome, pattern, lock, and results.
    """
    genome, pattern, lock, results = args
    m = len(pattern)
    n = len(genome)
    kmp_table = create_kmp_table(pattern)
    i = 0  # index for genome
    j = 0  # index for pattern
    while i < n:
        if pattern[j] == genome[i]:
            i += 1
            j += 1
        if j == m:
            with lock:
                results.append((pattern, i - j))
            j = kmp_table[j - 1]
        elif i < n and pattern[j] != genome[i]:
            if j != 0:
                j = kmp_table[j - 1]
            else:
                i += 1

def parse_arguments():
    """Parse command line arguments.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description='KMP Algorithm for Pattern Matching')
    parser.add_argument('genome_file', type=str, help='Path to the genome FASTA file')
    parser.add_argument('sequences_file', type=str, help='Path to the sequences FASTA file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    return parser.parse_args()

def load_genome(filename):
    """Load genome sequence from a FASTA file.

    Args:
        filename (str): Path to the genome FASTA.

    Returns:
        str: The genome sequence as an only string.
    """
    genome = ""
    try:
        for record in SeqIO.parse(filename, "fasta"):
            genome += str(record.seq)
        if not genome:
            raise ValueError(f"No sequences found in file {filename}")
        print(f"Loaded genome from {filename}")
    except IOError as e:
        print(f"Error reading file {filename}: {e}")
    except ValueError as e:
        print(e)
    return genome

def load_sequences(filename):
    """Load sequences for searching them into a FASTA file.

    Args:
        filename (str): Path to the FASTA sequences.

    Returns:
        list: of sequences.
    """
    sequences = []
    try:
        for record in SeqIO.parse(filename, "fasta"):
            sequences.append(str(record.seq))
        if not sequences:
            raise ValueError(f"No sequences found in file {filename}")
    except IOError as e:
        print(f"Error reading file {filename}: {e}")
    except ValueError as e:
        print(e)
    return sequences

def main():
    """Main function for the KMP search."""
    # Parse command line arguments
    args = parse_arguments()

    # Read the genome sequence
    genome = load_genome(args.genome_file)

    # Read the sequences to search for
    sequences = load_sequences(args.sequences_file)

    # Set up parallel processing
    manager = Manager()
    results = manager.list()
    lock = manager.Lock()

    # Create jobs for each sequence
    jobs = [(genome, seq, lock, results) for seq in sequences]

    # Process jobs in parallel
    with Pool() as pool:
        pool.map(kmp_search, jobs)

    # Convert results to DataFrame
    df = pd.DataFrame(list(results), columns=['Pattern', 'Position'])

    # Write results to the output file
    try:
        df.to_csv(args.output_file, index=False, header=True)
        print(f"Results written to {args.output_file}")
    except IOError as e:
        print(f"Error writing to file {args.output_file}: {e}")

if __name__ == "__main__":
    main()

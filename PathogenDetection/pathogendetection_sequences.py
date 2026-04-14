#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Caroline Scranton
Purpose: Extract actual k-mers mapped to specific taxIDs from Kraken2 outputs,
         reconstruct read sequences with Ns for unmatched bases,
         requiring >=3 k-mers per taxid,
         and match Kraken2/FASTQ pairs based on shared filename prefixes
         (e.g. AERO_123_un_merged.fastq <-> AERO_123_un_merged_output.txt).
"""

import argparse
import os
import re
from collections import defaultdict

def get_args():
    """Parameters for script"""
    parser = argparse.ArgumentParser(description='Extract k-mers matching taxID from Kraken2 outputs and FASTQ/A files in folders')
    parser.add_argument('kraken_folder', help='Folder containing Kraken2 output files')
    parser.add_argument('seq_folder', help='Folder containing corresponding FASTQ or FASTA files')
    parser.add_argument('-t', '--taxids', nargs='+', required=True, help='NCBI TaxIDs to extract')
    parser.add_argument('-k', '--kmer_length', type=int, default=35, help='K-mer length (default 35)')
    return parser.parse_args()

def detect_format(file_path):
    """Detects FASTA or FASTQ based on first character in first line"""
    with open(file_path) as f:
        first = f.readline()
        return 'fastq' if first.startswith('@') else 'fasta'

def normalize_read_id(header):
    """Extract cleaned read ID"""
    return header.strip()[1:].split()[0].split('|')[-1]

def parse_kraken_output(kraken_file, target_taxids):
    """Parse Kraken2 output: map read IDs to taxid-specific k-mer start positions"""
    taxid_set = set(target_taxids)
    read_tax_kmers = defaultdict(lambda: defaultdict(list))
    with open(kraken_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            _, raw_read_id, _, _, mapping_str = parts
            read_id = raw_read_id.split('|')[-1]
            kmers = mapping_str.strip().split()
            pos = 0
            for pair in kmers:
                taxid, count = pair.split(':')
                count = int(count)
                if taxid in taxid_set:
                    read_tax_kmers[read_id][taxid].append((pos, count))
                pos += count
    return read_tax_kmers

def get_base_prefix(filename):
    """
    Get the common prefix for matching FASTQ and Kraken2 output files.
    Example:
      AERO_123_un_merged.fastq         -> AERO_123_un_
      AERO_123_un_merged_output.txt    -> AERO_123_un_
    """
    name = filename
    name = re.sub(r'(_merged_output)?\.txt$', '', name)
    name = re.sub(r'\.f(ast)?q$', '', name)   # .fq or .fastq
    name = re.sub(r'\.fasta$', '', name)
    name = re.sub(r'_merged$', '', name)
    return name

def extract_reads_from_fasta(f, read_tax_kmers, kmer_length):
    """Extract reconstructed sequences from FASTA"""
    header, seq = '', ''
    for line in f:
        if line.startswith('>'):
            if header:
                read_id = normalize_read_id(header)
                yield from extract_kmers_for_read(read_id, seq, read_tax_kmers, kmer_length)
            header, seq = line, ''
        else:
            seq += line.strip()
    if header:
        read_id = normalize_read_id(header)
        yield from extract_kmers_for_read(read_id, seq, read_tax_kmers, kmer_length)

def extract_reads_from_fastq(f, read_tax_kmers, kmer_length):
    """Extract reconstructed sequences from FASTQ"""
    while True:
        header = f.readline()
        if not header:
            break
        seq = f.readline().strip()
        f.readline()  # +
        f.readline()  # quality
        read_id = normalize_read_id(header)
        yield from extract_kmers_for_read(read_id, seq, read_tax_kmers, kmer_length)

def extract_kmers_for_read(read_id, seq, read_tax_kmers, kmer_length):
    """Reconstruct a single contiguous sequence per read with Ns in between k-mers"""
    if read_id not in read_tax_kmers:
        return
    for taxid, starts in read_tax_kmers[read_id].items():
        if len(starts) < 3:
            continue
        # Create a single sequence for the read
        out_seq = ['N'] * len(seq)
        for pos, count in starts:
            for i in range(pos, min(pos + kmer_length, len(seq))):
                out_seq[i] = seq[i]
        # Yield **one sequence per read per taxid**
        yield taxid, read_id, ''.join(out_seq)

def write_reads_to_files(read_blocks, base_prefix, output_dir):
    """Write reconstructed sequences to taxid-specific files"""
    outputs = defaultdict(list)
    for taxid, read_id, seq in read_blocks:
        outputs[taxid].append((read_id, seq))

    for taxid, reads in outputs.items():
        if len(reads) < 3:
            print(f"Skipping output for taxID '{taxid}' due to insufficient reads (found {len(reads)}).")
            continue
        filename = os.path.join(output_dir, f"{base_prefix}_{taxid}.fasta")
        with open(filename, 'w') as f:
            for read_id, seq in sorted(reads):
                f.write(f">{read_id}\n{seq}\n")
        print(f"Wrote {len(reads)} sequences to {filename}")

def main():
    args = get_args()
    kraken_files = [f for f in os.listdir(args.kraken_folder) if f.endswith('.txt')]
    seq_files = [f for f in os.listdir(args.seq_folder) if f.endswith(('.fastq', '.fq', '.fasta'))]
    taxids = [str(t) for t in args.taxids]

    for kraken_file in kraken_files:
        kraken_path = os.path.join(args.kraken_folder, kraken_file)
        base_prefix = get_base_prefix(kraken_file)

        matched_seq_file = None
        for seq_file in seq_files:
            if get_base_prefix(seq_file) == base_prefix:
                matched_seq_file = seq_file
                break

        if not matched_seq_file:
            print(f"Warning: No matching sequence file found for {kraken_file}, skipping.")
            continue

        seq_path = os.path.join(args.seq_folder, matched_seq_file)
        print(f"Processing pair: Kraken2='{kraken_path}', Sequence='{seq_path}'")

        read_kmers = parse_kraken_output(kraken_path, taxids)
        file_format = detect_format(seq_path)
        if file_format not in ('fasta', 'fastq'):
            print(f"Warning: Unsupported file format for {seq_path}, skipping.")
            continue
        print(f"Format detected: {file_format.upper()}")

        with open(seq_path) as f:
            if file_format == 'fasta':
                read_iter = extract_reads_from_fasta(f, read_kmers, args.kmer_length)
            else:
                read_iter = extract_reads_from_fastq(f, read_kmers, args.kmer_length)

            output_dir = '.'
            write_reads_to_files(read_iter, base_prefix, output_dir)

if __name__ == '__main__':
    main()

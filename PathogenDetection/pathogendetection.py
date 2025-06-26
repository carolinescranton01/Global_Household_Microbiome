#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Caroline Scranton
Purpose: Pull read IDs matching specific taxID from Kraken2 output and corresponding FASTQ or FASTA files in folders,
         with permanent filtering requiring >=3 k-mers per taxid,
         and flexible matching of H and L numbers anywhere in filenames.
"""

import argparse
import os
import re
from collections import defaultdict

def get_args():
    """Parameters for script"""
    parser = argparse.ArgumentParser(description='Extract read IDs matching taxID from Kraken2 outputs and FASTQ/A files in folders')
    parser.add_argument('kraken_folder', help='Folder containing Kraken2 output files')
    parser.add_argument('seq_folder', help='Folder containing corresponding FASTQ or FASTA files')
    parser.add_argument('-t', '--taxids', nargs='+', required=True, help='NCBI TaxIDs to extract')
    parser.add_argument('-k', '--kmer_length', type=int, default=35, help='K-mer length (default 35')
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
    """Read the Kraken2 output to identify kmer locations for specified taxIDs"""
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

def extract_hl_numbers(filename):
    """Extract H and L numbers from anywhere in the filename."""
    h_match = re.search(r'H\d+', filename)
    l_match = re.search(r'L\d+', filename)
    h = h_match.group(0) if h_match else None
    l = l_match.group(0) if l_match else None
    return (h, l)

def extract_reads_from_fasta(f, read_tax_kmers):
    """Pull read IDs from a FASTA file"""
    header = ''
    for line in f:
        if line.startswith('>'):
            if header:
                read_id = normalize_read_id(header)
                yield from extract_kmers_for_read(read_id, read_tax_kmers)
            header = line
        else:
            continue  # sequence not used
    if header:
        read_id = normalize_read_id(header)
        yield from extract_kmers_for_read(read_id, read_tax_kmers)

def extract_reads_from_fastq(f, read_tax_kmers):
    """Pull read IDs from a FASTQ file"""
    while True:
        header = f.readline()
        if not header:
            break
        f.readline()       # +
        f.readline()       # quality
        read_id = normalize_read_id(header)
        yield from extract_kmers_for_read(read_id, read_tax_kmers)

def extract_kmers_for_read(read_id, read_tax_kmers):
    """Yield taxid and read_id only if read has >= 3 kmers for that taxid"""
    if read_id not in read_tax_kmers:
        return
    for taxid, starts in read_tax_kmers[read_id].items():
        if len(starts) < 3:
            continue
        yield taxid, read_id

def write_reads_to_files(read_blocks, base_prefix, output_dir):
    """Writes output files with only read IDs, grouped by taxID, if there are 3 or more reads."""
    outputs = defaultdict(list)  # Use list to collect reads for each taxID
    for taxid, read_id in read_blocks:
        outputs[taxid].append(read_id)
    
    for taxid, read_ids in outputs.items():
        if len(read_ids) < 3:  # Check if there are fewer than 3 reads
            print(f"Skipping output for taxID '{taxid}' due to insufficient reads (found {len(read_ids)}).")
            continue  # Skip writing this taxID's output if fewer than 3 reads
        filename = os.path.join(output_dir, f"{base_prefix}_{taxid}.txt")
        with open(filename, 'w') as f:
            for read_id in sorted(read_ids):
                f.write(f"{read_id}\n")
        print(f"Wrote {len(read_ids)} read IDs to {filename}")

def main():
    args = get_args()
    kraken_files = [f for f in os.listdir(args.kraken_folder) 
                    if os.path.isfile(os.path.join(args.kraken_folder, f)) and f.endswith('.txt')]
    seq_files = [f for f in os.listdir(args.seq_folder) 
                 if os.path.isfile(os.path.join(args.seq_folder, f)) and (f.endswith('.fastq') or f.endswith('.fasta'))]
    taxids = [str(t) for t in args.taxids]

    for kraken_file in kraken_files:
        base_name = os.path.splitext(kraken_file)[0]
        kraken_path = os.path.join(args.kraken_folder, kraken_file)

        # Extract H and L numbers anywhere in the Kraken filename
        h_l_kraken = extract_hl_numbers(base_name)

        matched_seq_file = None
        for seq_file in seq_files:
            h_l_seq = extract_hl_numbers(seq_file)
            if h_l_kraken == h_l_seq and None not in h_l_kraken:
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
                read_iter = extract_reads_from_fasta(f, read_kmers)
            else:
                read_iter = extract_reads_from_fastq(f, read_kmers)

            output_dir = '.'  # or specify an output folder
            write_reads_to_files(read_iter, base_name, output_dir)

if __name__ == '__main__':
    main()


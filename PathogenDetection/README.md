# Pathogen Detection from Kraken2 outputs

This folder contains a python script used to identify pathogens (by NCBI taxonomic IDs) from kraken2 outputs, and pull the read IDs from the raw sequence data. 

The script processes the input in the following way:
1. The sequence file type is detected (fastq or fasta)
2. The read ID is cleaned and extracted from each sequence
3. The kraken2 output is parsed - kraken2 outputs list info about the read followed by a long list of NCBI taxIDs:kmer positions. This step looks for NCBI taxIDs which match those specified in flag --taxids
4. If 3 or more kmers matching a specified taxID are located within a read, the readID is exported into a text file titled with the sequence/kraken matched name and the taxonomic ID, stored in the folder where the python script is ran from. 


To run the script, download the pathogendetection.py script and place it in a directory containing 2 folders - one with raw reads (in fastq or fasta format) and one with kraken2 **output** files. The raw reads and kraken2 outputs must have similar names - in this analysis, the names were H1_L1, H1_L2, ... H2_L4, etc (H for house number (1-10) and L for location (also 1-10)). This is specified in the code. To run this on samples with a different naming scheme, this segment of code will need to be edited:

```
def extract_hl_numbers(filename):
    """Extract H and L numbers from anywhere in the filename."""
    h_match = re.search(r'H\d+', filename)
    l_match = re.search(r'L\d+', filename)
    h = h_match.group(0) if h_match else None
    l = l_match.group(0) if l_match else None
    return (h, l)
```

Run the script using python3 with the arguments kraken_folder (required, change name to match your folder), seq_folder (required, change name to match your folder), --taxids (list of NCBI taxIDs to look for), and optional -k (kmer length if kmers are not 35 bases (Kraken2 defaults to 35-mers)). 

```
python3 pathogendetection.py kraken_folder seq_folder --taxids 123 456 789 1011 1213
```

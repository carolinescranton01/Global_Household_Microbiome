# Taxonomic, Antibiotic Resistance Gene, and Virulence Factor Identification from single-end reads on the UA HPC
The attached bash script (ARG_VF_and_Kraken_bashscript_CS.sh) can be used to run through all commands to analyze the ARGs (using DeepARG and abricate's CARD, Resfinder, and NCBI databases) and VFs (using abricate's VFDB) on the HPC, as well as assign taxonomy to the samples via kraken2 and generating a kraken-biom file. 

Upload the script to the HPC (alternatively, create a new bash script using the following command and copy and paste the contents of the script here on github into your script) into a folder titled with the project name. Within the folder, you should have a directory called raw_reads which contains all of the raw sequence data as .fastq files. This can be organized by sample name - so for example, if you had ten samples, there would be ten folders, labelled sample1 - sample10. Within these folders, there would be several .fastq files from the sequencer (pre-sorted by barcode, so all sample folders contain all of the sequence data for that sample from all sequencing runs). The first step in the script to create folders for all outputs, and the second step is to concatenate all of the sequencing files into one, and then rename it with the folder name. So if the sample1 folder has 200 .fastq files inside, the first loop in the script will concatenate those 200 files into one .fastq file called sample1, and copy that out into a separate directory called fastq_all which will only have the concatenated .fastq files within. If you have questions about this organization, email me (carolinescranton@arizona.edu). 

Additionally, if you are planning on analyzing your kraken2 taxonomic data in R, you will need to upload a metadata file into the kraken2 folder which is created when you run the script in order to make a .biom file (a type of file containing sequence information and metadata which can be analyzed in R). Metadata needs to be in a tab-separated values (.tsv, make sure there are no commas!) spreadsheet, and needs to have one column which corresponds to the sequencing sample IDs. An example metadata sheet is attached to this github repository, with the first column 'sampleid' matching the example sample names mentioned in the previous paragraph (sample1, sample2, etc) and then additional columns with types of metadata which may be important to include about the samples, like where they are from, what date they were collected, etc. The metadata will be specific to your project, and should be as detailed as possible - any way which you may want to group the data in future analysis should be included in the metadata. The file should be saved as metadata.tsv and uploaded to the KRAKEN2 folder after the script runs. See the end of this page for information on kraken-biom.

First, create the script either by downloading it from this repository and uploading it to your HPC folder, or by manually creating a file, copying and pasting the text in the script into the file, and saving it by using the following commands: 
```
nano analysis.sh
# copy and paste the contents into the GNU editor
# press control X, control Y, enter to save
```

Next, install the following packages to run the different ARG/VF detection pipelines and kraken. You may have some of these installed already (ie. Kraken2 in the Metagenomics environment), which you can try to install them again to check (which will just update them) or skip them. Documentation for installation is below:

**Deeparg**
**NOTE: as of 10-8-25, DeepARG install does not work. I am trying to figure it out but the database is not on the HPC and the command in the script will not run**
To run DeepARG, you should create a separate environment since it requires python 2.7 which is outdated. Then, you will install all of the packages as follows using one line of code which specifies which versions

```
conda create -n DeepARG python=2.7.18
conda activate DeepARG
conda install -c bioconda diamond==0.9.24 trimmomatic vsearch bedtools==2.29.2 bowtie2==2.3.5.1 samtools
conda deactivate
```

**Abricate**

To run the four abricate databases (NCBI, CARD, ResFinder, and VFDB), enter the following code. Abricate also has other databases which are not used in this tutorial, so check out the github page (https://github.com/tseemann/abricate) for more potential applications.

```
conda activate Antibiotics
conda install -c conda-forge -c bioconda -c defaults abricate
conda deactivate
```

**Kraken2**

To run kraken2 and to generate a kraken_biom file, install kraken2 and kraken_biom into the metagenomics environment

```
conda activate Metagenomics
conda install bioconda::kraken2
conda install bioconda::kraken-biom
conda deactivate
```

To run the script type the following (replace analysis with the script name, either ARG-VF-Kraken_script_CS.sh or whatever you named or renamed it as):

```
bash analysis.sh
```

If you are planning on analyzing the taxonomic data in R, you need to upload your metadata to the KRAKEN2 folder and run the following commands (make sure the metadata file is called metadata.tsv and that you have the Metagenomics environment active, which it should be after running the script):

```
cd KRAKEN2
kraken-biom sample_1_report.txt sample_2_report.txt sample_3_report.txt sample_4_report.txt -m metadata.tsv -o output.biom
```
Sample names must be listed out and match the sample_ID column (first column) of the metadata file (ie. sample_1_report.txt matches sample_1_report in the metadata file).
The output.biom file can be exported off the HPC and uploaded into R for analysis. 

After this script was run, all data was exported off the HPC in folders and compiled into excel spreadsheets for further analysis in R using the writeexcel.py and sortexcel.py python scripts, found in this github repository. 

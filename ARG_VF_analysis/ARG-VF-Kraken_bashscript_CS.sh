
mkdir fastq_all
mkdir KRAKEN2
mkdir DeepARG
mkdir VFDB
mkdir NCBI
mkdir CARD
mkdir resfinder

cd raw_reads
for d in */
do
cd "$d"
gunzip *.gz
cat *.fastq > "$(pwd).fastq"
cp "$(pwd).fastq" ../../fastq_all
cd ..
done

conda activate Antibiotics

cd ../fastq_all
for f in *.fastq
do
n=${f%%.fastq}
abricate --db vfdb --threads 94 --mincov 15 ${n}.fastq > ${n}_VFDB.tab
abricate --db resfinder --threads 94 --mincov 15 ${n}.fastq > ${n}_resfinder.tab
abricate --db ncbi --threads 94 --mincov 15 ${n}.fastq > ${n}_NCBI.tab
abricate --db card --threads 94 --mincov 15 ${n}.fastq > ${n}_CARD.tab
cp *CARD.tab ../CARD
cp *NCBI.tab ../NCBI
cp *VFDB.tab ../VFDB
cp *resfinder.tab ../resfinder
done

conda deactivate
conda activate DeepARG

for f in *.fastq
do
n=${f%%.fastq}
deeparg predict --model SS --type nucl --input ${n}.fastq --out ${n}_deeparg -d /xdisk/kcooper/carolinescranton/deeparg_database
cp *_deeparg.mapping.ARG ../DeepARG
done

conda deactivate
conda activate Metagenomics

for f in *.fastq
do
n=${f%%.fastq}
kraken2 --db /xdisk/kcooper/kcooper/Kraken_Special_DB --report ${n}_K2_report.txt --output ${n}_K2_output.txt ${n}.fastq
cp *.txt ../KRAKEN2
done




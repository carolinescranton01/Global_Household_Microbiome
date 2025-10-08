
mkdir KRAKEN2
mkdir DeepARG
mkdir VFDB
mkdir NCBI
mkdir CARD
mkdir resfinder

conda activate Antibiotics

for f in *.fastq
do
n=${f%%.fastq}
abricate --db vfdb --threads 94 --mincov 15 ${n}.fastq > ${n}_VFDB.tab
abricate --db resfinder --threads 94 --mincov 15 ${n}.fastq > ${n}_resfinder.tab
abricate --db ncbi --threads 94 --mincov 15 ${n}.fastq > ${n}_NCBI.tab
abricate --db card --threads 94 --mincov 15 ${n}.fastq > ${n}_CARD.tab
mv *CARD.tab CARD
mv *NCBI.tab NCBI
mv *VFDB.tab VFDB
mv *resfinder.tab resfinder
done

conda deactivate
conda activate DeepARG

for f in *.fastq
do
n=${f%%.fastq}
deeparg predict --model SS --type nucl --input ${n}.fastq --out ${n}_deeparg -d /xdisk/kcooper/carolinescranton/deeparg_database
mv *_deeparg* DeepARG
done

conda deactivate
conda activate Metagenomics

for f in *.fastq
do
n=${f%%.fastq}
kraken2 --db /xdisk/kcooper/kcooper/Kraken_Special_DB --report ${n}_K2_report.txt --output ${n}_K2_output.txt ${n}.fastq
mv *.txt KRAKEN2
done




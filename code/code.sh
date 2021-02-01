#!/usr/bin/env bash
# Códigos utilizados para responder as questões do arquivo QUESTION.txt
# Processo Seletivo Genomika/Einstein
# Janeiro 2021
# Autor: Antonio Victor Campos Coelho

# Indo da pasta code para a bioinfotest
cd ..

REF=data/hg19.fasta

# Questão 1

# Salvei as raízes dos nomes dos arquivos num arquivo de texto
for file in data/*.fastq.gz; do
  basename $file .fastq.gz >> data/reads.txt
done

# Usei o arquivo de texto criado na etapa anterior para resgatar o nome das reads e usá-los em um loop
while read line; do
  fastqc "data/${line}.fastq.gz" && unzip -o "data/${line}_fastqc.zip" -d data/
   
  cat "data/${line}_fastqc/fastqc_data.txt" | grep "Total Sequences" >> data/summaries.txt
done < data/reads.txt

cat data/summaries.txt

# Questão 2

grep ">" ${REF} | wc -l

# Questão 3

readarray -t READS < data/reads.txt # a partir do BASH versão 4.0

bwa mem ${REF} "data/${READS[0]}.fastq.gz" "data/${READS[1]}.fastq.gz" > output/alignment.sam

samtools sort output/alignment.sam -o output/alignment_sorted.bam
samtools index output/alignment_sorted.bam

echo -e "chr17\t41197694\t41197819" > data/regiao.bed

samtools view -L data/regiao.bed -c output/alignment_sorted.bam

# Questão 4

samtools view -c -f 4 output/alignment_sorted.bam

# Questão 5.1

freebayes -f ${REF} -t data/BRCA.list output/alignment_sorted.bam > output/variants.vcf

snpEff -v -stats output/variants.html hg19 output/variants.vcf > output/variants_annotated.vcf

bcftools query -f '%CHROM %POS %INFO/TYPE %REF %ALT[ %GT] %INFO/ANN\n' output/variants_annotated.vcf > output/genotypes.txt

echo -e "variants with HIGH impact" > output/questao51.txt

grep "HIGH" output/genotypes.txt | cut -d' ' -f 1,2,3 >> output/questao51.txt

echo -e "variants with MODERATE impact" >> output/questao51.txt

grep "MODERATE" output/genotypes.txt | cut -d' ' -f 1,2,3 >> output/questao51.txt

echo -e "variants with LOW impact" >> output/questao51.txt

grep "LOW" output/genotypes.txt | cut -d' ' -f 1,2,3 >> output/questao51.txt

cat output/questao51.txt

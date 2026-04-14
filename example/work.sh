#genomic data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/374/765/GCA_030374765.1_ASM3037476v1/GCA_030374765.1_ASM3037476v1_genomic.fna.gz

#WGS raw reads
parallel-fastq-dump --sra-id SRR23519128 --threads 20 --outdir SRR23519128_reads --split-files --gzip

#map the raw reads to the genome
bwa index Gracilaria_chilensis.genome.fa
bwa mem -t 12 Gracilaria_chilensis.genome.fa Gracilaria_chilensis_WGS_1.fq.gz Gracilaria_chilensis_WGS_2.fq.gz > Gchilensis.aln.sam  
samtools view -Sb --threads 12 -o Gchilensis.aln.bam Gchilensis.aln.sam  
samtools sort --threads 12 -o Gchilensis.aln.sorted.bam Gchilensis.aln.bam  
samtools index Gchilensis.aln.sorted.bam

#run pandepth
pandepth -i Gchilensis.aln.sorted.bam -w 500 -o depth

#run the python script
python gc-depth-plot.py Gracilaria_chilensis.genome.fa depth.win.stat.gz -w 500

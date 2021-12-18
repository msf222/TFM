# TFM

# Descarga de los datos 
### Creando una carpeta
```
mkdir Datos
```
### Descargando los datos
```
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2163030/SRR2163030.1
```
# Transformacion a fastaq

## Instalacion de la herramienta sratoolkit
### Creando una carpeta
```
mkdir sra-tools
```
### Descargando sratoolkit
```
curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz -o sratoolkit.2.9.6-ubuntu64.tar.gz
tar -xvzf sratoolkit.2.9.6-ubuntu64.tar.gz
export PATH=./sra-tools/sratoolkit.2.9.6-ubuntu64/bin:$PATH
```
### Usando fastq-dump para transformacion a fastaq
```
fastq-dump -gzip -split-files ./SRR2163030.1
```

# Herramienta seqQscorer (Control de calidad de los datos)
### Ref. https://github.com/salbrec/seqQscorer 

## Instalacion de docker
```
sudo apt-get install docker
sudo apt-get install docker.io
sudo docker login
sudo docker pull salbrec/seqqdocker
```

## Corriendo la imagen de docker
```
sudo docker run -i -t -v "/home/:/home/" salbrec/seqqdocker /bin/bash
```
## Descargando el indice genomico de referencia (humano en nuestro caso)
```
cd utils/genome_index
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip     # descarga
apt install unzip
unzip GRCh38_noalt_as.zip                                          # descompresion
```
## Corriendo la herramienta con los datos de origen
### Generando ficheros y fastqc
```
python deriveFeatureSets.py --fastq1 ./SRR2163030.1_1.fastq.gz --fastq2 ./SRR2163030.1_2.fastq.gz --cores 10 --btidx /var/seqQscorer/utils/genome_index/GRCh38_noalt_as/GRCh38_noalt_as --assembly GRCh38 --outdir ./Results_seQ/
```
### Generando estadistico de calidad
```
python seqQscorer.py --indir ./Results_seQ/ --species human --runtype paired-end --bestCalib --probOut ./the_probability.tsv --compOut ./comprehensive_output.txt --seed 42 --sampleID SRR2163030.1
```

# Cutadapt (Eliminacion de adaptadores y secuencias sobrerrepresentadas)
### Ref. https://cutadapt.readthedocs.io/en/stable/index.html

## Instalación
```
sudo apt install cutadapt
```
## Eliminacion del adaptador Nextera y el cebador sobrerrepresentado.
```
cutadapt -a CTGTCTCTTATA -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC -A CTGTCTCTTATA -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATG -j 8 -o cut_out.1.fastq.gz -p cut_out.2.fastq.gz ./SRR2163030.1_1.fastq.gz ./SRR2163030.1_2.fastq.gz 
```
# Herramienta seqQscorer (paso 2)

## Corriendo la herramienta con los datos tras la transformacion
### Generando ficheros y fastqc
```
python deriveFeatureSets.py --fastq1 ./cut_out.1.fastq.gz --fastq2 ./cut_out.2.fastq.gz --cores 10 --btidx ./utils/genome_index/GRCh38/GRCh38_noalt_as --assembly GRCh38 --outdir ./Results_seQ_2/
```
### Generando estadistico de calidad
```
python seqQscorer.py --indir ./Results_seQ_2/ --species human --runtype paired-end --bestCalib --probOut ./the_probability.tsv --compOut ./comprehensive_output.txt --seed 42 --sampleID cut_out.1.fastq
```

# Flash (Union de los reads Forward y Reverse)
### Ref. https://github.com/genome-vendor/FLASH/blob/master/MANUAL

## Instalación
```
wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-Linux-x86_64.tar.gz  # descarga
tar -xvf FLASH-1.2.11-Linux-x86_64.tar.gz                                # descompresion
```
## Merge
```
./flash ./cut_out.1.fastq.gz ./cut_out.2.fastq.gz -t 8 -z -M 101
```
# BwaMeme (Mapeo de referencia)
### Ref. https://github.com/kaist-ina/BWA-MEME

## Instalacion
```
git clone https://github.com/kaist-ina/BWA-MEME.git BWA-MEME
cd BWA-MEME
make clean
make -j<num_threads> arch=avx2
wget https://rustup.rs/
```
## Descarga genoma de referencia

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz (Chr21)

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/human_g1k_v37.fasta.gz 
gunzip human_g1k_v37.fasta.gz

## Construyendo el indice

./bwa-mem2 index -p in_21 chr21.fa.gz (Chr21)

./bwa-mem2 index -a meme -t 8 human_g1k_v37.fasta
./build_rmis_dna.sh human_g1k_v37.fasta

## Descargando ficheros necesarios

wget https://ina.kaist.ac.kr/~bwa-meme/human_g1k_v37.fasta.suffixarray_uint64_L1_PARAMETERS
wget https://ina.kaist.ac.kr/~bwa-meme/human_g1k_v37.fasta.suffixarray_uint64_L2_PARAMETERS

## Alineamiento

./bwa-mem2 mem -t 8 in_21 cut_out.1.fastq.gz > out_21.sam (Chr21)

./bwa-mem2 mem -Y -K 100000000 -t 8 -7 /home/martasf22/BWA-MEME/human_g1k_v37.fasta /home/martasf22/cut_out.1.fastq.gz -o output_meme.sam

# Samtools (Identificacion de pares duplicados)
Ref. http://www.htslib.org/

# Instalacion Samtools
```
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-1.3.1
make
sudo make prefix=/usr/local/bin install
```
## Ordenando el archivo y anadiendo marcadores
### Agrupando las lecturas del mismo nombre en grupos contiguos
```
samtools collate -l 5 ./BWA-MEME/out_21.sam.gz out_collate
```
### Completando las coordenadas de relación a la posición y agregando etiquetas ms que usara 'markdup' para seleccionar las mejores lecturas que seran guardadas.
```
samtools fixmate -m namecollate.bam fixmate.bam 
```
### Ordenando las alineaciones por coordenadas.
```
samtools sort -o positionsort.bam fixmate.bam
```
### marcando alineaciones duplicadas.
```
samtools markdup positionsort.bam markdup.bam
```


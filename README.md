# TFM

# Descarga de los datos 

## Descargando los datos
```
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2163030/SRR2163030.1
```
# Transformacion a fastaq

## Instalacion de la herramienta sratoolkit

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

# Herramienta Flash (Union de los reads Forward y Reverse)
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
# Herramienta BwaMeme (Mapeo de referencia)
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

## Construyendo el indice

./bwa-mem2 index -p in_21 chr21.fa.gz (Chr21)

## Alineamiento

./bwa-mem2 mem -t 8 in_21 cut_out.1.fastq.gz > out_21.sam (Chr21)

# Herramienta Samtools (Identificacion de pares duplicados)
### Ref. http://www.htslib.org/

## Instalacion Samtools
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
samtools fixmate -m out_collate.bam fixmate.bam 
```
### Ordenando las alineaciones por coordenadas.
```
samtools sort -o positionsort.bam fixmate.bam
```
### Marcando alineaciones duplicadas.
```
samtools markdup positionsort.bam markdup.bam
```
# Herramienta Deep-variant (Llamada a variantes)
### Ref. https://github.com/google/deepvariant
## Instalacion
```
sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:1.3.0
```
## Creando directorios de entrada y salida
```
OUTPUT_DIR="./quickstart-output"
mkdir -p "./quickstart-output"
INPUT_DIR="./quickstart-testdata"
mkdir -p "./quickstart-testdata"
```
## Generando archivo VCF
```
sudo docker run -v "./quickstart-testdata/":"/input" -v "./quickstart-output/":"/output" google/deepvariant:"1.3.0" /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/input/chr21.fa --reads=/input/markdup.bam --output_vcf=/output/output.vcf.gz --output_gvcf=/output/output.g.vcf.gz --num_shards=8
```

# Herramienta ForestQC (Refinamiento de llamada)
### Ref. https://github.com/avallonking/ForestQC
## Instalacion
```
conda install forestqc -c avallonking
conda update forestqc -c avallonking
```
## Dividiendo variantes por calidad
### Estableciendo umbrales para Outlier_GQ y Outlier_DP.
```
ForestQC set_outlier -m 1G -i ./quickstart-output/output.vcf.gz
```
### Generando tabla de contenido GC para nuestro genoma de referencia
```
ForestQC compute_gc -i ./chr21.fa -o GC_content
```
### Calcular las estadisticas del archivo VCF 
```
ForestQC stat -i ./quickstart-output/output.vcf.gz -o statics_Forest -c ./ForestQC/GC_content --gq 20 --dp 13
```
### Dividiendo el dataset por calidad
```
ForestQC split -i /./ForestQC/statics_Forest
```
###  Marcando posiciones de las variantes no clasificadas y filtrandolas
```
awk -F "\t" 'NR>1{print $2"\t"$3}' gray.statics_Forest > bad_variants_positions
vcftools --gzvcf output.vcf.gz --exclude-positions bad_variants_positions --recode --recode-INFO-all -c | gzip -c > output_vcf_good_variants
```
### Filtrando variantes no PASS
```
vcftools --gzvcf ./ForestQC/output_vcf_good_variants.vcf.gz --remove-filtered-all --recode --recode-INFO-all -c | gzip -c > output_vcf_filter_variants.vcf.gz
```
# Herramienta snpEff (Anotacion de variantes)
### Ref. http://pcingola.github.io/SnpEff/se_running/
## Instalacion
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```
## Descargando anotaciones del genoma hg38
```
java -jar snpEff.jar download -v hg38
```
## Anotando variantes
```
java -Xmx8g -jar snpEff.jar hg38 /./output_vcf_filter_variants.vcf.gz > test.chr21_filter.ann.vcf
```
# Código desarrollo de aplicación Shiny para la visualización y filtrado de variantes.
```
library(shiny)
library(ggplot2)
library(vcfR)
library(stringr)

input <- read.vcfR(
  "vcf.vcf",
  limit = 1e+07,
  nrows = -1,
  skip = 0,
  cols = NULL,
  convertNA = TRUE,
  checkFile = TRUE,
  check_keys = TRUE,
  verbose = TRUE
)

table <- (input@fix)[,-8]
table <- as.data.frame(table)
table$IMPACT <- sapply(str_split(input@fix[,8], "\\|"), function (x) (x[3]))
table$TYPE <- sapply(str_split(input@fix[,8], "\\|"), function (x) (x[6]))
table$ANOTATION <- sapply(str_split(input@fix[,8], "\\|"), function (x) (x[2]))
table$GEN <- sapply(str_split(input@fix[,8], "\\|"), function (x) (x[4]))
table <- table[,-3]



ui <- fluidPage(
  titlePanel("VisualizaciÃ³n y filtrado de variantes"),
  
  # Create a new Row in the UI for selectInputs
  fluidRow(
    column(4,
           selectInput("imp",
                       "Impacto:",
                       c("All",
                         unique(as.character(table$IMPACT)))
    ),
  ),
  
  column(4,
         selectInput("tip",
                     "AnotaciÃ³n:",
                     c("All",
                       unique(as.character(table$ANOTATION)))
         ),
  ),
  
  
  # Create a new row for the table.
  DT::dataTableOutput("table")
)
)


server <- function(input, output) {
  
  output$table <- DT::renderDataTable(DT::datatable({
    data <- table
    if (input$imp != "All") {
      data <- data[data$IMPACT == input$imp,]
    }
    if (input$tip != "All") {
      data <- data[data$ANOTATION == input$tip,]
    }
    data
  }))
  
}

shinyApp(ui = ui, server = server)
```

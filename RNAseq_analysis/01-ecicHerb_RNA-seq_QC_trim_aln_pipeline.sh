#!/usr/bin/env bash

#### SETTING THE PROJECT ####
# Define the path of the project
projectDir="/path/to/work/directory" ## SET CUSTOM PATH TO WORK DIRECTORY
cd $projectDir

# Generate all subfolders
mkdir -p $projectDir/src $projectDir/data/genome $projectDir/data/raw-seq \
$projectDir/results/fastqc-rawseq $projectDir/results/fastqc-trimseq  \
$projectDir/results/htseq-count  $projectDir/results/star \
$projectDir/results/trimmomatic $projectDir/results/featureCounts \
$projectDir/results/deseq2 $projectDir/results/topGO \


# Download the SRR_Acc_List.txt from bioproject PRJNA1129331 at SRA, 
# then download the accessions SRR29751793 to SRR29751828, 
# manually change their names to informative ones
cat SRR_Acc_List.txt | while read acc; do fasterq-dump --split-3 "$acc"; done

# Generate a list of file names from the project
ls -d $projectDir/data/raw-seq/*.fastq.gz | sed -e 's/\_[12].fastq.gz//g' -e "s#$projectDir/data/raw-seq/##g" | sort | uniq > $projectDir/src/ecic_basenames.txt

#### RUN QUALITY CHECK IN RAW SEQUENCES WITH FASTQC ####
# Run FastQC for all the project files
cat $projectDir/src/ecic_basenames.txt | while read line; do
	mkdir $projectDir/results/fastqc-rawseq/$line\_1.fastqc
	fastqc -o $projectDir/results/fastqc-rawseq/$line\_1.fastqc \
	$projectDir/data/raw-seq/$line\_1.fastq.gz;
	mkdir $projectDir/results/fastqc-rawseq/$line\_2.fastqc
	fastqc -o $projectDir/results/fastqc-rawseq/$line\_2.fastqc \
	$projectDir/data/raw-seq/$line\_2.fastq.gz;
done

# Generate a HTML quality chek report for raw reads with MultiQC
cd $projectDir/results/fastqc-rawseq
multiqc $projectDir/results/fastqc-rawseq/*/*.zip
cd $projectDir


#### FILTER READS AND NUCLEOTIDES BY QUALITY WITH TRIMMOMATIC ####
# Define a path for fasta file with Illumina adapters 
# Different adapters are used depending on type of sequencing machine
trimm="/home/usuario/Windows/Analyses/EcicEcaz_2022/Ecic_rnaseq_herbivory/src/adapters.fa" # Illumina TruSeq adapters
# Run Trimmomatic for all the fastq files
cat $projectDir/src/ecic_basenames.txt | while read line; do
  TrimmomaticPE \
  -phred33 \
  $projectDir/data/raw-seq/$line\_1.fastq.gz \
  $projectDir/data/raw-seq/$line\_2.fastq.gz \
  $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_1.trimmed.unpaired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_2.trimmed.unpaired.fastq.gz \
  ILLUMINACLIP:$trimm:2:30:10 \
  LEADING:28 TRAILING:28 \
  SLIDINGWINDOW:4:15 MINLEN:30;
  mkdir $projectDir/results/fastqc-trimseq/$line\_1.trimmed.paired.fastqc;
  fastqc -o $projectDir/results/fastqc-trimseq/$line\_1.trimmed.paired.fastqc \
  $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz;
  mkdir $projectDir/results/fastqc-trimseq/$line\_2.trimmed.paired.fastqc;
  fastqc -o $projectDir/results/fastqc-trimseq/$line\_2.trimmed.paired.fastqc \
  $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz;
done

# Generate a HTML quality chek report for trimmed reads with MultiQC
cd $projectDir/results/fastqc-trimseq
multiqc $projectDir/results/fastqc-trimseq/*/*.zip
cd $projectDir


#### EDIT GENOME FILES TO FIT THE PIPELINE ####
ln -s /home/usuario/Windows/Data/genomes/Ecic_genome_pacBio/MaSuRCA_assembly.fasta.PolcaCorrected.fa $projectDir/data/genome/ecic_genome.fa
ln -s /home/usuario/Windows/Data/genomes/Ecic_genome_pacBio/Gene_prediction_files/AUGUSTUS/CicuSeq2_augustus.gff3 $projectDir/data/genome/ecic_annotation.gff3
grep -v '#' ecic_annotation.gff3 > ecic_annotation.tmp.gff3


#### CONVERT GFF FILES INTO GTF FILES ####
# Use the Rscript gff2gtf_rtracklayer.R to convert the GFF into GTF
# NEED TO USE R version >= 3.5.0
/usr/local/bin/./R --slave --no-restore \
 --file=$projectDir/src/gff2gtf_rtracklayer.R \
 $projectDir/data/genome/ecic_annotation.tmp.gff3 \
 $projectDir/data/genome/ecic_annotation.gtf putEGin_gene_id
# Command above does not run.
# Used the modified Rscript gff2gtf_rtracklayer_windows.R in a Windows OS with Rstudio, copying and pasting the commands.
cd $projectDir


#### INDEX THE GENOME WITH STAR ####
mkdir -p $projectDir/data/genome/STAR-2.4.2a-R64_Index
STAR --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Index \
 --genomeFastaFiles $projectDir/data/genome/ecic_genome.fa \
 --sjdbGTFfile $projectDir/data/genome/ecic_annotation.gff3 \
 --sjdbOverhang 149


#### ALIGN THE READS TO THE REFERENCE GENOME WITH STAR ####
cat $projectDir/src/ecic_basenames.txt | while read line; do
mkdir $projectDir/results/star/$line.STAR;
STAR --runThreadN 12 \
 --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Index \
 --readFilesIn $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
 --readFilesCommand gunzip -c \
 --sjdbGTFfile $projectDir/data/genome/ecic_annotation.gff3 \
 --sjdbOverhang 149 \
 --limitGenomeGenerateRAM 5900000000 \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $projectDir/results/star/$line.STAR/$line.STAR \
 --quantMode GeneCounts \
 --twopassMode Basic;
 # Do not let STAR use more than 5.9 GB of memory
 # Output your alignment file in BAM format, sorted by coordinate
 # Give your output files a UNIQUE name/prefix (${line} also works)
 # Gene counting for generating gene count tables
 # Two-pass mode for refining splice junctions
done

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/star/
multiqc $projectDir/results/star/
cd $projectDir


#### GET MAPPED READ COUNTS WITH FEATURECOUNTS ####
mkdir -p $projectDir/results/featureCounts 
cat $projectDir/src/ecic_basenames.txt | while read line; do
mkdir $projectDir/results/featureCounts/$line.Counts;
featureCounts -T 1 -s 0 -g gene_id -t exon \
 -o $projectDir/results/featureCounts/$line.featCounts.txt \
 -a $projectDir/data/genome/ecic_annotation.gff3 \
 $projectDir/results/star/$line.STAR/$line.STARAligned.sortedByCoord.out.bam;
done

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/featureCounts/
multiqc $projectDir/results/featureCounts/
cd $projectDir

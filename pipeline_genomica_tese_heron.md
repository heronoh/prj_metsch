Yeast genome analysis pipeline\* used in *Heron Hilário* PhD Dissertation
================
Heron O. Hilário
15/03/2020

# Assembly
## reads quality assessment, trimming and genome assembly
### FastQC for visualization of reads quality

*M. australis* WGS reads from MiSeq and HiSeq platforms
``` bash
fastqc -f fastq -o . aus_miseq_001.fastq aus_miseq_002.fastq aus_hiseq_001.fastq aus_hiseq_002.fastq ;
```
Open html output on web browser.

### TRIMMOMATIC for removal of adpter and poor quality reads

Leaving only reads with phred score above 30.

#### MiSeq reads
``` bash
java -jar trimmomatic-0.38.jar PE -phred64 -threads 16 \
 aus_miseq_001.fastq aus_miseq_002.fastq \
 aus_miseq_001_P_trim.fastq aus_miseq_001_U_trim.fastq \
 aus_miseq_002_P_trim.fastq aus_miseq_002_U_trim.fastq \
 ILLUMINACLIP:$home/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:30 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50 ;
```
#### HiSeq reads
``` bash
java -jar trimmomatic-0.38.jar PE -phred64 -threads 16 \
 aus_hiseq_001.fastq aus_hiseq_002.fastq \
 aus_hiseq_001_P_trim.fastq aus_hiseq_001_U_trim.fastq \
 aus_hiseq_002_P_trim.fastq aus_hiseq_002_U_trim.fastq \
 ILLUMINACLIP:$home/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:30 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 ;
```
### FastQC for visualization of trimmed reads quality

``` bash
fastqc -f fastq -o . aus_hiseq_001_P_trim.fastq aus_hiseq_001_U_trim.fastq \
 aus_hiseq_002_P_trim.fastq aus_hiseq_002_U_trim.fastq \
 aus_miseq_001_P_trim.fastq aus_miseq_001_U_trim.fastq \
  aus_miseq_002_P_trim.fastq aus_miseq_002_U_trim.fastq ;
```
Open html output on web browser.

### Genome assembly with SPADES 3.9.0

``` bash
spades.py \
 -k 21,33,55,77 \
 --careful \
 -t 64 \
 --pe1-1 aus_miseq_001_P_trim.fastq \ 
 --pe1-2 aus_miseq_002_P_trim.fastq \
 --pe2-1 aus_hiseq_001_P_trim.fastq \
 --pe2-2 aus_hiseq_002_P_trim.fastq \
 --s1 aus_miseq_001_P_trim.fastq \
 --s2 aus_miseq_002_P_trim.fastq \
 --s3 aus_hiseq_001_U_trim.fastq \
 --s4 aus_hiseq_002_U_trim.fastq \
 -o ~/genome_assembly/aus/ ;
```

## Assembly quality assessment
### Scaffold and contigs statistics with scaffold_stats.pl

Contig metrics with scaffold_stats.pl by github.com/sujaikumar 
``` bash
scaffold_stats.pl -f aus_contigs.fasta -f 500 1000 -o aus_contigs_stats.txt;
```
### Genome completeness assessment with CEGMA 2.5

``` bash
cegma -g $genome.fasta -T 16 -o ~/.../cegma_out/$genome_cegma;
```

### Genome completeness assessment with BUSCO

Run BUSCO and then read the short summary file for completeness report.
``` bash
mkdir busco_out;
cd busco_out;
python ~/.../busco_v2/BUSCO.py -i $genome.fasta -l ~/.../busco/saccharomyceta_odb9/ -o $genome -m geno;
cat run_$genome/short_summary_$genome.txt
```







### Convert CEGMA gff output into HMM for maker

Conversion using SNAP scripts

``` bash
mkdir gff_cegma;
cd gff_cegma;
cegma2zff ~/.../cegma_out/$genome.cegma.gff ~/.../genomes/$genome.fasta; 
fathom genome.ann genome.dna -categorize 1000; 
fathom -export 1000 -plus uni.ann uni.dna; 
forge export.ann export.dna; 
hmm-assembler.pl ~/.../genomes/$genome.fasta . > \ 
 ~/.../gff_cegma/$genome.cegmasnap.hmm;
```


### Masking and counting repeats proprotion with RepeatMasker
``` bash
RepeatMasker -pa 16 -species fungi -s -noisy \
-dir ~/.../repeatMasker/ \
-html -small -lcambig -gccalc ~/.../genomes/$genome.fasta
```



## Including Plots

You can also embed plots, for
example:

![](pipeline_genomica_tese_heron_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

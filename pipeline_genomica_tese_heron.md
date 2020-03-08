Yeast genome analysis pipeline\* used in *Heron Hilário* PhD Dissertation
================
Heron O. Hilário
15/03/2020

# Assembly
## reads quality assessment and trimming
### FastQC for visualization of reads quality

*M. australis* WGS reads from MiSeq and HiSeq platforms
``` bash
fastqc -f fastq -o . aus_miseq_001.fastq aus_miseq_002.fastq aus_hiseq_001.fastq aus_hiseq_002.fastq
```
Open html output on web browser.

### TRIMMOMATIC for removal of adpter and poor quality reads
Leaving only reads with phred score above 30

MiSeq reads
``` bash
java -jar trimmomatic-0.38.jar PE -phred64 -threads 16 \
 aus_miseq_001.fastq aus_miseq_002.fastq \
 aus_miseq_001_P_trim.fastq aus_miseq_001_U_trim.fastq \
 aus_miseq_002_P_trim.fastq aus_miseq_002_U_trim.fastq \
 ILLUMINACLIP:$home/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:30 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50
```
HiSeq reads
``` bash
java -jar trimmomatic-0.38.jar PE -phred64 -threads 16 \
 aus_hiseq_001.fastq aus_hiseq_002.fastq \
 aus_hiseq_001_P_trim.fastq aus_hiseq_001_U_trim.fastq \
 aus_hiseq_002_P_trim.fastq aus_hiseq_002_U_trim.fastq \
 ILLUMINACLIP:$home/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:30 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```
### FastQC for visualization of trimmed reads quality

``` bash
fastqc -f fastq -o . aus_hiseq_001_P_trim.fastq aus_hiseq_001_U_trim.fastq \
 aus_hiseq_002_P_trim.fastq aus_hiseq_002_U_trim.fastq \
 aus_miseq_001_P_trim.fastq aus_miseq_001_U_trim.fastq \
  aus_miseq_002_P_trim.fastq aus_miseq_002_U_trim.fastq
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
 -o ~/genome_assembly/aus/
```

## Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for
example:

![](pipeline_genomica_tese_heron_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

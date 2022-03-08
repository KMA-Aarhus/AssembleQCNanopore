# AssembleQCNanopore
Input: 
- Raw nanopore fastq files as output from MinKnow. The directory must contain a folder titled "fastq_pass". This matches the default output format.   
- Samplesheet containing sample ids and barcodes. Barcodes must be in the format "NB01", "NB02" etc.  
Example:
![image](https://user-images.githubusercontent.com/90172976/157239868-b8989c11-0dce-4d4e-b886-5e89bc3bab1a.png)


Output:  
- QC report. 
- Taxonomic analysis. 
- Genome Assembly. 
- Genome annotation. 
## Installation. 
Requires conda, openpyxl and snakemake to run.  
https://docs.conda.io/en/latest/miniconda.html. 
https://openpyxl.readthedocs.io/en/stable/.  
https://snakemake.readthedocs.io/en/stable/. 

### Install with git
```
git clone https://github.com/KMA-Aarhus/AssembleQCNanopore.git
```
## How to run
Set up an alias for snakemake to run on slurm:
```
alias snakeslurm='mkdir -p logs/old; mv logs/*.{err,out} logs/old 2> /dev/null; snakemake --profile configs/slurm --use-conda --conda-frontend mamba'
```
Navigate to AssemblyQCNanopore directory where the snakefile is.  
The scripts does not allow overlaps in barcodes. If you have sequencing runs, the analysis must be started for each run.

### To run 
```
snakeslurm --config rundir="path_to_sequencing_folder" samplesheet="path_so_samplesheet/samplesheet"
Replace "path_to_sequencing_folder" with the path to the raw data. Replace "path_so_samplesheet/samplesheet" with the samplesheet location.

Additional optional options:
kraken2_db # Location of the kraken2 database. Default="/project/ClinicalMicrobio/faststorage/database/kraken2_20210423"
```
### Output
Output can be found in the specified output directory. Output contains:
* A multi-qc report for all samples in the run. The report contains read statistics, assembly statistics, taxonomic analysis, genome coverage and more.
* One folder per sample.
  * "sample_name"_consensus.fasta: assembled and polished genome
  * "sample_name"_consensus.gff: genome annotations
  * "sample_name"_consensus.gbk: annotated genome in genbank format
  * "sample_name"_kraken2_reads_report.txt: the output from kraken2. The same information can be found in the multi-qc report.
  * Addtional output for debugging can be found in:
    * trimmed
    * read_filtering
    * prokka
    * medaka
    * flye


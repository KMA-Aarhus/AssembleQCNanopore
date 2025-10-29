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
- Possible resistance genes
- Plasmids
## Installation. 
Requires conda, openpyxl and snakemake to run.  
https://docs.conda.io/en/latest/miniconda.html. 
https://openpyxl.readthedocs.io/en/stable/.  
https://snakemake.readthedocs.io/en/stable/. 

## DAG of workflow
![AssembleQCNanopore](https://user-images.githubusercontent.com/90172976/165507987-cd19bc60-4118-40e6-878d-4675c2476c13.png)


### Install with git
```
git clone https://github.com/KMA-Aarhus/AssembleQCNanopore.git
```
## How to run
Set up an alias for snakemake to run on slurm:
```
alias snakeslurm='mkdir -p logs/old; mv logs/*.{err,out} logs/old 2> /dev/null; snakemake --profile configs/slurm --use-conda --conda-frontend mamba --jobs 20'
```
Navigate to AssemblyQCNanopore directory where the snakefile is.  
The scripts does not allow overlaps in barcodes. If you have sequencing runs, the analysis must be started for each run.

The workflow installs required tools for each rule into environments via snakemake. YAML files generally containts specific versions to ensure identical outputs across new installs. Versions should however be updated periodically. Note that some of the environment files uses either bioconda or conda-forge due to "cannot solve" issues when snakemake is building the environments.

### To run 
```
snakeslurm --config rundir="path_to_sequencing_folder" samplesheet="path_so_samplesheet/samplesheet"
```
Replace "path_to_sequencing_folder" with the path to the raw data. Replace "path_so_samplesheet/samplesheet" with the samplesheet location.

Additional optional options: 
```
kraken2_db # Location of the kraken2 database. Default="/project/ClinicalMicrobio/faststorage/database/k2_standard_20250714"
amr_db: #Location of the amrfinder database. Default="/project/ClinicalMicrobio/faststorage/database/amrfinderplus/data/2025-07-16.1"
```
These must be updated if you are running the workflow outside of the ClinicalMicrobio projects on GenomeDK!


### Output
Output can be found in the specified output directory. Output contains:
* A multi-qc report for all samples in the run. The report contains read statistics, assembly statistics, taxonomic analysis, genome coverage and more.
* One folder per sample.
  * "sample_name"_consensus.fasta: assembled and polished genome
  * "sample_name"_consensus.gff: genome annotations
  * "sample_name"_consensus.gbk: annotated genome in genbank format
  * "sample_name"_kraken2_reads_report.txt: the output from kraken2. The same information can be found in the multi-qc report.
  * "sample_name"_amrfinder.txt: Resistance genes matching NCBI-AMR database.
  * "sample_name"_abricate_plasmid.txt: Detected plasmids.
  * Addtional output for debugging can be found in:
    * trimmed
    * read_filtering
    * prokka
    * medaka
    * flye

### Issues when running
* Snakemake version: Currently running on snakemake v7
* Issues building environments: Try building the environments from the yaml files outside the workflow.
* Logs: Logs are saved in the logs folder in the workflow location.
* Reproduce issues: Activate and run the individual commands outside snakemake.



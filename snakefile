# snakemake --profile configs/slurm


print("/*")
__author__ = "Tine Ebsen" # Please add your name here if you make changes.
__version__ = "0.1"

import sys
import os
from os import listdir
from os.path import isfile, isdir, join
import yaml
import pandas as pd
import numpy as np
from datetime import datetime
import glob
import time
import atexit
import datetime

# Actually given on the command line:
samplesheet = config["samplesheet"]
rundir = config["rundir"]
kraken2_db = config["kraken2_db"]


tab = "\t"
nl = "\n"


# Check that input was given.
if config["samplesheet"] == "NA":
    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config rundir=\"path/to/samplesheet/\" to the command line call.")
if config["rundir"] == "NA":
    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")
# TODO: Implement additional input validation, like checking that the objects given are file and dir respectively.

print("        /\\                          | |   | |     / __ \\ / ____| KMA,AUH")
print("       /  \\   ___ ___  ___ _ __ ___ | |__ | | ___| |  | | |      ")
print("      / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ |  | | |     ")
print("     / ____ \\__ \\__ \\  __/ | | | | | |_) | |  __/ |__| | |____  ")
print("    /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|\\___\\_\\_____| ")
print(" ")


print(f"These are the parameters given:")
print(f"  samplesheet: {samplesheet}")
print(f"  rundir: {rundir}")
print()


#########################
# Parse the samplesheet #
#########################

samplesheet_extension = samplesheet.split(".")[-1]
print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")

if samplesheet_extension == "xlsx":
    # Uses openpyxl
    df = pd.read_excel(samplesheet, dtype = str)

elif samplesheet_extension == "xls":
    df = pd.read_excel(samplesheet)

else:
    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")

# Clean up the spreadsheet
print("Cleaning sample sheet ...                              ", end = "", flush = True)
df.columns = map(str.lower, df.columns) # Lowercase
df.columns = map(str.strip, df.columns) # Remove edge-spaces
df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
df = df.dropna(subset = ["sample_id"])# remove rows not containing a sample ID
print("✓")

# Check that the spreadsheet complies
print("Checking that the necessary columns exist ...          ", end = "", flush = True)
for i in ["barcode", "sample_id"]:
    if not i in df.columns:
        raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
print("✓")

acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)]

print("Checking that the barcodes are correctly formatted and that the barcodes are unique... ", end = "", flush = True)
for i in df["barcode"]:
    #print("checking", i)
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode \'{i}\' is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
if not len(df["barcode"]) == len(set(df["barcode"])):
    raise Exception(f"One or more barcodes are duplicated. Each barcode may only be used once")
print("✓")


print("Checking that the sample id's are unique ...           ", end = "", flush = True)
if not len(df["sample_id"]) == len(set(df["sample_id"])):
    raise Exception(f"One or more sample_id's are duplicated. Each sample_id may only be used once")
print("✓")

print()
print("These are the samples from the samplesheet you have given:")
print(df.to_string())
print("//")
print()

###################
# Validate rundir #
###################

# Wait for the rundir to occur in the specified path. 
# If it doesn't occur after a specified waiting time, then stop the p

if rundir[-1] == "/":
    print("Removing trailing slash from rundir")
    rundir = rundir[0:-1]


print("Checking that the rundir exists ...                    ", end = "", flush = True)
if not os.path.isdir(rundir):
    raise Exception(f"The rundir does not exist.")
print("✓")

print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)

for i in range(200):
    print("  Looking ... ", end = "", flush = True)
    fastq_pass_bases = glob.glob(rundir + "/**/fastq_pass", recursive = True) # Find any occurrence of the wanted path
    if len(fastq_pass_bases) == 0:
        print("nothing found yet, waiting 10 secs ...")
        time.sleep(10) # Wait 10 seconds.
    elif(i == 10):
        print() # clean newline
        raise Exception("nothing found after 10 tries. Aborting.")
    else: 
        print(f"Found                                    ✓")
        break


if not len(fastq_pass_bases) == 1:
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found:{nl} {str(nl + ' ').join(fastq_pass_bases)}{nl}Please specify a more specific rundir.")


fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases

# base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"This is the batch base directory:{nl}  {base_dir}")



very_long_batch_id = base_dir.split("/")[-1]
print(f"This is the very long batch id:", very_long_batch_id)

date_parse, time_parse, minion_parse, flowcell_parse, arbhash_parse = very_long_batch_id.split("_")

print("date:    ", date_parse)
print("time:    ", time_parse)
print("minion:  ", minion_parse)
print("flowcell:", flowcell_parse)
print("arbhash: ", arbhash_parse)




batch_id = ".".join(very_long_batch_id.split("_")[0:2]) # The first two words (date, time), joined by a dot.
print(f"This is the parsed batch_id:", batch_id)



out_base = os.path.join(base_dir, "assembleQC_output") # out_base is the directory where the pipeline will write its output to.


# And here is the code from the rule wait_for_minknow
minutes_wait = 10
print("Checking that the sequencing_summary_*.txt-file has been written to disk ...")
while True:
    sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
    if len(sequencing_summary_file) == 0:
        print(f"  Still sequencing/basecalling; waiting {minutes_wait} minutes ...")
        time.sleep(60*minutes_wait)
    else:
        break



sequencing_summary_file = sequencing_summary_file[0]
print("  The sequencing summary has been found                ✓")
#print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")
print(f"  This is the sequencing_summary_*.txt-file (full): \"{sequencing_summary_file}\"")


sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
df.to_csv(sample_sheet_given_file, sep = "\t")
print("✓")

print()

disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})

disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])



print("Continuing with the following barcodes:")

# the workflow_table is the table that contains the records where the barcode could be found on the disk.
workflow_table = disk_barcodes_df.merge(df, how='left', on='barcode') # left join (merge) the present barcodes onto the df table.
workflow_table = workflow_table.dropna(subset = ["sample_id"])

#print(workflow_table[["barcode", "sample_id", "type"]].to_string(index = False))
print(workflow_table)
print("//")
print()


####################################################################
# Now we can run the analysis
####################################################################


rule all:
    input:
        expand(["{out_base}/{sample_id}/read_filtering/{sample_id}.fastq", \
                "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz", \
                "{out_base}/{sample_id}/{sample_id}_kraken2_reads_report.txt", \
                "{out_base}/{sample_id}/{sample_id}_kraken2_reads_report.txt", \
                "{out_base}/{sample_id}/{sample_id}_assembly-stats.tab", \
                "{out_base}/flye/{sample_id}/{sample_id}_assembly.fasta", \
                #"{out_base}/{sample_id}/{sample_id}_consensus.gff", \
                #"{out_base}/{sample_id}/{sample_id}_consensus.gbk", \
                #"{out_base}/multiqc_report.html" \
                ], \
                
               out_base = out_base, sample_id = df["sample_id"])





####################
# Setup for data analysis #
####################

rule merge_reads:
    input:
        barcode_dir = directory(lambda wildcards: workflow_table[workflow_table["sample_id"] == wildcards.sample_id]["barcode_path"].values[0])
    output: "{out_base}/{sample_id}/read_filtering/{sample_id}.fastq" #_barcode00.fastq
    shell: """
    cd {input.barcode_dir}

    cat * > {output}


    """

# Trim adapters
rule trim_adapt:
    input: 
        "{out_base}/{sample_id}/read_filtering/{sample_id}.fastq"
    output: 
        "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz"
    conda: "configs/conda.yaml"
    threads: 4
    shell: """
    mkdir -p {out_base}/{wildcards.sample_id}/trimmed
    porechop -i {input} --format fastq.gz -t 4 -o {output}           
    """

#Kraken2
# TODO: Should kraken use the full read set or the downsampled? ####
rule kraken2:
    input:
        "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz"
    output: 
        "{out_base}/{sample_id}/{sample_id}_kraken2_reads_report.txt"
    threads: 8
    conda: "configs/conda.yaml"
    shell: """
        kraken2 --db {kraken2_db} --report {output} --threads 8 {input}       
    """

rule assemble:
    input: 
        "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz"
    output: 
        contigs = "{out_base}/flye/{sample_id}/{sample_id}_assembly.fasta",
        assembly_stats = "{out_base}/{sample_id}/{sample_id}_assembly-stats.tab"
    conda: "configs/conda.yaml"
    threads: 4
    shell: """

        mkdir -p {out_base}/{wildcards.sample_id}/assembly

        flye -t 4 -i 2 -g 2.5m  --nano-hq trimmed/barcode01_trimmed.fastq.gz --asm-coverage 50 --out-dir {out_base}/{wildcards.sample_id}/flye
        cp {out_base}/{wildcards.sample_id}/flye/assembly.fasta {output.contigs}
        assembly-stats -t {output.contigs} > {output.assembly_stats}       
            
        """

rule minimap2:
    input:
        reads = "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fq.gz",
        contigs = "{out_base}/flye/{wildcards.sample_id}/{wildcards.sample_id}_assembly.fasta"
    output:
        "{out_base}/{sample_id}/mapped_reads/{sample_id}.bam"
    conda: "configs/conda.yaml"
    threads: 8
    shell: """


            """

# This rule polishes the assembly using the read mapping generated by minimap2    
rule racon:
    input:

    output:

    threads: 8
    shell: """
        racon -m 8 -x -6 -g -8 -w 500 -t 14 ~/workdir/basecall/basecall_trimmed.fastq.gz ~/workdir/nanopore_mapping/mapping.sam ~/workdir/assembly/assembly.contigs.fasta > racon/racon.fasta
            """

# This should further improve the assembly. Based on https://www.nature.com/articles/s41598-021-00178-w#Sec2 - maybe we only need racon or maybe we only need medaka?            
rule medaka:
    input:
        R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
        R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz",
        contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta"
    output:
        "{out_base}/{sample_id}/mapped_reads/{sample_id}.bam"
    conda: "configs/conda.yaml"
    threads: 8
    shell: """

            """


rule annotate_genes:
    input:
        "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
    output:
        gff = "{out_base}/{sample_id}/{sample_id}_consensus.gff",
        gbk = "{out_base}/{sample_id}/{sample_id}_consensus.gbk"
    conda: "configs/prokka.yaml"
    threads: 8
    shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/prokka
        prokka --outdir {out_base}/{wildcards.sample_id}/prokka --cpu 8 --force --prefix {wildcards.sample_id} {input}
        cp {out_base}/{wildcards.sample_id}/prokka/{wildcards.sample_id}.gff {output.gff}
        cp {out_base}/{wildcards.sample_id}/prokka/{wildcards.sample_id}.gbk {output.gbk}

        """

rule multiqc:
    input:
        expand("{out_base}/{sample_id}/{sample_id}.gff", out_base = out_base, sample_id = df["sample_id"])
    output:
        "{out_base}/multiqc_report.html"
    conda: "configs/conda.yaml"
    threads: 1
    shell: """
        multiqc -d {out_base} -o {out_base}

        """


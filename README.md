
# MEMDS output analysis, SLURM version

The MEMDS analysis pipeline is designed to run on distributed computing systems (clusters) managed by SLURM.
(https://slurm.schedmd.com/overview.html)

The pipeline steps need to be run in sequential order, unless marked otherwise.

## Prerequisite programs
<details>
<summary> A rundown of external programs utilized by the pipeline </summary>
<p></p>

The easiest way to organize and install all pipeline dependencies is through the Conda package manager.

#### Install miniconda: https://docs.conda.io/en/latest/miniconda.html
#### Install bioconda: https://bioconda.github.io

#### Use Conda to create environment with the following programs:
*conda create -n modules3  python=2.7  bwa cutadapt fastqc pear perl picard pysam biopython samtools seqtk trimmomatic*

#### Check that all the programs are present in the newly created environment (program version might be different):
biopython   1.73\
bwa         0.7.17\
cutadapt    1.18\
fastqc      0.11.9\
pear        0.9.6\
perl        5.32.1\
pysam       0.15.3\
python      2.7.15\
samtools    1.6\
seqtk       1.3\
trimmomatic 0.39

If the user does not wish to use a package manager, all the programs listed above can be manually downloaded and installed from their respective websites. However, in this case, program locations must be added to the $PATH environment variable to ensure that the pipeline scripts can detect them.
</details>

## Beginning the run
<details>
<summary> Quick start section </summary>
<p></p>

#### a) Prepare the parameter files needed to run the pipeline, as outlined in the **"Parameter_file_preparation"** guide.
#### b) Use the "MEMDS_pipeline_wrapper" script to run desired pipeline steps.
#### c) Alternatively, navigate to the pipeline scripts folder and run them directly. Note that the scripts might produce an error if activated from outside their folder.
</details>

## Pipeline steps
<details>
<summary> Step-by-step description of pipeline run </summary>
<p></p>

### 1) Joining partial \*.fastq files (if data is separated across several lanes):
a) **Prepare:** fastq_merging/samples_table_0.txt\
b) **Run:** *bash concatenate_partfiles.sh*\
c) **Check** output file "fastq_merging/samples_table_0.sh.concat.sh". Make sure it doesn't include comments indicating errors.\
d) **Run**: *srun bash fastq_merging/samples_table_0.sh.concat.sh*\
e) **Check**: concatenated output files to ensure that no data is missing.

### 2) Organizing parameter files:
a) **Prepare in the "config_files" folder**: params_1.sh; samples_table.txt; factors_table.txt\
b) In the **sequences directory** (as listed in **params_1.sh**), prepare reference fasta files. Each reference file should contain single reference sequence. All reference files should have the extension **.fa**.\
c) **Run**: \
   *bash setting_1-PE.sh* (Paired-end data) **or**\
   *bash setting_1-SE.sh* (Single-end data)\
d) **Check**: output file **config_files/samples_table.sh**. Make sure it contains all parameter data and doesn't include comments indicating errors.

### 3) Quality control and merging (for paired-end data):
a) **Run**: \
    *bash filter-PE4.sh $i* ($i is a number 1 to 5 (5 is optional)) for paired-end data **or**\
    *bash filter-SE4.sh $i* ($i is a number 1 to 4 (4 is optional)) for single-end data.\
   **Note**: All steps are running on Slurm. The steps should be run in a sequential manner; **don't start** a new step before completion of the previous one!\
b) In each step check that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** the output at each step:\
   **step 1:** Check Fastqc \*.html output for quality analysis of the raw data.\
   **step 2: For paired-end data only.** Check Pear results (merged fastq files). Check **\*.assembled.info** file for the percentage of reads that were merged.\
   **step 3:** Check quality-filtered fastq files - Cutadapt and Trimmomatic output - and log files. Each filtered fastq can be viewed with the tool "less" in Linux.\
   **step 4:** Check Fastqc \*.html output for quality analysis of the filtered data.\
   **step 5: Optional.** Check subsampled reads' data for "sanity checks" in the **"tests"** directory.

### 4) Separating barcode and gene sequences:
a) **Run:** *bash trim7.sh 1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** in the **filtered** directory trimming information files: **\*.trimmed.barcodes** files (correctly trimmed), **\*.wrongId.barcodes** files (unexpected identifiers), and **\*.trimmed.log** files (stating percentage of reads with correct identifiers).\
d) **Check** that correct barcode sequences are recognized by the script and that most reads in the data contain correct barcodes.

### 5) Sorting paralogous genes based on unique sequence signature:
a) **Run:** *bash sort2.sh 1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** sorted fastq files in the **"sorted"** directory. The **\*.others.fastq** files contain reads that couldn't be matched to either gene of origin.

### 6) Mapping filtered, trimmed and sorted reads to reference:
a) **Run:** *bash bwa9.sh 1*\
b) **Check** that dictionary files are generated for each **\*.fa** file in the **reference sequence directory**.\
c) **Run:** *bash bwa9.sh 2*\
d) **Check** that the jobs are completed on the grid using **"sacct"** command.\
e) **Check** that read-reference alignment files (**\*.sam** and **\*.bam**) are created in the **"mapping"** directory for each sorted fastq file.

### 7) Adjusting mapping results so they can be viewed with IGV viewer:
a) **Run:** *bash create_dummy_genome5.sh 1*\
b) Check the **\*.barcode.fasta** files are created in the **"mapping"** direcory.\
c) After this stage, all files for visual inspection of the alignments by IGV should be ready in the **"mapping"** directory: **\*.bam**, **\*.bam.bai**, **\*.barcode.fasta**.

### 8) Detecting potential mutations per primary-barcode family:
a) **Run:** *bash sam_to_mutation-table_5.2.sh 1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** the resulting output mutation tables (**\*.mutationFrequncyPerBarcode.txt**) and log files (**\*.mutationFrequncyPerBarcode.log.txt**) in the **"tables.BC3cutoff1"** directory. The log files list number of reads analyzed in each sorted fastq file for mutation presence.\
d) **Check** **"tag-BC3_missing"** output mutation tables for analysis of reads lacking secondary barcode, if the option is enabled in the script.

### 9) Detecting potential mutations and their quality per read:
a) **Run:** *bash sam_to_mutation-list-3.sh 1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** the resulting output mutation tables (**\*.mutations.txt**) and log files (**\*.log.txt**) in the **"mutations"** directory. The log files list number of reads analyzed in each sorted fastq file for mutation presence.\
   **Note**: The log files produced in this step should be identical to the log files produced in Step 8, since both scripts use same filters when analyzing reads.

### 10a) Collating mutation data per primary-barcode family:
a) **Run:** *bash consensus_15.1.sh 1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** the resulting summary tables in the **"tables_consensus.BC3cutoff1"** directory. In these tables all mutation information per primary-barcode family presented in a single row, as a semi-colon separated lists.  

### 10b) Generating consensus mutation tables:
a) **Run:** *bash consensus_15.1.sh 2*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** the resulting consensus tables in the **"tables_consensus.BC3cutoff1"** directory. Consensus tables list mutation profiles of each primary-bracode read family, after filtering them through a defined set of cutoff criteria.\
   Mutations passing the combined cutoff criteria are shown as is. Ambigous positions where both mutation and wild-type reads do not pass the combined cutoff criteria are shown as 'N'. Read families where all mutations are rejected, but wild-type reads pass the combined cutoff criteria at these positions are considered wild-type ('WT').
</details>

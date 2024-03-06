#### Activate module load environment ####
. /etc/profile.d/modules.sh

#### Activate Conda environment ####
. /data/apps/Easybuild/apps/Miniconda3/4.9.2/etc/profile.d/conda.sh
conda activate modules3
#############################

#Folder paths
params_adapters_1="$PWD/wildcard_adapters_1.fa" # File of adapter sequences to trim from the data
params_dir_out_1='/data/home/livnat/dmelamed/daniel/APOL1/Apl1_try/out1' # Pipeline output folder
params_dir_reference='/data/home/livnat/dmelamed/daniel/APOL1/Apl1_try/seqs' # Reference sequence folder

#Trimmomatic options
params_minimum_fastq_size_1=90 
qual_threshold=30
qual_window=3

is_SE=0

# Read sorting
offset=3

# Mutation calling
BC3_min=1 # Minimum 3' barcode count to call a mutation in mutation table
TSS=1102,-1
filter_pos=17

# Consensus cut-offs
min_freq="0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0"
min_count="0,1,2,3,4,5,6,7,8,9,10,25,50"
min_bc3="0,1,2,3,4,5"
bc3groupCountOK=0
bc3groupsWithReadsAbove1=0


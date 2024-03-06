# A wrapper script for paired-end raw ".fastq" file quality analysis, merging and trimming

# Help section
#if [[ "$1" != [1-4] ]]; then
#    echo "Use options 1-4 to run the actual jobs"
#    printf "\$1:\n1 - Fastqc raw data\n2 - Pear read merging\n3 - Filter adapters - Cutadapt and Trimmomatic\n4 - Fastqc trimmed files\n5 - subsample raw and trimmed data\n"
#fi
###########################################
# Gather pipeline parameter data and check that relevant parameters are non-empty strings
params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$params_1"
assert_
. "$params_2"
assert_
###########################################
# Check that required parameters and files exist
if [ ! -f "$params_adapters_1" ] || [ ! -n "${params_minimum_fastq_size_1+1}" ]; then
    exit "Error: not all parameters/files exit"
fi

# Create the output folders
if [ ! -d "$params_dir_out_1" ]; then
    mkdir "$params_dir_out_1"
fi

echo "Output dir: $params_dir_out_1"

for out_dir in "filtered" "fastqc" "tests"; do
    if [ ! -d "$params_dir_out_1/$out_dir" ]; then mkdir "$params_dir_out_1/$out_dir"; assert_; fi
    if [ ! -d "$params_dir_out_1/$out_dir/logs" ]; then mkdir "$params_dir_out_1/$out_dir/logs"; assert_; fi
done

runtime=$(date +"%y-%m-%dT%H%M")
logdir="logs_$runtime"

# Iterate over analyzed sample treatment (Cont/Exp)
for i in ${!title[@]}; do

    # Collect input data parameters
    echo $i
    title1="${title[i]}"_"idx$i"
    f1="${for1[i]}"
    
    echo "f1 = $f1"
    echo "title = $title1"
    
    out_prefix1="$params_dir_out_1/filtered/$title1"
    echo "filtered output = $out_prefix1"
    
    # Validate existence of the input files
    if [ -f "$f1" ]; then
        
        # Run FastQC to assess data quality
        if [ $1 -eq 1 ] || [ $1 -eq 10 ]; then
            log_out1="$params_dir_out_1/fastqc/logs/$logdir/$title1"
            if [ ! -d "$params_dir_out_1/fastqc/logs/$logdir" ]; then mkdir "$params_dir_out_1/fastqc/logs/$logdir"; assert_; fi
            echo "logout = $log_out1"        
            
            if [ $1 -eq 1 ]; then
                sbatch -N1 -n1 --ntasks-per-node=1 -o "$log_out1.f.out" \
                       -p hive1d,hive7d,hiveunlim,queen \
                       --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$f1\""
            fi       
        
        # Trim reads to remove contaminants and low quality data
        elif [ $1 -eq 2 ] || [ $1 -eq 20 ]; then
            log_out1="$params_dir_out_1/filtered/logs/$logdir/$title1"
            if [ ! -d "$params_dir_out_1/filtered/logs/$logdir" ]; then mkdir "$params_dir_out_1/filtered/logs/$logdir"; assert_; fi
            echo "logout = $log_out1"
            
            if [ $1 -eq 2 ]; then
                sbatch --output="$log_out1.filtered.out" --error="$log_out1.filtered.err" \
                       -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                       "$PWD/fastq-filter_job_3.sh" "SE" "$params_adapters_1" "$qual_threshold" "$qual_window" \
                       "$params_minimum_fastq_size_1" "$f1" "$out_prefix1.filtered.fastq"
            fi
            
        # Run FastQC to assess data quality after trimming
        elif [ $1 -eq 3 ] || [ $1 -eq 30 ]; then
            if [ -f "$out_prefix1.filtered.fastq" ] && [ -d "$params_dir_out_1/fastqc" ]; then
                log_out1="$params_dir_out_1/fastqc/logs/$logdir/$title1"
                if [ ! -d "$params_dir_out_1/fastqc/logs/$logdir" ]; then mkdir "$params_dir_out_1/fastqc/logs/$logdir"; assert_; fi
                echo "logout = $log_out1"
                
                if [ $1 -eq 3 ]; then
                    sbatch -N1 -n1 --ntasks-per-node=1 --output="$log_out1.fastqc.out" --error="log_out1.fastqc.err" \
                           -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$out_prefix1.filtered.fastq\""
                fi       
            fi
        
        # Create small subsamples of analyzed ".fastq" files for manual analyses, if needed
        elif [ $1 -eq 4 ] || [ $1 -eq 40 ]; then
            log_out1="$params_dir_out_1/tests/$logdir/$title1"
            if [ ! -d "$params_dir_out_1/tests/$logdir" ]; then mkdir "$params_dir_out_1/tests/$logdir"; assert_; fi
            echo "logout = $log_out1"
        
            # Subsample merged and trimmed output files
            if [ -f "$out_prefix1.filtered.fastq" ]; then
                if [ $1 -eq 4 ]; then
                    sbatch --output="$log_out1.subsample1000.out" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "head -n400 \"$out_prefix1.filtered.fastq\" > \"$params_dir_out_1/tests/$title1.filtered.head100.fastq\""
                        
                    sbatch --output="$log_out1.subsample1000.out" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "seqtk sample -s 100 \"$out_prefix1.filtered.fastq\" 10000 > \"$params_dir_out_1/tests/$title1.filtered.random-subsample1000.fastq\""
                fi    
            fi
            
            # Subsample input files
            if [ -f "$f1" ]; then
                if [ $1 -eq 4 ]; then
                    sbatch --output="$log_out1.f.subsample1000.out" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "head -n400 \"$f1\" > \"$params_dir_out_1/tests/$title1.head100.fastq\""
                        
                    sbatch --output="$log_out1.f.subsample1000.out" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "seqtk sample -s 100 \"$f1\" 10000 > \"$params_dir_out_1/tests/$title1.random-subsample1000.fastq\""
                fi    
            fi
        fi
    else
        echo 'At least one of the input files does not exist'
    fi
    echo '----------------'
done

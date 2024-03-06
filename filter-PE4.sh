# A wrapper script for paired-end raw ".fastq" file quality analysis, merging and trimming

# Help section
#if [[ "$1" != [1-5] ]]; then
#    echo "Use options 1-5 to run the actual jobs"
#    printf "\$1:\n1 - Fastqc raw data\n2 - Pear read merging\n3 - Filter adapters - Cutadapt and Trimmomatic\n4 - Fastqc trimmed files\n5 - Subsample raw and trimmed data\n"
#fi
###########################################
# Gather pipeline parameter data
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
    exit "Error: not all parameters/files exist"
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
    title1="${title[i]}"_"idx$i"
    f1="${for1[i]}"
    f2="${rev1[i]}"
    
    echo "forward reads = $f1"
    echo "reverse reads = $f2"
    echo "title = $title1"
    
    out_prefix1="$params_dir_out_1/filtered/$title1"
    echo "filtered output = $out_prefix1"
    
    # Validate existence of the input files
    if [ -f "$f1" ] && [ -f "$f2" ]; then
        
        # Run FastQC to assess data quality
        if [ $1 -eq 1 ] || [ $1 -eq 10 ]; then
            log_out1="$params_dir_out_1/fastqc/logs/$logdir/$title1"
            if [ ! -d "$params_dir_out_1/fastqc/logs/$logdir" ]; then mkdir "$params_dir_out_1/fastqc/logs/$logdir"; assert_; fi
            echo "logout = $log_out1"

            if [ $1 -eq 1 ]; then
                sbatch -N1 -n1 --ntasks-per-node=1 -o "$log_out1.f.out" -e "$log_out1.f.err" \
                       -p hive1d,hive7d,hiveunlim,queen \
                       --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$f1\""
                sbatch -N1 -n1 --ntasks-per-node=1 -o "$log_out1.r.out" -e "$log_out1.r.err" \
                       -p hive1d,hive7d,hiveunlim,queen \
                       --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$f2\""
            fi
        
        # Merge forward and reverse reads
        elif [ $1 -eq 2 ] || [ $1 -eq 20 ]; then
            log_out1="$params_dir_out_1/filtered/logs/$logdir/$title1"
            if [ ! -d "$params_dir_out_1/filtered/logs/$logdir" ]; then mkdir "$params_dir_out_1/filtered/logs/$logdir"; assert_; fi
            echo "logout = $log_out1"
             
            if [ $1 -eq 2 ]; then
                sbatch --output="$log_out1.assembled.out" --error="$log_out1.assembled.err" \
                       -N1 -n6 --ntasks-per-node=6 -p hive1d,hive7d,hiveunlim,queen \
                       --wrap "pear -f \"$f1\" -r \"$f2\" -o \"$out_prefix1\" -j 6 | tee \"$log_out1.assembled.info\""
            fi
        
        # Trim merged reads to remove contaminants and low quality data
        elif [ $1 -eq 3 ] || [ $1 -eq 30 ]; then
            if [ -f "$out_prefix1.assembled.fastq" ]; then
                log_out1="$params_dir_out_1/filtered/logs/$logdir/$title1"
                if [ ! -d "$params_dir_out_1/filtered/logs/$logdir" ]; then mkdir "$params_dir_out_1/filtered/logs/$logdir"; assert_; fi
                echo "logout = $log_out1"
                 
                if [ $1 -eq 3 ]; then
                    sbatch --output="$log_out1.assemb.filtered.out" --error="$log_out1.assemb.filtered.err" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           $PWD/fastq-filter_job_3.sh SE $params_adapters_1 $qual_threshold $qual_window \
                           $params_minimum_fastq_size_1 $out_prefix1.assembled.fastq \
                           $out_prefix1.assembled.filtered.fastq
                fi    
            fi
        
        # Run FastQC to assess data quality after merging and trimming
        elif [ $1 -eq 4 ] || [ $1 -eq 40 ]; then
            if [ -f "$out_prefix1.assembled.filtered.fastq" ]; then
                log_out1="$params_dir_out_1/fastqc/logs/$logdir/$title1"
                if [ ! -d "$params_dir_out_1/fastqc/logs/$logdir" ]; then mkdir "$params_dir_out_1/fastqc/logs/$logdir"; assert_; fi
                echo "logout = $log_out1"
                 
                if [ $1 -eq 4 ]; then
                    sbatch -N1 -n1 --ntasks-per-node=1 -o "$log_out1.assembled.fastqc.out" \
                           -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$out_prefix1.assembled.filtered.fastq\""
                fi    
            fi
        
        # Create small subsamples of analyzed ".fastq" files for manual analyses, if needed
        elif [ $1 -eq 5 ] || [ $1 -eq 50 ]; then
            log_out1="$params_dir_out_1/tests/$logdir/$title1"
            if [ ! -d "$params_dir_out_1/tests/$logdir" ]; then mkdir "$params_dir_out_1/tests/$logdir"; assert_; fi
            echo "logout = $log_out1"
            
            # Subsample merged and trimmed output files
            if [ -f "$out_prefix1.assembled.fastq" ]; then
                if [ $1 -eq 5 ]; then
                    sbatch --output="$log_out1.assembled.subsample1000.out" --error="$log_out1.assembled.subsample1000.err" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "seqtk sample \"$out_prefix1.assembled.fastq\" 10000 > \"$params_dir_out_1/tests/$title1.assembled.subsample1000.fastq\""
                fi    
            fi
            
            # Subsample input files
            if [ -f "$f1" ] && [ -f "$f2" ]; then
                if [ $1 -eq 5 ]; then
                    sbatch --output="$log_out1.f.subsample1000.out" --error="$log_out1.f.subsample1000.err" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "seqtk sample -s100 \"$f1\" 10000 > \"$params_dir_out_1/tests/$title1.f.subsample1000.fastq\""
                    
                    sbatch --output="$log_out1.r.subsample1000.out" --error="$log_out1.r.subsample1000.err" \
                           -N1 -n1 --ntasks-per-node=1 -p hive1d,hive7d,hiveunlim,queen \
                           --wrap "seqtk sample -s100 \"$f2\" 10000 > \"$params_dir_out_1/tests/$title1.r.subsample1000.fastq\""
                fi    
            fi
        fi
    else
        echo 'At least one of the input files does not exist'
    fi
    echo '----------------'
done

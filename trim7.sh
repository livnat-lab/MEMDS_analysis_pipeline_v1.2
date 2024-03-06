# A wrapper script for trimming reads (separating gene sequence from the 5' and 3' barcodes) and removing reads with bad barcodes

params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"

########################################################
# Gather pipeline parameter data and check that relevant parameters are non-empty strings
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$params_1"
assert_
. "$params_2"
assert_

if [ ! -d "$params_dir_out_1" ] || [ ! -d "$params_dir_out_1/filtered" ] || [ ! -n "${is_SE+1}" ]; then
    exit "Error: not all parameters/files exit"
fi

########################################################
# Iterate over analyzed sample treatment (Cont/Exp)
for i in ${!title[@]}; do
    
    # Collect input parameters
    title1="${title[i]}"_"idx$i"
    prefix1="$params_dir_out_1/filtered/$title1"
    
    # Define the log file directory
    runtime=$(date +"%y-%m-%dT%H%M")
    logdir="logs_$runtime"
    if [ ! -d "$params_dir_out_1/filtered/logs/$logdir" ]; then mkdir "$params_dir_out_1/filtered/logs/$logdir"; assert_; fi
    
    log_out1="$params_dir_out_1/filtered/logs/$logdir/$title1"
    
    # Define I/O prefixes
    if [ $is_SE == 0 ]; then
        f1=$prefix1.assembled.filtered.fastq
        f1_log=$log_out1.assembled.filtered.fastq
    else
        f1=$prefix1.filtered.fastq
        f1_log=$log_out1.filtered.fastq
    fi
    
    # Collect input parameters - cont.
    size_f1="${size_f[$i]}"; assert_
    size_r1="${size_r[$i]}"; assert_
    read_seq="${read_seq[$i]}"; assert_
    read_pos="${read_pos[$i]}"; assert_
    read_action="${read_action[$i]}"; assert_
    
    echo "fastq in = $f1"
    echo "logout = $f1_log"
    
    echo "bc_size_f = $size_f1"
    echo "bc_size_r = $size_r1"
    echo "bc_id_seq = $read_seq"
    echo "bc_id_pos  = $read_pos"
    
    # Check that input file exists
    if [ -f "$f1" ]; then
        echo 'Input files are found'
        
        # Run the jobs, if '1' is specified; otherwise run the wrapper w/o job execution, to test that input parameters are correct
        if [ $1 -eq 1 ]; then
            sbatch --output="$f1_log.trimmed.out" --error="$f1_log.trimmed.err" -N1 -n1 --ntasks-per-node=1 \
                    -p hive1d,hive7d,hiveunlim,queen \
                    --wrap "python trim7.py \"$f1\" \"$f1.trimmed\" \"$size_f1\" \"$size_r1\" \"$read_seq\" \"$read_pos\" \"$read_action\""
        fi
    else
        echo 'Not all files exist'
    fi
    
    echo "------------------"
done


# A wrapper script to run jobs sorting reads by their origin gene

params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"

########################################################
# Gather pipeline parameter data 

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$params_1"
assert_
. "$params_2"
assert_
########################################################

# Check that relevant parameters are properly defined
if [ ! -d "$params_dir_out_1" ] || [ ! -d "$params_dir_out_1/filtered" ] || [ ! -n "${is_SE+1}" ]; then
    exit "error: not all parameters/files exit"
fi

# Create ouptut directories
outdir="$params_dir_out_1/sorted"
if [ ! -d $outdir ]; then
    mkdir $outdir
    assert_
fi

runtime=$(date +"%y-%m-%dT%H%M")
logdir="$outdir/logs_$runtime"
if [ ! -d $logdir ]; then mkdir $logdir; assert_; fi

# Iterate over analyzed sample treatment (Cont/Exp)
for i in ${!title[@]}; do
    
    # Gather input data parameters 
    title1="${title[i]}"_"idx$i"
    assert_
    prefix1="$params_dir_out_1/filtered/$title1"
    
    if [ $is_SE == 0 ]; then
        f1=$prefix1.assembled.filtered.fastq.trimmed.fastq
    else
        f1=$prefix1.filtered.fastq.trimmed.fastq
    fi
    
    pos="${sort_pos[$i]}"
    nucl="${sort_nucl[$i]}"
    refs="${sort_refs[$i]}"
    ref="${sort_ref[$i]}"
    match_per="${sort_match[$i]}"
    
    echo "in_fastq = $f1"
    echo "refs = $refs;$ref"
    echo "sorting_pos = $pos"
    echo "sorting_nucl = $nucl"
    echo "sorting_match_% = $match_per"
    echo -e "match_offset = $offset\n"
    
    # Check that input file exists
    if [ -f "$f1" ]; then
        echo 'Input files are found'
        echo "logout = $logdir"
        
        # Run the jobs, if '1' is specified; otherwise run the wrapper w/o job execution, to test that input parameters are correct
        if [ $1 -eq 1 ]; then
            echo "Sorting $title1"
            assert_
            
            sbatch --output="$logdir/$title1.sorting.out" --error="$logdir/$title1.sorting.err" \
                   -N1 -n1 --ntasks-per-node=1 \
                   -p hive1d,hive7d,hiveunlim,queen \
                   --wrap "python sort2.py --fastq \"$f1\" --pos \"$pos\" --nucl \"$nucl\" \
                            --refs \"$refs\" --ref \"$ref\" --dirOut \"$outdir\" \
                            --title \"$title1\" --match_percentage \"$match_per\" --offset \"$offset\""
        fi
    else
        echo 'Not all files exist'
    fi
    
    echo '---------------'
done


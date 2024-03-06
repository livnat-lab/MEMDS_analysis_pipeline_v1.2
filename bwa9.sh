# A wrapper script to run BWA alignment jobs

# Check that valid input options were supplied to the script
if [[ $1 -ne 1 ]] && [ $1 -ne 2 ] && [ $1 -ne 10 ] && [ $1 -ne 20 ]; then
    printf "\$1:\n1 - Picard - create reference sequence dictionary\n2 - Map reads to the reference\n10 or 20 - Test run the above without executing the jobs"
    exit "Insert valid option\n"
fi
###########################################
# Gather pipeline parameter data
params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

# Check that all needed parameters and input data exist and well-defined
if [ ! -d "$params_dir_out_1" ] || [ ! -d "$params_dir_out_1/filtered" ] || [ ! -n "${is_SE+1}" ] || \
    [ ! -n "${params_dir_reference+1}" ] || [ ! -n "${sort_refs+1}" ] || [ ! -n "${sort_ref+1}" ] || [ ! -n "${reference_size+1}" ]; then
    exit "Error: not all parameters/files exit"
fi

# Create the output folders
outdir="$params_dir_out_1/mapping"
if [ ! -d $outdir ]; then mkdir $outdir; assert_; fi

runtime=$(date +"%y-%m-%dT%H%M")
logdir="$outdir/logs_$runtime"
if [ ! -d $logdir ]; then mkdir $logdir; assert_; fi

#########
# Create dictionary files from reference sequences for sequence alignment jobs
if [ $1 -eq 10 ] || [ $1 -eq 1 ]; then
    
    for ref1 in "$params_dir_reference/"*.fa; do
        echo $ref1
        if [ ! -f $ref1.indices.OK ] && [ $1 -eq 1 ] && [ -f "$ref1" ]; then
            bwa index $ref1
            assert_
            samtools faidx $ref1
            assert_
            #java -jar $params_picardDir/picard.jar CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=$ref1.dict 
            picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=$ref1.dict 
            assert_
            touch $ref1.indices.OK
            assert_
        fi
    done
fi
    
# Run read alignment jobs on the sorted reads
if [ $1 -eq 20 ] || [ $1 -eq 2 ]; then
    # Iterate over analyzed sample treatment (Cont/Exp)
    for i in ${!title[@]}; do

        # Gather input data parameters 
        title1="${title[i]}"_"idx$i"; assert_
        prefix1="$params_dir_out_1/sorted/$title1";
        refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
        
        echo "title = $title1"
        echo "prefix = $prefix1"        
        
        # Iterate over sorted read files to map the reads against their reference
        k=0
        for r1 in $refs1; do
            let k=$k+1
            echo $k
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # First item in loop refers to "others"
            
            # Define I/O variables
            f1="$prefix1.$r2.fastq"
            fout1="$params_dir_out_1/mapping/$title1.$r2.bwa"
            log_out1="$logdir/$title1.$r2.bwa"
            ref1="$params_dir_reference"/"$r1".fa
            ref1_length="${reference_size[$i]}"
            
            echo 'len = '$ref1_length
            echo 'ref = '$ref1
            echo 'f1 = '$f1
            echo 'fout = '$fout1
            echo 'logout = '$log_out1
            
            # Check that all required parameters are OK before executing the jobs
            if [ -f "$f1" ] && [ -d $outdir ] && [ -f "$ref1" ] && [ -f "$ref1.indices.OK" ] && [ $ref1_length != "" ]; then
                echo 'All input parameters and files are OK'
                
                # Run read alignment jobs
                if [ $1 -eq 2 ]; then
                    sbatch  -N1 -n20 --ntasks-per-node=20 \
                            -p hive1d,hive7d,hiveunlim,queen \
                            -o "$log_out1".out -e "$log_out1".err \
                            --export=ref1="$ref1",f1="$f1",fout1="$fout1",threads=20,ref1_length="$ref1_length",params_1="$params_1" \
                            bwa_job9.sh
                fi
            else
                echo 'Error: not all input files or parameters exist'
            fi
        done
        echo '--------------'
    done
fi


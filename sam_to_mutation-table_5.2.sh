# A wrapper script for identifying mutations in the aligned reads

params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"

run_on_slurm=1

#######################################
# Gather pipeline parameter data and check that read mapping step was performed
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

if [ ! -d "$params_dir_out_1/mapping" ]; then
    exit "Error: not all parameters/files exit"
fi

# Create log file capturing run information
run_log="$params_dir_out_1/Sam2mut_table_bash_log.txt"
if [ -f "$run_log" ]; then rm $run_log; fi

#########
# Iterate over analyzed samples to perform mutation identification
# BC3_min_cutoff variable allows to filter read groups sharing same BC3 within the read family that have less reads than the specified cutoff value
for BC3_min_cutoff in $BC3_min; do

    # Create the output folders
    outdir="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff"
    if [ ! -d "$outdir" ]; then mkdir "$outdir"; assert_; fi
    
    runtime=$(date +"%y-%m-%dT%H%M")
    logdir="$outdir/logs_$runtime"
    if [ ! -d "$logdir" ]; then mkdir "$logdir"; assert_; fi
    
    # Iterate over analyzed sample treatment (Cont/Exp)
    for i in ${!title[@]}; do
    
        # Collect input data parameters
        title1="${title[i]}"_"idx$i"
        refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
        limit_start="${limit_starts[$i]}"
        limit_end="${limit_ends[$i]}"
        seq_pos="${seq_pos[i]}"
        seq_mut="${seq_mut[i]}"
        seq_action="${seq_action[i]}"
        
        echo -e "******\n" >> $run_log
        echo "limits = $limit_start - $limit_end" >> $run_log
        echo "BC3_min_cutoff = $BC3_min_cutoff" >> $run_log

        # Iterate over sorted read files within each sample
        k=0
        for r1 in $refs1; do
            echo "r1 = $r1" >> $run_log
            let k=$k+1
            echo "k = $k" >> $run_log
            
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # First item in loop refers to "others"
            
            bam="$params_dir_out_1/mapping/$title1.$r2.bwa.sorted.bam"
            ref1="$params_dir_reference"/"$r1".fa
            
            echo "bam = $bam" >> $run_log
            echo "ref1 = $ref1" >> $run_log
            
            # Test that all parameters are defined and that input files exist before executing the jobs
            if [ -f "$bam" ] && [ -f $ref1 ]; then
                echo 'Alignment and reference sequence files are found' >> $run_log
                
                if [ ! -z "$limit_start" ] && [ ! -z "$limit_end" ] && [ ! -z "$seq_pos" ] && [ ! -z "$seq_mut" ] && [ ! -z "$seq_action" ] && [ ! -z "$title1" ] && [ ! -z "$bam" ] && [ ! -z "$ref1" ] && [ ! -z "$BC3_min_cutoff" ]; then
                    echo -e "Analysis parameters are found\n" >> $run_log
                
                    # Analyze BC3 with specified sequence(s); '' - analyze all BC3 sequences
                    for sam_BX_tag in ''; do # to include only specific sam BX tags
                        echo "sam_BX_tag = \"$sam_BX_tag\"" >> $run_log
                        echo "title = $title1" >> $run_log
                        echo -e "log name = ""$logdir"/"$title1.$r2.bwa.sorted.bam.$sam_BX_tag.out""\n" >> $run_log
    
                        # Run the jobs, if '1' is specified; otherwise run the wrapper w/o job execution, to test that input parameters are correct
                        if [ $1 -eq 1 ]; then
                            echo 'OK' >> $run_log
                            sleep 2
                            
                            sbatch --mem=20000  -N 1 -n 1 --ntasks-per-node=1 \
                            -p hive1d,hive7d,hiveunlim,queen \
                            -o "$logdir"/"$title1.$r2.bwa.sorted.bam.$sam_BX_tag.out" -e "$logdir"/"$title1.$r2.bwa.sorted.bam.$sam_BX_tag.err" \
                            --wrap "python sam_to_mutation-table_5.2.py \
                                        --sam_files_path=\"$bam\" \
                                        --out_dir=\"$outdir\" \
                                        --ref_fasta=\"$ref1\" \
                                        --limit_start=\"$limit_start\" \
                                        --limit_end=\"$limit_end\" \
                                        --queryAlnMustStartAt0=1 \
                                        --include_BX_tag=\"$sam_BX_tag\" \
                                        --BC3_min_cutoff=\"$BC3_min_cutoff\" \
                                        --checkMut_name=\"$seq_mut\" \
                                        --checkMut_pos=\"$seq_pos\" \
                                        --checkMut_action=\"$seq_action\""
                        fi
                    done
               
                    echo -e "----------" >> $run_log
                
                else
                    echo 'Error: not all parameters OK'
                    exit
                fi
            else
                echo 'Error: not all input files exist'
            fi
        done
        echo '---------'
    done
done


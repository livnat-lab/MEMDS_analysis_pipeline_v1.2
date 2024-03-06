# A wrapper script to run jobs summarizing mutation information per read family

params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"

###########################################
# Gather pipeline parameter data
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

###########################################
# Run mutation summary jobs
for BC3_min_cutoff in $BC3_min; do
    
    # Check existence of the input directory and
    if [ ! -d "$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff" ]; then
        exit "Error: not all parameters/files exit"
    fi
    
    # Create the output directory
    outdir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff"
    if [ ! -d $outdir ]; then mkdir $outdir; assert_; fi
    
    runtime=$(date +"%y-%m-%dT%H%M")
    logdir="$outdir/logs_$runtime"
    if [ ! -d $logdir ]; then mkdir $logdir; assert_; fi
    
    # Iterate over analyzed genes and sample treatments (Cont/Exp)
    for i in ${!title[@]}; do
        
        # Gather input data parameters 
        size_r1="${size_r[i]}"
        title1="${title[i]}"_"idx$i"
        refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
        limit_start="${limit_starts[$i]}"
        limit_end="${limit_ends[$i]}"
        seq_pos="${seq_pos[i]}"
        seq_mut="${seq_mut[i]}"
        seq_action="${seq_action[i]}"
        
        echo "size_r1 = $size_r1"
        echo "limits = $limit_start - $limit_end"
        echo "BC3_min_cutoff = $BC3_min_cutoff"
        echo "Filter ${seq_pos}${seq_mut}"
        echo "********"
        
        # Iterate over sorted read mutation data
        k=0
        for r1 in $refs1; do
            let k=$k+1
            echo "k = "$k
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # First item in loop refers to "others"
            
            # Analyzed file and reference names
            bam_name="$title1.$r2.bwa.sorted.bam"
            ref1="$params_dir_reference"/"$r1".fa
            
            echo "bam_name = $bam_name"
            echo "ref1 = $ref1"
            
            # Iterate over reads having specific BC3s
            for sam_BX_tag in ''; do 
                # Define BC3 tag variable added to the output file name
                echo "sam_BX_tag = \"$sam_BX_tag\""
                tag1=""
                if [ "$sam_BX_tag" != "" ]; then tag1="tag-$sam_BX_tag."; fi
                
                # Define I/O variables
                table_in1="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""mutationFrequncyPerBarcode.txt"
                table_out1="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""consensus.txt"
                log_out1="$logdir/$bam_name.""$tag1""consensus"
                
                table_out2_dir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""cutoffs-update" 
                table_out2="$table_out2_dir/$bam_name"
                
                echo "table_in1 = $table_in1"
                echo "table_out1 = $table_out1"
                echo "table_out2 = $table_out2"
                echo "logout = $log_out1"
                
                # Check that all required parameters are OK before executing the jobs
                if [ -f "$table_in1" ] && [ -d "$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff" ]; then
                    if [ ! -z "$table_in1" ] && [ ! -z "$table_out1" ] && [ ! -z "$table_out2" ] && [ ! -z "$size_r1" ] && [ ! -z "$title1" ] && [ ! -z "$BC3_min_cutoff" ]; then
                        echo 'All files and parameters are found'
                        
                        c=$(cat $table_in1 | wc -l)
                        echo "c="$c
                        if [ $c -gt 1 ]; then
                            echo 'Input table is non-empty'
                            
                            # Run jobs to create mutation profiles of read families
                            if [ $1 -eq 1 ]; then
                                sleep 2
                                
                                sbatch -o "$log_out1.out" -e "$log_out1.err" \
                                       -p hive1d,hive7d,hiveunlim,queen \
                                       --wrap "python consensus_15.py \"$table_in1\" \"$table_out1\" \"$TSS\" \"$filter_pos\" \"$BC3_min\""
                            
                            # Run jobs checking identified mutations against cutoff criteria sets to determine "true" mutations and ambigous positions for each cutoff criteria set
                            elif [ $1 -eq 2 ]; then
                                if [ -e "$table_out1" ]; then
                                    echo 'Input consensus table is found'
                                    
                                    # Separate job are initiated for each min_bc3 cutoff to speed up the analysis
                                    # Additional loops for other criteria sets can be added to further accelerate analysis rate
                                    for bc in $(echo "$min_bc3" | tr "," "\n" | sort -n); do
                                        if [ ! -d "$table_out2_dir" ]; then mkdir "$table_out2_dir"; assert_; fi
                                        
                                        sbatch -o "$log_out1.cutoffs-update.out" -e "$log_out1.cutoffs-update.err" --mem=128000 \
                                               -p hive1d,hive7d,hiveunlim,queen \
                                               --wrap "python consensus_cutoffs_15.py \"$table_out1\" \"$table_out2\" \"$size_r1\" \"$min_freq\" \"$min_count\" \"$bc\" \"$bc3groupCountOK\" \"$bc3groupsWithReadsAbove1\""
                                    done
                                fi
                            fi
                        fi
                    else
                         echo "Error: some of the required parameters are empty"
                        exit
                    fi
                else
                    echo "Error: not all I/O data exists"
                fi
            done
            echo '----------------'
        done
    done
done

conda deactivate

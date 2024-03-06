# A wrapper script for jobs generating files needed for manual inspection of read alignments via IGV (Integrative Genomics Viewer)

params_1="$PWD/config_files/params_1.sh"
params_2="$PWD/config_files/samples_table.sh"


# Gather pipeline parameter data and check that directory with read alignment files exists
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_
#####################################
# Check that input folder exists
if [ ! -d "$params_dir_out_1/mapping" ]; then
    exit "Error: aligned read directory is missing"
fi

# Iterate over analyzed sample treatment (Cont/Exp)
for i in ${!title[@]}; do
    
    # Analyzed file and reference names
    title1="${title[i]}"_"idx$i"; assert_
    refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
    
    # Iterate over sorted and aligned read files
    k=0
    for r1 in $refs1; do
        let k=$k+1
        echo $k
        r2="$r1"
        if [ $k -eq 1 ]; then r2="$r2.others"; fi # first item in loop refers to "others"
        
        # Define I/O variables
        bam1="$params_dir_out_1/mapping/$title1.$r2.bwa.sorted.bam"
        ref1="$params_dir_reference"/"$r1".fa
        header="$params_dir_out_1/mapping/$title1.$r2.bwa.sam.header"
        
        echo "ref1 = $ref1"
        echo "header = $header"
        echo "bam = $bam1"
        
        # Check that all required parameters are OK before executing the jobs
        if [ -f  "$header" ] && [ -f "$ref1" ] && [ -f "$bam1" ] && [ -f "$ref1.indices.OK" ]; then
            echo 'Input files are OK'
            
            # Run the jobs, if '1' is specified; otherwise run the wrapper w/o job execution, to test that input parameters are correct
            if [ $1 -eq 1 ]; then
                if [ ! -f $header.barcode ]; then
                    cat $header | perl -ne 'chomp; if($_ =~ m/SN:(\w+)/){printf("%s\n",$1)};' > $header.barcode
                    assert_
                    echo $header.barcode
                    python create_dummy_genome5.py "$ref1" "$header.barcode"
                    assert_
                else
                    echo 'Error: '$header.barcode' already exists; exiting'
                    exit
                fi
            fi
        fi
        echo '-----------------'
    done
done


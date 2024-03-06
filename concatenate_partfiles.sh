# A wrapper script to generate merging scripts for partial ".fastq" input files

in1=fastq_merging/samples_table_0.txt
out1=fastq_merging/samples_table_err.log

[ ! -d "../raw_concatenated" ] && mkdir "../raw_concatenated" # Create merged output directory

python mapper_write-SE-or-PE-list_11.py --input_table "$in1" --output_script "$out1" --concatenate_partfiles --output_dir "../raw_concatenated"


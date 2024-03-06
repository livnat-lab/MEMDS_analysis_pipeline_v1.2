# A wrapper script for a job organizing analyzed sample parameters in a bash file to be used by the pipeline. For SE input.

in1=config_files/samples_table.txt
in2=config_files/factors_table.txt
out1=config_files/samples_table.sh

python mapper_write-SE-or-PE-list_11.py --input_table $in1 --output_script $out1 --factors_table $in2 --SE

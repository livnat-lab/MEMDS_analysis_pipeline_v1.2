#!/bin/sh

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_

# Run BWA aligner
bwa mem -M -t $threads "$ref1" "$f1" > "$fout1.0.sam"
assert_

# Detect the unique 5' barcode sequence at the start of the query name and make it the reference name, whilst removing it from the query name
# Filter barcodes with Ns (unknown bases) and/or low MAPQ score
cat "$fout1.0.sam" \
  | perl -ne 'chomp; @x=split /\t/,$_; if($x[0] =~ m/^(\w+)_(.*)/){$x[2]=$1; $x[0]=$2; if((!($x[2] =~ m/.*N.*/))and($x[4]>0)){$j=join "\t",@x; printf("%s\n",$j)}}' \
  | perl -ne 'chomp; @x=split /\t/,$_; if($x[2] =~ m/^(.*?)x(.*)$/){$bc5=$1; $bc3=$2; if(length($bc5)<1){$bc5="BC5_missing"}; if(length($bc3)<1){$bc3="BC3_missing"}; $x[2]=$bc5; push @x,"XB:Z:$bc3"; $j=join "\t",@x; printf("%s\n",$j)}' \
  > "$fout1.temp.sam"
assert_

# Remove splited alignments and place them in a separate file
awk -F'\t' 'FNR==NR { a[$1]++; next } a[$1] > 1' "$fout1.temp.sam" "$fout1.temp.sam" > "$fout1.splited.sam"
awk -F'\t' 'FNR==NR { a[$1]++; next } a[$1] == 1' "$fout1.temp.sam" "$fout1.temp.sam" > "$fout1.sam"

# Prepare ".sam" file header and add it to the ".sam" file
cut -f3 "$fout1.sam" | sort -u | perl -ne 'chomp; printf("\@SQ\tSN:%s\tLN:%s\n",$_,'$ref1_length')' > "$fout1.sam.header"
assert_

cat "$fout1.sam.header" "$fout1.sam" > "$fout1.1.sam"
assert_

# Convert ".sam" file with the headers to ".bam"; sort and index the output ".bam"
samtools view -Sb  "$fout1.1.sam" > "$fout1.bam"
assert_

samtools sort -T "$fout1.sorted.temp" -o "$fout1.sorted.bam" -O bam "$fout1.bam"
assert_

samtools index "$fout1.sorted.bam"
assert_

# Remove unneeded files
rm "$fout1.1.sam"
assert_
rm "$fout1.bam"
assert_
rm "$fout1.temp.sam"
assert_

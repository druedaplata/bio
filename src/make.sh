# Each directory includes:
#   - meta.csv, a file with all the additional information for each sequence
#   - sequences.fasta, all sequences for that "month_year" in a single file
#   - seqs, a directory with all sequences for that "month_year" in individual files

year=$1
year_dir=$(ls $1)
header=0

for d in $year_dir
do
    current_dir=$year/$d
    # 1. Creates a combined meta where all information available is saved
    if [[ -d $current_dir  ]]
       then
        if [[ $header != '1' ]]
        then
            # get header from first file
            echo $(head -n 1 $current_dir/meta.csv) > $year/meta.csv
            header=1
        fi
        tail -n +2 $current_dir/meta.csv >> $year/meta.csv

    # 2. Creates a single fasta file where all sequences available are saved.
        cat $current_dir/sequences.fasta >> $year/sequences.fasta
    fi
done
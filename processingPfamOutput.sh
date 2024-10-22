#!/bin/bash

#Be sure that the pfam output file hmmer_pfam_output_tbl.txt is in th same folder as this script
#run this script 2nd after you have run  convertingTSVtoSED_unixVersion.sh
echo "Step 1.1: Delete all lines that begin with # and convert spaces to tabs"
grep -v '^#' hmmer_pfam_output_tbl.txt | tr ' ' '\t' | awk '{print $1, $2, $3}' > intermediate_file.txt1
grep -v '^$' intermediate_file.txt1 | sed 's/^ .*\n//' | sed '/^$/d' > intermediate_file.txt2

echo "Step 1.2: Merging overlapping domains and deleting subset domains the output of this process is hmmer_pfam_output_tbl_trimmed.txt"
sed '1d' intermediate_file.txt2 > intermediate_file.txt2.1
# Input and output files
input_file="intermediate_file.txt2.1"
output_file="hmmer_pfam_output_tbl_trimmed.txt"

# Temporary variables to hold previous line details
prev_protein=""
prev_start=0
prev_end=0

# Clear output file
> "$output_file"

# Read the input file line by line
while IFS=$' ' read -r protein start end; do
    # If we are still in the same protein
    if [[ "$protein" == "$prev_protein" ]]; then
        # Check if domains overlap or are adjacent
        if [[ $start -le $prev_end ]]; then
            # Merge the domains by extending the previous end
            if [[ $end -gt $prev_end ]]; then
                prev_end=$end
            fi
        else
            # No overlap, write the previous domain to output and reset to the new domain
            echo -e "$prev_protein\t$prev_start\t$prev_end" >> "$output_file"
            prev_start=$start
            prev_end=$end
        fi
    else
        # New protein, write the previous domain to output and reset the variables
        if [[ -n "$prev_protein" ]]; then
            echo -e "$prev_protein\t$prev_start\t$prev_end" >> "$output_file"
        fi
        prev_protein="$protein"
        prev_start=$start
        prev_end=$end
    fi
done < "$input_file"

# Write the last domain to the output file
if [[ -n "$prev_protein" ]]; then
    echo -e "$prev_protein\t$prev_start\t$prev_end" >> "$output_file"
fi

echo "Step 2: Calculate the difference between the third and second columns, and append it as a new column"
awk '{print $0 "\t" ($3 - $2)}' hmmer_pfam_output_tbl_trimmed.txt > intermediate_file.txt3

sed 's/^/>/g' < intermediate_file.txt3 > intermediate_file.txt3.tmp

echo "Step 3: Summing up the domain lengths within proteins that have multiple domains i.e. ---|||1|||-----|||2|||------ Domain 1+2"
awk '{arr[$1]+=$4} END {for (i in arr) print i"\t"arr[i]}' intermediate_file.txt3.tmp > hmmer_pfam_output_tbl_reduced.txt


echo "Step 4: Creating sed file from hmmer_pfam_output_tbl_reduced.txt<<<"
sed 's/$/£/g' < hmmer_pfam_output_tbl_reduced.txt > hmmer_pfam_output_tbl_reduced.txt.tmp2
sed 's/\(>.*\)@\t\([0-9]*\)£/s#\1@#\2£/g' < hmmer_pfam_output_tbl_reduced.txt.tmp2 > hmmer_pfam_output_tbl_reduced.txt.tmp3
tr -d '\n' < hmmer_pfam_output_tbl_reduced.txt.tmp3 > hmmer_pfam_output_tbl_reduced.txt.tmp4
sed 's/£/\n/g' < hmmer_pfam_output_tbl_reduced.txt.tmp4 > hmmer_pfam_output_tbl_reduced.txt.tmp5
sed 's/$/#/g' < hmmer_pfam_output_tbl_reduced.txt.tmp5 > hmmer_pfam_output_tbl_reduced.txt.tmp6
sed 's/^>.*//g' < hmmer_pfam_output_tbl_reduced.txt.tmp6 > sed_file_domainology.txt

echo "Step 5: Substituting protein names in ortholog_table.tsv with their respective domain lengths output is ortholog_table_domain_lengths.tsv"
sed 's/-//g' < ortholog_table.tsv > ortholog_table_for_hmmer_output_parsing.tsv
sed -f sed_file_domainology.txt < ortholog_table_for_hmmer_output_parsing.tsv > ortholog_table_domain_lengths.tsv



echo "Step 6: removing orthologs that did not have an hmmer output from ortholog_table_domain_lengths.tsv"
awk '{print $1}' < intermediate_file.txt3 > intermediate_file.txt3.tmp2
echo ">>>Step 6.1: converting ortholog_table_for_hmmer_output_parsing.tsv to list_of_proteins_from_ortholog_table.tsv"

# Input and output file names
input_file="ortholog_table_for_hmmer_output_parsing.tsv"
output_file="list_of_proteins_from_ortholog_table.tsv"

# Initialize the output file
> "$output_file"

# Process the input file
awk '
{
    for (i = 1; i <= NF; i++) {
        # Check if the field starts with ">" and ignore fields with "0"
        if ($i ~ /^>/) {
            # Split the field on commas and loop through each value
            n = split($i, proteins, ",")
            for (j = 1; j <= n; j++) {
                # Trim any leading/trailing spaces
                gsub(/^ +| +$/, "", proteins[j])
                # Only print non-empty lines to the output
                if (proteins[j] != "") {
                    print proteins[j] >> "'$output_file'"
                }
            }
        }
    }
}' "$input_file"

input_file="list_of_proteins_from_ortholog_table.tsv"
output_file="sed_commands.sh"
# Initialize counter to track the temp files
counter=1
> "output_file"
echo "#!/bin/bash" > "$output_file"
# Loop through each line of the input file
while IFS= read -r line; do
    # Extract just the identifier part without '>' and '##'
#    identifier=$(echo "$line" | sed 's/[>#]//g')

   # Check if it's the first line (to use ortholog_table_domain_lengths.tsv as input)
    if [ $counter -eq 1 ]; then
        echo "sed 's/$line//' < ortholog_table_domain_lengths.tsv > temp1" >> "$output_file"
    else
        previous_counter=$((counter - 1))
        echo "sed 's/$line//' < temp$previous_counter > temp$counter" >> "$output_file"
    fi

    # Increment the counter
    counter=$((counter + 1))

done < "$input_file"

final_counter=$((counter - 1))
echo "mv temp$final_counter final_sed_file.tsv" >> "$output_file"
chmod +x "$output_file"
./"$output_file"
find . -maxdepth 1 -name "temp*" -print0|xargs -0 rm

echo "Step 7: calculating average domain lengths of orthologs. output is ortholog_table_domain_lengths_averages.tsv"
awk '{gsub(/, ,/, ","); print}' final_sed_file.tsv > tmp1.tsv
awk '{gsub(/, \t|\t /, "\t"); print}' tmp1.tsv > tmp2.tsv
num_cols=$(head -1 tmp2.tsv | awk -F'\t' '{print NF}')
for i in $(seq 1 $num_cols); do
    cut -f$i tmp2.tsv > column$i.txt
done
for filename in ./column*.txt
do
awk -F',' 'NR==1 {print $0} NR>1 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > "$filename".average
#awk -F',' 'NR==1 {print $0} NR>1 && NF>0 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > temp && mv temp "$filename".average
done
paste column*.average > tmp3.tsv
awk '{gsub(/-nan/, "0"); print}' tmp3.tsv > ortholog_table_domain_lengths_averages.tsv

echo "Step 8: creating interdomain length file"
echo ">>>Step 8.1: creating interdomain sed file"
sed -f sed_file_domainology.txt < list_of_proteins_from_ortholog_table.tsv > temp_domain_length.txt
#Note that the sed_file.txt below is from the first script convertingTSVtoSED_unixVersion.sh
sed 's/-//g' < sed_file.txt > sed_file_for_hmmr_domain_length.txt 
sed -f sed_file_for_hmmr_domain_length.txt < list_of_proteins_from_ortholog_table.tsv > temp_full_length.txt
paste temp_domain_length.txt temp_full_length.txt list_of_proteins_from_ortholog_table.tsv > temp_concat1.txt
sed 's/^>.*//g' < temp_concat1.txt > temp_concat2.txt
awk '{print $0 "\t" ($2 - $1)}' temp_concat2.txt > temp_concat3.txt
awk '{print $3, $4}' < temp_concat3.txt > temp_concat4.txt

sed 's/$/£/g' < temp_concat4.txt > temp_concat5.txt
sed 's/\(>.*\) \([0-9]*\)£/s#\1#\2£/'g < temp_concat5.txt > temp_concat6.txt
tr -d '\n' < temp_concat6.txt > temp_concat7.txt
sed 's/£/\n/g' < temp_concat7.txt > temp_concat8.txt
sed 's/$/#/g' < temp_concat8.txt > sed_file_interdomain_length.txt

echo ">>>Step 8.2: Creating ortholog_table_interdomain_lengths.tsv"
sed -f sed_file_interdomain_length.txt < ortholog_table_for_hmmer_output_parsing.tsv > ortholog_table_interdomain_lengths.tsv
rm temp_concat*
echo "removing orthologs from ortholog_table_interdomain_lengths.tsv that did not have a hmmer pfam output"
input_file="list_of_proteins_from_ortholog_table.tsv"
output_file="sed_commands_interdomain_lengths.sh"
# Initialize counter to track the temp files
counter=1
> "output_file"
echo "#!/bin/bash" > "$output_file"
# Loop through each line of the input file
while IFS= read -r line; do
    # Extract just the identifier part without '>' and '##'
#    identifier=$(echo "$line" | sed 's/[>#]//g')

    # Check if it's the first line (to use ortholog_table_domain_lengths.tsv as input)
    if [ $counter -eq 1 ]; then
        echo "sed 's/$line//' < ortholog_table_interdomain_lengths.tsv > temp1" >> "$output_file"
    else
        previous_counter=$((counter - 1))
        echo "sed 's/$line//' < temp$previous_counter > temp$counter" >> "$output_file"
    fi

    # Increment the counter
    counter=$((counter + 1))

done < "$input_file"
final_counter=$((counter - 1))
echo "mv temp$final_counter final_sed_file_interdomain_lengths.tsv" >> "$output_file"
chmod +x "$output_file"
./"$output_file"
find . -maxdepth 1 -name "temp*" -print0|xargs -0 rm
echo "calculating average values of interdomain lengths"
awk '{gsub(/, ,/, ","); print}' final_sed_file_interdomain_lengths.tsv > tmp1.tsv
awk '{gsub(/, \t|\t /, "\t"); print}' tmp1.tsv > tmp2.tsv
num_cols=$(head -1 tmp2.tsv | awk -F'\t' '{print NF}')
for i in $(seq 1 $num_cols); do
    cut -f$i tmp2.tsv > column$i.txt
done
for filename in ./column*.txt
do
awk -F',' 'NR==1 {print $0} NR>1 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > "$filename".average
#awk -F',' 'NR==1 {print $0} NR>1 && NF>0 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > temp && mv temp "$filename".average
done
paste column*.average > tmp3.tsv
awk '{gsub(/-nan/, "0"); print}' tmp3.tsv > ortholog_table_interdomain_lengths_averages.tsv
rm tmp*.tsv


echo "Step 9: Looking up the lengths of proteins for which a pfam domain was predicted"
awk '{print $1}' < hmmer_pfam_output_tbl_trimmed.txt > tmp.txt
sed 's/^/>/g' < tmp.txt > tmp1.txt
sed 's/-//g' < sed_file.txt > sed_file_for_hmmr_full_protein_length.txt
sed -f sed_file_for_hmmr_full_protein_length.txt < tmp1.txt > tmp2.txt
paste tmp1.txt tmp2.txt > temp_concat4.txt
sed 's/$/£/g' < temp_concat4.txt > temp_concat5.txt
sed 's/\(>.*\)\t\([0-9]*\)£/s#\1#\2£/'g < temp_concat5.txt > temp_concat6.txt
tr -d '\n' < temp_concat6.txt > temp_concat7.txt
sed 's/£/\n/g' < temp_concat7.txt > temp_concat8.txt
sed 's/$/#/g' < temp_concat8.txt > sed_file_full_protein_length.txt

echo ">>>Step 9.2: Creating ortholog_table_full_protein_lengths.tsv"
sed -f sed_file_full_protein_length.txt < ortholog_table_for_hmmer_output_parsing.tsv > ortholog_table_full_protein_lengths.tsv
rm temp_concat*
echo "removing orthologs from ortholog_table_full_protein_lengths.tsv that did not have a hmmer pfam output"
input_file="list_of_proteins_from_ortholog_table.tsv"
output_file="sed_commands_full_protein_lengths.sh"
# Initialize counter to track the temp files
counter=1
> "output_file"
echo "#!/bin/bash" > "$output_file"
# Loop through each line of the input file
while IFS= read -r line; do
    # Extract just the identifier part without '>' and '##'
#    identifier=$(echo "$line" | sed 's/[>#]//g')

    # Check if it's the first line (to use ortholog_table_full_protein_lengths.tsv as input)
    if [ $counter -eq 1 ]; then
        echo "sed 's/$line//' < ortholog_table_full_protein_lengths.tsv > temp1" >> "$output_file"
    else
        previous_counter=$((counter - 1))
        echo "sed 's/$line//' < temp$previous_counter > temp$counter" >> "$output_file"
    fi

    # Increment the counter
    counter=$((counter + 1))

done < "$input_file"
final_counter=$((counter - 1))
echo "mv temp$final_counter final_sed_file_full_protein_lengths.tsv" >> "$output_file"
chmod +x "$output_file"
./"$output_file"
find . -maxdepth 1 -name "temp*" -print0|xargs -0 rm
echo "calculating average values of full protein lengths"
awk '{gsub(/, ,/, ","); print}' final_sed_file_full_protein_lengths.tsv > tmp1.tsv
awk '{gsub(/, \t|\t /, "\t"); print}' tmp1.tsv > tmp2.tsv
num_cols=$(head -1 tmp2.tsv | awk -F'\t' '{print NF}')
for i in $(seq 1 $num_cols); do
    cut -f$i tmp2.tsv > column$i.txt
done
for filename in ./column*.txt
do
awk -F',' 'NR==1 {print $0} NR>1 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > "$filename".average
#awk -F',' 'NR==1 {print $0} NR>1 && NF>0 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > temp && mv temp "$filename".average
done
paste column*.average > tmp3.tsv
awk '{gsub(/-nan/, "0"); print}' tmp3.tsv > ortholog_table_full_protein_lengths_averages.tsv
rm tmp*.tsv
rm column*.txt

#!/bin/sh
#run this script first
#this script lists the lengths of proteins in an orthologous family
#Note that you would have to manually create Homo.txt  and place it in  the same folder as this script before you  run it. It is just a list of your human proteins of interest formated this way >DPOD2_HUMAN@.
#written by Dominic Wiredu Boakye
for filename in ./*.tsv
do
echo ">>>EXTRACTING HOMOLOGS FOR "$filename"<<<"
	fgrep -f ../the_list.txt < "$filename" > "$filename".DDR
echo ">>>CONVERTING TSV FILES for "$filename" TO SED FILES<<<"
	sed 's/\t/#/g' < "$filename.DDR" > "$filename".1
	sed 's/OG......./s/' < "$filename".1 > "$filename".2
	sed "s/$(printf '\r')\$//" < "$filename".2 > "$filename".3
	sed 's/$/#/' < "$filename".3 > "$filename".4
	sed 's/#\([0-9A-Z]*_HUMAN\), .*_HUMAN\(#.*#\)/#\t\1\2/' < "$filename".4 > "$filename".5
	sed 's/^s#\([0-9A-Z]*_HUMAN.*\)/s#\t\1/' < "$filename".5 > "$filename".6
	sed 's/, /@, >/'g < "$filename".6 > "$filename".7
	sed 's/HUMAN#/HUMAN#>/'g < "$filename".7 > "$filename".8
	sed 's/#$/@#/'g < "$filename".8 > "$filename".9
echo ">>>CREATING ORTHOLOG SPREADSHEET<<<"
	echo ">>>ADDING ORTHOLOGS FROM "$filename" ORTHOLOG SPREADSHEET<<<"
	sed -f "$filename".9 < ../the_list.txt > "$filename".10
	awk 'BEGIN{print "'$filename'"}1' < "$filename".10 > "$filename".11
	sed "1s#\(.*__v__\)\(.*\).tsv#\2#" < "$filename".11 > "$filename".12
	sed 's/^\t.*/0/' < "$filename".12 > "$filename".13
	mv "$filename".13 "$(head -1 "$filename".13).txt"
	rm "$filename".*
done
echo ">>>Creating ortholog name table<<<"
paste Acanthamoeba.txt Albugo.txt Allomyces.txt Anopheles.txt Aphanomyces.txt Aspergillus.txt Babesia.txt Blastocystis.txt Blechomonas.txt Botrytis.txt Caenorhabditis.txt Carpediemonas.txt Chromera.txt Coprinopsis.txt Crithidia.txt Cryptococcus.txt Cryptosporidium.txt Cyclospora.txt Cytauxzoon.txt Danio.txt > temp1.txt
paste temp1.txt Drosophila.txt Eimeria.txt Encephalitozoon.txt Endotrypanum.txt Entamoeba.txt Enterospora.txt Eriocheir.txt Escherichia.txt Giardia.txt Glossina.txt Gregarina.txt Hyaloperonospora.txt Homo.txt Leishmania.txt Leptomonas.txt Magnaporthe.txt Mastigmoeba.txt Melampsora.txt Mitosporidium.txt Monocercomonoides.txt > temp2.txt
paste temp2.txt Naegleria.txt Nakaseomyces.txt Neospora.txt Nucleospora.txt Oncorhynchus.txt Oryza.txt Paramicrosporidium.txt Paramikrocytos.txt Phytophthora.txt Plasmodium.txt Pneumocystis.txt Pythium.txt Rhizopus.txt Rozella.txt Saccharomyces.txt Saprolegnia.txt Schizosaccharomyces.txt Spironucleus.txt Spizellomyces.txt > temp3.txt
paste temp3.txt Sporisorium.txt Theileria.txt Toxoplasma.txt Trichomonas.txt Trypanosoma.txt Ustilago.txt Vitrella.txt Xenopus.txt Zea.txt > ortholog_table.tsv
#echo ">>>joining all fasta protein files used in analyses from various species into one file<<<"
#cat *fasta > combined.fast 
#echo ">>>linearise combined.fasta<<<"
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < combined.fast > combined.fa
#sed 's/\(>.*\)/\1@/' < combined.fa > combined.fas1
#sed 's/:/_/'g < combined.fas1 > combined.fas
for filename in ./ortholog_table.tsv
do
echo ">>>deleting first line from "$filename"<<<"
sed '1d' < "$filename" > "$filename".1
echo ">>>replacing tabs with new line from "$filename"<<<"
tr '\t' '\n' < "$filename".1 > "$filename".2
tr ' ' '\n' < "$filename".2 > "$filename".3
echo ">>>deleting commas from "$filename"<<<"
sed 's/,//'g < "$filename".3 > "$filename".4 
echo ">>>deleting lines with 0s from "$filename"<<<"
sed 's/^0$//'g < "$filename".4 > "$filename".5
echo ">>>appending # to end of line in "$filename"<<<"
sed '/^$/d' < "$filename".5 > "$filename".6
sort "$filename".6 > "$filename".7
echo ">>>split grep file into several grep files! Like a lot of grep files!!"$filename"<<<"
split -l 100 "$filename".7 splitgrepfile.
done
echo ">>>THE NEXT STEPS COULD TAKE VERY LONG TO COMPLETE. SAY 8 HOURS IF YOUR ORIGINAL FASTA FILE WAS APPROX 500 MB<<<"
for filename in splitgrepfile.*
do
echo ">>>Grepping from tiny "$filename" grep file<<<"
fgrep -A 1 -f "$filename" < combined.fas > fa_"$filename"
done
echo ">>>concatenating fasta files from tiny grep files<<<"
cat fa_splitgrepfile.* > combinedDHHorthologs.fa
echo ">>>THE LONG GREPPING PROCESS IS NOW COMPLETE<<<"
echo ">>>getting protein lengths<<<"
sed 's/^--.*//'g < combinedDHHorthologs.fa > temp_length.txt
sed '/^$/d' < temp_length.txt > temp_length1.txt
awk '{print length($0);}' temp_length1.txt  > temp_length2.txt
paste temp_length1.txt temp_length2.txt > temp_length3.txt
echo ">>>creating sed file from temp_length3.txt<<<"
sed 's/^[A-Za-z]*\t\([0-9]*\)/\1£/' < temp_length3.txt > temp_length4.txt
sed 's/\(>.*\)@\t[0-9]*/s#\1@#/'g < temp_length4.txt > temp_length5.txt
tr -d '\n' < temp_length5.txt > temp_length6.txt
tr '£' '\n' < temp_length6.txt > temp_length7.txt
sed 's/$/#/'g < temp_length7.txt > sed_file.txt
echo ">>>Substituting protein names in ortholog_table.tsv with their respective lengths<<<"
sed -f sed_file.txt < ortholog_table.tsv > ortholog_table_protein_lengths.tsv
num_cols=$(head -1 ortholog_table_protein_lengths.tsv | awk -F'\t' '{print NF}')
for i in $(seq 1 $num_cols); do
    cut -f$i ortholog_table_protein_lengths.tsv > column$i.txt
done
for filename in ./column*.txt
do
awk -F',' 'NR==1 {print $0} NR>1 {sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum/NF}' < "$filename" > "$filename".average
done
paste column*.average > ortholog_table_protein_lengths_averages.tsv
#mkdir ExtractedFastaSeqs
#mv fa_* ExtractedFastaSeqs
#mkdir splitgrepfileDir
#mv splitgrepfile.* splitgrepfileDir
rm splitgrepfile.*
rm fa_splitgrepfile.*
rm temp_length*
rm combined.fast
rm combined.fa
rm column*.txt.average
rm temp*.txt
rm ortholog_table.tsv.*
rm column*.txt
sed 's/-//g' < combinedDHHorthologs.fa > combinedDHHorthologs.fas

#!/bin/bash

cut -d ' ' -f 1 < DDR_groupsV2.tsv > BER.output
cut -d ' ' -f 2 < DDR_groupsV2.tsv > MMR.output
cut -d ' ' -f 3 < DDR_groupsV2.tsv > NER.output
cut -d ' ' -f 4 < DDR_groupsV2.tsv > OTHER_SSR.output
cut -d ' ' -f 5 < DDR_groupsV2.tsv > FA.output
cut -d ' ' -f 6 < DDR_groupsV2.tsv > HR.output
cut -d ' ' -f 7 < DDR_groupsV2.tsv > NHEJ.output
cut -d ' ' -f 8 < DDR_groupsV2.tsv > OTHER_DSR.output
cut -d ' ' -f 9 < DDR_groupsV2.tsv > CHECKPOINT_FACTORS.output
cut -d ' ' -f 10 < DDR_groupsV2.tsv > CHROMATIN_REMODELLING.output
cut -d ' ' -f 11 < DDR_groupsV2.tsv > CHROMOSOME_SEGREGATION.output
cut -d ' ' -f 12 < DDR_groupsV2.tsv > DNA_REPLICATION.output
cut -d ' ' -f 13 < DDR_groupsV2.tsv > MODULATION_OF_NUCLEOTIDE_POOLS.output
cut -d ' ' -f 14 < DDR_groupsV2.tsv > P53_PATHWAY.output
cut -d ' ' -f 15 < DDR_groupsV2.tsv > TELOMERE_MAINTENANCE.output
cut -d ' ' -f 16 < DDR_groupsV2.tsv > TLS.output
cut -d ' ' -f 17 < DDR_groupsV2.tsv > TOPOISOMERASE_DAMAGE_REVERSAL.output
cut -d ' ' -f 18 < DDR_groupsV2.tsv > UBIQUITIN_RESPONSE.output
cut -d ' ' -f 19 < DDR_groupsV2.tsv > PROTEINS_WITH_PROBABLE_DDR_ROLE.output
cut -d ' ' -f 20 < DDR_groupsV2.tsv > DNA_REPLICATION.output
cut -d ' ' -f 21 < DDR_groupsV2.tsv > METABOLISM.output
cut -d ' ' -f 22 < DDR_groupsV2.tsv > HOUSEKEEPING.output
for filename in  *.output
do
sed 's/\r$//' "$filename" > "$filename".removed_windows_nonsense
sed '/^$/d' "$filename".removed_windows_nonsense > "$filename".clean
fgrep -f "$filename".clean < ortholog_table_protein_lengths_averages_for_grepping.tsv > "$filename".clean.grepped
sed '/^$/d' "$filename".clean.grepped > "$filename".clean.grepped1
tr -d '#' < "$filename".clean.grepped1 > "$filename".clean.grepped.tsv
mv "$filename".clean.grepped.tsv "$filename".tsv
done
rm *.clean
rm *.removed_windows_nonsense
rm *.clean.grepped
rm *.clean.grepped1
rm *.output

#!/usr/bin/env bash

DATADIR=/home/oppken001/PhD2019/EnrichmentAnalysis/GO

for file in $DATADIR/*; do
  UMLS_ID=$(echo $file | rev | cut -d'/' -f1 |rev);
  # echo $UMLS_ID;
  for ofile in $DATADIR/${UMLS_ID}/*.txt; do
   OUTFILE=$(readlink -e $ofile)
   # cut -f5,9 $OUTFILE | head -10
   awk -F'\t' -v awk_umls=$UMLS_ID 'FNR>1 match($2, /GO:[0-9]+/) \
   { print substr($2, RSTART, RLENGTH),$5,$8,$9, awk_umls;  OFS="\t"};' $OUTFILE \
   | awk -F"\t" '{ if ($2 < 0.05) { print } }' >> enriched_go_diseases.txt
  # | awk 'match($2, /GO:[0-9]+/){ print substr($2, RSTART, RLENGTH), $5, $8, $9 }' \
  # >> enriched_go_diseases.txt
 done

done

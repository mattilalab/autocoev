#!/bin/bash

# Calculate background for each protein pair (Prot 1 -> Prot 2) after our filtering
calc_back_final() {
local PROTCHI="${1}"
echo "$PROTCHI"

cat $PROTCHI | while read -r Seq1 Seq2 numPairs totalComp ; do
  back_calc=$(printf "%1.10f" `echo "($numPairs)/($totalComp)" |bc -l`)
  echo "$Seq1 $Seq2 $back_calc" >> $TMP/$RESULTS/chi/back_calc_final/${PROTCHI%.*}.txt
done

file_back=$(printf "%1.10f" `awk '{print $3}' $TMP/$RESULTS/chi/back_calc_final/${PROTCHI%.*}.txt | datamash mean 1`)
echo "${PROTCHI%.*} $file_back" >> $TMP/$RESULTS/chi/back_calc_final.tsv

cat $PROTCHI | while read -r Seq1 Seq2 numPairs totalComp ; do
  back_calc=$(printf "%1.10f" `echo "($numPairs)/($totalComp)" |bc -l`)
  exp=$(printf "%1.10f" `echo "($totalComp)*($file_back)" |bc -l`)
  if (( $(echo "$exp > 0" |bc -l) )); then
    echo "[COEVOL] $PROTCHI"
    chi=$(printf "%1.10f" `echo "($numPairs - $exp)^2/($exp)" |bc -l`)
    if (( $(echo "$numPairs > $exp" |bc -l) && $(echo "$numPairs >= 3.84146" |bc -l) )); then
      echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.5" >> $TMP/$RESULTS/chi/chi_final.tsv
      echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.5" >> $TMP/$RESULTS/chi/chi_test_final/${PROTCHI%.*}.tsv
    else
      echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.0" >> $TMP/$RESULTS/chi/chi_test_final/${PROTCHI%.*}.tsv
      echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.0" >> $TMP/$RESULTS/chi/chi_final.tsv
    fi
  else
    echo "[NOCOEV] $PROTCHI"
    echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp nan 0.0" >> $TMP/$RESULTS/chi/chi_test_final/${PROTCHI%.*}.tsv
  fi
done
}

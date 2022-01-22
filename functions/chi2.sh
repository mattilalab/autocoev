#!/bin/bash

# These are needed when the chi^2 is done
#coev_inter_collect() {
# while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
#   echo "${Seq1%.*} ${Seq2%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1%.*}.tsv
#   echo "${Seq1%.*}"
#   echo "${Seq2%.*} ${Seq1%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2%.*}.tsv
#   echo "${Seq2%.*}"
# done < $TMP/$RESULTS/coev_inter_nocoev.tsv
#
# # These are needed when the chi^2 is done after our filtering
# while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
#   echo "${Seq1%.*} ${Seq2%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteinsFinal/${Seq1%.*}.tsv
#   echo "${Seq1%.*}"
#   echo "${Seq2%.*} ${Seq1%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteinsFinal/${Seq2%.*}.tsv
#   echo "${Seq2%.*}"
# done < $TMP/$RESULTS/coev_inter_nocoev.tsv
#
# # These are needed when the chi^2 is done
# while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
#   echo "${Seq1%.*} ${Seq2%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1%.*}.tsv
#   echo "${Seq1%.*}"
#   echo "${Seq2%.*} ${Seq1%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2%.*}.tsv
#   echo "${Seq2%.*}"
# done < $TMP/$RESULTS/coev_inter_coev.tsv
# }

# # Calculate background for each protein pair (Prot 1 -> Prot 2)
# calc_back() {
# local PROTCHI="${1}"
# echo "$PROTCHI"
#
# cat $PROTCHI | while read -r Seq1 Seq2 numPairs totalComp ; do
#   back_calc=$(printf "%1.10f" `echo "($numPairs)/($totalComp)" |bc -l`)
#   echo "$Seq1 $Seq2 $back_calc" >> $TMP/$RESULTS/chi/back_calc/${PROTCHI%.*}.txt
# done
#
# file_back=$(printf "%1.10f" `awk '{print $3}' $TMP/$RESULTS/chi/back_calc/${PROTCHI%.*}.txt | datamash mean 1`)
# echo "${PROTCHI%.*} $file_back" >> $TMP/$RESULTS/chi/file_back.tsv
#
# cat $PROTCHI | while read -r Seq1 Seq2 numPairs totalComp ; do
#   back_calc=$(printf "%1.10f" `echo "($numPairs)/($totalComp)" |bc -l`)
#   exp=$(printf "%1.10f" `echo "($totalComp)*($file_back)" |bc -l`)
#   if (( $(echo "$exp > 0" |bc -l) )); then
#     echo "[COEVOL] $PROTCHI"
#     chi=$(printf "%1.10f" `echo "($numPairs - $exp)^2/($exp)" |bc -l`)
#     if (( $(echo "$numPairs > $exp" |bc -l) && $(echo "$numPairs >= 3.84146" |bc -l) )); then
#       echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.5" >> $TMP/$RESULTS/chi/chi.tsv
#       echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.5" >> $TMP/$RESULTS/chi/chi_test/${PROTCHI%.*}.tsv
#     else
#       echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.0" >> $TMP/$RESULTS/chi/chi_test/${PROTCHI%.*}.tsv
#       echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp $chi 0.0" >> $TMP/$RESULTS/chi/chi.tsv
#     fi
#   else
#     echo "[NOCOEV] $PROTCHI"
#     echo "${Seq1%.*} ${Seq2%.*} $back_calc $file_back $exp nan 0.0" >> $TMP/$RESULTS/chi/chi_test/${PROTCHI%.*}.tsv
#   fi
# done
# }

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

# # These two functions create summary files of the chi^2 (before our filtering). No longer used.
# coev_inter_chi_results() {
# while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
#   if [ -f "$TMP/$RESULTS/chi/chi_test/${Seq1%.*}.tsv" ] && [ -f "$TMP/$RESULTS/chi/chi_test/${Seq2%.*}.tsv" ] ; then
#     forward=$(grep "${Seq1%.*} ${Seq2%.*}" $TMP/$RESULTS/chi/chi_test/${Seq1%.*}.tsv | awk '{print $7}')
#     reverse=$(grep "${Seq2%.*} ${Seq1%.*}" $TMP/$RESULTS/chi/chi_test/${Seq2%.*}.tsv | awk '{print $7}')
#     chiboth=$(echo "$forward + $reverse" |bc -l)
#     echo "$Seq1 $Seq2 $numPairs $totalComp $CutOff $thresholdR $averageR $averageSigR $tree1length $tree2length $gapThreshold $bootCutOff $DistanceCoef $chiboth" >> $TMP/$RESULTS/coev_inter_chi.tsv
#     echo "${Seq1%.*} ${Seq2%.*} added after Chi test"
#   else
#     echo "${Seq1%.*} ${Seq2%.*} not found"
#   fi
# done < $TMP/$RESULTS/coev_inter_coev.tsv
# }
#
# coev_inter_chi_results_final() {
# while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
#   if [ -f "$TMP/$RESULTS/chi/chi_test_final/${Seq1%.*}.tsv" ] && [ -f "$TMP/$RESULTS/chi/chi_test_final/${Seq2%.*}.tsv" ] ; then
#     forward=$(grep "${Seq1%.*} ${Seq2%.*}" $TMP/$RESULTS/chi/chi_test_final/${Seq1%.*}.tsv | awk '{print $7}')
#     reverse=$(grep "${Seq2%.*} ${Seq1%.*}" $TMP/$RESULTS/chi/chi_test_final/${Seq2%.*}.tsv | awk '{print $7}')
#     chiboth=$(echo "$forward + $reverse" |bc -l)
#     echo "$Seq1 $Seq2 $numPairs $totalComp $CutOff $thresholdR $averageR $averageSigR $tree1length $tree2length $gapThreshold $bootCutOff $DistanceCoef $chiboth" >> $TMP/$RESULTS/coev_inter_chi_final.tsv
#     echo "${Seq1%.*} ${Seq2%.*} added after Chi test"
#   else
#     echo "${Seq1%.*} ${Seq2%.*} not found"
#   fi
# done < $TMP/$RESULTS/coev_inter_coev.tsv
# }

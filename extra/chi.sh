#!/bin/bash

# Export all
set -a

# What date and time is it?
DATESTAMP=$(date)

# Current work dir
CWD=$(pwd)

# Load settings first
. $CWD/settings.conf

# Load functions
. $CWD/functions/databases.sh
. $CWD/functions/retrieval.sh
. $CWD/functions/blast.sh
. $CWD/functions/msa.sh
. $CWD/functions/trees.sh
. $CWD/functions/pairing.sh
. $CWD/functions/caps.sh
. $CWD/functions/results.sh
. $CWD/functions/network.sh
. $CWD/functions/check.sh

mkdir -p $TMP/$RESULTS/chi/proteins

#while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
#  echo "${Seq1%.*} ${Seq2%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1%.*}.tsv
#  echo "${Seq2%.*} ${Seq1%.*} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2%.*}.tsv
#done < $TMP/$RESULTS/coev_inter.tsv

prot_chi() {
  local PROTCHI="${1}"
  declare -A file_back_arr
  cat $PROTCHI | while read -r Seq1 Seq2 numPairs totalComp ; do
  file_back=$(echo "($numPairs)/($totalComp)" |bc -l)
  file_back_arr+=($file_back)
}

mkdir -p $TMP/$RESULTS/chi/chi_pass
cd $TMP/$RESULTS/chi/proteins
protInt=$( ls ./ )
parallel $CORESCAPS prot_chi ::: "$protInt"
cd ..

#!/bin/bash

# Export all
set -a

# Current work dir
CWD=$(pwd)

# Specify the input file
echo -e "\n\e[96mEnter the full path to the pairs file:\e[39m"
read ETSV

# Define paths and files
echo -e "\n\e[96mEnter the name of the output file:\e[39m"
read EOUT

sed 1d "$ETSV" | while read -r Name1 msa1 Name2 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot corrT cCoev pvalA pvalB pMean pDiff corr1 corr2 bonferroni GapA GapB GapAB DivA DivB DivAB ; do
  echo ${msa1}-${realA}_${msa2}-${realB} $Name1 $msa1 $Name2 $msa2 $colA $realA $colB $realB $seq1 $seq2 $flt1 $flt2 $meanA $meanB $corr $boot $corrT $cCoev $pvalA $pvalB $pMean $pDiff $corr1 $corr2 $bonferroni $GapA $GapB $GapAB $DivA $DivB $DivAB >> $EOUT
  echo "Process ${msa1}-${realA} ${msa2}-${realB}"
done

sed -i "1i msa1-realA_msa2-realB Name1 msa1 Name2 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot corrT cCoev pvalA pvalB pMean pDiff corr1 corr2 bonferroni GapA GapB GapAB DivA DivB DivAB" $EOUT

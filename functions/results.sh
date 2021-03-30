#!/bin/bash

# Inspect CAPS results. Separate pair folders where CAPS run failed,
# skip pair folders where no coevolution was detected and copy those
# where co-evolving pairs were found. To do this, first check if row
# 2 in coev_inter.csv is empty (failed pairs), then check if its third
# column value (coevolving sites) is equal to zero (no coevolution) or
# larger than zero (coevolution detected).
caps_inspect() {
local pair="${1}"
  SUMMARY=$(sed -n '2p' $pair/coev_inter.csv)
  while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
    if [ -z "${SUMMARY}" ]; then
      echo -e "[\e[91mFAILED\e[39m] Copying pair where CAPS failed: $pair"
      cp -a $pair $TMP/$RESULTS/fail
    elif [ "$numPairs" -eq 0 ]; then
      echo -e "[\e[34mNOCOEV\e[39m] No co-evolving pairs found for: $pair"
      cp -a $pair $TMP/$RESULTS/nocoev
    elif [ "$numPairs" -gt 0 ]; then
      echo -e "[\e[92mCOEVOL\e[39m] Copying pair with co-evolution: $pair"
      mkdir -p $TMP/$RESULTS/coev/$folder
      cp -a $pair $TMP/$RESULTS/coev/$folder
    else
      echo -e "Something went wrong for $pair ... Check!"
    fi
  done <<< $(echo "$SUMMARY")
}

# Cleanup results to make them easily parsable. Leave protein names
# and their UniProt identifiers in coev_inter.tsv. For *.out files,
# leave only the info of residues positions of the co-evolving pairs.
results_cleanup() {
  local coevPair="${1}"
  SUMMARY=$(sed -n '2p' $coevPair/coev_inter.csv)
  while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
    cp $coevPair/${Seq1}_${Seq2}.out $coevPair/$coevPair.clean
    echo "Clean up $coevPair.out"
    sed -i -n '/Coevolving Pairs of amino acid sites/,/Overlapping groups of coevolving residues/p' $coevPair/$coevPair.clean
    sed -i '1,5d' $coevPair/$coevPair.clean
    sed -i 's/Overlapping groups of coevolving residues//' $coevPair/$coevPair.clean
    sed -i '/^$/d' $coevPair/$coevPair.clean
    sed -i "s:\t\t:\t:g" $coevPair/$coevPair.clean
    sed -i "s/(/ /g" $coevPair/$coevPair.clean
    sed -i "s/)/ /g" $coevPair/$coevPair.clean
    sed -i "s/	/ /g" $coevPair/$coevPair.clean
    sed -i "s/  / /g" $coevPair/$coevPair.clean
    sed -i "s/^/$Seq1 $Seq2 /" $coevPair/$coevPair.clean
    sed -i "s/\.fa/ /g" $coevPair/$coevPair.clean
    sed -i "1i msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2" $coevPair/$coevPair.clean
  done <<< $(echo "$SUMMARY")
}

adj_pVal(){
  local coevPair="${1}"
  sed 1d $coevPair/$coevPair.clean | while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 ; do
     corrdec=$(printf "%1.10f" $corr)
    corr1dec=$(printf "%1.10f" $corr1)
    corr2dec=$(printf "%1.10f" $corr2)
    pvalAdec=$(printf "%1.10f" $pvalA)
    pvalBdec=$(printf "%1.10f" $pvalB)
    pMeandec=$(printf "%1.10f" $pMean)
    if (( $(echo "$pMeandec <= $PVALUE" |bc -l) && \
    	  $(echo "$corr1dec > 0" |bc -l) && \
	  $(echo "$corr2dec > 0" |bc -l) )); then
      echo "Adding $pMeandec"
      echo "$msa1 $msa2 $colA $realA $colB $realB $meanA $meanB $corrdec $boot $pvalAdec $pvalBdec $pMeandec $corr1dec $corr2dec" >> $coevPair/$coevPair.$PVALUE.clean
    else
      echo "Skipping $pMeandec"
    fi
  done
    
  # Run Luqman's script
  if [ -s $coevPair/$coevPair.$PVALUE.clean ]; then
    sed -i "1i msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2" $coevPair/$coevPair.$PVALUE.clean
    echo -e "Running R for Bonferroni correction for $coevPair"
    Rscript $CWD/R/AdjPval.R $coevPair/$coevPair.$PVALUE.clean $coevPair/$coevPair.$PVALUE.corrected
  elif [ -z $coevPair/$coevPair.$PVALUE.clean ]; then
    echo -e "[SKIP] $coevPair.$PVALUE.clean not found. Skipping..."
  else
    echo -e "Check your input!"
  fi
}

# Extract columns for pairs that passed the Bonferroni correction test.
# Collect only pairs that exist in the reference organism (with "real" position > 0)
extract_columns(){
  local coevPair="${1}"
  sed 1d $coevPair/$coevPair.$PVALUE.corrected |
  while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 bonferroni ; do
    if (( $(echo "$bonferroni <= $BONFERRONI" |bc -l) && $(echo "$realA > 0" |bc -l) && $(echo "$realB > 0" |bc -l) )); then
      seq1=$( sed "${realA}q;d" $TMP/$MSA/$MSAMETHOD/${msa1}.fa.${ORGANISM}.ref | awk '{print $1$2}' )
      seq2=$( sed "${realB}q;d" $TMP/$MSA/$MSAMETHOD/${msa2}.fa.${ORGANISM}.ref | awk '{print $1$2}' )
      flt1=$( sed "${realA}q;d" $TMP/$MSA/$MSAMETHOD/${msa1}.fa.${ORGANISM}.ref | awk '{print $3}' )
      flt2=$( sed "${realB}q;d" $TMP/$MSA/$MSAMETHOD/${msa2}.fa.${ORGANISM}.ref | awk '{print $3}' )
      echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $flt1 $flt2 $meanA $meanB $corr $boot $pvalA $pvalB $pMean $corr1 $corr2 $bonferroni" >> $coevPair/$coevPair.$PVALUE.corrected-$BONFERRONI
      echo -e "Extracting $colA-$colB with corrected p-value=${bonferroni}"
      mkdir -p $coevPair/columnStats-$BONFERRONI/$colA-$colB/
      cat $coevPair/msa/$msa1.fa | seqkit subseq -r ${colA}:${colA} | seqkit sort -o $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colA-$msa1.txt
      cat $coevPair/msa/$msa2.fa | seqkit subseq -r ${colB}:${colB} | seqkit sort -o $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colB-$msa2.txt
     elif (( $(echo "$bonferroni > $BONFERRONI" |bc -l) )); then
       echo -e "Skipping $colA-$colB with corrected p-value=${bonferroni}"
     elif (( $(echo "$realA == 0" |bc -l) )); then
       echo -e "Skipping $colA-$colB with MSA $colA -> $realA in $ORGANISM"
     elif (( $(echo "$realB == 0" |bc -l) )); then
       echo -e "Skipping $colA-$colB with MSA $colB -> $realB in $ORGANISM"
     else
       echo "Check settings!"
     fi
  done
}

columns_stats(){
  local coevPair="${1}"
  if [ -s "$coevPair/${coevPair}.${PVALUE}.corrected-${BONFERRONI}" ]; then
    echo -e "Processing ${coevPair}.${PVALUE}.corrected-${BONFERRONI}"
    cat $coevPair/${coevPair}.${PVALUE}.corrected-${BONFERRONI} |
    while read -r msa1 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot pvalA pvalB pMean corr1 corr2 bonferroni ; do
      columnTotalA=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colA-$msa1.txt | wc -l)
       columnGapsA=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colA-$msa1.txt | grep "-"  | wc -l)
       percentageA=$(echo "1 - ${columnGapsA}/${columnTotalA}" | bc -l)
      roundPerGapA=$(printf "%1.5f" $percentageA)

      columnNoGapsA=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colA-$msa1.txt | grep -v "-" | wc -l)
      columnUniqueA=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colA-$msa1.txt | grep -v "-" | sort | uniq -c | sort -n -r | sed -n 1p | awk '{ print $1 }')
      diversityResA=$(echo "1 - ${columnUniqueA}/${columnNoGapsA}" | bc -l)
       roundDivResA=$(printf "%1.5f" $diversityResA)

      columnTotalB=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colB-$msa2.txt | wc -l)
       columnGapsB=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colB-$msa2.txt | grep "-" | wc -l)
       percentageB=$(echo "1 - ${columnGapsB}/${columnTotalB}" | bc -l)
      roundPerGapB=$(printf "%1.5f" $percentageB)

      columnNoGapsB=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colB-$msa2.txt | grep -v "-" | wc -l)
      columnUniqueB=$(sed '/^>/d' $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colB-$msa2.txt | grep -v "-" | sort | uniq -c | sort -n -r | sed -n 1p | awk '{ print $1 }')
      diversityResB=$(echo "1 - ${columnUniqueB}/${columnNoGapsB}" | bc -l)
       roundDivResB=$(printf "%1.5f" $diversityResB)

      # Calculate gaps and diversity    
      GapsAB=$(printf "${roundPerGapA}\n $roundPerGapB" | datamash mean 1 )
      DivsAB=$(printf "${roundDivResA}\n $roundDivResB" | datamash mean 1 )
	 
      # Calculate difference between raw p-values (absolute numbers)
      diffP=$(echo "${pvalA} - ${pvalB}" |bc -l | tr -d -)
      pDiff=$(printf "%1.10f" $diffP)
	
      # Include correlation threshold as well. Just in case, make sure value is not in scientific format
      thrsR=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $6}')
      corrT=$(printf "%1.10f" $thrsR)
	
      # Calculate coevolutionary correlation, normalized to threshold. For normalization:
      # Thanks to Dian Dimitrov.
      # https://www.mathworks.com/matlabcentral/answers/322438-normalize-data-with-a-threshold
      normC=$(echo "($corr - $corrT)/(1 - $corrT)" |bc -l)
      cCoev=$(printf "%1.10f" $normC)
      
      echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $flt1 $flt2 $meanA $meanB $corr $boot $corrT $cCoev $pvalA $pvalB $pMean $pDiff $corr1 $corr2 $bonferroni $roundPerGapA $roundPerGapB $GapsAB $roundDivResA $roundDivResB $DivsAB" >> $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv

    done
  else
    echo "Skipping $coevPair.$PVALUE.corrected"
  fi
}

summary_cleanup(){
  if [ -s $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv ]; then
    sed 1d $TMP/tsv/proteinsFound.tsv | while read -r idxml namexml attribute ; do
      sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv
    done
    sed -i \
    "1i Name1 msa1 Name2 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot corrT cCoev pvalA pvalB pMean pDiff corr1 corr2 bonferroni GapA GapB GapAB DivA DivB DivAB" \
    $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv
  else
    echo "No pairs-P${PVALUE}-B${BONFERRONI}.tsv !"
  fi
}

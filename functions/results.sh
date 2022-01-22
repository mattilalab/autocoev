#!/bin/bash

# Define paths and files
#ETSV="$TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv"
EOUT="$TMP/$RESULTS/filtered-P${PVALUE}-gaps${EGAPS}-pDiff${EPDIF}-Bonf${BONFERRONI}.tsv"

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
      
      # Export these, since they are needed for the Chi^2 test later
      mkdir -p $TMP/$RESULTS/chi/proteins/
      echo "${Seq1} ${Seq2} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1}.tsv
      echo "${Seq2} ${Seq1} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2}.tsv
      
      mkdir -p $TMP/$RESULTS/nocoev/$folder
      cp -a $pair $TMP/$RESULTS/nocoev/$folder
    elif [ "$numPairs" -gt 0 ]; then
      echo -e "[\e[92mCOEVOL\e[39m] Copying pair with co-evolution: $pair"
      mkdir -p $TMP/$RESULTS/coev/$folder
      cp -a $pair $TMP/$RESULTS/coev/$folder
     
      # Prepare the "reversed" MSA, so we run CAPS the other way round
      mv $TMP/$RESULTS/coev/$folder/$pair/coev_inter.csv $TMP/$RESULTS/coev/$folder/$pair/${Seq1}_${Seq2}-coev_inter.csv
      cd $TMP/$RESULTS/coev/$folder/$pair/msa
      firstMSA=$(ls | sed -n '1p')
     secondMSA=$(ls | sed -n '2p')
     mkdir -p ../msa-rev
     cp $firstMSA ../msa-rev/b_${firstMSA}
     cp $secondMSA ../msa-rev/a_${secondMSA}
    
    else
      echo -e "Something went wrong for $pair ... Check!"
    fi
  done <<< $(echo "$SUMMARY")
}

# Inspect the results from the "reversed" analysis
caps_reinspect() {
local pair="${1}"
  SUMMARY=$(sed -n '2p' $pair/coev_inter.csv)
  while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
    if [ -z "${SUMMARY}" ]; then
      echo -e "[\e[91mFAILED\e[39m] Copying pair where CAPS failed: $pair"
      mv $pair $TMP/$RESULTS/fail
    elif [ "$numPairs" -eq 0 ]; then
      echo -e "[\e[34mNOCOEV\e[39m] No co-evolving pairs found for: $pair"
      echo "$Seq1 $Seq2 $numPairs $totalComp $CutOff $thresholdR $averageR $averageSigR $tree1length $tree2length $gapThreshold $bootCutOff $DistanceCoef" >> $TMP/$RESULTS/coev_inter_nocoev.tsv
      mkdir -p $TMP/$RESULTS/nocoev/$folder
      mv $pair $TMP/$RESULTS/nocoev/$folder
      
      # These are needed for the Chi^2 test
      echo "${Seq1#*_} ${Seq2#*_} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1#*_}.tsv
      echo "${Seq2#*_} ${Seq1#*_} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2#*_}.tsv
      
     elif [ "$numPairs" -gt 0 ]; then
       echo -e "[\e[92mCOEVOL\e[39m]: $pair"
       # https://stackoverflow.com/a/15149278
       cd $pair
       mv coev_inter.csv ${Seq1#*_}_${Seq2#*_}-coev_inter.csv  
       mv ${Seq1}_${Seq2}.out ${Seq1#*_}_${Seq2#*_}.out
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
  cd $coevPair
  PROTEINONE=$(ls *.species | sed -n '1p')
  PROTEINTWO=$(ls *.species | sed -n '2p')
  
  # "Forward" Protein A vs Protein B
  cp ${PROTEINONE%.*}.fa_${PROTEINTWO%.*}.fa.out ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  echo "Clean up ${PROTEINONE%.*}_${PROTEINTWO%.*}.out"
  sed -i -n '/Coevolving Pairs of amino acid sites/,/Overlapping groups of coevolving residues/p' ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i '1,5d' ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i 's/Overlapping groups of coevolving residues//' ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i '/^$/d' ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s:\t\t:\t:g" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s/(/ /g" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s/)/ /g" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s/	/ /g" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s/  / /g" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s/^/${PROTEINONE%.*} ${PROTEINTWO%.*} /" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "s/\.fa/ /g" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  sed -i "1i msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  
  # "Reverse" Protein B vs Protein A
  cp ${PROTEINTWO%.*}.fa_${PROTEINONE%.*}.fa.out ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  echo "Clean up ${PROTEINTWO%.*}_${PROTEINONE%.*}.out"
  sed -i -n '/Coevolving Pairs of amino acid sites/,/Overlapping groups of coevolving residues/p' ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i '1,5d' ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i 's/Overlapping groups of coevolving residues//' ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i '/^$/d' ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s:\t\t:\t:g" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s/(/ /g" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s/)/ /g" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s/	/ /g" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s/  / /g" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s/^/${PROTEINTWO%.*} ${PROTEINONE%.*} /" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "s/\.fa/ /g" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  sed -i "1i msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  
  sed 1d ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean | while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 ; do
    if rowmatch=$(grep "$msa2 $msa1 $colB $realB $colA $realA" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean) ; then
      echo "Found a match bothways: $rowmatch"
    
      # Do decimal values for fwd. This is from the while read loop
       fcorrdec=$(printf "%1.10f" $corr)
      fcorr1dec=$(printf "%1.10f" $corr1)
      fcorr2dec=$(printf "%1.10f" $corr2)
      #fpvalAdec=$(printf "%1.10f" $pvalA)
      #fpvalBdec=$(printf "%1.10f" $pvalB)
      fpMeandec=$(printf "%1.10f" $pMean)
    
      # Define rev. This comes from the grep step
       rcorrdec=$(printf "%1.10f" `echo $rowmatch | awk '{print $9}'`)
      rcorr1dec=$(printf "%1.10f" `echo $rowmatch | awk '{print $14}'`)
      rcorr2dec=$(printf "%1.10f" `echo $rowmatch | awk '{print $15}'`)
      #rpvalAdec=$(printf "%1.10f" `echo $rowmatch | awk '{print $11}'`)
      #rpvalBdec=$(printf "%1.10f" `echo $rowmatch | awk '{print $12}'`)
      rpMeandec=$(printf "%1.10f" `echo $rowmatch | awk '{print $13}'`)
  
      # Make sure we use reliable values. And we do not need all this in the output (comment out for now)
      if (( $(echo "$fpMeandec <= $PVALUE" |bc -l) && $(echo "$rpMeandec <= $PVALUE" |bc -l) && \
    	  $(echo "$fcorr1dec > 0" |bc -l) && $(echo "$rcorr1dec > 0" |bc -l) && \
	  $(echo "$rcorr2dec > 0" |bc -l) && $(echo "$rcorr2dec > 0" |bc -l) )); then

          # Save the two correlation s (FWD and REV) mean value. The correlation is, in turn,
	  # estimated in two directions for each FWD and REV (e.g. $fcorr1dec & $fcorr2dec),
	  # but we do not want to output a crazy amount of stuff, do we? Same goes for the
	  # p-values (e.g. fpvalAdec & fpvalBdec).	  
	  corrdec=$(printf "${fcorrdec}\n $rcorrdec" | datamash mean 1)
          #corr1dec=$(printf "${fcorr1dec}\n $rcorr1dec" | datamash mean 1)
          #corr2dec=$(printf "${fcorr2dec}\n $rcorr2dec" | datamash mean 1)
          #pvalAdec=$(printf "${fpvalAdec}\n $rpvalAdec" | datamash mean 1)
          #pvalBdec=$(printf "${fpvalBdec}\n $rpvalBdec" | datamash mean 1)
          pMeandec=$(printf "${fpMeandec}\n $rpMeandec" | datamash mean 1)
      
        echo "Adding $pMeandec"
        
	#echo "$msa1 $msa2 $colA $realA $colB $realB $meanA $meanB $corrdec $boot $pvalAdec $pvalBdec $pMeandec $corr1dec $corr2dec" >> bothWays.tsv
        echo "$msa1 $msa2 $colA $realA $colB $realB $corrdec $boot $pMeandec" >> bothWays.tsv
	
       elif (( $(echo "$fpMeandec <= $PVALUE" |bc -l) || $(echo "$rpMeandec <= $PVALUE" |bc -l) || \
    	  $(echo "$fcorr1dec > 0" |bc -l) || $(echo "$rcorr1dec > 0" |bc -l) || \
	  $(echo "$rcorr2dec > 0" |bc -l) || $(echo "$rcorr2dec > 0" |bc -l) )); then
	  echo "Skipping $pMeandec"
       else
         echo "Something went wrong at $coevPair/${PROTEINONE%.*}_${PROTEINTWO%.*}.clean"
       fi
    else
      echo "Did not find a match bothways: $rowmatch"
    fi
    done
    
  #sed -i "1i msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2" bothWays.tsv
  sed -i "1i msa1 msa2 colA realA colB realB corr boot p_value" bothWays.tsv
  
  # How many pairs do we have left? This is needed for Chi^2 and this is
  # how CAPS2 counts pairs. Also get the totalcomp and export it.
  pairsNumber=$(sed 1d bothWays.tsv | awk '{print $3}' | datamash count 1)
  totCompares=$(sed -n 2p ${PROTEINONE%.*}.fa_${PROTEINTWO%.*}.fa-coev_inter.csv | awk '{print $4}')
  echo "${PROTEINONE%.*}.fa ${PROTEINTWO%.*}.fa $pairsNumber $totCompares" >> $TMP/$RESULTS/chi/proteins/${PROTEINONE%.*}.fa.tsv
  echo "${PROTEINTWO%.*}.fa ${PROTEINONE%.*}.fa $pairsNumber $totCompares" >> $TMP/$RESULTS/chi/proteins/${PROTEINTWO%.*}.fa.tsv
  
  # Run Luqman's script for Bonferroni & co correction of p-values
  echo -e "Running R for Bonferroni correction for $coevPair"
  Rscript $CWD/R/AdjPval.R bothWays.tsv bothWays-corrected.tsv
  cd ..
}

# get the amino acid from the reference organism's sequence and calculate
# column statistics
extract_columns_stats(){
  local coevPair="${1}"
  cd $coevPair
  sed 1d bothWays-corrected.tsv | while read -r msa1 msa2 colA realA colB realB corr boot p_value bonferroni holm bh hochberg hommel by fdr; do
    
    # make sure the amino acid position exists in the reference organism
    if (( $(echo "$realA > 0" |bc -l) && $(echo "$realB > 0" |bc -l) )); then
      
      # get amino acid and Gblocks score
      seq1=$( sed "${realA}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa1}.fa.${ORGANISM}.ref | awk '{print $1$2}' )
      seq2=$( sed "${realB}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa2}.fa.${ORGANISM}.ref | awk '{print $1$2}' )
      flt1=$( sed "${realA}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa1}.fa.${ORGANISM}.ref | awk '{print $3}' )
      flt2=$( sed "${realB}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa2}.fa.${ORGANISM}.ref | awk '{print $3}' )

      # Convert Gblocks scores to numbers (. = 0 ; # = 1)
      if [ "$flt1" = "#" ]; then
        gblscore1="1"
      elif [ "$flt1" = "." ]; then
        gblscore1="0"
      else
        echo "Check Gblocks scores!"
      fi
      if [ "$flt2" = "#" ]; then
        gblscore2="1"
      elif [ "$flt2" = "." ]; then
        gblscore2="0"
      else
        echo "Check Gblocks scores!"
      fi
      
      ### Calculate mean Gblocks score. Do we really need this?
      #gblscore=$(printf "${gblscore1}\n $gblscore2" | datamash mean 1)
      
      # Extract MSA columns for the co-evolving amino acids  
      echo -e "Extracting $colA-$colB"
      cat ./msa/$msa1.fa | seqkit subseq -r ${colA}:${colA} | seqkit sort -o columnStats/$msa1-$colA.txt
      cat ./msa/$msa2.fa | seqkit subseq -r ${colB}:${colB} | seqkit sort -o columnStats/$msa2-$colB.txt

      columnTotalA=$(sed '/^>/d' columnStats/$msa1-$colA.txt | wc -l)
       columnGapsA=$(sed '/^>/d' columnStats/$msa1-$colA.txt | grep "-"  | wc -l)
       roundPerGapA=$(printf "%1.5f" `echo "1 - ${columnGapsA}/${columnTotalA}" | bc -l`)

      columnNoGapsA=$(sed '/^>/d' columnStats/$msa1-$colA.txt | grep -v "-" | wc -l)
      columnUniqueA=$(sed '/^>/d' columnStats/$msa1-$colA.txt | grep -v "-" | sort | uniq -c | sort -n -r | sed -n 1p | awk '{ print $1 }')
      roundDivResA=$(printf "%1.5f" `echo "1 - ${columnUniqueA}/${columnNoGapsA}" | bc -l`)

      columnTotalB=$(sed '/^>/d' columnStats/$msa2-$colB.txt | wc -l)
       columnGapsB=$(sed '/^>/d' columnStats/$msa2-$colB.txt | grep "-" | wc -l)
       roundPerGapB=$(printf "%1.5f" `echo "1 - ${columnGapsB}/${columnTotalB}" | bc -l`)

      columnNoGapsB=$(sed '/^>/d' columnStats/$msa2-$colB.txt | grep -v "-" | wc -l)
      columnUniqueB=$(sed '/^>/d' columnStats/$msa2-$colB.txt | grep -v "-" | sort | uniq -c | sort -n -r | sed -n 1p | awk '{ print $1 }')
      roundDivResB=$(printf "%1.5f" `echo "1 - ${columnUniqueB}/${columnNoGapsB}" | bc -l`)

      # Calculate gaps and diversity    
      GapsAB=$(printf "${roundPerGapA}\n $roundPerGapB" | datamash mean 1 )
      DivsAB=$(printf "${roundDivResA}\n $roundDivResB" | datamash mean 1 )
      
      # Include correlation threshold as well. Just in case, make sure value is not in scientific format
      corrT1=$(printf "%1.10f" `sed -n '2p' ${msa1}.fa_${msa2}.fa-coev_inter.csv | awk '{print $6}'`)
      corrT2=$(printf "%1.10f" `sed -n '2p' ${msa2}.fa_${msa1}.fa-coev_inter.csv | awk '{print $6}'`)
       corrT=$(printf "${corrT1}\n $corrT2" | datamash mean 1)
      	
      # Calculate coevolutionary correlation, normalized to threshold. For normalization:
      # https://www.mathworks.com/matlabcentral/answers/322438-normalize-data-with-a-threshold
      # Thanks to Dian Dimitrov.
      normC=$(printf "%1.10f" `echo "($corr - $corrT)/(1 - $corrT)" |bc -l`)
      else
      # Make it echo sth else...
        echo "No amino acid pairs passed the Bonferroni correction for $msa1 $msa2!"
      fi
      echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $gblscore1 $gblscore2 $GapsAB $DivsAB $corrT $corr $normC $boot $p_value $bonferroni $holm $bh $hochberg $hommel $by $fdr" >> bothWays-corrected-columns.tsv
  done
  sed -i "1i msa1 msa2 colA realA colB realB seq1 seq2 gblscore1 gblscore2 GapsAB DivsAB corrT corr normC boot p_value bonferroni holm bh hochberg hommel by fdr" bothWays-corrected-columns.tsv
cd ..
}

# Export column statistics
protein_pairs_stats() {
  local coevPair="${1}"
  cd $coevPair

  # Define UniProt numbers of the two proteins
  msa_1=$(sed -n '2p' bothWays-corrected-columns.tsv | awk '{print $1}')
  msa_2=$(sed -n '2p' bothWays-corrected-columns.tsv | awk '{print $2}')
  
  ### Count the number of co-evolving sites from each protein. We have headers, so skip them...
  sitesCountA=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $4}' | datamash countunique 1)
  sitesCountB=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $6}' | datamash countunique 1) 
 
  ### Report total comparisons (MSA1*MSA2). They are the same FWD and REV.
  totCompar=$(sed -n '2p' ${msa_1}.fa_${msa_2}.fa-coev_inter.csv | awk '{print $4}')

  ### Report correlation threshold
  coevThr=$(sed -n '2p' bothWays-corrected-columns.tsv | awk '{print $13}')

  ### Include average correlation
  avgRa=$(printf "%1.10f" `sed -n '2p' ${msa_1}.fa_${msa_2}.fa-coev_inter.csv | awk '{print $7}'`)
  avgRb=$(printf "%1.10f" `sed -n '2p' ${msa_2}.fa_${msa_1}.fa-coev_inter.csv | awk '{print $7}'`)
  averR=$(printf "${avgRa}\n $avgRb" | datamash mean 1)
      
  ### Include average significant correlation
  avgSigRa=$(printf "%1.10f" `sed -n '2p' ${msa_1}.fa_${msa_2}.fa-coev_inter.csv | awk '{print $8}'`)
  avgSigRb=$(printf "%1.10f" `sed -n '2p' ${msa_2}.fa_${msa_1}.fa-coev_inter.csv | awk '{print $8}'`)
  averSigR=$(printf "${avgSigRa}\n $avgSigRb" | datamash mean 1)
  
  #coevNumber=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash count 1)

        cCoevMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $15}' | datamash min 1)
        cCoevMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $15}' | datamash max 1)
       cCoevMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $15}' | datamash mean 1)
       #cCoevSDEV=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $15}' | datamash sstdev 1)

         bootMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash min 1)
         bootMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash max 1)
        bootMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash mean 1)
	#bootSDEV=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash sstdev 1)
       
        pMeanMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash min 1)
        pMeanMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash max 1)
       pMeanMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash mean 1)
       #pMeanSDEV=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash sstdev 1)
       
   BonferroniMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $20}' | datamash min 1)
   BonferroniMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $20}' | datamash max 1)
  BonferroniMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $20}' | datamash mean 1)
  #BonferroniSDEV=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $20}' | datamash sstdev 1)

#       gblocksMIN=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
#       gblocksMAX=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
#      gblocksMEAN=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
#      gblocksSDEV=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)
     
#          coevMIN=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
#          coevMAX=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
#         coevMEAN=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
#         coevSDEV=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)


#       # Add score after Chi squared based on non-filtered results
#       if [ -s "$TMP/$RESULTS/chi/chi_test/${msa_1}.tsv" ] && [ -s "$TMP/$RESULTS/chi/chi_test/${msa_2}.tsv" ]; then
#         forward=$(grep "${msa_1} ${msa_2}" $TMP/$RESULTS/chi/chi_test/${msa_1}.tsv | awk '{print $7}')
#         reverse=$(grep "${msa_2} ${msa_1}" $TMP/$RESULTS/chi/chi_test/${msa_2}.tsv | awk '{print $7}')
#         chiboth=$(echo "$forward + $reverse" |bc -l)
#       else
#         echo "${msa_1} or ${msa_2} missing" >> $TMP/$RESULTS/miss.non-filtered
#       fi
      
        # Add score after Chi squared based on filtered results
	forward_fin=$(grep "${msa_1} ${msa_2}" $TMP/$RESULTS/chi/chi_test_final/${msa_1}.fa.tsv | awk '{print $7}')
	reverse_fin=$(grep "${msa_2} ${msa_1}" $TMP/$RESULTS/chi/chi_test_final/${msa_2}.fa.tsv | awk '{print $7}')
	chiboth_fin=$(echo "$forward_fin + $reverse_fin" |bc -l)
      
      # Collect data
      #echo "$msa_1 $msa_2 $coevThr $avgCor $avgSignCor $totCompar $coevNumAll $coevMIN $coevMAX $coevMEAN $coevSDEV $coevNumber $cCoevMIN $cCoevMAX $cCoevMEAN $cCoevSDEV $bootMIN $bootMAX $bootMEAN $bootSDEV $pMeanMIN $pMeanMAX $pMeanMEAN $pMeanSDEV $BonferroniMIN $BonferroniMAX $BonferroniMEAN $BonferroniSDEV $gblocksMIN $gblocksMAX $gblocksMEAN $gblocksSDEV $chiboth $chiboth_fin" >> $EOUT
      echo "$msa_1 $msa_2 $coevThr $averR $averSigR $totCompar $sitesCountA $sitesCountB $cCoevMIN $cCoevMAX $cCoevMEAN $cCoevSDEV $bootMIN $bootMAX $bootMEAN $pMeanMIN $pMeanMAX $pMeanMEAN $BonferroniMIN $BonferroniMAX $BonferroniMEAN $BonferroniSDEV $chiboth_fin" >> $EOUT
      echo "${msa_1} ${msa_2} added"
      
cd ..
}

# # Save residue pairs with p-values below threshold and positive correlations
# # for both directions.
# adj_pVal(){
#   local coevPair="${1}"
#   sed 1d $coevPair/$coevPair.clean | while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 ; do
#      corrdec=$(printf "%1.10f" $corr)
#     corr1dec=$(printf "%1.10f" $corr1)
#     corr2dec=$(printf "%1.10f" $corr2)
#     pvalAdec=$(printf "%1.10f" $pvalA)
#     pvalBdec=$(printf "%1.10f" $pvalB)
#     pMeandec=$(printf "%1.10f" $pMean)
#     if (( $(echo "$pMeandec <= $PVALUE" |bc -l) && \
#     	  $(echo "$corr1dec > 0" |bc -l) && \
# 	  $(echo "$corr2dec > 0" |bc -l) )); then
#       echo "Adding $pMeandec"
#       echo "$msa1 $msa2 $colA $realA $colB $realB $meanA $meanB $corrdec $boot $pvalAdec $pvalBdec $pMeandec $corr1dec $corr2dec" >> $coevPair/$coevPair.$PVALUE.clean
#     else
#       echo "Skipping $pMeandec"
#     fi
#   done
#     
#   # Run Luqman's script for Bonferroni correction of p-values
#   if [ -s $coevPair/$coevPair.$PVALUE.clean ]; then
#     sed -i "1i msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2" $coevPair/$coevPair.$PVALUE.clean
#     echo -e "Running R for Bonferroni correction for $coevPair"
#     Rscript $CWD/R/AdjPval.R $coevPair/$coevPair.$PVALUE.clean $coevPair/$coevPair.$PVALUE.corrected
#   elif [ -z $coevPair/$coevPair.$PVALUE.clean ]; then
#     echo -e "[SKIP] $coevPair.$PVALUE.clean is empty. Skipping..."
#   else
#     echo -e "[SKIP] $coevPair.$PVALUE.clean does not exist. Skipping..."
#   fi
# }

# # Extract columns for pairs that passed the Bonferroni correction test.
# # Collect only pairs that exist in the reference organism (with "real" position > 0)
# extract_columns(){
#   local coevPair="${1}"
#   if [ -s "$coevPair/${coevPair}.${PVALUE}.corrected" ]; then
#     sed 1d $coevPair/$coevPair.$PVALUE.corrected |
#     while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 bonferroni ; do
#       if (( $(echo "$bonferroni <= $BONFERRONI" |bc -l) && $(echo "$realA > 0" |bc -l) && $(echo "$realB > 0" |bc -l) )); then
#         seq1=$( sed "${realA}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa1}.fa.${ORGANISM}.ref | awk '{print $1$2}' )
#         seq2=$( sed "${realB}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa2}.fa.${ORGANISM}.ref | awk '{print $1$2}' )
#         flt1=$( sed "${realA}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa1}.fa.${ORGANISM}.ref | awk '{print $3}' )
#         flt2=$( sed "${realB}q;d" $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/${msa2}.fa.${ORGANISM}.ref | awk '{print $3}' )
#         echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $flt1 $flt2 $meanA $meanB $corr $boot $pvalA $pvalB $pMean $corr1 $corr2 $bonferroni" >> $coevPair/$coevPair.$PVALUE.corrected-$BONFERRONI
#         echo -e "Extracting $colA-$colB with corrected p-value=${bonferroni}"
#         mkdir -p $coevPair/columnStats-$BONFERRONI/$colA-$colB/
#         cat $coevPair/msa/$msa1.fa | seqkit subseq -r ${colA}:${colA} | seqkit sort -o $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colA-$msa1.txt
#         cat $coevPair/msa/$msa2.fa | seqkit subseq -r ${colB}:${colB} | seqkit sort -o $coevPair/columnStats-$BONFERRONI/$colA-$colB/$colB-$msa2.txt
#        elif (( $(echo "$bonferroni > $BONFERRONI" |bc -l) )); then
#          echo -e "Skipping $colA-$colB with corrected p-value=${bonferroni}"
#        elif (( $(echo "$realA == 0" |bc -l) )); then
#          echo -e "Skipping $colA-$colB with MSA $colA -> $realA in $ORGANISM"
#        elif (( $(echo "$realB == 0" |bc -l) )); then
#          echo -e "Skipping $colA-$colB with MSA $colB -> $realB in $ORGANISM"
#        else
#          echo "Check settings!"
#        fi
#     done
#   else
#     echo "Skipping $coevPair.$PVALUE.corrected"
#   fi
# }

# Calculate statistics on columns for the protein pairs with co-evolving
# residues that passed the Bonferroni correction. First, check if file
# exists and has non-zero size "if [ -s ... ]; then"
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


      # Report total number of pairs
      totNumP=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $3}')
      
      # Report total comparisons (MSA1*MSA2)
      totComp=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $4}')
	
      # Include correlation threshold as well. Just in case, make sure value is not in scientific format
      thrsR=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $6}')
      corrT=$(printf "%1.10f" $thrsR)
      
      # Include average correlation
      averageR=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $7}')
      averR=$(printf "%1.10f" $averageR)
      
      # Include average significant correlation
      averageSigR=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $8}')
      averSigR=$(printf "%1.10f" $averageSigR)
	
      # Calculate coevolutionary correlation, normalized to threshold. For normalization:
      # Thanks to Dian Dimitrov.
      # https://www.mathworks.com/matlabcentral/answers/322438-normalize-data-with-a-threshold
      normC=$(echo "($corr - $corrT)/(1 - $corrT)" |bc -l)
      cCoev=$(printf "%1.10f" $normC)
      
      # Collect info of each amino acid pair into a singe file
      echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $flt1 $flt2 $meanA $meanB $corr $boot $corrT $averR $averSigR $cCoev $pvalA $pvalB $pMean $pDiff $corr1 $corr2 $bonferroni $roundPerGapA $roundPerGapB $GapsAB $roundDivResA $roundDivResB $DivsAB $totNumP $totComp" >> $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv

      # Convert Gblocks data (. or #) into numerical score. Collect data regardless of the Gblocks score.
      if (( $(echo "$roundPerGapA > $EGAPS" |bc -l) && $(echo "$roundPerGapB > $EGAPS" |bc -l) && $(echo "$pDiff < $EPDIF" |bc -l) )); then
        if [ "$flt1" = "#" ] && [ "$flt2" = "#" ]; then
          echo -e "$msa1 <-> $msa2 with score 1"
          echo "$msa1 $msa2 $cCoev $pMean $bonferroni $boot 1 $corr $corrT $averR $averSigR $totNumP $totComp" >> $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv
        elif [ "$flt1" = "#" ] && [ "$flt2" = "." ]; then
          echo -e "$msa1 <-> $msa2 with score 0.5"
          echo "$msa1 $msa2 $cCoev $pMean $bonferroni $boot 0.5 $corr $corrT $averR $averSigR $totNumP $totComp" >> $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv
        elif [ "$flt1" = "." ] && [ "$flt2" = "#" ]; then
          echo -e "$msa1 <-> $msa2 with score 0.5"
          echo "$msa1 $msa2 $cCoev $pMean $bonferroni $boot 0.5 $corr $corrT $averR $averSigR $totNumP $totComp" >> $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv
        elif [ "$flt1" = "." ] && [ "$flt2" = "." ]; then
          echo -e "$msa1 <-> $msa2 with score 0"
          echo "$msa1 $msa2 $cCoev $pMean $bonferroni $boot 0 $corr $corrT $averR $averSigR $totNumP $totComp" >> $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv
        else
          echo "Something went wrong... CHECK!"
        fi
      else
        # Make it echo sth else...
        echo "No amino acid pairs passed the Bonferroni correction for $msa1 $msa2!"
      fi
  done

  # If no amino acid pairs passed the Bonferroni correction, prepare results
  # for Chi^2 test as a protein pair with 0 coevolving amino acids.
  elif [ ! -s "$coevPair/${coevPair}.${PVALUE}.corrected-${BONFERRONI}" ]; then
    emsa_1=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $1}')
    emsa_2=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $2}')
    etotComp=$(sed -n '2p' $coevPair/coev_inter.csv | awk '{print $4}')
    echo "${emsa_1%.*} ${emsa_2%.*} 0 $etotComp" >> $TMP/$RESULTS/chi/proteinsFinal/${emsa_1%.*}.tsv
    echo "${emsa_2%.*} ${emsa_1%.*} 0 $etotComp" >> $TMP/$RESULTS/chi/proteinsFinal/${emsa_2%.*}.tsv
    echo "[NOPAIR] Prepare results for Chi_2 $coevPair.$PVALUE.corrected-${BONFERRONI}"
  else
    echo "[PASSED] Will analyze later: Chi_2 $coevPair.$PVALUE.corrected-${BONFERRONI}"
  fi

  # When column statistics are done, prepare results for Chi^2 test
  if [ -s "$coevPair/summary-${PVALUE}-${BONFERRONI}.tsv" ]; then
    msa_1=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $1}')
    msa_2=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $2}')
    coevNumber=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash count 1)
    totCompar=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $13}')

    echo "${msa_1} ${msa_2} $coevNumber $totCompar" >> $TMP/$RESULTS/chi/proteinsFinal/${msa_1}.tsv
    echo "${msa_2} ${msa_1} $coevNumber $totCompar" >> $TMP/$RESULTS/chi/proteinsFinal/${msa_2}.tsv
    echo "PICKING: $coevPair"
  else
    echo "SKIPPING: $coevPair"
  fi
}

# Export column statistics
exp_column_stats() {
  local coevPair="${1}"     
  if [ -s "$coevPair/summary-${PVALUE}-${BONFERRONI}.tsv" ]; then
           msa_1=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $1}')
           msa_2=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $2}')
  
      coevNumber=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash count 1)

        cCoevMIN=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
        cCoevMAX=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
       cCoevMEAN=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
       cCoevSDEV=$(awk '{print $3}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)
       
        pMeanMIN=$(awk '{print $4}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
        pMeanMAX=$(awk '{print $4}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
       pMeanMEAN=$(awk '{print $4}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
       pMeanSDEV=$(awk '{print $4}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)
       
   BonferroniMIN=$(awk '{print $5}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
   BonferroniMAX=$(awk '{print $5}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
  BonferroniMEAN=$(awk '{print $5}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
  BonferroniSDEV=$(awk '{print $5}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)

         bootMIN=$(awk '{print $6}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
         bootMAX=$(awk '{print $6}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
        bootMEAN=$(awk '{print $6}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
	bootSDEV=$(awk '{print $6}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)
  
      gblocksMIN=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
      gblocksMAX=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
     gblocksMEAN=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
     gblocksSDEV=$(awk '{print $7}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)
     
         coevMIN=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash min 1)
         coevMAX=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash max 1)
        coevMEAN=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash mean 1)
        coevSDEV=$(awk '{print $8}' $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | datamash sstdev 1)
     
         coevThr=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $9}')
          avgCor=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $10}')
      avgSignCor=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $11}')
      coevNumAll=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $12}')
       totCompar=$(head -n 1 $coevPair/summary-${PVALUE}-${BONFERRONI}.tsv | awk '{print $13}')

      # Add score after Chi squared based on non-filtered results
      if [ -s "$TMP/$RESULTS/chi/chi_test/${msa_1}.tsv" ] && [ -s "$TMP/$RESULTS/chi/chi_test/${msa_2}.tsv" ]; then
        forward=$(grep "${msa_1} ${msa_2}" $TMP/$RESULTS/chi/chi_test/${msa_1}.tsv | awk '{print $7}')
        reverse=$(grep "${msa_2} ${msa_1}" $TMP/$RESULTS/chi/chi_test/${msa_2}.tsv | awk '{print $7}')
        chiboth=$(echo "$forward + $reverse" |bc -l)
      else
        echo "${msa_1} or ${msa_2} missing" >> $TMP/$RESULTS/miss.non-filtered
      fi
      
      # Add score after Chi squared based on filtered results
      if [ -s "$TMP/$RESULTS/chi/chi_test_final/${msa_1}.tsv" ] && [ -s "$TMP/$RESULTS/chi/chi_test_final/${msa_2}.tsv" ]; then
	forward_fin=$(grep "${msa_1} ${msa_2}" $TMP/$RESULTS/chi/chi_test_final/${msa_1}.tsv | awk '{print $7}')
	reverse_fin=$(grep "${msa_2} ${msa_1}" $TMP/$RESULTS/chi/chi_test_final/${msa_2}.tsv | awk '{print $7}')
	chiboth_fin=$(echo "$forward_fin + $reverse_fin" |bc -l)
      else
        echo "${msa_1} or ${msa_2} missing" >> $TMP/$RESULTS/miss.filtered
      fi
      
      # Collect data
      echo "$msa_1 $msa_2 $coevThr $avgCor $avgSignCor $totCompar $coevNumAll $coevMIN $coevMAX $coevMEAN $coevSDEV $coevNumber $cCoevMIN $cCoevMAX $cCoevMEAN $cCoevSDEV $bootMIN $bootMAX $bootMEAN $bootSDEV $pMeanMIN $pMeanMAX $pMeanMEAN $pMeanSDEV $BonferroniMIN $BonferroniMAX $BonferroniMEAN $BonferroniSDEV $gblocksMIN $gblocksMAX $gblocksMEAN $gblocksSDEV $chiboth $chiboth_fin" >> $EOUT
      echo "${msa_1} ${msa_2} added"
      
    else
      echo "SKIPPING THE WHOLE PAIR: $coevPair"
    fi
}

# Clean up the results, add column headers
summary_cleanup(){
  if [ -s $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv ] && [ -s $TMP/$RESULTS/filtered-P${PVALUE}-gaps${EGAPS}-pDiff${EPDIF}-Bonf${BONFERRONI}.tsv ] ; then
    sed 1d $TMP/tsv/proteinsFound.tsv | while read -r idxml namexml attribute ; do
      sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv
      sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/filtered-P${PVALUE}-gaps${EGAPS}-pDiff${EPDIF}-Bonf${BONFERRONI}.tsv
    done
    sed -i \
    "1i Name1 msa1 Name2 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot corrT averR averSigR cCoev pvalA pvalB pMean pDiff corr1 corr2 bonferroni GapA GapB GapAB DivA DivB DivAB totNumP totComp" \
    $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv
    sed -i \
    "1i Name1 msa1 Name2 msa2 coevThr avgCor avgSignCor totCompar coevNumAll coevMIN coevMAX coevMEAN coevSDEV coevNumber cCoevMIN cCoevMAX cCoevMEAN cCoevSDEV bootMIN bootMAX bootMEAN bootSDEV pMeanMIN pMeanMAX pMeanMEAN pMeanSDEV BonferroniMIN BonferroniMAX BonferroniMEAN BonferroniSDEV gblocksMIN gblocksMAX gblocksMEAN gblocksSDEV chiboth chiboth_fin" \
    $TMP/$RESULTS/filtered-P${PVALUE}-gaps${EGAPS}-pDiff${EPDIF}-Bonf${BONFERRONI}.tsv
  else
    echo "No pairs-P${PVALUE}-B${BONFERRONI}.tsv or $TMP/$RESULTS/filtered-P${PVALUE}-gaps${EGAPS}-pDiff${EPDIF}-Bonf${BONFERRONI}.tsv!"
  fi
}


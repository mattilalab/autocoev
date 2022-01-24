#!/bin/bash

# Define paths and files
EOUT="$TMP/$RESULTS/proteinPairs.tsv"

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
    
    # Did the run fail? Is there any output at all on row 2 of coev_inter.csv?
    if [ -z "${SUMMARY}" ]; then
      echo -e "[\e[91mFAILED RUN\e[39m] $pair"
      cp -a $pair $TMP/$RESULTS/fail
      
    # Is there coevolution detected at all? If we have an output, is the number of coevolving amino acids 0?
    elif [ "$numPairs" -eq 0 ]; then
      echo -e "[\e[34mNON-COEVOL\e[39m] $pair"
      
      ### Chi^2 ### Export these, since they are needed for the Chi^2 test later (background calculation)
      mkdir -p $TMP/$RESULTS/chi/proteins/
      echo "${Seq1} ${Seq2} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1}.tsv
      echo "${Seq2} ${Seq1} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2}.tsv
      
      # Copy the non-coevolving pairs to a separate folder of results
      mkdir -p $TMP/$RESULTS/nocoev/$folder
      cp -a $pair $TMP/$RESULTS/nocoev/$folder

    # If coevolution was detected, copy to a separate folder of results
    elif [ "$numPairs" -gt 0 ]; then
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
    
      # If we have PhyML generated trees, do the same for them
      if [ "$TREESCAPS" = "phyml" ]; then
        echo -e "[REVPREPARE] $pair (PhyML trees)"
        mkdir -p ../tre-rev
	cd ../tre
	cp ${firstMSA%.*}.tre ../tre-rev/b_${firstMSA%.*}.tre
	cp ${secondMSA%.*}.tre ../tre-rev/a_${secondMSA%.*}.tre
      elif [ "$TREESCAPS" = "auto" ]; then
        echo -e "[REVPREPARE] $pair (auto trees)"
      else
        echo "Something went wronng at trees REV step!"
      fi

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
      echo -e "[\e[91mFAILED RUN\e[39m]: $pair"
      mv $pair $TMP/$RESULTS/fail
    elif [ "$numPairs" -eq 0 ]; then
      echo -e "[\e[34mNON COEVOL\e[39m]: $pair"
      mkdir -p $TMP/$RESULTS/nocoev/$folder
      mv $pair $TMP/$RESULTS/nocoev/$folder
      
      ### Chi^2 ### These are needed for the Chi^2 test
      echo "${Seq1#*_} ${Seq2#*_} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq1#*_}.tsv
      echo "${Seq2#*_} ${Seq1#*_} $numPairs $totalComp" >> $TMP/$RESULTS/chi/proteins/${Seq2#*_}.tsv
      
     elif [ "$numPairs" -gt 0 ]; then
       echo -e "[COEVOL REV] $pair"
       # https://stackoverflow.com/a/15149278
       cd $pair
       mv coev_inter.csv ${Seq1#*_}_${Seq2#*_}-coev_inter.csv  
       mv ${Seq1}_${Seq2}.out ${Seq1#*_}_${Seq2#*_}.out
    else
      echo -e "Something went wrong for $folder/$pair ... Check!"
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
  sed -i "1i msaA msaB colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corrA corrB" ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean
  
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
  sed -i "1i msaA msaB colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corrA corrB" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean
  
  sed 1d ${PROTEINONE%.*}_${PROTEINTWO%.*}.clean | while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 ; do
    if rowmatch=$(LANG=C grep -F -w "$msa2 $msa1 $colB $realB $colA $realA" ${PROTEINTWO%.*}_${PROTEINONE%.*}.clean) ; then
      #echo "Found a match bothways: $rowmatch"
    
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
	    $(echo "$fcorr2dec > 0" |bc -l) && $(echo "$rcorr2dec > 0" |bc -l) )); then

          # Save the two correlation s (FWD and REV) mean value. The correlation is, in turn,
	  # estimated in two directions for each FWD and REV (e.g. $fcorr1dec & $fcorr2dec),
	  # but we do not want to output a crazy amount of stuff, do we? Same goes for the
	  # p-values (e.g. fpvalAdec & fpvalBdec).	  
	  corrdec=$(printf "${fcorrdec}\n ${rcorrdec}" | datamash mean 1)
          #corr1dec=$(printf "${fcorr1dec}\n $rcorr1dec" | datamash mean 1)
          #corr2dec=$(printf "${fcorr2dec}\n $rcorr2dec" | datamash mean 1)
          #pvalAdec=$(printf "${fpvalAdec}\n $rpvalAdec" | datamash mean 1)
          #pvalBdec=$(printf "${fpvalBdec}\n $rpvalBdec" | datamash mean 1)
          pMeandec=$(printf "${fpMeandec}\n ${rpMeandec}" | datamash mean 1)
      
        echo "$msa1 $msa2 $colA $realA $colB $realB $corrdec $boot $pMeandec" >> bothWays.tsv
        echo -e "[COEV RESID] $coevPair"
	
       elif (( $(echo "$fpMeandec >= $PVALUE" |bc -l) || $(echo "$rpMeandec >= $PVALUE" |bc -l) || \
    	  $(echo "$fcorr1dec <= 0" |bc -l) || $(echo "$rcorr1dec <= 0" |bc -l) || \
	  $(echo "$rcorr2dec <= 0" |bc -l) || $(echo "$rcorr2dec <= 0" |bc -l) )); then
	  echo -e "[SKIP RESID] $coevPair"
       else
         echo "Something went wrong at $folder/$coevPair/${PROTEINONE%.*}_${PROTEINTWO%.*}.clean"
       fi
    else
      echo -e "[NO REVERSE] $coevPair"
    fi
    done

  # Check if we have co-evolving pairs detected bidirectionally
  if [ -f bothWays.tsv ]; then
    echo "[RWD-REV EQ] $coevPair"
    sed -i "1i msa1 msa2 colA realA colB realB corr boot p_value" bothWays.tsv
  
    ### Chi^2 ### How many pairs do we have left? This is how CAPS2
    # counts pairs. Also get the totalcomp and export it.
    pairsNumber=$(sed 1d bothWays.tsv | awk '{print $3}' | datamash count 1)
    totCompares=$(sed -n 2p ${PROTEINONE%.*}.fa_${PROTEINTWO%.*}.fa-coev_inter.csv | awk '{print $4}')
    echo "${PROTEINONE%.*}.fa ${PROTEINTWO%.*}.fa $pairsNumber $totCompares" >> $TMP/$RESULTS/chi/proteins/${PROTEINONE%.*}.fa.tsv
    echo "${PROTEINTWO%.*}.fa ${PROTEINONE%.*}.fa $pairsNumber $totCompares" >> $TMP/$RESULTS/chi/proteins/${PROTEINTWO%.*}.fa.tsv
    cd ..
  elif [ ! -f bothWays.tsv ]; then
    echo -e "[FW-RE DIFF] $coevPair"
    cd ..
    mkdir -p $TMP/$RESULTS/noBothWays
    mv $coevPair $TMP/$RESULTS/noBothWays
  else
    echo "Something went wrong at $folder/$coevPair"
  fi
}

# get the amino acid from the reference organism's sequence and calculate
# column statistics
extract_columns_stats(){
  local coevPair="${1}"
  
  if [ -d "$coevPair" ]; then
    cd $coevPair
  
    # Run Luqman's script for Bonferroni & co correction of p-values
    echo -e "[MULTI HYPO] $coevPair"
    Rscript $CWD/functions/AdjPval.R bothWays.tsv bothWays-corrected.tsv
  
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
      
        # Calculate mean Gblocks score. Do we really need this?
        gblscore=$(printf "${gblscore1}\n $gblscore2" | datamash mean 1)
      
        # Extract MSA columns for the co-evolving amino acids  
        echo -e "[PROPERTIES] $coevPair/$colA-$colB"
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
        echo -e "[NOT IN REF] $coevPair/$colA-$colB not in $ORGANISM!"
      fi
      echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $gblscore1 $gblscore2 $gblscore $GapsAB $DivsAB $corrT $corr $normC $boot $p_value $bonferroni $holm $bh $hochberg $hommel $by $fdr" >> bothWays-corrected-columns.tsv
      echo "$msa1 $msa2 $colA $realA $colB $realB $seq1 $seq2 $gblscore $GapsAB $DivsAB $corrT $corr $normC $boot $p_value $bonferroni $holm $bh $hochberg $hommel $by $fdr" >> $TMP/$RESULTS/allResidues.tsv
    done
    sed -i "1i msaA msaB colA realA colB realB seqA seqB gblAB GapsAB DivsAB corrT corr normC boot p_value bonferroni holm bh hochberg hommel by fdr" bothWays-corrected-columns.tsv
  elif [ ! -d "$coevPair" ]; then
    echo "No protein pairs"
  else
    echo "Something went wrong!"
  fi
cd ..
}

# Generate and export column statistics
protein_pairs_stats() {
  local coevPair="${1}"

  # Do we have protein pairs left?
  if [ -d "$coevPair" ]; then
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
  
    ### Gblocks scores
    gblocksMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $11}' | datamash min 1)
    gblocksMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $11}' | datamash max 1)
    gblocksMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $11}' | datamash mean 1)

    # MSA column gaps
    GapsMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $12}' | datamash min 1)
    GapsMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $12}' | datamash max 1)
    GapsMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $12}' | datamash mean 1)

    # MSA column diversity  
    DivsMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $13}' | datamash min 1)
    DivsMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $13}' | datamash max 1)
    DivsMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $13}' | datamash mean 1)

    # Normalized coevolution, corrected to threshold
    cCoevMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $16}' | datamash min 1)
    cCoevMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $16}' | datamash max 1)
    cCoevMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $16}' | datamash mean 1)

    # Bootstrap of CAPS
    bootMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $17}' | datamash min 1)
    bootMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $17}' | datamash max 1)
    bootMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $17}' | datamash mean 1)

    # Mean p-value of both FWD and REV coevolution run       
    pMeanMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash min 1)
    pMeanMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash max 1)
    pMeanMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $18}' | datamash mean 1)

    # Bonferroni corrected p-values, same as above       
    BonferroniMIN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash min 1)
    BonferroniMAX=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash max 1)
    BonferroniMEAN=$(sed 1d bothWays-corrected-columns.tsv | awk '{print $19}' | datamash mean 1)

    # Add score after Chi squared based on filtered results
    forward_fin=$(grep "${msa_1} ${msa_2}" $TMP/$RESULTS/chi/chi_test_final/${msa_1}.fa.tsv | awk '{print $7}')
    reverse_fin=$(grep "${msa_2} ${msa_1}" $TMP/$RESULTS/chi/chi_test_final/${msa_2}.fa.tsv | awk '{print $7}')
    chiboth_fin=$(echo "$forward_fin + $reverse_fin" |bc -l)
      
    # Collect data in a single file, which can be imported in Cytoscape
    echo "$msa_1 $msa_2 $coevThr $averR $averSigR $totCompar $sitesCountA $sitesCountB $gblocksMIN $gblocksMAX $gblocksMEAN $GapsMIN $GapsMAX $GapsMEAN $DivsMIN $DivsMAX $DivsMEAN $cCoevMIN $cCoevMAX $cCoevMEAN $bootMIN $bootMAX $bootMEAN $pMeanMIN $pMeanMAX $pMeanMEAN $BonferroniMIN $BonferroniMAX $BonferroniMEAN $chiboth_fin" >> $EOUT
    echo "${msa_1} ${msa_2} added"
    
    cd ..

  # Skip if $folder was left without protein pairs  
  elif [ ! -d "$coevPair" ]; then
    echo "No protein pairs"
  else
    echo "Something went wrong!"
  fi
}

# Clean up the results, add column headers
summary_cleanup(){
  if [ -s $TMP/$RESULTS/allResidues.tsv ] && [ -s $EOUT ] ; then
    sed 1d $TMP/tsv/proteinsFound.tsv | while read -r idxml namexml attribute ; do
      sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/allResidues.tsv
      sed -i "s/$idxml/$namexml $idxml/g" $EOUT
    done
    sed -i "1i NameA msaA msaB NameB colA realA colB realB seqA seqB GblAB GapsAB DivsAB corrT corr normC boot p_value bonferroni holm bh hochberg hommel by fdr" $TMP/$RESULTS/allResidues.tsv
    sed -i \
    "1i NameA msaA NameB msaB coevThr averR averSigR totCompar sitesCountA sitesCountB gblocksMIN gblocksMAX gblocksMEAN GapsMIN GapsMAX GapsMEAN DivsMIN DivsMAX DivsMEAN cCoevMIN cCoevMAX cCoevMEAN bootMIN bootMAX bootMEAN p_valueMIN p_valueMAX p_valueMEAN BonferroniMIN BonferroniMAX BonferroniMEAN chiboth_fin" \
    $EOUT
  else
    echo "No Protein pairs!"
  fi
}

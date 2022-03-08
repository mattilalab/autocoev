#!/bin/bash

# Export all
set -a

# What date and time is it?
DATESTAMP=$(date)

# Current work dir
CWD=$(pwd)

NAME="AutoCoEv" # script name
VER="0.2.7beta" # version

# Load settings first
. $CWD/settings.conf

# Load functions
. $CWD/functions/databases.sh
. $CWD/functions/retrieval.sh
. $CWD/functions/blast.sh
. $CWD/functions/guidance.sh
. $CWD/functions/msa.sh
. $CWD/functions/trees.sh
. $CWD/functions/pairing.sh
. $CWD/functions/caps.sh
. $CWD/functions/results.sh
. $CWD/functions/check.sh
. $CWD/functions/chi2.sh

# Make sure we have the right decimal separator...
export LC_ALL=C

# User provided input
echo -e "\nWelcome to \e[46m${NAME} version ${VER}!\e[49m $DATESTAMP\n"

# Check for binaries
echo -e "Checking for executables..."
check_bin
echo ""

# Let user change work dir here as well
read -p "Change working dir or press ENTER to accept: " -e -i $TMP TMP
echo -e "\n\e[96m$TMP\e[39m\n"
mkdir -p $TMP

# Show databases menu
# Help from:
# https://askubuntu.com/questions/1705/how-can-i-create-a-select-menu-in-a-shell-script
# https://unix.stackexchange.com/questions/293340/bash-how-can-i-re-display-selection-menu-after-a-selection-is-chosen-and-perfo
echo -e "Prepare databases or skip [11]:"
PS3="Your choice: "
PREPOPT=( "Download $GENEXREFALL"
	 "Download $OG2GENESALL"
	 "Download $ALLFASTA"
	 "Check MD5sum of databases"
         "Extract $GENEXREFALL"
	 "Extract $OG2GENESALL"
	 "Extract $ALLFASTA"
         "Index $ALLFASTA"
         "Trim $GENEXREF"
	 "Trim $OG2GENES"
         "[DONE AND CONTINUE]" )

select prep in "${PREPOPT[@]}" ; do
case $prep in

 "Download $GENEXREFALL")
 download_db $GENEXREFALL
 echo -e "\nDone with 1)"
 ;;
  
 "Download $OG2GENESALL")
 download_db $OG2GENESALL
 echo -e "\nDone with 2)"
 ;;
 
 "Download $ALLFASTA")
 download_db $ALLFASTA
 echo -e "\nDone with 3)"
 ;;
 
 "Check MD5sum of databases")
 md5sum_check "$GENEXREFALL" "$GENEXREFALLM5"
 md5sum_check "$OG2GENESALL" "$OG2GENESALLM5"
 md5sum_check "$ALLFASTA" "$ALLFASTAM5"
 echo -e "\nDone with 4"
 ;;
 
 "Extract $GENEXREFALL")
 extract_db $GENEXREFALL
 echo -e "\nDone with 5)"
 ;;
 
 "Extract $OG2GENESALL")
 extract_db $OG2GENESALL
 echo -e "\nDone with 6)"
 ;;
 
 "Extract $ALLFASTA")
 extract_db $ALLFASTA
 echo -e "\nDone with 7)"
 ;;
 
 "Index $ALLFASTA")
 index_all_fasta
 echo -e "\nDone with 8)"
 ;;
 
 "Trim $GENEXREF")
 trim_db "$GENEXREFALL" "\b${ORGANISM}_" "${ORGANISM}"
 echo -e "\nDone with 9)"
 ;;
 
 "Trim $OG2GENES")
 trim_db "$OG2GENESALL" "at${LEVEL}\b" ${LEVEL}
 echo -e "\nDone with 10)"
 ;;
 
"[DONE AND CONTINUE]")
echo -e "Continue...\n"
break 2
;;

*)
echo "Invalid option"
;;

esac
REPLY=
done

cd $CWD

# Preview databases in their location
preview_tab $DTB/$GENEXREF.tab
preview_tab $DTB/$OG2GENES.tab
preview_tab $DTB/$ALLFASTA.tab
preview_tab $DTB/$ALLFASTA.tab.index

# Let user change the proteins list and show preview
read -p "Change proteins folder or press ENTER to accept: " -e -i $PROTEIN PROTEIN
ls $PROTEIN/

# Let user change the species list and show preview
read -p "Change species list or press ENTER to continue: " -e -i $SPECIES SPECIES
preview_tab $SPECIES

#read -p "How many threads do we use?: " -e -i $THREADS THREADS
echo -e "CPU threads to be used are: \e[92m$THREADS\e[39m\n"
echo -e "We achieve parallelization via \e[92mGNU/Parallel\e[39m"
echo -e "Run 'parallel --citation' once to silence its citation notice.\n\n"

# Output a summary of settings
echo -e "Work directory: \e[92m${TMP}\e[39m\n"
echo -e "UniProt reference sequences are from:.........\e[92m$ORGANISM\e[39m"
echo -e "Reverse BLAST identity (%) cutoff:............\e[92m$PIDENT\e[39m"
echo -e "Reverse BLAST gaps (%) cutoff:................\e[92m$PGAPS\e[39m"
echo -e "MSAs for orthologues assessment created by....\e[92m$GUIDANCEMSA\e[39m"
echo -e "Guidance orthologues cutoff...................\e[92m$GUIDANCECUT\e[39m"
echo -e "MSAs for PhyML and/or CAPS to be created by:..\e[92m$MSAMETHOD\e[39m"
echo -e "Use rooted PhyML trees (if selected)?.........\e[92m$TREESROOT\e[39m"
echo -e "Use guide tree for PhyML (if selected)?.......\e[92m$PHYMLGUIDE\e[39m"
echo -e "Trees to be used with CAPS:...................\e[92m$TREESCAPS\e[39m"
echo -e "Minimum number of common species in a pair:...\e[92m$MINCOMMONSPCS\e[39m"
echo -e "CAPS alpha-value cutoff at runtime:...........\e[92m$ALPHA\e[39m"
echo -e "CAPS bootstrap value at runtime:..............\e[92m$BOOT\e[39m"
echo -e "Postrun P-value correlation cutoff:...........\e[92m$PVALUE\e[39m"
echo -e "\n"

cd $TMP

# Main menu
echo -e "Select a step:"
PS3="Your choice: "
SEQOPT=( "Pair UniProt <-> OrthoDB <-> OGuniqueID"
         "Get orthologues (level: $LEVEL)"
         "Download sequences from UniProt (organism: $ORGANISM)"
         "BLAST orthologues against UniProt sequence ($ORGANISM, detailed: $DETBLAST)"
         "Get FASTA sequences of the best hits (identity: $PIDENT; gaps: $PGAPS)"
	 "[MSA] Exclude too divergent sequences"
         "[MSA] Create MSA with selected method ($MSAMETHOD)"
         "[TRE] Prepare trees ($TREESCAPS, $MSAMETHOD, $PHYMLGUIDE, $TREESROOT)"
         "[RUN] Create pairs ($PAIRINGMANNER)"
         "[RUN] CAPS run (alpha: $ALPHA, $MSAMETHOD, $TREESCAPS)"
         "[RUN] Inspect results and re-run CAPS (REV)"
         "[RES] Generate columns stats"
         "[Exit script]" )

select opt in "${SEQOPT[@]}" ; do
case $opt in

"Pair UniProt <-> OrthoDB <-> OGuniqueID")

  report_tsv
  
  # Define the protein categories
  PROTLST=$(ls $CWD/$PROTEIN/)
  
  for prlst in ${PROTLST[@]} ; do
    mkdir -p $TMP/$ORTHO
    cd $TMP/$ORTHO
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    parallel $CORESCAPS pair_uniprot_vs_orthodb ::: $UNIPROTID
    echo ""
  done

  for prlst in ${PROTLST[@]} ; do
    cd $TMP/$ORTHO
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    echo -e "\nExtracting OGuniqueID for proteins in $prlst"
    parallel $CORESCAPS extract_oguniqueid ::: $UNIPROTID
    echo -e "\nGenerating summary report for proteins in $prlst"
    parallel $CORESCAPS report_gen ::: $UNIPROTID
  done

  for prlst in $CWD/$PROTEIN/* ; do
    summary_clarify $prlst
  done

  report_duplicates
  headers_insert
echo -e "\nDone with 1)"
;;

"Get orthologues (level: $LEVEL)")
if [ -d "$TMP/$ORTHO" ]; then
  PROTLST=$(ls $CWD/$PROTEIN/)
  cd $TMP/$ORTHO/

  for prlst in ${PROTLST[@]} ; do
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    sed "/$ORGANISM/d" $CWD/$SPECIES | \
    while read -r taxidNCBI latinName; do
      parallel $CORESCAPS prepare_orthologues_list ::: $UNIPROTID
    done
  done
   
  echo -e "Checking for proteins with no homologues..."
  parallel $CORESCAPS no_homologues_check ::: $UNIPROTID
  
  for prlst in ${PROTLST[@]} ; do
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    parallel $CORESCAPS get_ortho_fasta ::: $UNIPROTID
  done
else
  echo "NOT FOUND: $TMP/$ORTHO"
fi  
echo -e "\nDone with 2)"
;;

"Download sequences from UniProt (organism: $ORGANISM)")
if [ -d "$TMP/$ORTHO" ]; then
  cd $TMP/$ORTHO/
  uniprot_download
  cd -
else
 echo "NOT FOUND: $TMP/$ORTHO"
fi
echo -e "\nDone with 3)"
;;

"BLAST orthologues against UniProt sequence ($ORGANISM, detailed: $DETBLAST)")
if [ -d "$TMP/$ORTHO" ]; then
  PROTLST=$(ls $CWD/$PROTEIN/)
  cd $TMP/$ORTHO/
  for prlst in ${PROTLST[@]} ; do
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    parallel $CORESCAPS blast_db_prep ::: $UNIPROTID
    parallel $CORESCAPS reciprocal_blast ::: $UNIPROTID
  done
else
 echo "NOT FOUND: $TMP/$ORTHO"
fi
echo -e "\nDone with 4)"
;;

"Get FASTA sequences of the best hits (identity: $PIDENT; gaps: $PGAPS)")
if [ -d "$TMP/$ORTHO" ]; then
  PROTLST=$(ls $CWD/$PROTEIN/)
  cd $TMP/$ORTHO/
  mkdir -p $TMP/$GETFA
  for prlst in ${PROTLST[@]} ; do
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    parallel $CORESCAPS best_hits ::: $UNIPROTID
    parallel $CORESCAPS species_names ::: $UNIPROTID
  done

  clarify_blast
else
 echo "NOT FOUND: $TMP/$ORTHO"
fi
echo -e "\nDone with 5)"
;;

"[MSA] Exclude too divergent sequences")
if [ -d "$TMP/$GETFA" ]; then
  PROTLST=$(ls $CWD/$PROTEIN/)
  cd $TMP/$GETFA
  rm -rf $TMP/tsv/excluded-$GUIDANCEMSA-$GUIDANCECUT.tsv
  rm -rf $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT
  mkdir -p $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT/
  rm -rf $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT
  mkdir -p $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/
  for prlst in ${PROTLST[@]} ; do
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
    parallel $CORESCAPS run_guidance ::: $UNIPROTID
  done
else
 echo "NOT FOUND: $TMP/$GETFA"
fi
echo -e "\nDone with 6)"
;;

"[MSA] Create MSA with selected method ($MSAMETHOD)")
if [ -d "$TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT" ]; then
  PROTLST=$(ls $CWD/$PROTEIN/)
  cd $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT
  mkdir -p $TMP/$MSA/${GUIDANCEMSA}-${GUIDANCECUT}-${MSAMETHOD}/
  for prlst in ${PROTLST[@]} ; do
    UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
  
    if [ "$MSAMETHOD" = "muscle" ]; then
      parallel $CORESMUSCLE musclefn ::: $UNIPROTID
    elif [ "$MSAMETHOD" = "prank" ]; then
      prankprep
      parallel $CORESPRANK prankfn ::: $UNIPROTID
    elif [ "$MSAMETHOD" = "mafft" ] || [ "$MSAMETHOD" = "mafft-linsi" ] || [ "$MSAMETHOD" = "mafft-ginsi" ] || [ "$MSAMETHOD" = "mafft-einsi" ] || [ "$MSAMETHOD" = "mafft-fftns" ] || [ "$MSAMETHOD" = "mafft-fftnsi" ]; then
      parallel $CORESCAPS mafftfn ::: $UNIPROTID
    else
      echo "[ERROR] Specify MSA method properly!"
    fi
  done

  cd $TMP/$MSA/${GUIDANCEMSA}-${GUIDANCECUT}-${MSAMETHOD}/
  msa_process
  position_reference
else
  echo "NOT FOUND: $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT"
fi
echo -e "\nDone with 7)"
;;

"[TRE] Prepare trees ($TREESCAPS, $MSAMETHOD, $PHYMLGUIDE, $TREESROOT)")
if [ -d "$TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD" ]; then
  if [ "$TREESCAPS" = "auto" ]; then
    echo -e "CAPS will generate its own trees. Skipping..."
  elif [ "$TREESCAPS" = "phyml" ] && [ "$PHYMLGUIDE" = "exguide" ]; then
    phyml_ext
    phyml_guide
    phyml_prep
  
    PROTLST=$(ls $CWD/$PROTEIN/)
    for prlst in ${PROTLST[@]} ; do
      UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
      parallel $CORESCAPS phymlfn ::: $UNIPROTID
    done

    phyml_process
  
  elif [ "$TREESCAPS" = "phyml" ] && [ "$PHYMLGUIDE" = "noguide" ]; then
    phyml_prep
    
    PROTLST=$(ls $CWD/$PROTEIN/)
    for prlst in ${PROTLST[@]} ; do
      UNIPROTID=$( awk '{print $1}' $CWD/$PROTEIN/$prlst | datamash transpose )
      parallel $CORESCAPS phymlfn ::: $UNIPROTID
    done

    phyml_process
  else
    echo -e "Check your tree settings!"
  fi
else
  echo "NOT FOUND: $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD"
fi
echo -e "\nDone with 8)"
;;

"[RUN] Create pairs ($PAIRINGMANNER)")
if [ -d "$TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD" ]; then
  if [ "$PAIRINGMANNER" = "all" ]; then
    pair_msa
    pair_tree
  elif [ "$PAIRINGMANNER" = "defined" ]; then
    pair_defined_msa
    pair_tree
  else
    echo -e "Check your pairing settings!"
  fi
else
  echo "NOT FOUND: $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD"
fi
echo -e "\nDone with 9"
;;

"[RUN] CAPS run (alpha: $ALPHA, $MSAMETHOD, $TREESCAPS)")
if [ -d "$TMP/$PAIRM" ]; then
 split_dirs
 caps_set
  cd $TMP/$CAPSM
  for folder in * ; do
    cd $folder
    echo -e "Processing \e[92m${folder}\e[39m"
    echo "$DATESTAMP ${folder}" >> $TMP/progress-$ALPHA-$MSAMETHOD-$TREESCAPS.txt
    PAIRLIST=$(ls)
    msa=msa
    parallel $CORESCAPS --progress capsfn ::: "$PAIRLIST"
    echo -e "Done in \e[92m${folder}\e[39m"
    echo "" >> $TMP/progress-$ALPHA-$MSAMETHOD-$TREESCAPS.txt
    cd ..
  done
else
  echo "NOT FOUND: $TMP/$PAIRM"
fi
echo -e "\nDone with 10"
;;

"[RUN] Inspect results and re-run CAPS (REV)")
if [ -d "$TMP/$CAPSM" ]; then
  mkdir -p $TMP/$RESULTS/{fail,nocoev,coev}
  cd $TMP/$CAPSM/
  for folder in * ; do
    echo -e "Processing $folder"
    cd $folder
    PAIRLIST=$( ls ./ )
    parallel $CORESCAPS caps_inspect ::: "$PAIRLIST"
    cd ..
  done
  
  # Prepare for REV CAPS run
  caps_set_rev

  cd $TMP/$RESULTS/coev
  for folder in * ; do
    cd $folder
    echo -e "Processing \e[92m${folder}\e[39m"
    echo "$DATESTAMP ${folder}" >> $TMP/progress-$ALPHA-$MSAMETHOD-$TREESCAPS.txt
    PAIRLIST=$(ls)
    msa=msa-rev
    parallel $CORESCAPS --progress capsfn ::: "$PAIRLIST"
    echo -e "Done in \e[92m${folder}\e[39m"
    echo "" >> $TMP/progress-$ALPHA-$MSAMETHOD-$TREESCAPS.txt
    cd ..
  done
else
  echo "NOT FOUND: $TMP/$CAPSM"
fi

echo -e "\nDone with 11"
;;

"[RES] Generate columns stats")
if [ -d "$TMP/$RESULTS/coev" ]; then
  cd $TMP/$RESULTS/coev
  for folder in * ; do
  echo -e "Processing $folder"
    cd $folder
    PAIRLIST=$( ls ./ )
    parallel $CORESCAPS caps_reinspect ::: "$PAIRLIST"
    cd ..
  done
  echo -e "\n\e[92mResults inspections done!\e[39m\n"

  cd $TMP/$RESULTS/coev/
  for resfold in * ; do
    echo -e "[PROCESSING] $resfold"
    cd $resfold
    SUBFOLD=$( ls ./ )
    parallel $CORESCAPS results_cleanup ::: "$SUBFOLD"
    cd ..
  done

  for resfold in * ; do
    echo -e "[PROCESSING] $resfold"
    cd $resfold
    SUBFOLD=$( ls ./ )
    parallel $CORESCAPS extract_columns_stats ::: "$SUBFOLD"
    cd ..
  done

  mkdir -p $TMP/$RESULTS/chi/{back_calc_final,chi_test_final}
  cd $TMP/$RESULTS/chi/proteins
  protInt=$( ls ./ )
  parallel $CORESCAPS calc_back_final ::: "$protInt"

  cd $TMP/$RESULTS/coev/
  for resfold in * ; do
    echo -e "[PROCESSING] $resfold"
    cd $resfold
    SUBFOLD=$( ls ./ )
    parallel $CORESCAPS protein_pairs_stats ::: "$SUBFOLD"
    cd ..
  done

  summary_cleanup

else
  echo "NOT FOUND: $TMP/$RESULTS/coev"
fi

echo -e "\nDone with 12"
;;

"[Exit script]")
echo "Quitting..."
exit
;;

*)
echo "Invalid option"
;;

esac
REPLY=
done

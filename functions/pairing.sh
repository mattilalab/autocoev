#!/bin/bash

MSASUFF="fa" # define MSA extension

# Construct the PAIRS folder
if [ "$TREESCAPS" = "phyml" ]; then
  PAIRM="Pairs-${PAIRINGMANNER}/${MSAMETHOD}..PhyML_${PHYMLGUIDE}-${TREESROOT}"
  PAIRL="Pairs-${PAIRINGMANNER}-excluded/${MSAMETHOD}..PhyML_${PHYMLGUIDE}${TREESROOT}"
  CAPSM="CAPS-${PAIRINGMANNER}/${MSAMETHOD}..PhyML_${PHYMLGUIDE}-${TREESROOT}/Alpha${ALPHA}"
  RESULTS="Results-${PAIRINGMANNER}/${MSAMETHOD}..PhyML_${PHYMLGUIDE}-${TREESROOT}/Alpha${ALPHA}"
elif [ "$TREESCAPS" = "auto" ]; then
  PAIRM="Pairs-${PAIRINGMANNER}/${MSAMETHOD}..TREE_auto"
  PAIRL="Pairs-${PAIRINGMANNER}-excluded/${MSAMETHOD}..TREE_auto"
  CAPSM="CAPS-${PAIRINGMANNER}/${MSAMETHOD}..TREE_auto/Alpha${ALPHA}"
  RESULTS="Results-${PAIRINGMANNER}/${MSAMETHOD}..TREE_auto/Alpha${ALPHA}"
else
  echo -e "Check your CAPS trees settings ($TREESCAPS)!"
fi

# Loop through the list defined above twice in order to create sequence pairs.
# We can use an inequality to handle this; only display the filename if the
# first file comes before the second file alphabetically. This will ensure only
# one of each matches. Help is from: https://unix.stackexchange.com/a/490660
# Also, pick the corresponding tree from the "tree pool" defined by $TRIMTRE.
# This is the only way CAPS will use the user-supplied trees.
# Output only sequences in common!
pair_msa() {
  mkdir -p $TMP/$PAIRM

  echo -e "Protein1\tProtein2\tCommon" > $TMP/tsv/pairsGood.tsv
  echo -e "Protein1\tProtein2\tCommon" > $TMP/tsv/pairsBad.tsv

  cd $TMP/$MSA/$MSAMETHOD
  MSALIST=$(ls *.${MSASUFF})

  for i in ${MSALIST[@]} ; do
    nami=$(basename $i .${MSASUFF})
    for j in ${MSALIST[@]} ; do
      namj=$(basename $j .${MSASUFF})
      if [ "$i" \< "$j" ]; then
        COMMONSPECIES=$(seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$i $TMP/$MSA/$MSAMETHOD/$j | seqkit --threads ${THREADS} seq -n | wc -l)
        if [ "$COMMONSPECIES" -ge "$MINCOMMONSPCS" ]; then
          mkdir -p $TMP/$PAIRM/${nami}_vs_${namj}/msa
          seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$i $TMP/$MSA/$MSAMETHOD/$j -o $TMP/$PAIRM/${nami}_vs_${namj}/msa/$nami.fa
          seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$j $TMP/$MSA/$MSAMETHOD/$i -o $TMP/$PAIRM/${nami}_vs_${namj}/msa/$namj.fa
          seqkit --threads ${THREADS} seq $TMP/$PAIRM/${nami}_vs_${namj}/msa/$nami.fa -n > $TMP/$PAIRM/${nami}_vs_${namj}/$nami.species
          seqkit --threads ${THREADS} seq $TMP/$PAIRM/${nami}_vs_${namj}/msa/$namj.fa -n > $TMP/$PAIRM/${nami}_vs_${namj}/$namj.species
          echo -e "[\e[92mMORE\e[39m] Paired MSA ($MSAMETHOD): ${nami} and ${namj}"
          echo -e "$nami\t$namj\t$COMMONSPECIES" >> $TMP/tsv/pairsGood.tsv
        elif [ "$COMMONSPECIES" -lt "$MINCOMMONSPCS" ]; then
          mkdir -p $TMP/$PAIRL/${nami}_vs_${namj}/msa
          seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$i $TMP/$MSA/$MSAMETHOD/$j -o $TMP/$PAIRL/${nami}_vs_${namj}/msa/$nami.fa
          seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$j $TMP/$MSA/$MSAMETHOD/$i -o $TMP/$PAIRL/${nami}_vs_${namj}/msa/$namj.fa
          seqkit --threads ${THREADS} seq $TMP/$PAIRL/${nami}_vs_${namj}/msa/$nami.fa -n > $TMP/$PAIRL/${nami}_vs_${namj}/$nami.species
          seqkit --threads ${THREADS} seq $TMP/$PAIRL/${nami}_vs_${namj}/msa/$namj.fa -n > $TMP/$PAIRL/${nami}_vs_${namj}/$namj.species
          echo -e "[\e[37mLESS\e[39m] Paired MSA ($MSAMETHOD): ${nami} and ${namj}"
          echo -e "$nami\t$namj\t$COMMONSPECIES" >> $TMP/tsv/pairsBad.tsv
        else
          echo -e "\e[91mSomething went wrong in pairing. Check!\e[39m"
        fi
      fi
    done
  done
}

pair_defined_msa() {
  mkdir -p $TMP/$PAIRM

  echo -e "Protein1\tProtein2\tCommon" > $TMP/tsv/pairsGood.tsv
  echo -e "Protein1\tProtein2\tCommon" > $TMP/tsv/pairsBad.tsv

  cd $TMP/$MSA/$MSAMETHOD
  while read -r nami namj ; do
    i="${nami}.${MSASUFF}"
    j="${namj}.${MSASUFF}"
    if [ -f "$i" ] && [ -f "$j" ]; then
      COMMONSPECIES=$(seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$i $TMP/$MSA/$MSAMETHOD/$j | seqkit --threads ${THREADS} seq -n | wc -l)
      if [ "$COMMONSPECIES" -ge "$MINCOMMONSPCS" ]; then
        mkdir -p $TMP/$PAIRM/${nami}_vs_${namj}/msa
        seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$i $TMP/$MSA/$MSAMETHOD/$j -o $TMP/$PAIRM/${nami}_vs_${namj}/msa/$nami.fa
        seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$j $TMP/$MSA/$MSAMETHOD/$i -o $TMP/$PAIRM/${nami}_vs_${namj}/msa/$namj.fa
        seqkit --threads ${THREADS} seq $TMP/$PAIRM/${nami}_vs_${namj}/msa/$nami.fa -n > $TMP/$PAIRM/${nami}_vs_${namj}/$nami.species
        seqkit --threads ${THREADS} seq $TMP/$PAIRM/${nami}_vs_${namj}/msa/$namj.fa -n > $TMP/$PAIRM/${nami}_vs_${namj}/$namj.species
        echo -e "[\e[92mMORE\e[39m] Paired MSA ($MSAMETHOD): ${nami} and ${namj}"
        echo -e "$nami\t$namj\t$COMMONSPECIES" >> $TMP/tsv/pairsGood.tsv
      elif [ "$COMMONSPECIES" -lt "$MINCOMMONSPCS" ]; then
        mkdir -p $TMP/$PAIRL/${nami}_vs_${namj}/msa
        seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$i $TMP/$MSA/$MSAMETHOD/$j -o $TMP/$PAIRL/${nami}_vs_${namj}/msa/$nami.fa
        seqkit --threads ${THREADS} common $TMP/$MSA/$MSAMETHOD/$j $TMP/$MSA/$MSAMETHOD/$i -o $TMP/$PAIRL/${nami}_vs_${namj}/msa/$namj.fa
        seqkit --threads ${THREADS} seq $TMP/$PAIRL/${nami}_vs_${namj}/msa/$nami.fa -n > $TMP/$PAIRL/${nami}_vs_${namj}/$nami.species
        seqkit --threads ${THREADS} seq $TMP/$PAIRL/${nami}_vs_${namj}/msa/$namj.fa -n > $TMP/$PAIRL/${nami}_vs_${namj}/$namj.species
        echo -e "[\e[37mLESS\e[39m] Paired MSA ($MSAMETHOD): ${nami} and ${namj}"
        echo -e "$nami\t$namj\t$COMMONSPECIES" >> $TMP/tsv/pairsBad.tsv
      else
        echo -e "\e[91mSomething went wrong with pairing. Check!\e[39m"
      fi
    else
      echo -e "[SKIP] $nami and $namj from $PAIRLST not found among your proteins."
    fi
  done < $CWD/$PAIRLST
}

pair_tree() {

  if [ "$TREESCAPS" = "auto" ]; then
    echo -e "\nCAPS will generate its own trees.\n"
    return 1
  elif [ "$TREESCAPS" = "phyml" ] && [ "$TREESROOT" = "noroot" ]; then
    TREPATH="$TMP/$TREES/phyml/$MSATRETAG/noroot"
  elif [ "$TREESCAPS" = "phyml" ] && [ "$TREESROOT" = "rooted" ]; then
    TREPATH="$TMP/$TREES/phyml/$MSATRETAG/rooted"
  elif [ "$TREESCAPS" = "external" ] && [ "$ETREEROOT" = "noroot" ]; then
    TREPATH="$TMP/$TREES/external/noroot"
  elif [ "$TREESCAPS" = "external" ] && [ "$ETREEROOT" = "rooted" ]; then
    TREPATH="$TMP/$TREES/external/rooted"
  else
    echo -e "\nSomething went wrong with pairing.\n"
  fi

  sed 1d $TMP/tsv/pairsGood.tsv | \
  while read -r Protein1 Protein2 Common ; do
    mkdir -p $TMP/$PAIRM/${Protein1}_vs_${Protein2}/tre
    treebest subtree $TREPATH/${Protein1}.tre $TMP/$PAIRM/${Protein1}_vs_${Protein2}/$Protein1.species > $TMP/$PAIRM/${Protein1}_vs_${Protein2}/tre/$Protein1.tre
    treebest subtree $TREPATH/${Protein2}.tre $TMP/$PAIRM/${Protein1}_vs_${Protein2}/$Protein2.species > $TMP/$PAIRM/${Protein1}_vs_${Protein2}/tre/$Protein2.tre
    echo -e "Paired trees ($TREESCAPS): ${Protein1} and ${Protein2}"
  done
  
}

# Determine how many pairs there are
# https://linuxize.com/post/bash-increment-decrement-variable/
# http://en.tldp.org/HOWTO/Bash-Prompt-HOWTO/x700.html
# https://apple.stackexchange.com/questions/336427/copy-the-first-n-files-from-one-directory-to-another
split_dirs(){
  cd $TMP/$PAIRM
  PNUM=$(ls -1 | wc -l) # Count the number of folders
  FIRST=0

  until [ $FIRST -gt $PNUM ] ; do
    echo "Making $CAPSM/$FIRST"
    mkdir -p $TMP/$CAPSM/$FIRST
    COPY=$(ls | tail -n +${FIRST} | head -$INCR | xargs)
    cp -a $COPY $TMP/$CAPSM/$FIRST/
    (( FIRST=FIRST+$INCR ))
  done

  cd $TMP
  echo -e "\nSplit pairs in groups of $INCR are placed in:\n\e[96m$TMP/$CAPSM\e[39m"
}

#!/bin/bash

TREES="Trees"

# Construct the folder where trees will be calculated
MSATRETAG="$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD-${PHYMLGUIDE}"
PHYMLCORES="--j ${THREADS}" # "12"

# Trim the external tree following the species in each MSA. These can be
# used as guide trees for PhyML or as 'external' trees in the CAPS run.
# Use with care! This usually results in a very high number of detected
# co-evolving pairs and does not seem to be the recommended way!
ext_trees(){
  mkdir -p $TMP/$TREES/external/{noroot,rooted}
  cp $CWD/$EXTTREE $TMP/$TREES/external
  cat $CWD/$SPECIES | while read -r a b ; do
    sed -i "s/$b/$a/g" $TMP/$TREES/external/$EXTTREE
  done
  echo -e "Species names in the external tree changed to TAXID!"
  
  cd $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD
  SUBSPEC=$( ls *.species )
  for spcs in ${SUBSPEC[@]} ; do
    treebest subtree $TMP/$TREES/external/$EXTTREE $spcs > $TMP/$TREES/external/noroot/${spcs%.*}.tre
    treebest root $TMP/$TREES/external/noroot/${spcs%.*}.tre > $TMP/$TREES/external/rooted/${spcs%.*}.tre
    echo -e "Trimmed external tree for $spcs"
  done
}

phyml_ext(){
  mkdir -p $TMP/$TREES/phyml/ext
  cp $CWD/$EXTTREE $TMP/$TREES/phyml
  cat $CWD/$SPECIES | while read -r a b ; do
    sed -i "s/$b/$a/g" $TMP/$TREES/phyml/$EXTTREE
  done
  echo -e "Species names in the external tree changed to TAXID!"
}

phyml_guide(){
  cd $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD
  SUBSPEC=$( ls *.species )
  for spcs in ${SUBSPEC[@]} ; do
    treebest subtree $TMP/$TREES/phyml/$EXTTREE $spcs > $TMP/$TREES/phyml/ext/${spcs%.*}.tre
    echo -e "Trim external tree for $spcs"
  done
}

# Copy the necessary MSA files, converting them to PHYLIPI format.
phyml_prep(){
  mkdir -p $TMP/$TREES/phyml/$MSATRETAG/phy/
  cd $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD
  VANMSA=$(ls *.fa)
  for vmsa in ${VANMSA[@]} ; do
    squizz $vmsa -c PHYLIPI > $TMP/$TREES/phyml/$MSATRETAG/phy/${vmsa%.*}.phy
  done
}

phymlfn() {
  local PTRE="${1}"
  cd $TMP/$TREES/phyml/$MSATRETAG/phy
  if [ -s $PTRE.phy ] ; then
    if [ "$PHYMLGUIDE" = "exguide" ]; then
      local ETRE="--inputtree $TMP/$TREES/phyml/ext/${PTRE%.*}.tre"
      echo -e "PhyML will use $TMP/$TREES/phyml/ext/${PTRE%.*}.tre!"
    elif [ "$PHYMLGUIDE" = "noguide" ]; then
      ETRE=""
      echo -e "PhyML will not use a guide tree!"
    else
      echo -e "Check you input tree settings!"
    fi

  phyml -d aa \
    $PHYMLOPTIONS $ETRE \
    --no_memory_check \
    --leave_duplicates \
    -i $PTRE.phy
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${PTRE}"
  fi 
}

phyml_process() {
  mkdir -p $TMP/$TREES/phyml/$MSATRETAG/{rooted,noroot}
  cd $TMP/$TREES/phyml/$MSATRETAG/phy
  TREESOUT=$( ls *.phy_phyml_tree.txt )
  for t in ${TREESOUT[@]} ; do
    tname=$(basename $t .phy_phyml_tree.txt)
    cp $tname.phy_phyml_tree.txt $TMP/$TREES/phyml/$MSATRETAG/noroot/$tname.tre
    echo -e "$tname: Copy non-rooted"
    treebest root $tname.phy_phyml_tree.txt > $TMP/$TREES/phyml/$MSATRETAG/rooted/$tname.tre
    echo -e "$tname: Root and copy"
  done
}

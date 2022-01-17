#!/bin/bash

MSA="MSA" # MSA dir

CORESPRANK="--j ${THREADS}"
CORESMUSCLE="--j ${THREADS}"
THREADMAFFT="--thread ${THREADS}"

musclefn() {
  local SEQMSA="${1}"
  if [ -s $SEQMSA.fa ] ; then
    $MSAMETHOD $MUSCLEOPTIONS > $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/$SEQMSA.fa < $SEQMSA.fa
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${SEQMSA}"
  fi 
}

prankprep() {
  if [ "$PRANKGUIDE" = "exguide" ]; then
    mkdir -p $TMP/$MSA
    cp $CWD/$EXTTREE $TMP/$MSA
    cat $CWD/$SPECIES | while read -r a b ; do
    sed -i "s/$b/$a/g" $TMP/$MSA/$EXTTREE
    done
    echo -e "Species names in the external tree changed to TAXID!"
    PRGUID="-t=$TMP/$MSA/$EXTTREE -prunetree"
  else
    PRGUID=""
  fi
}

prankfn() {
  local SEQMSA="${1}"
  if [ -s $SEQMSA.fa ] ; then
    $MSAMETHOD $PRGUID $PRANKOPTIONS \
      -o=$TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/$SEQMSA.fa \
      -d=$SEQMSA.fa
    mv $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/$SEQMSA.fa.best.fas $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/$SEQMSA.fa
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${SEQMSA}"
  fi 
}

# MAFFT has a nice multithreading option, which we take advantage of
mafftfn() {
  local SEQMSA="${1}"
  if [ -s $SEQMSA.fa ] ; then
    $MSAMETHOD $MAFFTOPTIONS --anysymbol $SEQMSA.fa > $TMP/$MSA/$GUIDANCEMSA-$GUIDANCECUT-$MSAMETHOD/$SEQMSA.fa
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${SEQMSA}"
  fi 
}

# Run Gblocks and create a list of species of each MSA. Some help from:
# https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
# Also, make a list of species for each MSA
msa_process(){
  #cd $TMP/$MSA/$MSAMETHOD/$subcat
  PROCESSMSA=$(ls *.fa)
  for van in ${PROCESSMSA[@]} ; do
    Gblocks $van -t=p $GBLOCKSOPT -p=t -e=".gbl"
    sed -i "s: ::g" $van.gbl
    seqkit --threads ${THREADS} seq $van -n > ${van%.*}.species
  done
}

# Extract positions between vanilla and Gblocks filtered sequences from
# $ORGANISM. Number the first "vanilla" column, excluding gaps ("-"),
# then number the third "gblocks-filtered" column ("#"). Help from here:
# https://stackoverflow.com/questions/8696751/add-space-between-every-letter
# https://stackoverflow.com/questions/16750911/count-line-lengths-in-file-using-command-line-tools
# https://www.gnu.org/software/gawk/manual/html_node/Using-Shell-Variables.html
# https://www.shortcutfoo.com/app/dojos/awk/cheatsheet
position_reference(){
  #cd $TMP/$MSA/$MSAMETHOD/$subcat
  PROCESSMSA=$(ls *.fa)
  for van in ${PROCESSMSA[@]} ; do
    sed -n -e "/$ORGANISM/p" $van.gbl.txt | tr -d '\n' | tr -d ' ' | tr -d "$ORGANISM" | sed -e '$a\' > $van.$ORGANISM
    sed -n -e '/Gblocks  /p' $van.gbl.txt | tr -d '\n' | tr -d ' ' | tr -d 'Gblocks' >> $van.$ORGANISM
    PRLNGT=$(head -n1 $van.$ORGANISM | awk '{print length}')
    sed -i 's/./& /g' $van.$ORGANISM
    #cat $van.$ORGANISM | datamash -W transpose | awk '$1 !~ /-/' | awk -v BIGNUM="$PRLNGT" '{print ((NR-1)%BIGNUM)+1, $0}' | awk '$3 ~ /#/' | awk -v BIGNUM="$PRLNGT" '{print ((NR-1)%BIGNUM)+1, $0}' > $van.$ORGANISM.col
    cat $van.$ORGANISM | datamash -W transpose | awk '$1 !~ /-/' | awk -v BIGNUM="3000" '{print ((NR-1)%BIGNUM)+1, $0}' > $van.$ORGANISM.ref
    rm $van.$ORGANISM
    #sed -i "s: :\t:g" $van.$ORGANISM.col
    sed -i "s: :\t:g" $van.$ORGANISM.ref
  done
}

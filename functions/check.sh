#!/bin/bash

EXELIST=( blastp blast_formatter $CAPS datamash fastafetch \
          fastaindex Gblocks mafft makeblastdb muscle parallel phyml \
	  prank seqkit squizz treebest )

check_bin(){
  for e in ${EXELIST[@]} ; do
    BINPATH=$(which $e)
    if [ -x "$BINPATH" ] ; then
      echo -e "[\e[92mIN PATH\e[39m] $BINPATH"
    else
      echo -e "[\e[91mMISSING\e[39m] $e not found!"
    fi
  done
}

preview_tab(){
  local PRVW="${1}"
  echo -e "Preview of \e[96m${PRVW}\e[39m"
  head -n 5 $PRVW
  echo ""
}

report_tsv(){
  mkdir -p $TMP/tsv
  touch $TMP/tsv/{OrthoDB_Missing.tsv,proteinsFound.tsv,Summary.tsv,duplicates_OrthoDB.tsv,duplicates_UniProt.tsv,duplicates_OrthoGroup.tsv,blastBestFasta.tsv,blastBestExclude.tsv}
}

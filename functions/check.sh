#!/bin/bash

EXELIST=( blastp blast_formatter $CAPS datamash extractalign fastafetch \
          fastaindex Gblocks mafft makeblastdb muscle parallel phyml \
	  prank seqkit squizz treebest )

check_bin(){
  for e in ${EXELIST[@]} ; do
    BINPATH=$(which $e)
    if [ -x "$BINPATH" ] ; then
      echo -e "[\e[92mFOUND\e[39m] $BINPATH"
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

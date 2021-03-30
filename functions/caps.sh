#!/bin/bash

# CAPS executable
#CAPS=caps
CAPS=vCAPS

# Number of threads to use for CAPS
CORESCAPS="--j ${THREADS}"

# Determine what trees to use
caps_set() {
  if [ "$TREESCAPS" = "phyml" ]; then
    TREMSA="-T ./tre"
    echo -e "\nCAPS will run with PhyML trees!"
  elif [ "$TREESCAPS" = "external" ]; then
    TREMSA="-T ./tre"
    echo -e "\nCAPS will run with external trees!"
  elif [ "$TREESCAPS" = "auto" ]; then
    TREMSA=""
    echo -e "\nCAPS will estimate its own trees!"
  else
  echo -e "\nSometing is wrong with your trees!"
  fi 
}

# Set a CAPS function with user-defined variables to be executed via GNU
# Parallel for each individual folder of MSA pairs.
capsfn() {
  local PAIR="${1}"
  cd $PAIR
  $CAPS \
    --inter \
    -a ${ALPHA} \
    -b ${BOOT} \
    ${ALNS} \
    ${CONV} \
    ${GAPS} \
    ${REFER} \
    ${CAPSOPTIONS} \
    ${TREMSA} \
    -F ./msa 2> log.txt
  cd ..
  echo -e Processed "$PAIR"
  echo "$PAIR" >> $TMP/progress-$ALPHA-$MSAMETHOD-$TREESCAPS.txt
}

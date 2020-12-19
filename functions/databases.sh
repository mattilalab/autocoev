#!/bin/bash

# Define trimmed databases
GENEXREF="${GENEXREFALL}.${ORGANISM}" # Genes from a certain organism (e.g. ORGANISM="10090" for mouse)
OG2GENES="${OG2GENESALL}.${LEVEL}"    # Orthologues groups db at specified level (e.g. LEVEL="32523")

# Download the databases
download_db(){
  local DBDW="${1}"
  mkdir -p $DTB
  cd $DTB
  if [ -s "${DBDW}.tab.gz" ]; then
    echo -e "[\e[93mPRESENT\e[39m] ${DBDW}.tab.gz is already in $DTB! Skipping..."
  elif [ ! -f ${DBDW}.tab.gz ] ; then
    echo -e "Downloading ${DBDW}.tab.gz in \e[96m$DTB\e[39m"
    wget -c "https://$ORTHODBVER.orthodb.org/download/${DBDW}.tab.gz"
  else
    echo -e "[\e[91mERROR\e[39m] Check $DTB and your settings!"
  fi  
}

# Check MD5SUMs of databases
md5sum_check(){
  local DBCH="${1}"
  local DB5S="${2}"
  mkdir -p $DTB
  cd $DTB
  if [ -f "$DBCH.tab.gz" ]; then
    m5s=$(md5sum $DBCH.tab.gz | awk '{print $1}' )
    if [ "$m5s" = "$DB5S" ]; then
      echo -e "[\e[92mMATCH\e[39m] $m5s = $DB5S $DBCH.tab.gz"
    elif [ "$m5s" != "$DB5S" ]; then
      echo -e "[\e[91mWRONG\e[39m] $m5s != $DB5S $DBCH.tab.gz"
    else
      echo -e "Check your settings"
    fi
  elif [ ! -f "$DBCH.tab.gz" ]; then
    echo -e "[\e[91mMISSING\e[39m] $DBCH.tab.gz not found! Download first!"
  else
    echo -e "[\e[91mERROR\e[39m] Check your settings!"
  fi
}

# Extract databases
extract_db(){
  local DBEX="${1}"
  mkdir -p $DTB
  cd $DTB
  if [ -s "${DBEX}.tab" ]; then
    echo -e "[\e[93mPRESENT\e[39m] ${DBEX}.tab is already in $DTB! Skipping..."
  elif [ ! -s "${DBEX}.tab.gz" ]; then
    echo -e "[\e[91mMISSING\e[39m] Archive ${DBEX}.tab.gz not found! Download first!"
  elif [ -s "${DBEX}.tab.gz" ]; then
    echo -e "Extracting ${DBEX}.tab.gz in \e[96m$DTB\e[39m"
    gunzip -c ${DBEX}.tab.gz > ${DBEX}.tab
    echo -e "[\e[92mEXTRACT\e[39m] Extracted ${DBEX}.tab! You may delete the archive."
  else
    echo -e "[\e[91mERROR\e[39m] Check your setings!"
  fi
}

# Index FASTA database
index_all_fasta(){
  mkdir -p $DTB
  cd $DTB
  if [ -s "$ALLFASTA.tab.index" ]; then
    echo -e "Indexing $ALLFASTA.tab in \e[96m$DTB\e[39m"
    fastaindex \
      $ALLFASTA.tab \
      $ALLFASTA.tab.index
    echo -e "[\e[92mINDEXED\e[39m] Indexing $ALLFASTA.tab complete!"
  elif [ -s "$ALLFASTA.tab.index" ]; then
    echo -e "[\e[93mPRESENT\e[39m] Indexed $ALLFASTA.tab already present in $DTB! Skipping..."
  elif [ ! -s "$ALLFASTA.tab" ]; then
    echo -e "[\e[91mMISSING\e[39m] No $ALLFASTA.tab found in $DTB! Download/extract first!"
  else
    echo -e "Check your settings!"
  fi
}

# Trim databases following organism or level
trim_db(){
  local DBTR="${1}"
  local IDTR="${2}"
  local SBTR="${3}"
  mkdir -p $DTB
  cd $DTB
  if [ -s "${DBTR}.${SBTR}.tab" ]; then
    echo -e "[\e[93mPRESENT\e[39m] ${DBTR}.${SBTR}.tab is already in $DTB! Skipping..."
  elif [ ! -s "${DBTR}.tab" ]; then
    echo -e "[\e[91mMISSING\e[39m] No ${DBTR}.tab found in $DTB! Download/extract first!"
  elif [ -s "${DBTR}.tab" ]; then
    echo -e "Trimming ${DBTR}.tab in \e[96m$DTB\e[39m"
    grep "$IDTR" ${DBTR}.tab > ${DBTR}.${SBTR}.tab
    echo -e "[\e[92mTRIMMED\e[39m] Trimmed ${DBTR}.tab for ${SBTR}!"
  else
    echo -e "Check your setings!"
  fi
}

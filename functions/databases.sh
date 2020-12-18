#!/bin/bash

# Define trimmed databases
GENEXREF="${GENEXREFALL}.${ORGANISM}" # Genes from a certain organism (e.g. ORGANISM="10090" for mouse)
OG2GENES="${OG2GENESALL}.${LEVEL}"    # Orthologues groups db at specified level (e.g. LEVEL="32523")

# Download the databases
download_db(){
  local DBDW="${1}"
  mkdir -p $DTB
  cd $DTB
  if [ -z "${DBDW}.tab.gz" ]; then
    echo -e "Downloading ${DBDW}.tab.gz in \e[96m$DTB\e[39m"
    wget -c "https://$ORTHODBVER.orthodb.org/download/${DBDW}.tab.gz"
  elif [ -f "${DBDW}.tab.gz" ]; then
    echo -e "${DBDW}.tab.gz is already in $DTB! Skipping..."
  else
    echo -e "Check your settings!"
  fi  
}

# Check if databases are downloaded
db_exist(){
  DBLIST=( "${GENEXREFALL}" "${OG2GENESALL}" "${ALLFASTA}" )
  mkdir -p $DTB
  cd $DTB
  for d in ${DBLIST[@]} ; do
    if [ -z "$d.tab.gz" ]; then
      echo -e "\e[91mMISSING\e[39m] $d.tab.gz not found! Download first!"
      return 1
    elif [ -f "$d.tab.gz" ]; then
      echo -e "[\e[92mFOUND\e[39m] $DTB/$d.tar.gz"
    else
      echo -e "Check your settings!"
      return 1
    fi
  done
}

# Check MD5SUMs of databases
md5sum_check(){
  local DBCH="${1}"
  local DB5S="${2}"
  mkdir -p $DTB
  cd $DTB
  m5s=$(md5sum $DBCH.tab.gz | awk '{print $1}' )
  if [ "$m5s" = "$DB5S" ]; then
    echo -e "[\e[92mMATCH\e[39m] $m5s = $DB5S $DBCH.tab.gz"
  elif [ "$m5s" != "$DB5S" ]; then
    echo -e "\e[91mWRONG\e[39m] $m5s != $DB5S $DBCH.tab.gz"
  else
    echo -e "Check your settings"
  fi
}

# Extract databases
extract_db(){
  local DBEX="${1}"
  mkdir -p $DTB
  cd $DTB
  if [ -z "${DBEX}.tab" ]; then
    echo -e "Extracting ${DBEX}.tab.gz in \e[96m$DTB\e[39m"
    gunzip -v -c ${DBEX}.tab.gz > ${DBEX}.tab
    echo -e "Extracted ${DBEX}.tab! You may delete the archive."
  elif [ -f "${DBEX}.tab" ]; then
    echo -e "${DBEX}.tab is already in $DTB! Skipping..."
  elif [ -z "${DBEX}.tab.gz" ]; then
    echo -e "Archive ${DBEX}.tab.gz not found! Download first!"
  else
    echo -e "Check your setings!"
  fi
}

# Index FASTA database
index_all_fasta(){
  mkdir -p $DTB
  cd $DTB
  if [ -z "$ALLFASTA.tab.index" ]; then
    echo -e "Indexing $ALLFASTA.tab in \e[96m$DTB\e[39m"
    fastaindex \
      $ALLFASTA.tab \
      $ALLFASTA.tab.index
    echo -e "Indexing $ALLFASTA.tab complete!"
  elif [ -f "$ALLFASTA.tab.index" ]; then
    echo -e "Indexed $ALLFASTA.tab already present in $DTB! Skipping..."
  elif [ -z "$ALLFASTA.tab" ]; then
    echo -e "No $ALLFASTA.tab found in $DTB! Download/extract first!"
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
  if [ -z "${DBTR}.${SBTR}.tab" ]; then
    echo -e "Trimming ${DBTR}.tab in \e[96m$DTB\e[39m"
    grep "$IDTR" ${DBTR}.${SBTR}.tab > ${DBTR}.${SBTR}.tab
    echo -e "Trimmed ${DBTR}.tab for ${SBTR}!"
  elif [ -f "${DBTR}.${SBTR}.tab" ]; then
    echo -e "${DBTR}.${SBTR}.tab is already in $DTB! Skipping..."
  elif [ -z "${DBTR}.tab" ]; then
    echo -e "No ${DBTR}.tab found in $DTB! Download/extract first!"
  else
    echo -e "Check your setings!"
  fi
}

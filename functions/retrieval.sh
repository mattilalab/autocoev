#!/bin/bash

# Where do we collect the orthologues?
ORTHO="Orthologues"

report_tsv(){
  mkdir -p $TMP/tsv
  touch $TMP/tsv/{OrthoDB_Missing.tsv,proteinsFound.tsv,Summary.tsv,duplicates_OrthoDB.tsv,duplicates_UniProt.tsv,duplicates_OrthoGroup.tsv,blastBestFasta.tsv,blastBestExclude.tsv}
}

# Pair UniProt ID to the corresponding OrthoDB ID. Help with defining word borders with grep:
# https://www.linuxquestions.org/questions/linux-newbie-8/how-to-grep-for-an-exact-word-4175439257/
pair_uniprot_vs_orthodb(){
  local UPOD="${1}"
  mkdir -p $UPOD
  if LANG=C grep -F -w "$UPOD" $DTB/$GENEXREF.tab | grep "\b${ORGANISM}_" >> $UPOD/UniProt_OrthoDB.tsv ; then
    echo -e "[\e[92mMATCHED\e[39m] ${UPOD}\t${prlst}"
  else
    echo -e "[\e[91mMISSING\e[39m] ${UPOD}\t${prlst}"
    echo -e "$UPOD" >> $TMP/tsv/OrthoDB_Missing.tsv
  fi
}

# Extract OGuniqueID. Use grep with LANG=C and -F option to speed up the process (no strings accepted).
extract_oguniqueid() {
  local UPOD="${1}"
  if [ -s $UPOD/UniProt_OrthoDB.tsv ]; then
    while read -r OrthoDBgeneID UniProtID externlDBname ; do
      if LANG=C grep -F -w "$OrthoDBgeneID" $DTB/$OG2GENES.tab | grep "at${LEVEL}\b" >> $UPOD/OrthoGroup.tsv ; then
        echo -e "[\e[92mEXTRACT\e[39m] ${UPOD}\t${prlst}\t${OrthoDBgeneID}"
      else
        echo -e "[\e[91mMISSING\e[39m] ${UPOD}\t${prlst}\t${OrthoDBgeneID}"
      fi
    done < $UPOD/UniProt_OrthoDB.tsv
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UPOD}\t${prlst}\t${OrthoDBgeneID}"
  fi
}

# Generate summary report. Some help with awk from:
# https://unix.stackexchange.com/questions/222121/how-to-remove-a-column-or-multiple-columns-from-file-using-shell-command
report_gen(){
  local UPOD="${1}"
  if [ -s $UPOD/OrthoGroup.tsv ]; then
    cd $UPOD
    paste OrthoGroup.tsv UniProt_OrthoDB.tsv | awk '!($3=$5="")'>> $TMP/tsv/Summary.tsv
    echo -e "[\e[92mSUMMARY\e[39m] ${UPOD}\t${prlst}"
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UPOD}\t${prlst}"
    cd ..
  fi
}

summary_clarify(){
  local PSLS="${1}"
  while read -r UniProtID geneName assign ; do
    sed -i "s:$UniProtID:$UniProtID $geneName $assign:g" $TMP/tsv/Summary.tsv
    sed -i "s:  : :g" $TMP/tsv/Summary.tsv
    sed -i "s: :	:g" $TMP/tsv/Summary.tsv
  done < $PSLS
}

# Deal with duplicate entries
# https://askubuntu.com/questions/434545/identify-duplicate-lines-in-a-file-without-deleting-them
# https://unix.stackexchange.com/questions/224433/grep-first-column-uniq-values/224434
# https://superuser.com/questions/1092282/bash-sort-by-not-first-character
report_duplicates(){

  # Output duplicated Ortho Groups entries
  for dup in $(awk '{ print $1 }' $TMP/tsv/Summary.tsv | sort | uniq -d) ; do
    grep "$dup" $TMP/tsv/Summary.tsv >> $TMP/tsv/duplicates_OrthoGroup.tsv
  done
  
  # Output duplicated OrthoDB entries
  for dup in $(awk '{ print $2 }' $TMP/tsv/Summary.tsv | sort | uniq -d) ; do
    grep "$dup" $TMP/tsv/Summary.tsv >> $TMP/tsv/duplicates_OrthoDB.tsv
  done

  # Output duplicated UniProt entries
  for dup in $(awk '{ print $3 }' $TMP/tsv/Summary.tsv | sort |  uniq -d) ; do
    grep "$dup" $TMP/tsv/Summary.tsv >> $TMP/tsv/duplicates_UniProt.tsv
  done
  
  # Output found proteins
  sort --key 3 -u $TMP/tsv/Summary.tsv | sort --key 4 | awk '{ print $3, $4, $5 }' >> $TMP/tsv/proteinsFound.tsv
  
  echo -e "Reporting duplicates..."
}

# Output dir and column headers. Do this in the end, so headers do not mess with the while loops
headers_insert(){
  sed "1i UniProt\tName\tAssign" -i $TMP/tsv/OrthoDB_Missing.tsv
  sed "1i UniProt\tName\tAssign" -i $TMP/tsv/proteinsFound.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/Summary.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/duplicates_OrthoDB.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/duplicates_UniProt.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/duplicates_OrthoGroup.tsv
}

# Prepare the list of orthologues, collect all identifiers for all isoforms for each protein.
# Skip reference organism (e.g. mouse), as we will get the sequences from UniProt anyway.
prepare_orthologues_list() {
  local UniProt="${1}"
  if [ -s $UniProt/OrthoGroup.tsv ]; then
    while read -r OGuniqueID OrthoDBgeneID ; do
      if LANG=C grep -F -w "$OGuniqueID" $DTB/$OG2GENES.tab | grep "\<${taxidNCBI}_\B" >> $UniProt/orthologues.tsv ; then
        echo -e "[\e[92mMATCHED\e[39m] ${UniProt}\t${prlst}\t$OGuniqueID\t$latinName"
        echo -e "$taxidNCBI\t$latinName" >> $UniProt/${OrthoDBgeneID}.speciesFound.tsv
      else
        echo -e "[\e[91mMISSING\e[39m] ${UniProt}\t${prlst}\t$OGuniqueID\t$latinName"
        echo -e "$taxidNCBI\t$latinName" >> $UniProt/${OrthoDBgeneID}.speciesMissing.tsv
      fi
    done < $UniProt/OrthoGroup.tsv
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}\t${prlst}"
  fi
}

# Check if no homologues were found. These will be removed
no_homologues_check(){
  find $TMP/$ORTHO/*/orthologues.tsv -size 0 -print >> $TMP/tsv/homologues_None.tsv \
  -exec rm -rf $UniProt {} \;
}

# Get the fasta sequences of orthologues all isoforms. Remove any duplicates.
# https://stackoverflow.com/questions/2664740/extract-file-basename-without-path-and-extension-in-bash
get_ortho_fasta() {
  local UniProt="${1}"
  if [ -s $UniProt/orthologues.tsv ]; then
    cat $UniProt/orthologues.tsv | sort | uniq | \
    while read -r CLIDatCLADE OrthoDBgeneID ; do
      mkdir -p $UniProt/FASTA
      TAXID=$(basename "$OrthoDBgeneID" | cut -d_ -f1)
      fastafetch -s \
        -f $DTB/$ALLFASTA.tab \
        -i $DTB/$ALLFASTA.tab.index \
        -q $OrthoDBgeneID >> $UniProt/FASTA/$TAXID.fa
      echo -e "[\e[92mFAFETCH\e[39m] $UniProt\t${prlst}\t$OrthoDBgeneID\t$TAXID"
    done
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}\t${prlst}"
  fi
}

# Retrieve the mouse protein sequences for each gene from UniProt and store them separately.
# Retry up to 50 times (or set "-t inf" for infinite retries). Check for empty files if
# download was unsuccessful. I can parallelize this, but the server may block the download.
uniprot_download() {
  for UniProt in * ; do
    if [ -s $UniProt/orthologues.tsv ]; then
      mkdir -p $UniProt/$ORGANISM
      wget \
        -t 50 \
        -c "https://www.uniprot.org/uniprot/${UniProt}.fasta" \
        -O $UniProt/$ORGANISM/$UniProt.fa
      find $UniProt/$ORGANISM/ -size 0 -print >> $TMP/tsv/UniProt.failed
    else
      echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}"
    fi
  done
}

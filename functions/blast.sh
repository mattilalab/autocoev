#!/bin/bash

GETFA="BestBLASTfasta" # Best pBLAST hits

# Create BLAST database for each mouse sequence, against which we will run the sequences of
# collected orthologues in order to determine best matching isoform. Use only the found proteins.
blast_db_prep() {
  local UniProt=${1}
  if [ -s $UniProt/orthologues.tsv ] ; then
    makeblastdb \
      -in $UniProt/$ORGANISM/$UniProt.fa \
      -parse_seqids \
      -blastdb_version 5 \
      -title "$UniProt" \
      -dbtype prot
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}"
  fi
}

# Run reciprocal BLASTP against mouse protein and sort by bit score. Then sort by "bit score"
# (12th column, --key 12). Output format info: http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
reciprocal_blast() {
  local UniProt="${1}"
  if [ -s $UniProt/orthologues.tsv ] ; then
    cd $UniProt/FASTA/
    for blfas in * ; do
      mkdir -p ../BLAST/${blfas%.*}
      blastp -num_threads $BLASTCORES \
        -query $blfas \
        -db ../$ORGANISM/$UniProt.fa \
        -out ../BLAST/${blfas%.*}/${blfas%.*}.out \
        -outfmt="6 qseqid sacc pident ppos length mismatch gapopen gaps qstart qend sstart send evalue bitscore"
      sort --key 14 --numeric-sort --reverse ../BLAST/${blfas%.*}/${blfas%.*}.out | head -n 1 >> ../blastBest.tsv
      echo -e "BLAST ${blfas%.*} against $UniProt short."

      # Shall we run BLAST again to generate detailed results?
      if [ "$DETBLAST" = "yes" ]; then
        blastp -num_threads $BLASTCORES \
          -query $blfas \
          -db ../$ORGANISM/$UniProt.fa \
          -out ../BLAST/${blfas%.*}/${blfas%.*}.default \
          -outfmt="0"
         echo -e "BLAST ${blfas%.*} against $UniProt detailed."
       elif [ "$DETBLAST" = "no" ]; then
         echo -e "BLAST ${blfas%.*} against $UniProt no details."
      else
        "Check your BLAST settings!"
      fi
    done
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}"
  fi
  echo ""
}

# Get the fasta sequences of the best hits with certain identity, defined by user. Help from:
# https://stackoverflow.com/questions/8654051/how-to-compare-two-floating-point-numbers-in-bash
best_hits() {
  local UniProt=${1}
  if [ -s $UniProt/orthologues.tsv ] ; then
    cat $UniProt/blastBest.tsv | \
    while read -r qseqid sacc pident ppos length mismatch gapopen gaps qstart qend sstart send evalue bitscore ; do
      PERCENTGAPS=$(echo "${gaps}/${length}*100" |bc -l)
      if (( $(echo "$pident > $PIDENT" |bc -l) && $(echo "$PERCENTGAPS < $PGAPS" |bc -l) )); then
        echo -e "[\e[92mGOOD\e[39m] Identity and gaps: $sacc <- $qseqid: $pident and $PERCENTGAPS"
        fastafetch \
          -f $DTB/$ALLFASTA.tab \
          -i $DTB/$ALLFASTA.tab.index \
          -q $qseqid >> $TMP/$GETFA/$UniProt.fa
        echo -e "$qseqid\t$sacc\t$pident\t$ppos\t$length\t$mismatch\t$gapopen\t$gaps\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore\t$PERCENTGAPS" >> $TMP/tsv/blastBestFasta.tsv
      else
        echo -e "[\e[91mPOOR\e[39m] Identity and gaps: $sacc <- $qseqid: $pident and $PERCENTGAPS"
        echo -e "$qseqid\t$sacc\t$pident\t$ppos\t$length\t$mismatch\t$gapopen\t$gaps\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore\t$PERCENTGAPS" >> $TMP/tsv/blastBestExclude.tsv
      fi
    done
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}"
  fi
}

# Deal with names
species_names() {
  local UniProt=${1}
  if [ -s $UniProt/orthologues.tsv ] ; then
    sed -i "s/_.*//" $TMP/$GETFA/$UniProt.fa
    echo ">${ORGANISM}" >> $TMP/$GETFA/$UniProt.fa
    sed 1d $UniProt/$ORGANISM/$UniProt.fa >> $TMP/$GETFA/$UniProt.fa
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${UniProt}"
  fi
}

clarify_blast(){
  while read -r UniProtID geneName assign ; do
    sed -i "s/$UniProtID/$UniProtID\t$geneName/g" $TMP/tsv/blastBestFasta.tsv
    sed -i "s/$UniProtID/$UniProtID\t$geneName/g" $TMP/tsv/blastBestExclude.tsv
  done < $TMP/tsv/proteinsFound.tsv
  
  while read -r taxidID speciesName ; do
    sed -i "s/${taxidID}_/$speciesName\t${taxidID}\t${taxidID}_/g" $TMP/tsv/blastBestFasta.tsv
    sed -i "s/${taxidID}_/$speciesName\t${taxidID}\t${taxidID}_/g" $TMP/tsv/blastBestExclude.tsv
  done < $CWD/$SPECIES
  
  sed -i "1i Species\ttaxid\tqseqid\tsacc\tname\tpident\tppos\tlength\tmismatch\tgapopen\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tgaps" $TMP/tsv/blastBestFasta.tsv
  sed -i "1i Species\ttaxid\tqseqid\tsacc\tname\tpident\tppos\tlength\tmismatch\tgapopen\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tgaps" $TMP/tsv/blastBestExclude.tsv
}

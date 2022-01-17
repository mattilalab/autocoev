#!/bin/bash

GUIDANCE="Guidance"
GUIDEDFA="BestGuidedFasta"

# Determine the MSA program to use with Guidance
if [ "$GUIDANCEMSA" = "muscle" ]; then
  MSACAPSLOCK="MUSCLE"
  MSABINBATH="--muscle"
elif [ "$GUIDANCEMSA" = "prank" ]; then
  MSACAPSLOCK="PRANK"
  MSABINBATH="--prank"
elif [ "$GUIDANCEMSA" = "mafft" ] || [ "$GUIDANCEMSA" = "mafft-linsi" ] || [ "$GUIDANCEMSA" = "mafft-ginsi" ] || [ "$GUIDANCEMSA" = "mafft-einsi" ] || [ "$GUIDANCEMSA" = "mafft-fftns" ] || [ "$GUIDANCEMSA" = "mafft-fftnsi" ]; then
  MSACAPSLOCK="MAFFT"
  MSABINPATH="--mafft"
else
  echo "Check your Guidance settings!"
fi


run_guidance(){
  local SEQMSA="${1}"
  if [ -s $SEQMSA.fa ] ; then
  perl /usr/lib/guidance/www/Guidance/guidance.pl \
    --seqFile $SEQMSA.fa \
    --program GUIDANCE \
    --msaProgram $MSACAPSLOCK \
    $MSABINPATH /usr/bin/$GUIDANCEMSA \
    --seqType aa \
    --seqCutoff $GUIDANCECUT \
    --colCutoff 0 \
    --outDir $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA
  else
    echo -e "[\e[37mNOTINDB\e[39m] ${SEQMSA}"
  fi
  mkdir -p $TMP/$GUIDEDFA/
  
  # Save the filtered sequences separately
  if [ -s $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ] ; then
    cp $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names \
    $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.fa
  else
    echo -e "[\e[91mEMPTY\e[39m] ${SEQMSA} - no sequences above $GUIDANCECUT"
  fi
  
  # Save the excluded sequences separately for reference
  if [ -s $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA/Seqs.Orig.fas.FIXED.Removed_Seq.With_Names ] ; then
    cp $TMP/$GUIDANCE/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA/Seqs.Orig.fas.FIXED.Removed_Seq.With_Names \
    $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.removed.fa
    EXCLUDEDSEQS=$(seqkit seq $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.removed.fa -n | datamash transpose)
    echo -e "${SEQMSA}\t${EXCLUDEDSEQS}" >> $TMP/tsv/excluded-$GUIDANCEMSA-$GUIDANCECUT.tsv
  else
    echo -e "${SEQMSA} - no sequences below $GUIDANCECUT"
  fi
  
  # Check if the reference organism happens to be among the excludes
  if [ -s $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.fa ] && [ -s $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.removed.fa ] ; then
    if LANG=C grep -F -w "$ORGANISM" $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.removed.fa ; then
      mv $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.fa $TMP/$GUIDEDFA/$GUIDANCEMSA-$GUIDANCECUT/$SEQMSA.fa.no-$ORGANISM
      echo -e "[\e[91mNOREFF\e[39m] reference organism $ORGANISM sequence excluded from $SEQMSA!"
    fi
  fi
}

#!/bin/bash

# Export all
set -a

# What date and time is it?
DATESTAMP=$(date)

# Current work dir
CWD=$(pwd)

# Load settings first
. $CWD/settings.conf

# Load functions
. $CWD/functions/databases.sh
. $CWD/functions/retrieval.sh
. $CWD/functions/blast.sh
. $CWD/functions/msa.sh
. $CWD/functions/trees.sh
. $CWD/functions/pairing.sh
. $CWD/functions/caps.sh
. $CWD/functions/results.sh
. $CWD/functions/network.sh
. $CWD/functions/check.sh

EGAPS="0.8" # Gaps score cut
EPDIF="0.0011"

ETSV="$TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv"
EOUT="filtered-P${PVALUE}-pDiff${EPDIF}-Bonf${BONFERRONI}.tsv"
ENET="FILTERED-P${PVALUE}-pDiff${EPDIF}-Bonf${BONFERRONI}.xml"

# Let user change work dir here as well
read -p "Change working dir or press ENTER to accept: " -e -i $TMP TMP
mkdir -p $TMP

# Specify the input file
echo -e "\n\e[96mEnter the full path to the pairs-P${PVALUE}-B${BONFERRONI}.tsv file:\e[39m"
read -e -i $ETSV ETSV

rm -rf $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}
mkdir -p $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}

sed 1d "$ETSV" | while read -r Name1 msa1 Name2 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot corrT cCoev pvalA pvalB pMean pDiff corr1 corr2 bonferroni GapA GapB GapAB DivA DivB DivAB ; do
  if [ "$flt1" = "#" ] && [ "$flt2" = "#" ]; then
    echo "Pair from a good quality MSA"
    if (( $(echo "$GapA > $EGAPS" |bc -l) && $(echo "$GapB > $EGAPS" |bc -l) && $(echo "$pDiff < $EPDIF" |bc -l) )); then
      echo "$Name1 $msa1 $Name2 $msa2 $cCoev $pvalA $pvalB $pMean $bonferroni $boot" >> $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}/$msa1-$msa2.tsv
      echo "Consider $Name1 $msa1 <-> $Name2 $msa2"
    else
      echo "Skipping $Name1 $msa1 <-> $Name2 $msa2"
    fi
  else
    echo "pair from a bad quality MSA"
  fi
done

cd $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}
  for pair in * ; do
        namesIDs=$(head -n 1 $pair | awk '{print $1,$2,$3,$4}')
      coevNumber=$(awk '{print $5}' $pair | datamash count 1)
      
        cCoevMAX=$(awk '{print $5}' $pair | datamash max 1)
       cCoevMEAN=$(awk '{print $5}' $pair | datamash mean 1)
        cCoevMIN=$(awk '{print $5}' $pair | datamash min 1)
       cCoevSDEV=$(awk '{print $5}' $pair | datamash sstdev 1)
       
         bootMAX=$(awk '{print $10}' $pair | datamash max 1)
        bootMEAN=$(awk '{print $10}' $pair | datamash mean 1)
         bootMIN=$(awk '{print $10}' $pair | datamash min 1)
	bootSDEV=$(awk '{print $10}' $pair | datamash sstdev 1)
       
        pMeanMIN=$(awk '{print $8}' $pair | datamash min 1)
       pMeanMEAN=$(awk '{print $8}' $pair | datamash mean 1)
        pMeanMAX=$(awk '{print $8}' $pair | datamash max 1)
       pMeanSDEV=$(awk '{print $8}' $pair | datamash sstdev 1)
       
   BonferroniMIN=$(awk '{print $9}' $pair | datamash min 1)
  BonferroniMEAN=$(awk '{print $9}' $pair | datamash mean 1)
   BonferroniMAX=$(awk '{print $8}' $pair | datamash max 1)
  BonferroniSDEV=$(awk '{print $8}' $pair | datamash sstdev 1)

  echo "$namesIDs $coevNumber $cCoevMIN $cCoevMAX $cCoevMEAN $cCoevSDEV $bootMIN $bootMEAN $bootMAX $bootSDEV $pMeanMIN $pMeanMAX $pMeanMEAN $pMeanSDEV $BonferroniMIN $BonferroniMAX $BonferroniMEAN $BonferroniSDEV" >> $TMP/$RESULTS/$EOUT
  echo "Processed $pair"
done

sed -i "1i Name1 msa1 Name2 msa2 coevNumber cCoevMIN cCoevMAX cCoevMEAN cCoevSDEV bootMIN bootMEAN bootMAX bootSDEV pMeanMIN pMeanMAX pMeanMEAN pMeanSDEV BonferroniMIN BonferroniMAX BonferroniMEAN BonferroniSDEV" $TMP/$RESULTS/$EOUT

# Prepare the XML header
cat << EOF >> $TMP/$RESULTS/$ENET
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<graph label="AutoCoEv" directed="0" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">
  <att name="networkData">
    <rdf:RDF>
      <rdf:Description rdf:about="http://www.cytoscape.org/">
        <dc:type>Protein-Protein Interaction</dc:type>
        <dc:description>N/A</dc:description>
        <dc:identifier>N/A</dc:identifier>
        <dc:date>$DATE</dc:date>
        <dc:title>CAPS</dc:title>
        <dc:source>http://www.cytoscape.org/</dc:source>
        <dc:format>Cytoscape-XGMML</dc:format>
      </rdf:Description>
    </rdf:RDF>
  </att>
  <att name="name" value="${MSAMETHOD}; ${GBLOCKS}; ${TREESCAPS}; ${PHYMLGUIDE}; ${TREESROOT}" type="string" cy:type="String"/>
  <att name="MSA method" value="${MSAMETHOD}" type="string" cy:type="String"/>
  <att name="Tree calculation" value="${TREESCAPS}" type="string" cy:type="String"/>
  <att name="Guide tree used?" value="${PHYMLGUIDE}" type="string" cy:type="String"/>
  <att name="Tree rooted?" value="${TREESROOT}" type="string" cy:type="String"/>
EOF

# Add nodes.
sed 1d $TMP/tsv/proteinsFound.tsv | while read -r idxml namexml attribute; do
echo -e "Adding node $idxml $namexml $attribute in $TMP/$RESULTS/$ENET"
cat << EOF >> $TMP/$RESULTS/$ENET
  <node id="$idxml" label="$namexml">
    <att name="shared name" value="$idxml" type="string" cy:type="String"/>
    <att name="name" value="$namexml" type="string" cy:type="String"/>
    <att name="attribute" value="$attribute" type="string" cy:type="String"/>
    <att name="URL" value="https://www.uniprot.org/uniprot/$idxml" type="string" cy:type="String"/>
  </node>
EOF
done

# Add edges
sed 1d $TMP/$RESULTS/$EOUT | \
while read -r Name1 msa1 Name2 msa2 coevNumber cCoevMIN cCoevMAX cCoevMEAN cCoevSDEV bootMIN bootMEAN bootMAX bootSDEV pMeanMIN pMeanMAX pMeanMEAN pMeanSDEV BonferroniMIN BonferroniMAX BonferroniMEAN BonferroniSDEV ; do
echo -e "Adding edge: $msa1-$msa2 in $TMP/$RESULTS/$ENET"
cat << EOF >> $TMP/$RESULTS/$ENET
  <edge id="$Name1-$Name2" label="$Name1 - $Name2" source="$msa1" target="$msa2" cy:directed="0">
    <att name="shared name" value="$msa1 - $msa2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="$Name1 - $Name2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="Number of coevolving sites" value="$coevNumber" type="real" cy:type="Double"/>
    <att name="Normalized correlation (max)" value="$cCoevMAX" type="real" cy:type="Double"/>
    <att name="Normalized correlation (mean)" value="$cCoevMEAN" type="real" cy:type="Double"/>
    <att name="Bootstrap value (max)" value="$bootMAX" type="real" cy:type="Double"/>
    <att name="Bootstrap value (mean)" value="$bootMEAN" type="real" cy:type="Double"/>
    <att name="P-value non-corrected (min)" value="$pMeanMIN" type="real" cy:type="Double"/>
    <att name="P-value non-corrected (mean)" value="$pMeanMEAN" type="real" cy:type="Double"/>
    <att name="P-value Bonferroni (min)" value="$BonferroniMIN" type="real" cy:type="Double"/>
    <att name="P-value Bonferroni (mean)" value="$BonferroniMEAN" type="real" cy:type="Double"/>
  </edge>
EOF
done
echo "</graph>" >> $TMP/$RESULTS/$ENET

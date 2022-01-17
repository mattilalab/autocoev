#!/bin/bash
ENET="$TMP/$RESULTS/FILTERED-P${PVALUE}-gaps${EGAPS}-pDiff${EPDIF}-Bonf${BONFERRONI}.xml"

write_xml() {
if [ -s "$TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv" ]; then

# Prepare the XML header
cat << EOF >> $TMP/$RESULTS/CAPS.P${PVALUE}-B${BONFERRONI}.xml
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
echo -e "Adding node $idxml $namexml $attribute in $TMP/$RESULTS/CAPS.P${PVALUE}-B${BONFERRONI}.xml"
cat << EOF >> $TMP/$RESULTS/CAPS.P${PVALUE}-B${BONFERRONI}.xml
  <node id="$idxml" label="$namexml">
    <att name="shared name" value="$idxml" type="string" cy:type="String"/>
    <att name="name" value="$namexml" type="string" cy:type="String"/>
    <att name="attribute" value="$attribute" type="string" cy:type="String"/>
    <att name="URL" value="https://www.uniprot.org/uniprot/$idxml" type="string" cy:type="String"/>
  </node>
EOF
done

# Add edges
sed 1d $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv | \
while read -r Name1 msa1 Name2 msa2 colA realA colB realB seq1 seq2 flt1 flt2 meanA meanB corr boot corrT averR averSigR cCoev pvalA pvalB pMean pDiff corr1 corr2 bonferroni GapA GapB GapAB DivA DivB DivAB totNumP totComp ; do
echo -e "Adding edge: $msa1-$msa2 in $TMP/$RESULTS/CAPS.P${PVALUE}-B${BONFERRONI}.xml"
cat << EOF >> $TMP/$RESULTS/CAPS.P${PVALUE}-B${BONFERRONI}.xml
  <edge id="$Name1-$Name2" label="$Name1 - $Name2" source="$msa1" target="$msa2" cy:directed="0">
    <att name="shared name" value="$msa1 - $msa2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="$Name1 - $Name2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="Residues" value="$seq1 - $seq2" type="string" cy:type="String"/>
    <att name="Gblocks" value="${flt1}${flt2}" type="string" cy:type="String"/>
    <att name="Gaps" value="$GapA - $GapB" type="string" cy:type="String"/>
    <att name="Variability" value="$DivA - $DivB" type="string" cy:type="String"/>
    <att name="P-value (Non-corrected)" value="$pMean" type="real" cy:type="Double"/>
    <att name="P-values difference" value="$pDiff" type="real" cy:type="Double"/>
    <att name="P-value (Bonferroni corrected)" value="$bonferroni" type="real" cy:type="Double"/>
    <att name="cCoev (normalized correlation)" value="$cCoev" type="real" cy:type="Double"/>
    <att name="Bootstrap" value="$boot" type="real" cy:type="Double"/>
    <att name="Mean gaps" value="$GapAB" type="real" cy:type="Double"/>
    <att name="Mean variability" value="$DivAB" type="real" cy:type="Double"/>
  </edge>
EOF
done
echo "</graph>" >> $TMP/$RESULTS/CAPS.P${PVALUE}-B${BONFERRONI}.xml
else
  echo -e "No codependence found"
fi
}

# View results in R Shiny
view_shiny(){
  if [ -s "$TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv" ]; then
    Rscript $CWD/R/howToRun.R \
      $TMP/$RESULTS/pairs-P${PVALUE}-B${BONFERRONI}.tsv \
      $CWD/R/autoCoEvShinyApp.R
  else
    echo -e "No codependence found"
  fi
}

write_filtered_xml(){

# Prepare the XML header
cat << EOF >> $ENET
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
  <att name="name" value="${MSAMETHOD}; ${TREESCAPS}; ${PHYMLGUIDE}; ${TREESROOT}" type="string" cy:type="String"/>
  <att name="MSA method" value="${MSAMETHOD}" type="string" cy:type="String"/>
  <att name="Tree calculation" value="${TREESCAPS}" type="string" cy:type="String"/>
  <att name="Guide tree used?" value="${PHYMLGUIDE}" type="string" cy:type="String"/>
  <att name="Tree rooted?" value="${TREESROOT}" type="string" cy:type="String"/>
EOF

# Add nodes.
sed 1d $TMP/tsv/proteinsFound.tsv | while read -r idxml namexml attribute; do
echo -e "Adding node $idxml $namexml $attribute in $ENET"
cat << EOF >> $ENET
  <node id="$idxml" label="$namexml">
    <att name="shared name" value="$idxml" type="string" cy:type="String"/>
    <att name="name" value="$namexml" type="string" cy:type="String"/>
    <att name="attribute" value="$attribute" type="string" cy:type="String"/>
    <att name="URL" value="https://www.uniprot.org/uniprot/$idxml" type="string" cy:type="String"/>
  </node>
EOF
done

# Add edges
sed 1d $TMPRES/$EOUT | \
while read -r Name1 msa1 Name2 msa2 coevThr avgCor avgSignCor totCompar coevNumAll coevMIN coevMAX coevMEAN coevSDEV coevNumber cCoevMIN cCoevMAX cCoevMEAN cCoevSDEV bootMIN bootMAX bootMEAN bootSDEV pMeanMIN pMeanMAX pMeanMEAN pMeanSDEV BonferroniMIN BonferroniMAX BonferroniMEAN BonferroniSDEV gblocksMIN gblocksMAX gblocksMEAN gblocksSDEV chiboth ; do
echo -e "Adding edge: $msa1-$msa2 in $ENET"
cat << EOF >> $ENET
  <edge id="$Name1-$Name2" label="$Name1 - $Name2" source="$msa1" target="$msa2" cy:directed="0">
    <att name="shared name" value="$msa1 - $msa2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="$Name1 - $Name2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="Number of coevolving sites" value="$coevNumber" type="real" cy:type="Double"/>
    <att name="Normalized correlation (min)" value="$cCoevMIN" type="real" cy:type="Double"/>
    <att name="Normalized correlation (max)" value="$cCoevMAX" type="real" cy:type="Double"/>
    <att name="Normalized correlation (mean)" value="$cCoevMEAN" type="real" cy:type="Double"/>
    <att name="Bootstrap value (min)" value="$bootMIN" type="real" cy:type="Double"/>
    <att name="Bootstrap value (max)" value="$bootMAX" type="real" cy:type="Double"/>
    <att name="Bootstrap value (mean)" value="$bootMEAN" type="real" cy:type="Double"/>
    <att name="P-value non-corrected (min)" value="$pMeanMIN" type="real" cy:type="Double"/>
    <att name="P-value non-corrected (max)" value="$pMeanMAX" type="real" cy:type="Double"/>
    <att name="P-value non-corrected (mean)" value="$pMeanMEAN" type="real" cy:type="Double"/>
    <att name="P-value Bonferroni (min)" value="$BonferroniMIN" type="real" cy:type="Double"/>
    <att name="P-value Bonferroni (max)" value="$BonferroniMAX" type="real" cy:type="Double"/>
    <att name="P-value Bonferroni (mean)" value="$BonferroniMEAN" type="real" cy:type="Double"/>
    <att name="Gblocks score (min)" value="$gblocksMIN" type="real" cy:type="Double"/>
    <att name="Gblocks score (max)" value="$gblocksMAX" type="real" cy:type="Double"/>
    <att name="Gblocks score (mean)" value="$gblocksMEAN" type="real" cy:type="Double"/>
  </edge>
EOF
done
echo "</graph>" >> $ENET
}

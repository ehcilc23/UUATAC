toolDir=/software/path/prefix/tools/

## 1 process cCRE: .bdg -> .bw 
# Tcell
grep -Ff mm10_chrom.txt Tcell/Mouse_Tcell_ATAC_treat_pileup.bdg > Tcell/tmp.bdg
grep -Fv "random" Tcell/tmp.bdg > Tcell/Mouse_Tcell_ATAC_treat_pileup_NormChrom.bdg
${toolDir}/bedGraphToBigWig Tcell/Mouse_Tcell_ATAC_treat_pileup_NormChrom.bdg ${toolDir}/mm10.normal.chrom.sizes Tcell/Tcell_cCRE.bw && rm Tcell/tmp.bdg && rm Tcell/Mouse_Tcell_ATAC_treat_pileup_NormChrom.bdg 

# NK cell
grep -Ff mm10_chrom.txt NK/Mouse_NK_ATAC_treat_pileup.bdg > NK/tmp.bdg
grep -Fv "random" NK/tmp.bdg > NK/Mouse_NK_ATAC_treat_pileup_NormChrom.bdg
${toolDir}/bedGraphToBigWig NK/Mouse_NK_ATAC_treat_pileup_NormChrom.bdg ${toolDir}/mm10.normal.chrom.sizes NK/NK_cCRE.bw && rm NK/tmp.bdg && rm NK/Mouse_NK_ATAC_treat_pileup_NormChrom.bdg

# Monocyte
grep -Ff mm10_chrom.txt Monocyte/Mouse_Monocyte_ATAC_treat_pileup.bdg > Monocyte/tmp.bdg
grep -Fv "random" Monocyte/tmp.bdg > Monocyte/Mouse_Monocyte_ATAC_treat_pileup_NormChrom.bdg
${toolDir}/bedGraphToBigWig Monocyte/Mouse_Monocyte_ATAC_treat_pileup_NormChrom.bdg ${toolDir}/mm10.normal.chrom.sizes Monocyte/Monocyte_cCRE.bw && rm Monocyte/tmp.bdg && rm Monocyte/Mouse_Monocyte_ATAC_treat_pileup_NormChrom.bdg


## 2 process conventional peak: .bdg -> .bw
# Tcell
sort -k1,1 -k2,2n Tcell/Tcell_conventional_peak_unbinary.bdg > Tcell/Peak.sort.bdg
${toolDir}/bedGraphToBigWig Tcell/Peak.sort.bdg ${toolDir}/mm10.chrom.sizes Tcell/Tcell_conventional_peak_unbinary.bw && rm Tcell/Peak.sort.bdg

# NK cell
sort -k1,1 -k2,2n NK/NK_conventional_peak_unbinary.bdg > NK/Peak.sort.bdg
${toolDir}/bedGraphToBigWig NK/Peak.sort.bdg ${toolDir}/mm10.chrom.sizes NK/NK_conventional_peak_unbinary.bw && rm NK/Peak.sort.bdg

# Monocyte
sort -k1,1 -k2,2n Monocyte/Monocyte_conventional_peak_unbinary.bdg > Monocyte/Peak.sort.bdg
${toolDir}/bedGraphToBigWig Monocyte/Peak.sort.bdg ${toolDir}/mm10.chrom.sizes Monocyte/Monocyte_conventional_peak_unbinary.bw && rm Monocyte/Peak.sort.bdg

## 3 download DNase-seq data from ENCODE project
## https://www.encodeproject.org/
# B cell: ENCFF392YEX, ENCFF937VYZ, ENCFF213OBH
# T cell: ENCFF191BSS, ENCFF215LYP, ENCFF788RUT
# NK cell: ENCFF049UXF, ENCFF140CJI, ENCFF142QWS
# Monocyte: ENCFF405TPZ, ENCFF502MGP, ENCFF293VUC


## 4 calculate Spearman correlation using deeptools
multiBigwigSummary bins -b Bcell/DNase/ENCFF392YEX_Bcell_rep3.bigWig Bcell/DNase/ENCFF937VYZ_Bcell_rep4.bigWig Bcell/Bcell_cCRE_macs3.bw Bcell/Bcell_conventional_peak_unbinary.bw \
    Tcell/DNase/ENCFF191BSS_Tcell_Male_rep1.bigWig Tcell/DNase/ENCFF215LYP_Tcell_Male_rep4.bigWig Tcell/Tcell_cCRE.bw Tcell/Tcell_conventional_peak_unbinary.bw \
    NK/DNase/ENCFF049UXF_NK_Male_rep2.bigWig NK/DNase/ENCFF140CJI_NK_Female_rep3.bigWig NK/NK_cCRE.bw NK/NK_conventional_peak_unbinary_Female.bw \
    Monocyte/DNase/ENCFF405TPZ_Monocyte_Female_rep1.bigWig Monocyte/DNase/ENCFF502MGP_Monocyte_Male_rep1.bigWig Monocyte/Monocyte_cCRE.bw Monocyte/Monocyte_conventional_peak_unbinary.bw \
    -o results.npz

plotCorrelation \
    -in results.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts (DNase,cCRE,Conventional Peak)" \
    --whatToPlot heatmap --colorMap coolwarm --plotNumbers \
    -o heatmap_SpearmanCorr_readCounts_DNase_cCRE_Peak.pdf   \
    --outFileCorMatrix SpearmanCorr_readCounts_DNase_cCRE_Peak.tab


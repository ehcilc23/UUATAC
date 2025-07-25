#BSUB -q normal
#BSUB -J cCRE
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=9]"
#BSUB -n 18



## 1 filter fragments by barcode for each tissue 
## (take mouse as an example, chicken, gecko, axolotl and zebrafish are processed by the same procedure)
bcDir=/file/path/prefix/UU/Mouse
outDir=/file/path/prefix/UU/Mouse/Mouse_Fragments_bcFiltered

# Mouse (Male)
for tissue in Bladder BoneMarrow BrainStem Cerebellum Heart Hypothalamus Kidney LargeIntestine Liver \
    Lung LymphNode Muscle Pancreas Pituitary Prostate SmallIntestine SpinalCord Spleen Stomach Telencephalon Testis Thymus; do

    zcat ${tissue}.sort.bed.gz | grep -Ff $bcDir/bc_Global_MouseMale.txt >> $outDir/MouseMale_Global.bed # 'bc_Global_MouseMale.txt' contains all barcodes used for Mouse (Male)
done

# Mouse (Female)
for tissue in AdrenalGland Bladder BoneMarrow BrainStem Cerebellum Heart Hypothalamus Kidney LargeIntestine Liver \
    Lung LymphNode Muscle Ovary Pancreas Pituitary SmallIntestine SpinalCord Spleen Stomach; do

    zcat ${tissue}.sort.bed.gz | grep -Ff $bcDir/bc_Global_MouseFemale.txt >> $outDir/MouseFemale_Global.bed # 'bc_Global_MouseFemale.txt' contains all barcodes used for Mouse (Female)
done

# merge into one
cat $outDir/MouseMale_Global.bed $outDir/MouseFemale_Global.bed >> $outDir/Mouse_Global.bed

# get each lineage compartment
for lineage in Endothelial Epithelial Erythroid Hepatocyte Immune Muscle Neuronal Reproductive Secretory Stromal; do
    grep -Ff $bcDir/bc_${lineage}_Mouse.txt $outDir/Mouse_Global.bed > $outDir/Mouse_${lineage}.bed
done


## 2 sort and compress
sort -k1,1V -k2,2n -k3,3n $outDir/Mouse_Global.bed > $outDir/Mouse_Global.sort.bed
bgzip $outDir/Mouse_Global.sort.bed


## 3 call cCREs using macs3
cd $outDir
macs3 callpeak -t Mouse_Global.sort.bed.gz \
-g 2.7e9 -f BEDPE -q 0.01 --keep-dup all \
-n MouseMale_Global_ATAC --outdir ./macs3_q0.01_wholegenome_keepdup/ -B --call-summits

macs3 callpeak -t Mouse_Global.sort.bed.gz \
-g 2.7e9 -f BEDPE -q 0.005 --keep-dup all \
-n MouseMale_Global_ATAC --outdir ./macs3_q0.005_wholegenome_keepdup/ -B --call-summits

## run for each lineage (the detailed code is omitted here for concision)

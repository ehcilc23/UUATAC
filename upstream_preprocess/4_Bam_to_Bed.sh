#BSUB -q normal
#BSUB -J bam2bed
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=9]"
#BSUB -n 18



sample=$(basename `pwd`)
# merge all bam files for one tissue into one bam
bam_to_merge=$(find -name *bam)
samtools merge -n $sample.shift.all.bam $bam_to_merge

# bam2bed and filter
samtools view -b -q 30 -F 4 -F 256 \
    $sample.shift.all.bam | \
    bedtools bamtobed -bedpe | \
    awk '{if ($1 == $4 && $9 != $10){print $0}}'> $sample.bed

# filter out fragments shorter than 10bp and longer than 1000bp
python3 /software/path/prefix/tools/UU_barcode/Fragments_filter1.py \
    $sample.bed $sample.clean.bed && rm $sample.bed

# deduplicate uniquely mapped fragments (UM) to unique fragments (UQ) and select top 50,000 barcodes
python3 /software/path/prefix/tools/UU_barcode/Fragments_filter2.py \
    $sample.clean.bed filter.bed 50000

# sort
cat filter.bed | sort -k1,1V -k2,2n -k3,3n > final.sort.bed && rm filter.bed
bgzip final.sort.bed

#BSUB -q normal
#BSUB -J Tagbc
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 18



outdir=bwa_out
tmpdir=tmp
shift=shift

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

if [ ! -d $tmpdir ]; then
    mkdir $tmpdir
fi

if [ ! -d $shift ]; then
    mkdir $shift
fi

chipID=E100000000
sample_name=$(basename `pwd`)
col=$(basename `pwd`)
dropseq_root=software/path/prefix/tools/Drop-seq_tools-2.5.1
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

# fastq --> bam
java -jar ${picard_jar} FastqToSam F1=$chipID\_L01_$col\_1.fq.gz F2=$chipID\_L01_$col\_2.fq.gz  O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir && rm $chipID\_L01_$col\_1.fq.gz && rm $chipID\_L01_$col\_2.fq.gz

# -----  Cell Barcode ------
#tag barcodes_BS
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=unaligned_tagged_Cellular.bam_summary_R1.txt \
   BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=True DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=2 \
 	INPUT=H.bam OUTPUT=$tmpdir/unaligned_tagged_Cell_R1.bam COMPRESSION_LEVEL=1
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=unaligned_tagged_Cellular.bam_summary_R2.txt \
   BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=2 \
 	INPUT=$tmpdir/unaligned_tagged_Cell_R1.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=1 && rm $tmpdir/unaligned_tagged_Cell_R1.bam

#add RT barcode
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=unaligned_tagged_RT_summary_R1.txt \
   BASE_RANGE=1-20 BASE_QUALITY=10 BARCODED_READ=2 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=RT NUM_BASES_BELOW_QUALITY=4 \
   INPUT=$tmpdir/unaligned_tagged_Cell.bam OUTPUT=$tmpdir/unaligned_tagged_RT_R1.bam COMPRESSION_LEVEL=1 && rm $tmpdir/unaligned_tagged_Cell.bam
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=unaligned_tagged_RT_summary_R2.txt \
   BASE_RANGE=1-20 BASE_QUALITY=10 BARCODED_READ=2 TAG_BARCODED_READ=True DISCARD_READ=false TAG_NAME=RT NUM_BASES_BELOW_QUALITY=4 \
   INPUT=$tmpdir/unaligned_tagged_RT_R1.bam OUTPUT=$tmpdir/unaligned_tagged_RT.bam COMPRESSION_LEVEL=1 && rm $tmpdir/unaligned_tagged_RT_R1.bam

# FilterBAM
${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=$tmpdir/unaligned_tagged_RT.bam OUTPUT=$tmpdir/unaligned_tagged_filtered.bam && rm $tmpdir/unaligned_tagged_RT.bam

# corrected bam for one mismatch
barcodepath=/software/path/prefix/tools/UU_barcode/

col=$(basename `pwd`)
python3 ${barcodepath}/UU_correct_barcode.py $barcodepath \
   $tmpdir/unaligned_tagged_filtered.bam \
   col$col\_ $tmpdir/filtered.bam && rm $tmpdir/unaligned_tagged_filtered.bam

echo "correct sam files done"

bamtools split -in $tmpdir/filtered.bam -tag TN
sh /software/path/prefix/tools/UU_barcode/Change_Tn5_names.sh

# mv
mkdir merge1
merge1=merge1
split=split

for i in Pituitary AdrenalGland LymphNode BoneMarrow SmallIntestine Telencephalon Cerebellum Hypothalamus Ovary BrainStem SpinalCord; do
   mkdir $merge1/$i
done

for i in {97..104}; do
   mv $split/$i/$i.bam $merge1/Pituitary
done
for i in {105..112}; do
   mv $split/$i/$i.bam $merge1/AdrenalGland
done
for i in {113..120} {189..192}; do
   mv $split/$i/$i.bam $merge1/LymphNode
done
for i in {121..128} {185..188}; do
   mv $split/$i/$i.bam $merge1/BoneMarrow
done
for i in {129..136}; do
   mv $split/$i/$i.bam $merge1/SmallIntestine
done
for i in {137..144}; do
   mv $split/$i/$i.bam $merge1/Telencephalon
done
for i in {145..152}; do
   mv $split/$i/$i.bam $merge1/Cerebellum
done
for i in {153..160}; do
   mv $split/$i/$i.bam $merge1/Hypothalamus
done
for i in {161..168}; do
   mv $split/$i/$i.bam $merge1/Ovary
done
for i in {169..176}; do
   mv $split/$i/$i.bam $merge1/BrainStem
done
for i in {177..184}; do
   mv $split/$i/$i.bam $merge1/SpinalCord
done

# # merge & align
# cp ../jobs2.sh .
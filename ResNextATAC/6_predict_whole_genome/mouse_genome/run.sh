#!/bin/bash 
#SBATCH --job-name=UU_pred_mouse
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=%J.out
#SBATCH --error=%J.err

python pred_Mouse_chr.py --region Region.*.slide.Mouse.input.txt --anno Female_mus_anno_5w.txt

species_chrom_size=mm10.chrom.sizes

#ct=`tail -n+2 ./Female_mus_anno_5w.txt |cut -f 5|sort|uniq`
#for i in $ct;do
#cat Predict_track/test_predict_"$i"_*.txt >Predict_track/merged_"$i".txt;
#done

# get bedgraph
for x in `ls ./Predict_track/merged_*`; do 
cat $x | perl -ne '@a=($_ =~ /(.*):(\d+)-(\d+)\t(.*)/); print "".join("\t", $a[0],  $a[1]+250-50, $a[1]+250+50, $a[3]), "\n";' > `echo ${x}|sed 's/\.txt//g'`.bedGraph;
done

# get bw
cat ${species_chrom_size} |awk '{print $1"\t"0"\t"$2}' >chrom.size

for x in `ls ./Predict_track/merged_*.bedGraph`; do 
subanno_bedgraph=`echo ${x} | sed 's#merged_##g' `
subanno_id=`echo ${subanno_bedgraph} | sed 's#.bedGraph##g' `
echo $subanno_id
bedtools intersect -a $x -b chrom.size >tmp
sort -k1,1V -k2,2n -k3,3n tmp >$subanno_bedgraph 
bedGraphToBigWig $subanno_bedgraph ${species_chrom_size} ${subanno_id}.bw
rm tmp;
done

mkdir bw
mv ./Predict_track/*.bw ./bw/
rm -rf Predict_track/

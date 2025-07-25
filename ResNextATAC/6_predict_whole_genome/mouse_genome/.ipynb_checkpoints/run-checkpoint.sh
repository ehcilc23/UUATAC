#!/bin/bash 
#SBATCH --job-name=UU_pred_mouse
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=%J.out
#SBATCH --error=%J.err


python pred_Mouse_chr.py --region Region.*.slide.Mouse.input.txt --anno Female_mus_anno_5w.txt
ct=`tail -n+2 ./Female_mus_anno_5w.txt |cut -f 5|sort|uniq`
for i in $ct;do
cat Predict_track/test_predict_"$i"_*.txt >Predict_track/merged_"$i".txt;
done

#get bedgraph
for x in `ls ./Predict_track/merged_*`; do 
cat $x | perl -ne '@a=($_ =~ /(.*):(\d+)-(\d+)\t(.*)/); print "".join("\t", $a[0],  $a[1]+250-50, $a[1]+250+50, $a[3]), "\n";' > `echo ${x}|sed 's/\.txt//g'`.bedGraph;
done

#get bw
for x in `ls ./Predict_track/merged_*.bedGraph`; do 
tmp=`echo ${x} | sed 's#merged_##g' `
sort -k1,1V -k2,2n -k3,3n $x  >tmp
bedClip tmp chrom.size $tmp
tmp1=`echo ${tmp} | sed 's#.bedGraph##g' `
bedGraphToBigWig $tmp chrom.size ${tmp1}.bw; 
done

mkdir bw
mv ./Predict_track/*.bw ./bw/
rm -rf Predict_track/

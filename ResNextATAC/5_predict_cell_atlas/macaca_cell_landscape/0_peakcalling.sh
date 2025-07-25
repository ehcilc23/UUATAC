#BSUB -q normal
#BSUB -J macs3
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=9]"
#BSUB -n 36
cd /file/path/prefix/macaca_bedgraph_sub/
for i in `ls *.bedGraph`;do 
lineage=`echo ${i} | sed 's#.bedGraph##g' `;
macs3 bdgpeakcall -i $i --cutoff 0.05 --o-prefix $lineage --outdir OUTPUT
 done

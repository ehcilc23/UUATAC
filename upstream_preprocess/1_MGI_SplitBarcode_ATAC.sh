#BSUB -q normal
#BSUB -J splitBarcode
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 36



/software/path/prefix/tools/splitBarcode-master/V0.1.6_release/linux/splitBarcode \
/software/path/prefix/tools/splitBarcode-master/index_HUADA.txt E100000000_L01_read_1.fq.gz -2 E100000000_L01_read_2.fq.gz \
-o OUTPUT -b 211 10 2 -r

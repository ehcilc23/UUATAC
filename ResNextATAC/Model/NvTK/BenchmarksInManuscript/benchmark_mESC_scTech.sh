BASEDIR=~/JiaqiLi/NvTK/BenchmarksInManuscript

# 'BulkSeq', 'CELseq', 'DropSeq', 'MARSseq', 'MicrowellSeq','SCRBseq', 'SmartSeq'
# 'MAGIC-CELseq', 'MAGIC-DropSeq', 'MAGIC-MARSseq', 'MAGIC-MicrowellSeq', 'MAGIC-SCRBseq', 'MAGIC-SmartSeq'
arr_Tech=("BulkSeq" "CELseq" "DropSeq" "MARSseq" "MicrowellSeq" "SCRBseq" "SmartSeq" "SmartSeq2" \
"MAGIC-CELseq" "MAGIC-DropSeq" "MAGIC-MARSseq" "MAGIC-MicrowellSeq" "MAGIC-SCRBseq" "MAGIC-SmartSeq" "MAGIC-SmartSeq2")

for Tech in ${arr_Tech[@]}
do
    mkdir ${BASEDIR}/Benchmark_scTech_mESC/mESC_$Tech
    cd ${BASEDIR}/Benchmark_scTech_mESC/mESC_$Tech && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
        --gpu-device 1 \
        --subset_task $Tech --subset_task_by Cluster \
        --lr 1e-5 --patience 10 --batch_size 128
    # cd Motif/ && tomtom -oc tomtom_JASPAR meme_conv1_thres9.txt ${BASEDIR}/Dataset/JASPAR2022_CORE_non-redundant_pfms_meme.txt
done

# mess up all single cell technologies to simulate batch-effect
mkdir ${BASEDIR}/Benchmark_scTech_mESC/mESC_BatchEffect
cd ${BASEDIR}/Benchmark_scTech_mESC/mESC_BatchEffect && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
    --gpu-device 1 \
    --subset_task Mouse --subset_task_by Species \
    --lr 1e-5 --patience 10 --batch_size 128
# cd Motif/ && tomtom -oc tomtom_JASPAR meme_conv1_thres9.txt ${BASEDIR}/Dataset/JASPAR2022_CORE_non-redundant_pfms_meme.txt

# compare the experimental batch-effect
arr_Tech=("CELseq2A" "CELseq2B" "DropSeqA" "DropSeqB" "MARSseqA" "MARSseqB" \
            "SCRBseqA" "SCRBseqB" "SmartSeq2A" "SmartSeq2B" "SmartSeqA" "SmartSeqB")

for Tech in ${arr_Tech[@]}
do
    mkdir ${BASEDIR}/Benchmark_scTech_mESC/mESC_$Tech
    cd ${BASEDIR}/Benchmark_scTech_mESC/mESC_$Tech && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
        --gpu-device 1 \
        --subset_task $Tech --subset_task_by Celltype \
        --lr 1e-5 --patience 10 --batch_size 128
    # cd Motif/ && tomtom -oc tomtom_JASPAR meme_conv1_thres9.txt ${BASEDIR}/Dataset/JASPAR2022_CORE_non-redundant_pfms_meme.txt
done

# for i in `ls .. |grep mESC_`;do
# echo ../$i/runs/Mar*
# cp -r ../$i/runs/Mar* $i
# done
# tensorboard --logdir=runs --bind_all

#!/bin/bash 
#SBATCH --job-name=UU_chicken
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=~/WHY/AI/CodeTest/20230901_UU_chicken_scATAC_Lineage/5wCells_16wPeak_8wNega
cd ${BASEDIR}/
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Chicken_5wCells_16wPeaks_8wnegative.shuffled.noimpute.500bp.20230901.h5 \
    --use_ResNet 34 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 0


# get cell embedding 
#python GetCellEmbedding.py \
#    --anno ${BASEDIR}/cell_info_1wCells_hq.csv


# # plot PT heatmap
# python PT_heatmap.py \
#     --data ${BASEDIR}/Zebrafish_1wCells_hq_30wPeaks.noimpute.500bp.20230714.h5 \
#     --anno ${BASEDIR}/cell_info_1wCells_hq.csv

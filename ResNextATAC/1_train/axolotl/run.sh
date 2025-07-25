#!/bin/bash 
#SBATCH --job-name=UU_ax
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=~/WHY/AI/CodeTest/20230902_UU_axolotl_scATAC_Lineage/5wCells_8wPeaks_4wNega
cd ${BASEDIR}/
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Axolotl_5wCells_8wPeaks_4wnegative.shuffled.noimpute.500bp.20230904.h5 \
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

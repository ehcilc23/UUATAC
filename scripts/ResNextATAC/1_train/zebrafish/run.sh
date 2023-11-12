#!/bin/bash 
#SBATCH --job-name=UU_zf
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=/file/path/prefix/20230915_UU_zebrafish_scATAC_Lineage/5wCells_7wPeaks_3wNega
cd ${BASEDIR}/
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Zebrafish_5wCells_7wPeaks_3wnegative.shuffled.noimpute.500bp.20230916.h5 \
    --use_ResNet 34 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 0
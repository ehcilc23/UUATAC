#!/bin/bash 
#SBATCH --job-name=UU_gecko
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=/file/path/prefix/20230901_UU_gecko_scATAC_Lineage/5wCells_13wPeak_6wNega
cd ${BASEDIR}/
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Gecko_5wCells_13wPeaks_6wnegative.shuffled.noimpute.500bp.20230901.h5 \
    --use_ResNet 34 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 0
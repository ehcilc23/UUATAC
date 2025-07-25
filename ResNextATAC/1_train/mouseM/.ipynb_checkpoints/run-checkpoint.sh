#!/bin/bash 
#SBATCH --job-name=UU_malemus
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=/file/path/prefix/20230831_UU_MusM_scATAC_Lineage/5wCells_15wPeak_7wNega
cd ${BASEDIR}/
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/malemus_5wCells_15wPeaks_7wnegative.shuffled.noimpute.500bp.20230901.final.h5 \
    --use_ResNet 34 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 0
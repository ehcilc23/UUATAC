#!/bin/bash 
#SBATCH --job-name=UU_femalemus
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=/file/path/prefix/20230511_UU_MusF_scATAC_Lineage/5wCells_17wPeak_8wNega

cd ${BASEDIR}/
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Femalemus_5wCells_17wPeaks_8wnegative.shuffled.noimpute.500bp.20230822.h5 \
    --use_ResNet 34 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 0



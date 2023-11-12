#!/bin/bash 
#SBATCH --job-name=benchmark
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err


BASEDIR=/file/path/prefix/
cd $BASEDIR/e=resnext
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Femalemus_5kCells_17wPeaks_8wnegative.shuffled.noimpute.500bp.20231021.h5 \
    --use_ResNet 34 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 10 \
    --gpu-device 0


#!/bin/bash 

BASEDIR=~/JiaqiLi/NvTK/BenchmarksInManuscript

# DeepCNN
mkdir ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_DeepCNN
cd ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_DeepCNN && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --use_DeepCNN 4 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 2000 \
    --gpu-device 3 


# CBAM
mkdir ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_CBAM
cd ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_CBAM && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --use_CBAM True \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 2000 \
    --gpu-device 3 


# Transformer
mkdir ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_Transformer
cd ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_Transformer && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --use_Transformer \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 3 


# ResNet
mkdir ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_ResNet18
cd ${BASEDIR}/Benchmark_Arch_scATAC/scATAC_ResNet18 && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --use_ResNet 18 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 100 \
    --gpu-device 3 

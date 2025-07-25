#!/bin/bash 

BASEDIR=~/JiaqiLi/NvTK/BenchmarksInManuscript

# DeepCNN
mkdir ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_DeepCNN
cd ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_DeepCNN && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
    --use_DeepCNN 4 \
    --subset_task SmartSeq2 --subset_task_by Cluster \
    --lr 1e-5 --patience 100 --batch_size 16 \
    --gpu-device 1


# CBAM
mkdir ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_CBAM
cd ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_CBAM && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
    --use_CBAM True \
    --subset_task SmartSeq2 --subset_task_by Cluster \
    --lr 1e-5 --patience 100 --batch_size 16 \
    --gpu-device 1


# Transformer
mkdir ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_Transformer
cd ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_Transformer && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
    --use_Transformer \
    --subset_task SmartSeq2 --subset_task_by Cluster \
    --lr 1e-5 --patience 10 --batch_size 4 \
    --gpu-device 1


# ResNet
mkdir ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_ResNet18
cd ${BASEDIR}/Benchmark_Arch_scRNA/scRNA_ResNet18 && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
    --use_ResNet 18 \
    --subset_task SmartSeq2 --subset_task_by Cluster \
    --lr 1e-5 --patience 10 --batch_size 4 \
    --gpu-device 1

BASEDIR=~/JiaqiLi/NvTK/BenchmarksInManuscript

# baseline-model
N=128; L=15; P=15; PT=avgpool; AT=ReLU
mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}
cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
    --tasktype classification \
    --lr 1e-5 --patience 10 --batch_size 1000 \
    --gpu-device 2


# Length
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for L in 5 25
do
    mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
        --tasktype classification \
        --filterLenConv1 $L \
        --lr 1e-5 --patience 10 --batch_size 1000 \
        --gpu-device 2
done


# Filter_Number
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for N in 8 32 512
do
    mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
        --tasktype classification \
        --numFiltersConv1 $N \
        --lr 1e-5 --patience 10 --batch_size 1000 \
        --gpu-device 2
done


# Pooling size
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for P in 5 25
do
    mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
        --tasktype classification \
        --Pool1 $P \
        --lr 1e-5 --patience 10 --batch_size 1000 \
        --gpu-device 2
done


# Pooling Type
N=128; L=15; P=15; PT=maxpool; AT=ReLU
mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}
cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
    --tasktype classification \
    --pooltype $PT \
    --lr 1e-5 --patience 10 --batch_size 1000 \
    --gpu-device 2


# activation type
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for AT in Sigmoid Tanh Exp LeakyReLU
do
    mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
        --tasktype classification \
        --activation $AT \
        --lr 1e-5 --patience 10 --batch_size 1000 \
        --gpu-device 2
done


# BatchNorm
N=128; L=15; P=15; PT=avgpool; AT=ReLU
mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}BN
cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}BN && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
    --tasktype classification \
    --use_BN True \
    --lr 1e-5 --patience 10 --batch_size 1000 \
    --gpu-device 2


# activation type
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for AT in Sigmoid Tanh Exp
do
    mkdir ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}BN
    cd ${BASEDIR}/Benchmark_cnn_synthetic/Synthetic_N${N}L${L}P${PT:0:1}P${P}${AT}BN && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/synthetic_dataset_simple.h5 \
        --tasktype classification \
        --activation $AT \
        --use_BN True \
        --lr 1e-5 --patience 10 --batch_size 1000 \
        --gpu-device 2
done

# for dir in `ls -lt $BASEDIR/Benchmark_cnn_synthetic | awk '/^d/ {print $NF}'`
# do
#     cd $BASEDIR/Benchmark_cnn_synthetic/${dir}/Motif/ && echo `pwd`
#     tomtom -oc tomtom_JASPAR meme_conv1_thres9.txt ${BASEDIR}/Dataset_compare/JASPAR2022_CORE_non-redundant_pfms_meme.txt
# done

# for i in `ls .. |grep Synthetic_`;do
# echo ../$i/runs/Mar*/
# cp -r ../$i/runs/Mar*/ $i
# done
# tensorboard --logdir=runs --bind_all

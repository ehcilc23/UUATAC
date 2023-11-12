BASEDIR=~/JiaqiLi/NvTK/BenchmarksInManuscript

# baseline-model
N=128; L=15; P=15; PT=avgpool; AT=ReLU
mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}
cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --tasktype binary_classification \
    --lr 1e-3 --patience 10 \
    --batch_size 5000 \
    --gpu-device 3 


# Length
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for L in 5 25
do
    mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
        --tasktype binary_classification \
        --filterLenConv1 $L \
        --lr 1e-3 --patience 10 \
        --batch_size 5000 \
        --gpu-device 2
done


# Filter_Number
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for N in 8 32 512
do
    mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
        --tasktype binary_classification \
        --numFiltersConv1 $N \
        --lr 1e-3 --patience 10 \
        --batch_size 5000 \
        --gpu-device 2
done


# Pooling size
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for P in 5 25
do
    mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
        --tasktype binary_classification \
        --Pool1 $P \
        --lr 1e-3 --patience 10 \
        --batch_size 5000 \
        --gpu-device 2
done


# Pooling Type
N=128; L=15; P=15; PT=maxpool; AT=ReLU
mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}
cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --tasktype binary_classification \
    --pooltype $PT \
    --lr 1e-3 --patience 10 \
    --batch_size 5000 \
    --gpu-device 2


# activation type
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for AT in Sigmoid Tanh Exp LeakyReLU
do
    mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}
    cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT} && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
        --tasktype binary_classification \
        --activation $AT \
        --lr 1e-3 --patience 10 \
        --batch_size 5000 \
        --gpu-device 2
done


# BatchNorm
N=128; L=15; P=15; PT=avgpool; AT=ReLU
mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}BN
cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}BN && echo `pwd`
python ${BASEDIR}/nvtk-benchmark.py \
    ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
    --tasktype binary_classification \
    --use_BN True \
    --lr 1e-3 --patience 10 \
    --batch_size 5000 \
    --gpu-device 2


# activation type
N=128; L=15; P=15; PT=avgpool; AT=ReLU
for AT in Sigmoid Tanh Exp
do
    mkdir ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}BN
    cd ${BASEDIR}/Benchmark_cnn_scATAC/scATAC_N${N}L${L}P${PT:0:1}P${P}${AT}BN && echo `pwd`
    python ${BASEDIR}/nvtk-benchmark.py \
        ${BASEDIR}/Dataset/Dataset.sciATAC1_train_test.h5 \
        --tasktype binary_classification \
        --activation $AT \
        --use_BN True \
        --lr 1e-3 --patience 10 \
        --batch_size 5000 \
        --gpu-device 2
done

# for dir in `ls -lt $BASEDIR/Benchmark_cnn_scATAC | awk '/^d/ {print $NF}'`
# do
#     cd $BASEDIR/Benchmark_cnn_scATAC/${dir}/Motif/ && echo `pwd`
#     tomtom -oc tomtom_JASPAR meme_conv1_thres9.txt ${BASEDIR}/Dataset_compare/JASPAR2022_CORE_non-redundant_pfms_meme.txt
# done

# for i in `ls .. |grep scATAC_`;do
# echo ../$i/runs/Mar*/
# cp -r ../$i/runs/Mar*/ $i
# done
# tensorboard --logdir=runs --bind_all

# for i in `ls . |grep scATAC_`;do
# echo $i
# mkdir results/$i results/$i/Test results/$i/Motif
# mv $i/Test/test* results/$i/Test/.
# mv $i/Motif/seqlets.fasta results/$i/Motif/.
# mv $i/Motif/W.npy results/$i/Motif/.
# done
# tar zvcf Benchmark_Arch_scATAC.tar.gz scATAC_*


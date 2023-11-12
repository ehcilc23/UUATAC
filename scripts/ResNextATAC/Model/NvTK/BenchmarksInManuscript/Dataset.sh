BASEDIR=/media/ggj/Files/mount/Nvwa_Benchmark

# generate mESC data set
cd $BASEDIR/Dataset_compare
python ./2_propare_datasets.py leave_chrom \
    /media/ggj/Files/NvWA/PreprocGenome/mouse_updown10k_official.onehot.p \
    ./scTech_datasets/adata_MAGIC_All.p ./scTech_datasets/Mouse_sctechnologies_MAGIC_All.annotation.tsv \
    ./Dataset.Mouse_sctechnologies_leave_chrom8_mESC_MAGIC.h5 \
    /media/ggj/Files/NvWA/PreprocGenome/Mouse.gtf_annotation.csv 8
mv log.leave_chrom.txt scTech_datasets/log.leave_chrom8_mESC.20220314.txt


# generate Synthetic data set
# run 0_Generate_synthetic_dataset.ipynb

# Human GRCh38
wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens//Homo_sapiens.GRCh38.102.gtf.gz

python 0_proc_prom_region_seq_official.py database_genome/Human_Homo_sapiens/Homo_sapiens.GRCh38.fa database_genome/Human_Homo_sapiens/Homo_sapiens.GRCh38.gtf 10000 10000 >Human_updown10k_official.fa 2>log.Human_updown10k_official.fa
python ./0_onehot_geome.py Human_updown10k_official.fa Human_updown10k_official.onehot.p

Onehot=../onehot/Human_updown10k.rm_official.onehot.p
GTF=../onehot/Human.gtf_annotation.csv
Label=./HCL_MAGIC_merge_breast_testis_500gene_009.p
Annotation=../HCL_microwell_twotissue_preEmbryo.cellatlas.annotation.txt

python ../../2_propare_datasets.py leave_chrom $Onehot $Label $Annotation Dataset.Human_Chrom8_train_test.h5 $GTF 8
python ../../2_propare_datasets.py leave_chrom_CV $Onehot $Label $Annotation ./ $GTF 8

python ./2_propare_datasets.py leave_chrom \
    /media/ggj/Files/NvWA/PreprocGenome/Human_updown10k_official.onehot.p \
    ./HPSC.regression.p ./HPSC.regression.anno.tsv \
    ./Dataset.HPSC_regression_leaveChrom8.h5 \
    /media/ggj/Files/NvWA/PreprocGenome/Human.gtf_annotation.csv 8

# Mouse GRCm38
wget ftp://ftp.ensembl.org/pub/release-88/fasta/mus_musculus//dna/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-88/gtf/mus_musculus//Mus_musculus.GRCm38.102.gtf.gz

python 0_proc_prom_region_seq_official.py database_genome/Human_Homo_sapiens/Homo_sapiens.GRCh38.fa database_genome/Human_Homo_sapiens/Homo_sapiens.GRCh38.gtf 10000 10000 >Human_updown10k_official.fa 2>log.Human_updown10k_official.fa
python ./0_onehot_geome.py Human_updown10k_official.fa Human_updown10k_official.onehot.p

Onehot=../onehot/Human_updown10k.rm_official.onehot.p
GTF=../onehot/Human.gtf_annotation.csv
Label=./HCL_MAGIC_merge_breast_testis_500gene_009.p
Annotation=../HCL_microwell_twotissue_preEmbryo.cellatlas.annotation.txt

python ../../2_propare_datasets.py leave_chrom $Onehot $Label $Annotation Dataset.Human_Chrom8_train_test.h5 $GTF 8
python ../../2_propare_datasets.py leave_chrom_CV $Onehot $Label $Annotation ./ $GTF 8

# DeepSEA 919 epigenetic features
cd /share/home/guoguoji/JiaqiLi/nvwa-official/pretrain_deepsea/data
# /share/home/guoguoji/JiaqiLi/nvwa/predict_epigenetic/data
python prepare_dataset_DeepSEA.py

# sci-ATAC
cd /share/home/guoguoji/NvWA/sciATAC
python prepare_dataset_sciATAC.py

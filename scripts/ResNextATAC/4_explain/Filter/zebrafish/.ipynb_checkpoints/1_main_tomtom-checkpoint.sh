#JASPAR
MOTIFDB=/media/ggj/Guo-4T-AI/FYT/UU/TFBS/download_meme_motifDB/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
tomtom -oc tomtom_conv1_JASPAR -thresh 0.1 -dist allr -no-ssc ./Motif/meme_conv1.txt $MOTIFDB &>log.tomtom_conv1_JASPAR &


#CISBP
MOTIFDB=/media/ggj/Guo-4T-AI/FYT/UU/TFBS/download_meme_motifDB/motif_databases/CIS-BP_2.00/Danio_rerio.meme
tomtom -oc tomtom_conv1_CISBP -thresh 0.1 -dist allr -no-ssc ./Motif/meme_conv1.txt $MOTIFDB &>log.tomtom_conv1_CISBP &

#CISTARGET
MOTIFDB=/media/ggj/Guo-4T-AI/FYT/UU/TFBS/Cistarget/out.meme
tomtom -oc tomtom_conv1_Cistarget -thresh 0.1 -dist allr -no-ssc ./Motif/meme_conv1.txt $MOTIFDB &>log.tomtom_conv1_Cistarget &

#MOTIFDB=/media/ggj/Guo-4T-AI/FYT/UU/TFBS/Cistarget/out.meme
#tomtom -oc test -thresh 0.5 /media/ggj/Guo-4T-AI/FYT/UU/TFBS/meme_conv1_thres9.txt $MOTIFDB &>log.tomtom_nvwa_conv1_Cistarget &
#tomtom -oc test1 -thresh 0.1 -dist allr -no-ssc /media/ggj/Guo-4T-AI/FYT/UU/TFBS/meme_conv1_thres9.txt $MOTIFDB &>log.tomtom_conv1 &
#CISTARGET
MOTIFDB=/file/path/prefix/Cistarget/out.meme
tomtom -oc tomtom_conv1_Cistarget -thresh 0.1 -dist allr -no-ssc ./Motif/meme_conv1.txt $MOTIFDB &>log.tomtom_conv1_Cistarget &

##mouseF with other species
#mouseF vs mouseF
MOTIFDB=/file/path/prefix/Cross_species/mouseF_conv1.meme
tomtom -oc tomtom_conv1_mouseF_mouseF -thresh 0.1 -dist allr -no-ssc  mouseF_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

#mouseF vs chicken
MOTIFDB=/file/path/prefix/Cross_species/mouseF_conv1.meme
tomtom -oc tomtom_conv1_mouseF_chicken -thresh 0.1 -dist allr -no-ssc  chicken_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

MOTIFDB=/file/path/prefix/Cross_species/chicken_conv1.meme
tomtom -oc tomtom_conv1_chicken_mouseF -thresh 0.1 -dist allr -no-ssc  mouseF_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

#mouseF vs gecko
MOTIFDB=/file/path/prefix/Cross_species/mouseF_conv1.meme
tomtom -oc tomtom_conv1_mouseF_gecko -thresh 0.1 -dist allr -no-ssc  gecko_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

MOTIFDB=/file/path/prefix/Cross_species/gecko_conv1.meme
tomtom -oc tomtom_conv1_gecko_mouseF -thresh 0.1 -dist allr -no-ssc  mouseF_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

#mouseF vs axolotl
MOTIFDB=/file/path/prefix/Cross_species/mouseF_conv1.meme
tomtom -oc tomtom_conv1_mouseF_axolotl -thresh 0.1 -dist allr -no-ssc  axolotl_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

MOTIFDB=/file/path/prefix/Cross_species/axolotl_conv1.meme
tomtom -oc tomtom_conv1_axolotl_mouseF -thresh 0.1 -dist allr -no-ssc  mouseF_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

#mouseF vs zebrafish
MOTIFDB=/file/path/prefix/Cross_species/mouseF_conv1.meme
tomtom -oc tomtom_conv1_mouseF_zebrafish -thresh 0.1 -dist allr -no-ssc  zebrafish_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 

MOTIFDB=/file/path/prefix/Cross_species/zebrafish_conv1.meme
tomtom -oc tomtom_conv1_zebrafish_mouseF -thresh 0.1 -dist allr -no-ssc  mouseF_conv1.meme $MOTIFDB &>log.tomtom_conv1_Cistarget 



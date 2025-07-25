###### set working dir ######
parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)

###### select mouse filters ######
df=read.delim('../female_mus/tomtom_conv1_Cistarget/tomtom.tsv')
df=df[df$E.value<=0.05,]
df=df[!is.na(df$Query_ID),]
mouseF_motif=df$Query_ID%>%unique()
mouseF_motif=paste0('mouseF_',mouseF_motif)


###### compare mouse filters with other vetebrates'######
#mouseF with chicken
library(data.table)
df=read.delim('./tomtom_conv1_mouseF_chicken/tomtom.tsv')

df=df[,c("Query_ID","Target_ID", "E.value")]

df=df[df$E.value<=0.1,]

df=df[!is.na(df$Query_ID),]
mouseF_chicken_motif=paste0('mouseF_',df$Target_ID%>%unique())

#mouseF with gecko
df=read.delim('./tomtom_conv1_mouseF_gecko/tomtom.tsv')

df=df[,c("Query_ID","Target_ID", "E.value")]

df=df[df$E.value<=0.1,]

df=df[!is.na(df$Query_ID),]
mouseF_gecko_motif=paste0('mouseF_',df$Target_ID%>%unique())


#mouseF with axolotl
df=read.delim('./tomtom_conv1_mouseF_axolotl/tomtom.tsv')

df=df[,c("Query_ID","Target_ID", "E.value")]

df=df[df$E.value<=0.1,]

df=df[!is.na(df$Query_ID),]
mouseF_axolotl_motif=paste0('mouseF_',df$Target_ID%>%unique())

#mouseF with zebrafish
df=read.delim('./tomtom_conv1_mouseF_zebrafish/tomtom.tsv')

df=df[,c("Query_ID","Target_ID", "E.value")]

df=df[df$E.value<=0.1,]

df=df[!is.na(df$Query_ID),]
mouseF_zebrafish_motif=paste0('mouseF_',df$Target_ID%>%unique())

mouseF_chicken_motif1=mouseF_chicken_motif[mouseF_chicken_motif%in%mouseF_motif]
mouseF_gecko_motif1=mouseF_gecko_motif[mouseF_gecko_motif%in%mouseF_motif]
mouseF_axolotl_motif1=mouseF_axolotl_motif[mouseF_axolotl_motif%in%mouseF_motif]
mouseF_zebrafish_motif1=mouseF_zebrafish_motif[mouseF_zebrafish_motif%in%mouseF_motif]

tmp=c(mouseF_chicken_motif1,mouseF_gecko_motif1,mouseF_axolotl_motif1,mouseF_zebrafish_motif1)
table(tmp)%>%as.data.frame()->tmp_df

tmp_df$Freq=tmp_df$Freq+1

motif_mouseF_only=setdiff(mouseF_motif,tmp_df$tmp)
tmp_df1=data.frame(tmp=motif_mouseF_only,Freq=1)

all_df=rbind(tmp_df,tmp_df1)
colnames(all_df)=c('Motif','Level')
all_df$Level=paste0('Level_',all_df$Level)

###### pie plot ######
type_df=table(all_df$Level)%>%as.data.frame()
type_df$ratio=paste0(round(type_df$Freq/sum(type_df$Freq)*100,2),'%')
cols=rev(c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854"))
library(scales)
show_col(cols)

label=as.character(type_df$Var1)
pdf('./Cross_species/mouse_motif_level.pdf',height = 10,width = 10)
pie(type_df$Freq,labels=type_df$ratio, main="Mouse motif",col = cols)
legend("topright",label, cex=0.8,fill= cols)
dev.off()

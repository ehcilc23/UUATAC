###### set working dir ######
parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)

###### merge peaks ######
file.ls=list.files()

all.use=c()
for (i in file.ls) {
  df=fread(i)
  df$summit=df$V2+df$V10
  df$chr=df[,1]
  use=paste0(df$chr,'_',df$summit)
  all.use=c(all.use,use)
}

all.use[1:5]
all.use=unique(all.use)
chr=reshape2::colsplit(all.use,pattern = '_',names = c('c1','c2'))$c1
summit=reshape2::colsplit(all.use,pattern = '_',names = c('c1','c2'))$c2

df=data.frame('chr'=chr,'summit'=summit)
head(df)

df$start=df$summit-250
df$end=df$summit+250
#df$chr=df$`track type=narrowPeak name="human" description="human" nextItemButton=on`
df=df[,c('chr','start','end')]
table(df$start<0)


###### save ######
head(df,20)
df$start=as.integer(df$start)
df$end=as.integer(df$end)
table(df$chr)
fwrite(df,file = './pred_peak_monkey_sub.bed',quote = F,row.names = F,col.names = F,sep = '\t')

#run in terminal
#bedClip pred_peak_monkey_sub.bed chrom.size pred_peak_monkey_sub_filter.bed
#bedtools sort -i pred_peak_monkey_sub_filter.bed >pred_peak_monkey_sub_filter_sort.bed &
#bedRemoveOverlap pred_peak_monkey_sub_filter_sort.bed pred_peak_monkey_sub_filter_sort_rmdup.bed
#GENOME=/file/path/prefix/fasc6/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa
#bedtools getfasta -tab -fi $GENOME -bed pred_peak_monkey_sub_filter_sort_rmdup.bed -fo Region.pred.139w.monkey.input.txt
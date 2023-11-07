import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("inBed")
parser.add_argument("outBed")
parser.add_argument("top",default=50000)

args = parser.parse_args()


bed = pd.read_csv(args.inBed,sep='\t',header=None)
bed.columns = ['chr','start','end','barcode']

d = {'UM':[bed.shape[0]]}
UM_UQ = pd.DataFrame(data=d)

bed = bed.drop_duplicates()
UM_UQ['UQ'] = bed.shape[0]
UM_UQ['UM/UQ'] = UM_UQ['UM'] / UM_UQ['UQ']
UM_UQ.to_csv('UM_UQ.txt',sep='\t')


## 1.extract the top x barcodes
# sort
stat = pd.DataFrame(pd.value_counts(bed.iloc[:,3]))
stat.columns = ['UQ']
stat = stat.sort_values(by='UQ',axis=0,ascending=False)

# top x barcodes
stat_top = stat.iloc[range(int(args.top)),]

## 2.filter barcodes that UQ<10 
stat_pass = stat_top.loc[stat['UQ']>10,:]
stat_pass

## top
bc_pass = stat_pass.index
bed_pass = bed.loc[bed['barcode'].isin(bc_pass),:]

bed_pass.to_csv(args.outBed,sep='\t',header=None,index=None)

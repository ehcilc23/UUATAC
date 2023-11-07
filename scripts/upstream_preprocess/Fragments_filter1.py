import sys
import pickle
import pandas as pd

def bedClean(in_bed,out_bed):
	# read in bamfile and get header
	bf = open(in_bed,"r")
	
	with open(out_bed,"w") as outf:
		lines = bf.readlines()
		for f in lines:
			line = f.split('\t')
			chr = line[0]
			bc = line[6].split(':')[0]
			if line[8] == "+":
				start = int(line[1]) + 1
				end = int(line[5]) + 1
			elif line[8] == "-":
				start = int(line[4]) + 1
				end = int(line[2]) + 1
			if end - start > 10 and end - start < 1000:
				line_out = chr + "\t" + str(start) +"\t"+ str(end) +"\t"+ bc + "\n"
				outf.write(line_out)
		outf.close()

if __name__ == "__main__":
	in_bed = sys.argv[1]
	out_bed = sys.argv[2]

	bedClean(in_bed,out_bed)
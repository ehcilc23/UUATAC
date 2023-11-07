import sys
import pickle
import pysam


def correct_barcodes(barcodepath, in_bam, col, out_bam):
	# load the barcode_list
	print("Load barcode whitelist dictionary...")
	bc1 = open(barcodepath+"/BS_BC_61_76.pickle2","rb")
	barcode_list1 = pickle.load(bc1)
	bc1.close()
	RT_bc2 = open(barcodepath+"/AS_BC_1_768.pickle2","rb")
	barcode_list2 = pickle.load(RT_bc2)
	RT_bc2.close()
	RT_bc3 = open(barcodepath+"/Tn5_1_192.pickle2","rb")
	barcode_list3 = pickle.load(RT_bc3)
	RT_bc3.close()

	print("Start to correct barcodes...")

	# read in bamfile and get header
	bf = pysam.AlignmentFile(in_bam, 'rb', check_sq=False)
	bf_head_dict = dict(bf.header)
	# count total and filtered lines
	total_line = 0
	filtered_line = 0
	with pysam.AlignmentFile(out_bam, "wb", header=bf_head_dict) as outf:
		for r in bf:
			total_line += 1
			
			barcode1 = r.get_tag('CB')
			barcode = r.get_tag('RT')
			barcode2 = barcode[:10]
			barcode3 = barcode[10:]
			if barcode1 in barcode_list1 and barcode2 in barcode_list2 and barcode3 in barcode_list3:
				filtered_line += 1
				rtbarcode_match1 = barcode_list1[barcode1]
				rtbarcode_match2 = barcode_list2[barcode2]
				rtbarcode_match3 = barcode_list3[barcode3]
				barcode_corrected = col + rtbarcode_match1 + rtbarcode_match2 + rtbarcode_match3

				r.set_tag('BC', barcode_corrected)
				r.set_tag('TN', rtbarcode_match3)
				outf.write(r)
		
		outf.close()

	# message filter rate
	print("total line: %f, filtered line: %f, filter rate: %f"%(total_line, filtered_line, float(filtered_line) / float(total_line)))
	com_message = '''~~~~~~~~~~~~~~~correct barcodes done~~~~~~~~~~~~~~~~~~'''
	print(com_message)


if __name__ == "__main__":
	barcodepath = sys.argv[1]
	in_bam = sys.argv[2]
	col = sys.argv[3]
	out_bam = sys.argv[4]

	correct_barcodes(barcodepath, in_bam, col, out_bam)

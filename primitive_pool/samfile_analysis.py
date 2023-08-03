# import pysam
import pandas as pd
import csv
import re
import matplotlib.pyplot as plt
import numpy as np
import sys
from cigar import Cigar

mapped_reads_dict = {}
unmapped_reads_dict = {}
chimeric_reads_2048_dict = {}
chimeric_reads_2064_dict = {}

unique_S255s1_reads = {}
unique_S236_reads = {}
unique_S155L_reads = {}

S255s1_reads_counter = 0
S155L_reads_counter = 0
S236s_reads_counter = 0
Ecoli_reads_counter = 0


duplicate_list = []
line_counter = 0

chimeric_seq_len_list = []
SA_len_list = []

csv.field_size_limit(sys.maxsize)

#/Users/macbookpro/Documents/uni/Msc/Thesis Msc/sample_minimap_catref_LS36concatenated.sam (first try)
#/Users/macbookpro/Documents/uni/Msc/minimap_catref_LS36concatenated_canucorrected.sam (Michael)
#/Users/macbookpro/Documents/uni/Msc/minimap_catref_LS36concatenated_canucorrectedparameter.sam (Michael)
#/Users/macbookpro/Documents/uni/Msc/minimap_LS36_assembly_allreads_regularsettingsfile_sept2021f_versus_cat_references_Ecolibl21DE_NCBIEcoliDE3_LE392_S155L_S236s_S255s1.sam"

with open("/Users/macbookpro/Documents/p1_S155l.sam",'r') as infile:
	fieldnames = ['qname','flag','rname','pos','mapq','cigar','rnext','pnext','tlen','seq','qual','nm','ms','AS','nn','tp','cm','s1','s2','de','SA','rl']
	
	reader = csv.DictReader(infile,delimiter='\t',fieldnames=fieldnames)
	for line in reader:
		if "@" in line:
			continue
		else:
			#fill dict with unmapped reads based on flag "4"
			line_counter +=1
			if line['flag'] == "4":
				dict_id_name = line['qname']
				# print(dict_id_name)
				value_data = {}
				value_data['qname'] = line['qname']
				value_data['flag'] = line['flag']
				value_data['rname'] = line['rname']
				value_data['pos'] = line['pos']
				value_data['mapq'] = line['mapq']
				value_data['cigar'] = line['cigar']
				value_data['seq'] = line['seq']
				unmapped_reads_dict[dict_id_name] = value_data


			# elif line['flag'] != "4":
			# 	dict_id_name = line['qname']
			# 	# print(dict_id_name)
			# 	value_data = {}
			# 	value_data['qname'] = line['qname']
			# 	value_data['flag'] = line['flag']
			# 	value_data['rname'] = line['rname']
			# 	value_data['pos'] = line['pos']
			# 	value_data['mapq'] = line['mapq']
			# 	value_data['SA'] = line['SA']
			# 	mapped_reads_dict[dict_id_name] = value_data

			#fill chimeric reads dict with chimeric reads based on flag 2048 (supplementary alignment)
			elif line['flag'] == "2048":
				dict_id_name = line['qname']
				value_data = {}
				value_data['qname'] = line['qname']
				value_data['flag'] = line['flag']
				value_data['rname'] = line['rname']
				value_data['pos'] = line['pos']
				value_data['mapq'] = line['mapq']
				value_data['cigar'] = line['cigar']
				value_data['seq'] = line['seq']
				value_data['SA'] = line['SA']
				chimeric_reads_2048_dict[dict_id_name] = value_data
				seq_len = len(line['seq'])
				chimeric_seq_len_list.append(seq_len)

				SA_string = line['SA']
				# print("SA_string",SA_string)
				if 'SA' in SA_string:
					multiple_value_list = re.split(";",SA_string)
					multiple_value_list_filtered = list(filter(None, multiple_value_list))
					# print(multiple_value_list_filtered)
					multiple_value_list_filtered[0] = re.sub("(SA:Z:)","",multiple_value_list[0])
					# print("multiple_value_list_filtered:",multiple_value_list_filtered)
					
					if len(multiple_value_list_filtered) > 1:
						# print("multiple_value_list_filtered:",multiple_value_list_filtered)
						for string in multiple_value_list_filtered:
							string_list = re.split(",",string)
							# print("stringlist",string_list)
							cigar_string = string_list[3]
							# print("cigar_string",cigar_string)
							c = Cigar(cigar_string)
							# print(c)
							SA_cigar_len = len(c)
							SA_len_list.append(SA_cigar_len)

					elif len(multiple_value_list_filtered) == 1:
						for string in multiple_value_list_filtered:
							string_list = re.split(",",string)
							# print("stringlist",string_list)
							cigar_string = string_list[3]
							# print("cigar_string",cigar_string)
							c = Cigar(cigar_string)
							# print(c)
							SA_cigar_len = len(c)
							SA_len_list.append(SA_cigar_len)

			#fill chimeric reads dict with chimeric reads based on flag 2064 (supplementary alignment reverse strand)
			elif line['flag'] == "2064":
				dict_id_name = line['qname']
				value_data = {}
				value_data['qname'] = line['qname']
				value_data['flag'] = line['flag']
				value_data['rname'] = line['rname']
				value_data['pos'] = line['pos']
				value_data['mapq'] = line['mapq']
				value_data['cigar'] = line['cigar']
				value_data['seq'] = line['seq']
				value_data['SA'] = line['SA']
				chimeric_reads_2048_dict[dict_id_name] = value_data
				#NOTE: HAVE CHANGED ABOVE TO 2048 TO COLLECT ALL SECONDARY READS FOR PLOTTING
				
				SA_string = line['SA']
				# print("SA_string",SA_string)
				if 'SA' in SA_string:
					multiple_value_list = re.split(";",SA_string)
					multiple_value_list_filtered = list(filter(None, multiple_value_list))
					# print(multiple_value_list_filtered)
					multiple_value_list_filtered[0] = re.sub("(SA:Z:)","",multiple_value_list[0])
					# print("multiple_value_list_filtered:",multiple_value_list_filtered)
					
					if len(multiple_value_list_filtered) > 1:
						# print("multiple_value_list_filtered:",multiple_value_list_filtered)
						for string in multiple_value_list_filtered:
							string_list = re.split(",",string)
							# print("stringlist",string_list)
							cigar_string = string_list[3]
							# print("cigar_string",cigar_string)
							c = Cigar(cigar_string)
							# print(c)
							SA_cigar_len = len(c)
							SA_len_list.append(SA_cigar_len)

					elif len(multiple_value_list_filtered) == 1:
						for string in multiple_value_list_filtered:
							string_list = re.split(",",string)
							# print("stringlist",string_list)
							cigar_string = string_list[3]
							# print("cigar_string",cigar_string)
							c = Cigar(cigar_string)
							# print(c)
							SA_cigar_len = len(c)
							SA_len_list.append(SA_cigar_len)

print("chimeric_reads_2048_dict",len(chimeric_reads_2048_dict.keys()))
print("chimeric_reads_2064",len(chimeric_reads_2064_dict.keys()))
print("chimeric_seq_len_list",len(chimeric_seq_len_list))
print("SA_len_list",len(SA_len_list))
print("line_counter",line_counter)
print(chimeric_seq_len_list)

# plt.title('SA_len_list_total_redo_normal')
# binlist = range(0,10000,10)
# plt.hist(chimeric_seq_len_list,bins=binlist)
# plt.hist(SA_len_list,bins=binlist)
# plt.show()

# primary alignments filtered into dictionaries, now filter for prim.a to genome X, SA alignment to genome Y/Z/etc.
def iterdict(read_dictionary):
	for seq_id,value_data in read_dictionary.items():
		# print(seq_id,value_data)
		SA_string = value_data['SA']
		multiple_value_list = re.split(";",SA_string)
		multiple_value_list_filtered = list(filter(None, multiple_value_list))
		multiple_value_list_filtered[0] = re.sub("(SA:Z:)","",multiple_value_list[0])
		# print("multiple_value_list_filtered:",multiple_value_list_filtered)
		SA_dict = {}
		for string in multiple_value_list_filtered:
			string_list = re.split(",",string)
			# print(string_list)
			genome_name_dict_key = string_list[0]
			genome_value_dict = {}
			genome_value_dict[genome_name_dict_key] = string_list[1:]
			# print(genome_value_dict)
			SA_dict[genome_name_dict_key] = string_list
		# print(SA_dict)

		read_dictionary[seq_id]['SA'] = SA_dict
	# return(read_dictionary)

	mapping_dict = {}
	
	prim_S255s1_sec_S155L = 0
	prim_S255s1_sec_S255s1 = 0
	prim_S255s1_sec_S236s = 0
	prim_S255s1_sec_Ecoli = 0
	prim_S255s1_sec_NCBIEcoli = 0
	prim_S255s1_sec_LE392 = 0
	prim_S255s1_sec_double = 0


	prim_S155L_sec_S255s1 = 0
	prim_S155L_sec_S155L = 0
	prim_S155L_sec_S236s = 0
	prim_S155L_sec_Ecoli = 0
	prim_S155L_sec_NCBIEcoli = 0
	prim_S155L_sec_LE392 = 0
	prim_S155L_sec_double = 0

	prim_S236s_sec_S155L = 0
	prim_S236s_sec_S236s = 0
	prim_S236s_sec_S255s1 = 0
	prim_S236s_sec_Ecoli = 0
	prim_S236s_sec_NCBIEcoli = 0
	prim_S236s_sec_LE392 = 0
	prim_S236s_sec_double = 0

	prim_Ecoli_sec_S155L = 0
	prim_Ecoli_sec_S236s = 0
	prim_Ecoli_sec_S255s1 = 0
	prim_Ecoli_sec_Ecoli = 0
	prim_Ecoli_sec_NCBIEcoli = 0
	prim_Ecoli_sec_LE392 = 0
	prim_Ecoli_sec_double = 0

	prim_NCBIEcoli_sec_S155L = 0
	prim_NCBIEcoli_sec_S236s = 0
	prim_NCBIEcoli_sec_S255s1 = 0
	prim_NCBIEcoli_sec_Ecoli = 0
	prim_NCBIEcoli_sec_NCBIEcoli = 0
	prim_NCBIEcoli_sec_LE392 = 0
	prim_NCBIEcoli_sec_double = 0
	
	prim_LE392_sec_S155L = 0
	prim_LE392_sec_S236s = 0
	prim_LE392_sec_S255s1 = 0
	prim_LE392_sec_Ecoli = 0
	prim_LE392_sec_NCBIEcoli = 0
	prim_LE392_sec_LE392 = 0
	prim_LE392_sec_double = 0
	for seq_id,value_data in read_dictionary.items():
		# print(seq_id)
		# print(value_data)
		primary_mapping_genome_name = value_data['rname']
		SA_dictionary = value_data['SA']
		if value_data['rname'] == 'S2_55s1':
			SA_dictionary_keys_list = SA_dictionary.keys()
			# print(SA_dictionary_keys_list)
			if len(SA_dictionary_keys_list) >1:
				prim_S255s1_sec_double +=1
			elif len(SA_dictionary_keys_list) == 1:
				if "S2_55s1" in SA_dictionary_keys_list:
					prim_S255s1_sec_S255s1 +=1
				elif "S1_55L" in SA_dictionary_keys_list:
					print(seq_id,value_data)
					prim_S255s1_sec_S155L +=1
				elif "S2_36s" in SA_dictionary_keys_list:
					prim_S255s1_sec_S236s +=1
				elif "E_coli_bl21_DE3_polished" in SA_dictionary_keys_list:
					prim_S255s1_sec_Ecoli +=1
				elif "CP053602.1" in SA_dictionary_keys_list:
					prim_S255s1_sec_NCBIEcoli +=1
				elif "LE392" in SA_dictionary_keys_list:
					prim_S255s1_sec_LE392 +=1

		elif value_data['rname'] == 'S1_55L':
			SA_dictionary_keys_list = SA_dictionary.keys()
			# print(SA_dictionary_keys_list)
			if len(SA_dictionary_keys_list) >1:
				prim_S155L_sec_double +=1
			elif len(SA_dictionary_keys_list) == 1:
				if "S2_55s1" in SA_dictionary_keys_list:
					prim_S155L_sec_S255s1 +=1
				elif "S1_55L" in SA_dictionary_keys_list:
					prim_S155L_sec_S155L +=1
				elif "S2_36s" in SA_dictionary_keys_list:
					prim_S155L_sec_S236s +=1
				elif "E_coli_bl21_DE3_polished" in SA_dictionary_keys_list:
					prim_S155L_sec_Ecoli +=1
				elif "CP053602.1" in SA_dictionary_keys_list:
					prim_S155L_sec_NCBIEcoli +=1
				elif "LE392" in SA_dictionary_keys_list:
					prim_S155L_sec_LE392 +=1

		elif value_data['rname'] == 'S2_36s':
			SA_dictionary_keys_list = SA_dictionary.keys()
			# print(SA_dictionary_keys_list)
			if len(SA_dictionary_keys_list) >1:
				prim_S236s_sec_double +=1
			elif len(SA_dictionary_keys_list) == 1:
				if "S2_55s1" in SA_dictionary_keys_list:
					prim_S236s_sec_S255s1 +=1
				elif "S1_55L" in SA_dictionary_keys_list:
					prim_S236s_sec_S155L +=1
				elif "S2_36s" in SA_dictionary_keys_list:
					prim_S236s_sec_S236s +=1
				elif "E_coli_bl21_DE3_polished" in SA_dictionary_keys_list:
					prim_S236s_sec_Ecoli +=1
				elif "CP053602.1" in SA_dictionary_keys_list:
					prim_S236s_sec_NCBIEcoli +=1
				elif "LE392" in SA_dictionary_keys_list:
					prim_S236s_sec_LE392 +=1

		elif value_data['rname'] == 'E_coli_bl21_DE3_polished':
			SA_dictionary_keys_list = SA_dictionary.keys()
			# print(SA_dictionary_keys_list)
			if len(SA_dictionary_keys_list) >1:
				prim_Ecoli_sec_double +=1
			elif len(SA_dictionary_keys_list) == 1:
				if "S2_55s1" in SA_dictionary_keys_list:
					prim_Ecoli_sec_S255s1 +=1
				elif "S1_55L" in SA_dictionary_keys_list:
					prim_Ecoli_sec_S155L +=1
				elif "S2_36s" in SA_dictionary_keys_list:
					prim_Ecoli_sec_S236s +=1
				elif "E_coli_bl21_DE3_polished" in SA_dictionary_keys_list:
					prim_Ecoli_sec_Ecoli +=1
				elif "CP053602.1" in SA_dictionary_keys_list:
					prim_Ecoli_sec_NCBIEcoli +=1
				elif "LE392" in SA_dictionary_keys_list:
					prim_Ecoli_sec_LE392 +=1

		elif value_data['rname'] == 'CP053602.1':
			SA_dictionary_keys_list = SA_dictionary.keys()
			# print(SA_dictionary_keys_list)
			if len(SA_dictionary_keys_list) >1:
				prim_Ecoli_sec_double +=1
			elif len(SA_dictionary_keys_list) == 1:
				if "S2_55s1" in SA_dictionary_keys_list:
					prim_NCBIEcoli_sec_S255s1 +=1
				elif "S1_55L" in SA_dictionary_keys_list:
					prim_NCBIEcoli_sec_S155L +=1
				elif "S2_36s" in SA_dictionary_keys_list:
					prim_NCBIEcoli_sec_S236s +=1
				elif "E_coli_bl21_DE3_polished" in SA_dictionary_keys_list:
					prim_NCBIEcoli_sec_Ecoli +=1
				elif "CP053602.1" in SA_dictionary_keys_list:
					prim_NCBIEcoli_sec_NCBIEcoli +=1
				elif "LE392" in SA_dictionary_keys_list:
					prim_NCBIEcoli_sec_LE392 +=1

		elif value_data['rname'] == 'LE392':
			SA_dictionary_keys_list = SA_dictionary.keys()
			# print(SA_dictionary_keys_list)
			if len(SA_dictionary_keys_list) >1:
				prim_Ecoli_sec_double +=1
			elif len(SA_dictionary_keys_list) == 1:
				if "S2_55s1" in SA_dictionary_keys_list:
					prim_LE392_sec_S255s1 +=1
				elif "S1_55L" in SA_dictionary_keys_list:
					prim_LE392_sec_S155L +=1
				elif "S2_36s" in SA_dictionary_keys_list:
					prim_LE392_sec_S236s +=1
				elif "E_coli_bl21_DE3_polished" in SA_dictionary_keys_list:
					prim_LE392_sec_Ecoli +=1
				elif "CP053602.1" in SA_dictionary_keys_list:
					prim_LE392_sec_NCBIEcoli +=1
				elif "LE392" in SA_dictionary_keys_list:
					prim_LE392_sec_LE392 +=1										

	print("prim_S255s1_sec_S155L",prim_S255s1_sec_S155L)
	print("prim_S255s1_sec_S255s1",prim_S255s1_sec_S255s1)
	print("prim_S255s1_sec_S236s",prim_S255s1_sec_S236s)
	print("prim_S255s1_sec_Ecoli",prim_S255s1_sec_Ecoli)
	print("prim_S255s1_sec_LE392",prim_S255s1_sec_LE392)
	print("prim_S255s1_sec_NCBIEcoli",prim_S255s1_sec_NCBIEcoli)
	print("prim_S255s1_sec_double",prim_S255s1_sec_double)


	print("prim_S155L_sec_S255s1",prim_S155L_sec_S255s1)
	print("prim_S155L_sec_S155L",prim_S155L_sec_S155L)
	print("prim_S155L_sec_S236s",prim_S155L_sec_S236s)
	print("prim_S155L_sec_Ecoli",prim_S155L_sec_Ecoli)
	print("prim_S155L_sec_LE392",prim_S155L_sec_LE392)
	print("prim_S155L_sec_NCBIEcoli",prim_S155L_sec_NCBIEcoli)
	print("prim_S155L_sec_double",prim_S155L_sec_double)

	print("prim_S236s_sec_S155L",prim_S236s_sec_S155L)
	print("prim_S236s_sec_S236s",prim_S236s_sec_S236s)
	print("prim_S236s_sec_S255s1",prim_S236s_sec_S255s1)
	print("prim_S236s_sec_Ecoli",prim_S236s_sec_Ecoli)
	print("prim_S236s_sec_LE392",prim_S236s_sec_LE392)
	print("prim_S236s_sec_NCBIEcoli",prim_S236s_sec_NCBIEcoli)
	print("prim_S236s_sec_double",prim_S236s_sec_double)

	print("prim_Ecoli_sec_S155L",prim_Ecoli_sec_S155L)
	print("prim_Ecoli_sec_S236s",prim_Ecoli_sec_S236s)
	print("prim_Ecoli_sec_S255s1",prim_Ecoli_sec_S255s1)
	print("prim_Ecoli_sec_Ecoli",prim_Ecoli_sec_Ecoli)
	print("prim_Ecoli_sec_LE392",prim_Ecoli_sec_LE392)
	print("prim_Ecoli_sec_NCBIEcoli",prim_Ecoli_sec_NCBIEcoli)
	print("prim_Ecoli_sec_double",prim_Ecoli_sec_double)

	print("prim_NCBIEcoli_sec_S155L",prim_NCBIEcoli_sec_S155L)
	print("prim_NCBIEcoli_sec_S255s1",prim_NCBIEcoli_sec_S255s1)
	print("prim_NCBIEcoli_sec_S236s",prim_NCBIEcoli_sec_S236s)
	print("prim_NCBIEcoli_sec_Ecoli",prim_NCBIEcoli_sec_Ecoli)
	print("prim_NCBIEcoli_sec_LE392",prim_NCBIEcoli_sec_LE392)
	print("prim_NCBIEcoli_sec_NCBIEcoli",prim_NCBIEcoli_sec_NCBIEcoli)
	print("prim_NCBIEcoli_sec_double",prim_NCBIEcoli_sec_double)

	print("prim_LE392_sec_S155L",prim_LE392_sec_S155L)
	print("prim_LE392_sec_S255s1",prim_LE392_sec_S255s1)
	print("prim_LE392_sec_S236s",prim_LE392_sec_S236s)
	print("prim_LE392_sec_Ecoli",prim_LE392_sec_Ecoli)
	print("prim_LE392_sec_LE392",prim_LE392_sec_LE392)
	print("prim_LE392_sec_NCBIEcoli",prim_LE392_sec_NCBIEcoli)
	print("prim_LE392_sec_double",prim_LE392_sec_double)

	S255s1_bar_height = int(prim_S255s1_sec_S155L+prim_S255s1_sec_S255s1+prim_S255s1_sec_S236s+prim_S255s1_sec_Ecoli+prim_S255s1_sec_double+prim_S255s1_sec_LE392+prim_S255s1_sec_NCBIEcoli)
	S155L_bar_height = int(prim_S155L_sec_S155L+prim_S155L_sec_S255s1+prim_S155L_sec_S236s+prim_S155L_sec_Ecoli+prim_S155L_sec_double+prim_S155L_sec_LE392+prim_S155L_sec_NCBIEcoli)
	S236s_bar_height = int(prim_S236s_sec_S155L+prim_S236s_sec_S255s1+prim_S236s_sec_S236s+prim_S236s_sec_Ecoli+prim_S236s_sec_double+prim_S155L_sec_LE392+prim_S155L_sec_NCBIEcoli)
	Ecoli_bar_height = int(prim_Ecoli_sec_S155L+prim_Ecoli_sec_S255s1+prim_Ecoli_sec_S236s+prim_Ecoli_sec_Ecoli+prim_Ecoli_sec_double+prim_S155L_sec_LE392+prim_S155L_sec_NCBIEcoli)
	NCBIEcoli_bar_height = int(prim_NCBIEcoli_sec_S155L+prim_NCBIEcoli_sec_S255s1+prim_NCBIEcoli_sec_S236s+prim_NCBIEcoli_sec_Ecoli+prim_NCBIEcoli_sec_double+prim_NCBIEcoli_sec_LE392+prim_NCBIEcoli_sec_NCBIEcoli)
	LE392_bar_height = int(prim_LE392_sec_S155L+prim_LE392_sec_S255s1+prim_LE392_sec_S236s+prim_LE392_sec_Ecoli+prim_LE392_sec_double+prim_LE392_sec_LE392+prim_LE392_sec_NCBIEcoli)

	print("s255s1_bar_height",S255s1_bar_height)
	print("s155L_bar_height",S155L_bar_height)
	print("S236s_bar_height",S236s_bar_height)
	print("Ecoli_bar_height",Ecoli_bar_height)
	print("NCBIEcoli_bar_height",NCBIEcoli_bar_height)
	print("LE392_bar_height",LE392_bar_height)

	# PLOT number of reads with secondary mappings per genome
	genomes=["S1_55L","S2_55s1","S2_36s","Ecoli","NCBI_Ecoli","LE392"]
	sec_mapping_nr = {
		"sec_S155L":[prim_S155L_sec_S155L,prim_S255s1_sec_S155L,prim_S236s_sec_S155L,prim_Ecoli_sec_S155L,prim_NCBIEcoli_sec_S155L,prim_LE392_sec_S155L],
		"sec_S255s1":[prim_S155L_sec_S255s1,prim_S155L_sec_S255s1,prim_S155L_sec_S255s1,prim_Ecoli_sec_S255s1,prim_NCBIEcoli_sec_S255s1,prim_LE392_sec_S255s1],
		"sec_S236s":[prim_S155L_sec_S236s,prim_S255s1_sec_S236s,prim_S236s_sec_S236s,prim_Ecoli_sec_S236s,prim_NCBIEcoli_sec_S236s,prim_LE392_sec_S236s],
		"sec_Ecoli":[prim_S155L_sec_Ecoli,prim_S255s1_sec_Ecoli,prim_S236s_sec_Ecoli,prim_Ecoli_sec_Ecoli,prim_NCBIEcoli_sec_Ecoli,prim_LE392_sec_Ecoli],
		"sec_NCBIEcoli":[prim_S155L_sec_NCBIEcoli,prim_S255s1_sec_NCBIEcoli,prim_S236s_sec_NCBIEcoli,prim_Ecoli_sec_NCBIEcoli,prim_NCBIEcoli_sec_NCBIEcoli,prim_LE392_sec_NCBIEcoli],
		"sec_LE392":[prim_S155L_sec_LE392,prim_S255s1_sec_LE392,prim_S236s_sec_LE392,prim_Ecoli_sec_LE392,prim_NCBIEcoli_sec_LE392,prim_LE392_sec_LE392],
		"sec_multiple":[prim_S155L_sec_double,prim_S255s1_sec_double,prim_S236s_sec_double,prim_Ecoli_sec_double,prim_NCBIEcoli_sec_double,prim_LE392_sec_double]
	}

	df = pd.DataFrame(sec_mapping_nr,index=genomes)
	# df.plot(kind="bar",stacked=True,figsize=(10,8),rot=0)
	# plt.legend(loc="upper right")
	# plt.title('Number of reads with secondary mappings from corrected primitive pool')

	genomes=["S1_55L","S2_55s1","S2_36s","Ecoli","NCBI_Ecoli","LE392","unmapped"]
	mapping_nr={'Number of assemblies':[48,328,37,0,0,0,25]}
	df = pd.DataFrame(mapping_nr,index=genomes)
	df.plot(kind="bar",stacked=True,figsize=(10,8),rot=0)
	plt.title('Number of assemblies from primitive pool mapping to references')
	plt.show()	

print(iterdict(chimeric_reads_2048_dict))
# print(iterdict(chimeric_reads_2064_dict))

infile.close()

with open("/Users/macbookpro/Documents/uni/Msc/minimap_aug_LS36corrected_S155L_S236s_S255s1_Ecolibl21DE_LE392.sam",'r') as infile:
	fieldnames = ['qname','flag','rname','pos','mapq','cigar','rnext','pnext','tlen','seq','qual','nm','ms','AS','nn','tp','cm','s1','s2','de','SA','rl']
	
	reader = csv.DictReader(infile,delimiter='\t',fieldnames=fieldnames)
	for line in reader:
		if "@" in line:
			continue
		else:
			if line['rname'] == "S2_36s":
				S236s_reads_counter +=1
			elif line['rname'] == "S1_55L":
				S155L_reads_counter +=1
			elif line['rname'] == "S2_55s1":
				S255s1_reads_counter +=1
			elif line['rname'] == "E_coli_bl21_DE3_polished":
				Ecoli_reads_counter +=1
print("S255s1_reads_counter",S255s1_reads_counter)
print("S155L_reads_counter",S155L_reads_counter)
print("S236s_reads_counter",S236s_reads_counter)
print("Ecoli_reads_counter",Ecoli_reads_counter)

# S255s1_single_map_reads = int(S255s1_reads_counter) - int(prim_S255s1_sec_S155L) - int(prim_S255s1_sec_S255s1) - int(prim_S255s1_sec_S236s) - int(prim_S255s1_sec_Ecoli) - int(prim_S255s1_sec_double)
# print(prim_S255s1_sec_S155L)

# print("prim_S255s1_sec_S155L",prim_S255s1_sec_S155L)
# 	print("prim_S255s1_sec_S255s1",prim_S255s1_sec_S255s1)
# 	print("prim_S255s1_sec_S236s",prim_S255s1_sec_S236s)
# 	print("prim_S255s1_sec_Ecoli",prim_S255s1_sec_Ecoli)
# 	print("prim_S255s1_sec_double",prim_S255s1_sec_double)
# 	print("prim_S155L_sec_S255s1",prim_S155L_sec_S255s1)
# 	print("prim_S155L_sec_S155L",prim_S155L_sec_S155L)
# 	print("prim_S155L_sec_S236s",prim_S155L_sec_S236s)
# 	print("prim_S155L_sec_Ecoli",prim_S155L_sec_Ecoli)
# 	print("prim_S155L_sec_double",prim_S155L_sec_double)
# 	print("prim_S236s_sec_S155L",prim_S236s_sec_S155L)
# 	print("prim_S236s_sec_S236s",prim_S236s_sec_S236s)
# 	print("prim_S236s_sec_S255s1",prim_S236s_sec_S255s1)
# 	print("prim_S236s_sec_Ecoli",prim_S236s_sec_Ecoli)
# 	print("prim_S236s_sec_double",prim_S236s_sec_double)
# 	print("prim_Ecoli_sec_S155L",prim_Ecoli_sec_S155L)
# 	print("prim_Ecoli_sec_S236s",prim_Ecoli_sec_S236s)
# 	print("prim_Ecoli_sec_S255s1",prim_Ecoli_sec_S255s1)
# 	print("prim_Ecoli_sec_Ecoli",prim_Ecoli_sec_Ecoli)
# 	print("prim_Ecoli_sec_double",prim_Ecoli_sec_double)

# counter = 0 
# set1 = set(chimeric_reads_2064_dict)
# set2 = set(chimeric_reads_2048_dict)
# for name in set1.intersection(set2):
#     print(name, chimeric_reads_2048_dict[name])
#     counter +=1

# print(counter# print('pauze')
# print(chimeric_reads_2048_dict)
# print(chimeric_reads_2048_dict['d47c59bd-5df2-44d3-9114-c11909dd470b'])
# print(chimeric_reads_2064_dict['d47c59bd-5df2-44d3-9114-c11909dd470b'])


# def Repeat(x): 
#     _size = len(x) 
#     repeated = [] 
#     for i in range(_size): 
#         k = i + 1
#         for j in range(k, _size): 
#             if x[i] == x[j] and x[i] not in repeated: 
#                 repeated.append(x[i]) 
#     return repeated 

infile.close()

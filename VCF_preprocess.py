


def Preprocessing_original_vcf_file(ori_vcf:str,new1_vcf:str, new2_vcf:str):
	"""
	This function splits all multiallele variants into two normal variants and convert
	them to one of the three categories substitution, insertion or deletion.
	There are changes in the Position because of this and maybe warnings, if the data is erroneous.
	:param ori_vcf: The ingoing file including path and name.
	:param new1_vcf: The first multiallele entry and all normal variants.
	:param new2_vcf:  The second multiallele variant entry and all normal variants.
	:return: Nothing.
	"""
	def formatMultiallelVariant(splitline:list)-> str:
		"""
		This inner function only converts the multiallele-entry in the ALT column.
		:param splitline: One vcf-data-row devided with split and it contains only the first or second ALT-Entry.
		:return: New vcf-data-row as a string.
		"""
		reflen = len(splitline[3])
		altlen = len(splitline[4])
		if reflen == 1 and altlen == 1:
			# sub,changes nothing
			return "\t".join(splitline)
		elif reflen == 1 :
			#insert, no need to change anything
			return "\t".join(splitline)
		elif altlen == 1:
			#deletion, no need to change anything
			return "\t".join(splitline)
		elif reflen > altlen:
			#Chr1	12997	.	CTT	C,CT
			#here:
			# Chr1	12997	.	CTT	CT
			# transform to
			# Chr1	12998	.	TT	T
			new_entry = [str(splitline[0])]#chr
			new_entry.append(str(int(splitline[1]) + altlen-1))#pos
			new_entry.append(str(splitline[2]))# .
			new_entry.append(str(splitline[3][altlen-1:])) # new ref
			new_entry.append(str(splitline[4][altlen-1])) # new alt
			new_entry.append("".join(splitline[5:]))
			return("\t".join(new_entry))
		elif reflen == altlen:
			# Chr1	12997	.	CT	C,CT (made up data)
			# here:
			# Chr1	12997	.	CT	CT
			# transform to
			# Chr1	12998	.	T	T
			new_entry = [str(splitline[0])] # chr
			new_entry.append(str(int(splitline[1]) + reflen - 1)) # pos
			new_entry.append(str(splitline[2])) # .
			new_entry.append(str(splitline[3][reflen -1])) # new ref
			new_entry.append(str(splitline[4][reflen -1])) # new alt
			new_entry.append("".join(splitline[5:])) # anything after
			return ("\t".join(new_entry))
		elif altlen > reflen:
			# here:
			# Chr1	993202	.	CA	CAA
			# transform to:
			# Chr1	993203	.	A	AA
			new_entry = [str(splitline[0])]  # chr
			new_entry.append(str(int(splitline[1]) + reflen-1))  # pos
			new_entry.append(str(splitline[2]))  # .
			new_entry.append(str(splitline[3][reflen - 1]))  # new ref
			new_entry.append(str(splitline[4][reflen - 1:]))  # new alt
			new_entry.append("".join(splitline[5:]))  # anything after
			return ("\t".join(new_entry))
		else:
			print("Warning: Preprocessing_original_vcf_file")
			print("\t".join(splitline))
			return("#error: " + "\t".join(splitline))


	ori_vcf_file = open(ori_vcf,"r")
	print("Preprocessing")
	line = ori_vcf_file.readline()
	lines1 = []
	lines2 = []
	while line:
		splitline = line.split("\t")
		if len(splitline) == 1:
			#no tabs or the empty line in the end
			line = ori_vcf_file.readline()
			continue

		elif line.startswith("#"):
			lines1.append(line)
			lines2.append(line)
		elif "," in splitline[4]:
			split_allele = splitline[4].split(",")
			#first allele
			new_entry = splitline[0:4]
			new_entry.append(split_allele[0])
			new_entry.append("\t".join(splitline[5:]))
			lines1.append(formatMultiallelVariant(new_entry))
			#second allele
			new_entry = splitline[0:4]
			new_entry.append(split_allele[1])
			new_entry.append("\t".join(splitline[5:]))
			lines2.append(formatMultiallelVariant(new_entry))
		else:
			lines1.append(line)
			lines2.append(line)
		line = ori_vcf_file.readline()
	ori_vcf_file.close()

	# sorting, first after chromosome, then after position
	# necessary, because the standardization of triallelvariation could change the order of the lines
	lines1 = sorted(lines1,key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
	lines2 = sorted(lines2, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
	#check if there are lines, which have the same position
	warning = True
	for i,line in enumerate(lines1):
		spline = line.split("\t")
		if i == 0 or len(spline) == 1:
			continue
		if spline[1] == lines1[i-1].split("\t")[1]:
			if warning:
				print("VCF-Preprocessing Warning for File 1:")
				warning = False
			print(str(lines1[i - 1]))
			print("Removed: " + str(line))
			lines1.remove(line)

	warning = True
	for i, line in enumerate(lines2):
		spline = line.split("\t")
		if i == 0 or len(spline) == 1:
			continue
		if spline[1] == lines2[i - 1].split("\t")[1]:
			if warning:
				print("VCF-Preprocessing Warning for File 2:")
				warning = False
			print(str(lines2[i - 1]))
			print("Removed:" + str(line))
			lines2.remove(line)



	nvcf_1 = open (new1_vcf,"w")
	nvcf_1.write("".join(lines1))
	nvcf_1.close()

	nvcf_2 = open (new2_vcf, "w")
	nvcf_2.write("".join(lines2))
	nvcf_2.close()
	print("preprocessing: done")


def vcf_preprocessing(invcf,outpath):
	"""
	Main function inside this module.
	There is no need to address any other function inside this module.
	:param invcf: Path and name of the original vcf file.
	:param outpath: Path to a existing folder, in which the new files will be created.
	:return: None.
	"""
	new1_vcf = outpath + "first.vcf"
	new2_vcf = outpath + "second.vcf"
	### split vcf file data rows into two new vcf files, because of the multiallele variants
	Preprocessing_original_vcf_file(invcf, new1_vcf, new2_vcf)

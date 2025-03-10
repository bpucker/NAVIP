__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"

from datetime import datetime

VCF_preproc_log = []
VCF_preproc_log_short = []

def preprocess_vcf_file(original_vcf:str, new_vcf1:str, new_vcf2:str, outpath:str):
	"""
	This function splits all multiallele variants into two normal variants and converts
	them to one of the three categories: substitution, insertion, or deletion.
	There are changes in the position because of this and maybe warnings, if the data is erroneous.
	:param original_vcf: The input file including path and name.
	:param new_vcf1: The first multiallele entry and all normal variants.
	:param new_vcf2: The second multiallele variant entry and all normal variants.
	:param outpath: Output path for the logfile.
	:return: Nothing.
	"""

	def format_multiallele_variant(splitline:list)-> str:
		"""
		This inner function only converts the multiallele entry in the ALT column.
		:param splitline: One VCF data row divided with split and it contains only the first or second ALT entry.
		:return: New VCF data row as a string.
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
		elif reflen != altlen:
			#Chr1	12997	.	CTT	C,CT
			#here:
			# Chr1	12997	.	CTT	CT
			# transform to
			# Chr1	12998	.	TT	T
			nref = splitline[3]
			nalt = splitline[4]
			skipped = 0
			removed = 0
			# remove all starting bases, if they are identical
			# chr1	222991	.	TG	TGAGAGAGAGAGAGAGAG,T
			# chr1	584140	.	TAAAA	TAAAAA,T
			# chr1	332796	.	C TCTCTTTCTCTCATTA	ATCTCTTTCTCTCATTA,C
			# chr1	332812	.	A	C

			#if splitline[0] == "chr1" and splitline[1] == "222991":
			#	print("bugsearch 222991")
			#if splitline[0] == "chr1" and splitline[1] == "584140":
			#	print("bugsearch 584140")
			#if splitline[0] == "chr1" and splitline[1] == "332796":
			#	print("bugsearch 332796")
			while nref[0] == nalt[0] and nref[1] == nalt[1]:
				nref = nref[1:]
				nalt = nalt[1:]
				skipped += 1
				if len(nref) == 1 or len(nalt) == 1:
					break
			if len(nref) > 1 and len(nalt) > 1:
				# remove end bases, if they are identical
				removed = 0
				#print(nref)
				#print(nalt)
				while nref[len(nref) - 1] == nalt[len(nalt) - 1]:
					nref = nref[0:len(nref) - 1]
					nalt = nalt[0:len(nalt) - 1]
					removed += 1
					if len(nref) == 1 or len(nalt) == 1:
						break
				if len(nref) > 1 and len(nalt) > 1:
					# something like: Chr1	12998	.	ATT	ACCC
					# but this (hopefully) never happens
					VCF_preproc_log.append("Warning: Preprocessing_original_vcf_file\n.")
					VCF_preproc_log.append(str("\t".join(splitline)))
					return "#error: " + "\t".join(splitline)
			#else:
				#ref, alt => one or both have length 2
				#skipped += 1
			new_entry = [str(splitline[0]),  # chr
						 str(int(splitline[1]) + skipped),  # pos
						 str(splitline[2]),  # .
						 str(splitline[3][skipped:len(splitline[3]) - removed]),  # new ref
						 str(splitline[4][skipped:len(splitline[4]) - removed]),  # new alt
						 "\t".join(splitline[5:])]  # anything after
			return "\t".join(new_entry)
		elif reflen == altlen:
			# Chr1	12997	.	CT	C,CT (made up data)
			# here:
			# Chr1	12997	.	CT	CT
			# transform to
			# Chr1	12998	.	T	T

			#chr1	332796	.	CTCTCTTTCTCTCATTA	ATCTCTTTCTCTCATTA,C
			#chr1	332796	.	CTCTCTTTCTCTCATTA	ATCTCTTTCTCTCATTA

			nref = splitline[3]
			nalt = splitline[4]
			skipped = 0
			# remove all starting bases, if they are identical
			while  nref[0] == nalt[0]:
				nref = nref[1:]
				nalt = nalt[1:]
				skipped += 1
				if len(nref) == 1 and len(nalt) ==1:
					if nref[0] == nalt[0]:
						VCF_preproc_log.append("Warning: Bases are somehow still identical: "
											   + str(nref)
											   + "\t"
											   + str(nalt)
											   + "\n")
						VCF_preproc_log.append("Original:\t" + str("\t".join(splitline)) + "\n")
						return "#error: " + "\t".join(splitline)
					break

			removed  = 0
			while  nref[len(nref)-1] == nalt[len(nalt)-1]:
				nref = nref[0:len(nref)-1]
				nalt = nalt[0:len(nalt)-1]
				removed += 1
				if len(nref) == 1 and len(nalt) ==1:
					if nref[0] == nalt[0]:
						VCF_preproc_log.append("Warning: Bases are somehow still identical: "
											   + str(nref)
											   + "\t"
											   + str(nalt)
											   + "\n")
						VCF_preproc_log.append("Original:\t" + str("\t".join(splitline)) + "\n")
						return "#error: " + "\t".join(splitline)
					break
			new_entry = [str(splitline[0]),  # chr
						 str(int(splitline[1]) + skipped),  # pos
						 str(splitline[2]),  # .
						 str(splitline[3][skipped:len(splitline[3]) - removed]),  # new ref
						 str(splitline[4][skipped:len(splitline[4]) - removed]),  # new alt
						 "".join(splitline[5:])]  # anything after
			return "\t".join(new_entry)
		elif altlen > reflen:
			print("mhhhhhh")
			# todo remove this case, if everything worked fine
			# here:
			# Chr1	993202	.	CA	CAA
			# transform to:
			# Chr1	993203	.	A	AA
			"""
			new_entry = [str(splitline[0])]  # chr
			new_entry.append(str(int(splitline[1]) + reflen-1))  # pos
			new_entry.append(str(splitline[2]))  # .
			new_entry.append(str(splitline[3][reflen - 1]))  # new ref
			new_entry.append(str(splitline[4][reflen - 1:]))  # new alt
			new_entry.append("".join(splitline[5:]))  # anything after
			return ("\t".join(new_entry))
			"""
			nref = splitline[3]
			nalt = splitline[4]
			skipped = 0
			# remove all starting bases, if they are identical
			while nref[0] == nalt[0]:
				nref = nref[1:]
				nalt = nalt[1:]
				skipped += 1
				if len(nref) == 1 or len(nalt) == 1:
					if nref[0] == nalt[0]:
						VCF_preproc_log.append("Warning: Bases are somehow still identical: "
											   + str(nref)
											   + "\t"
											   + str(nalt)
											   + "\n")
						VCF_preproc_log.append("Original:\t" + str("\t".join(splitline)) + "\n")
						return "#error: " + "\t".join(splitline)
					break

			if len(nref) > 1 and len(nalt) > 1:
				# something like: Chr1	12998	.	TT	CC
				# but this (hopefully) never happens
				VCF_preproc_log.append("Warning: Preprocessing_original_vcf_file\n.")
				VCF_preproc_log.append(str("\t".join(splitline)))
				return "#error: " + "\t".join(splitline)
			new_entry = [str(splitline[0]),  # chr
						 str(int(splitline[1]) + skipped),  # pos
						 str(splitline[2]),  # .
						 str(splitline[3][skipped:]),  # new ref
						 str(splitline[4][skipped:]),  # new alt
						 "".join(splitline[5:])]  # anything after
			return "\t".join(new_entry)
		else:
			print("Warning (critical?): Preprocessing_original_vcf_file (logfile)")
			VCF_preproc_log.append("Warning: Preprocessing_original_vcf_file\n.")
			VCF_preproc_log.append(str("\t".join(splitline)))
			return "#error: " + "\t".join(splitline)


	vcf_file = open(original_vcf, "r")
	print("Preprocessing")
	starttime = datetime.now()

	line = vcf_file.readline()
	infos = []
	lines1 = []
	lines2 = []
	format_multiallele_variant_warning = True
	while line:
		if line.startswith("#"):
			infos.append(line)
			line = vcf_file.readline()
			continue

		splitline = line.split("\t")
		if len(splitline) == 1:
			#no tabs or the empty line in the end
			print("No tabs or the empty line in the end.")
			line = vcf_file.readline()
			continue
		elif "," in splitline[4]:
			split_allele = splitline[4].split(",")

			#first allele
			new_entry = splitline[0:4]
			new_entry.append(split_allele[0])
			new_entry.append("\t".join(splitline[5:]))
			erg = format_multiallele_variant(new_entry)
			if not erg.startswith("#"):
				lines1.append(erg)
			elif format_multiallele_variant_warning:
				format_multiallele_variant_warning = False
				print("Warning: Preprocessing_original_vcf_file (logfile)")

			#second allele
			new_entry = splitline[0:4]
			new_entry.append(split_allele[1])
			new_entry.append("\t".join(splitline[5:]))
			erg = format_multiallele_variant(new_entry)
			if not erg.startswith("#"):
				lines2.append(erg)
			elif format_multiallele_variant_warning:
				format_multiallele_variant_warning = False
				print("Warning: Preprocessing_original_vcf_file (logfile)")
		else:
			if len(splitline[3]) > 1 and len(splitline[4]) > 1:
				new_entry = splitline[0:5]
				new_entry.append("\t".join(splitline[5:]))
				line = format_multiallele_variant(new_entry)
			lines1.append(line)
			lines2.append(line)
		line = vcf_file.readline()
	vcf_file.close()

	# sorting, first after chromosome, then after position
	# necessary, because the standardization of variations could change the order of the lines
	lines1 = sorted(lines1, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
	lines2 = sorted(lines2, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
	#check if there are lines, which have the same position

	print("1/3 done: " + str(datetime.now()- starttime))
	warning = True
	i = 0
	max = len(lines1)-1
	while i <= max:
		spline = lines1[i].split("\t")
		if i == 0 or len(spline) == 1:
			i += 1
			continue
		if spline[1] == lines1[i-1].split("\t")[1]:
			if warning:
				print("VCF-Preprocessing Warning for File 1 in logfile.")
				VCF_preproc_log_short.append("VCF-Preprocessing Warning for File 1:\n")
				warning = False
			VCF_preproc_log_short.append("Take:\t" + str("\t".join(lines1[i - 1].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append(str(lines1[i - 1]))
			VCF_preproc_log_short.append("Rem.:\t" + str("\t".join(lines1[i].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append("Removed: " + str(lines1[i]))
			max -= 1
			del lines1[i]
			continue
		i += 1

	print("2/3 done: " + str(datetime.now() - starttime))
	warning = True
	i = 0
	max = len(lines2) - 1
	while i <= max:
		spline = lines2[i].split("\t")
		if i == 0 or len(spline) == 1:
			i += 1
			continue
		if spline[1] == lines2[i - 1].split("\t")[1]:
			if warning:
				print("VCF-Preprocessing Warning for File 2 in logfile.")
				VCF_preproc_log_short.append("VCF-Preprocessing Warning for File 2:\n")
				warning = False
			VCF_preproc_log_short.append("Take:\t" + str("\t".join(lines2[i - 1].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append(str(lines2[i - 1]))
			VCF_preproc_log_short.append("Rem.:\t" + str("\t".join(lines2[i].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append("Removed: " + str(lines2[i]))
			max -= 1
			del lines2[i]
			continue
		i += 1

	def remove_data(lines:list) -> list:
		i = 0
		removed_entries = 0
		pos = 0
		oldline = ""
		lines_new = []
		for line in lines:
			spline = line.split("\t")
			if int(spline[1]) > pos:
				oldline = line
				pos = int(spline[1])
				pos += -1 + len(spline[3])  # + length of ref, for deletions
				i += 1
				lines_new.append(line)
			elif spline[0] == oldline.split("\t")[0]:
				# print("Removed:\n" + line + "Because of:\n" + oldline)
				VCF_preproc_log_short.append("Rem.:\t" + str("\t".join(lines[i].split("\t")[0:5])) + "\n")
				VCF_preproc_log.append("Removed: " + str(lines[i]))
				VCF_preproc_log_short.append("Because:\t" + str("\t".join(oldline.split("\t")[0:5])) + "\n")
				VCF_preproc_log.append(oldline)
				#del lines[i]
				removed_entries += 1
				i +=1
			else:
				lines_new.append(line)
				oldline = line
				pos = int(spline[1])
				pos += -1 + len(spline[3])  # + length of ref, for deletions
				i += 1
		print(str(removed_entries) + " lines removed, because of conflicting data (see log).")
		return lines_new

	lines1 = sorted(lines1, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
	lines2 = sorted(lines2, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))

	VCF_preproc_log_short.append("Data Conflicts in File 1:\n")
	lines1 = remove_data(lines1)
	print("Remove data 1 done: "+ str(datetime.now() - starttime))
	VCF_preproc_log_short.append("Data Conflicts in File 2:\n")
	lines2 = remove_data(lines2)
	print("Remove data 2 done: " + str(datetime.now() - starttime))

	nvcf1 = open (new_vcf1, "w")
	nvcf1.write("".join(infos))
	nvcf1.write("".join(lines1))
	nvcf1.close()

	nvcf2 = open (new_vcf2, "w")
	nvcf2.write("".join(infos))
	nvcf2.write("".join(lines2))
	nvcf2.close()
	print("preprocessing done: " + str(datetime.now() - starttime))

	print("Start writing log.")
	logfile = open(outpath + "logfile.txt", "w")
	logfile.write("".join(VCF_preproc_log_short) + "\n" + "".join(VCF_preproc_log))
	logfile.close()
	print("Colliding variants (half of them removed): " + str(len(VCF_preproc_log_short) - 2))
	print("Writing log done: " + str(datetime.now() - starttime))


def vcf_preprocessing(invcf, outpath):
	"""
	Main function inside this module.
	There is no need to address any other function inside this module.
	:param invcf: Path and name of the original vcf file.
	:param outpath: Path to an existing folder, in which the new files will be created.
	:return: None.
	"""
	new_vcf1 = outpath + "first.vcf"
	new_vcf2 = outpath + "second.vcf"
	### split vcf file data rows into two new vcf files, because of the multiallele variants
	preprocess_vcf_file(invcf, new_vcf1, new_vcf2, outpath)

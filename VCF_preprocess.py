__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"

from datetime import datetime

VCF_preproc_log = []
VCF_preproc_log_in_short = []

def Preprocessing_original_vcf_file(ori_vcf:str,new1_vcf:str, new2_vcf:str, outpath:str):
	"""
	This function splits all multiallele variants into two normal variants and convert
	them to one of the three categories substitution, insertion or deletion.
	There are changes in the Position because of this and maybe warnings, if the data is erroneous.
	:param ori_vcf: The ingoing file including path and name.
	:param new1_vcf: The first multiallele entry and all normal variants.
	:param new2_vcf:  The second multiallele variant entry and all normal variants.
	:param outpath: Outpath for the logfile.
	:return: Nothing.
	"""

	def formatMultiallelVariant(splitline:list)-> str:
		"""
		This inner function only converts the multiallele-entry in the ALT column.
		:param splitline: One vcf-data-row divided with split and it contains only the first or second ALT-Entry.
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
		elif reflen > altlen or reflen < altlen:
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
					return ("#error: " + "\t".join(splitline))
			#else:
				#ref, alt => one or both have length 2
				#skipped += 1
			new_entry = [str(splitline[0])]  # chr
			new_entry.append(str(int(splitline[1]) + skipped))  # pos
			new_entry.append(str(splitline[2]))  # .
			new_entry.append(str(splitline[3][skipped:len(splitline[3]) - removed]))  # new ref
			new_entry.append(str(splitline[4][skipped:len(splitline[4]) - removed]))  # new alt
			new_entry.append("".join(splitline[5:]))  # anything after
			return ("\t".join(new_entry))
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
						return ("#error: " + "\t".join(splitline))
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
						return ("#error: " + "\t".join(splitline))
					break
			new_entry = [str(splitline[0])]  # chr
			new_entry.append(str(int(splitline[1]) + skipped))  # pos
			new_entry.append(str(splitline[2]))  # .
			new_entry.append(str(splitline[3][skipped:len(splitline[3]) - removed]))  # new ref
			new_entry.append(str(splitline[4][skipped:len(splitline[4]) - removed]))  # new alt
			new_entry.append("".join(splitline[5:]))  # anything after
			return ("\t".join(new_entry))
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
						return ("#error: " + "\t".join(splitline))
					break

			if len(nref) > 1 and len(nalt) > 1:
				# something like: Chr1	12998	.	TT	CC
				# but this (hopefully) never happens
				VCF_preproc_log.append("Warning: Preprocessing_original_vcf_file\n.")
				VCF_preproc_log.append(str("\t".join(splitline)))
				return ("#error: " + "\t".join(splitline))
			new_entry = [str(splitline[0])]  # chr
			new_entry.append(str(int(splitline[1]) + skipped))  # pos
			new_entry.append(str(splitline[2]))  # .
			new_entry.append(str(splitline[3][skipped:]))  # new ref
			new_entry.append(str(splitline[4][skipped:]))  # new alt
			new_entry.append("".join(splitline[5:]))  # anything after
			return ("\t".join(new_entry))
		else:
			print("Warning (critical?): Preprocessing_original_vcf_file (logfile)")
			VCF_preproc_log.append("Warning: Preprocessing_original_vcf_file\n.")
			VCF_preproc_log.append(str("\t".join(splitline)))
			return("#error: " + "\t".join(splitline))


	ori_vcf_file = open(ori_vcf,"r")
	print("Preprocessing")
	starttime = datetime.now()

	line = ori_vcf_file.readline()
	lines1 = []
	lines2 = []
	formatMultiallelVariantWarning = True
	while line:
		if line.startswith("#"):
			line = ori_vcf_file.readline()
			continue

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
			erg = formatMultiallelVariant(new_entry)
			if not erg.startswith("#"):
				lines1.append(erg)
			elif formatMultiallelVariantWarning:
				formatMultiallelVariantWarning = False
				print("Warning: Preprocessing_original_vcf_file (logfile)")

			#second allele
			new_entry = splitline[0:4]
			new_entry.append(split_allele[1])
			new_entry.append("\t".join(splitline[5:]))
			erg = formatMultiallelVariant(new_entry)
			if not erg.startswith("#"):
				lines2.append(erg)
			elif formatMultiallelVariantWarning:
				formatMultiallelVariantWarning = False
				print("Warning: Preprocessing_original_vcf_file (logfile)")
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


	""" old and slow
	warning = True
	for i,line in enumerate(lines1):
		spline = line.split("\t")
		if i == 0 or len(spline) == 1:
			continue
		if spline[1] == lines1[i-1].split("\t")[1]:
			if warning:
				print("VCF-Preprocessing Warning for File 1 in logfile.")
				VCF_preproc_log_in_short.append("VCF-Preprocessing Warning for File 1:\n")
				warning = False
			#print(str(lines1[i - 1]))
			#print("Removed: " + str(line))
			VCF_preproc_log_in_short.append("Take:\t" + str("\t".join(lines1[i - 1].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append(str(lines1[i - 1]))
			VCF_preproc_log_in_short.append("Rem.:\t" + str("\t".join(line.split("\t")[0:5])) + "\n")
			VCF_preproc_log.append("Removed: " + str(line))

			lines1.remove(line)
	"""
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
				VCF_preproc_log_in_short.append("VCF-Preprocessing Warning for File 1:\n")
				warning = False
			VCF_preproc_log_in_short.append("Take:\t" + str("\t".join(lines1[i - 1].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append(str(lines1[i - 1]))
			VCF_preproc_log_in_short.append("Rem.:\t" + str("\t".join(lines1[i].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append("Removed: " + str(lines1[i]))
			max -= 1
			del lines1[i]
			continue
		i += 1

	""" old and slow
	warning = True
	for i, line in enumerate(lines2):
		spline = line.split("\t")
		if i == 0 or len(spline) == 1:
			continue
		if spline[1] == lines2[i - 1].split("\t")[1]:
			if warning:
				print("VCF-Preprocessing Warning for File 2 in logfile.")
				VCF_preproc_log_in_short.append("VCF-Preprocessing Warning for File 2:\n")
				warning = False
			#print(str(lines2[i - 1]))
			#print("Removed:" + str(line))
			VCF_preproc_log_in_short.append("Take:\t" + str("\t".join(lines2[i - 1].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append(str(lines2[i - 1]))
			VCF_preproc_log_in_short.append("Rem.:\t" + str("\t".join(line.split("\t")[0:5])) + "\n")
			VCF_preproc_log.append("Removed:" + str(line))
			lines2.remove(line)
	"""
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
				VCF_preproc_log_in_short.append("VCF-Preprocessing Warning for File 2:\n")
				warning = False
			VCF_preproc_log_in_short.append("Take:\t" + str("\t".join(lines2[i - 1].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append(str(lines2[i - 1]))
			VCF_preproc_log_in_short.append("Rem.:\t" + str("\t".join(lines2[i].split("\t")[0:5])) + "\n")
			VCF_preproc_log.append("Removed: " + str(lines2[i]))
			max -= 1
			del lines2[i]
			continue
		i += 1

	def remove_shitty_data(lines:list)->list:
		i = 0
		removedEntries = 0
		max = len(lines) - 1
		pos = 0
		oldline = ""
		while i <= max:
			line = lines[i]
			spline = line.split("\t")
			if int(spline[1]) > pos:
				oldline = line
				pos = int(spline[1])
				pos += -1 + len(spline[3]) #+ length of ref, for deletions
				i +=1
			elif spline[0] == oldline.split("\t")[0]:
				#print("Removed:\n" + line + "Because of:\n" + oldline)

				VCF_preproc_log_in_short.append("Rem.:\t" + str("\t".join(lines[i].split("\t")[0:5])) + "\n")
				VCF_preproc_log.append("Removed: " + str(lines[i]))
				VCF_preproc_log_in_short.append("Because:\t" + str("\t".join(oldline.split("\t")[0:5])) + "\n")
				VCF_preproc_log.append(oldline)
				del lines[i]
				removedEntries +=1
				max -= 1
			else:
				oldline = line
				pos = int(spline[1])
				pos += -1 + len(spline[3])  # + length of ref, for deletions
				i += 1
		print(str(removedEntries) + " lines removed, because of conflicting data (see log).")
		return lines

	VCF_preproc_log_in_short.append("Data Conflicts in File 1:\n")
	lines1 = remove_shitty_data(lines1)
	print("Remove Shitty data 1 done: "+ str(datetime.now() - starttime))
	VCF_preproc_log_in_short.append("Data Conflicts in File 2:\n")
	lines2 = remove_shitty_data(lines2)
	print("Remove Shitty data 2 done: " + str(datetime.now() - starttime))

	nvcf_1 = open (new1_vcf,"w")
	nvcf_1.write("".join(lines1))
	nvcf_1.close()

	nvcf_2 = open (new2_vcf, "w")
	nvcf_2.write("".join(lines2))
	nvcf_2.close()
	print("preprocessing done: " + str(datetime.now() - starttime))

	print("Start writing log.")
	logfile = open(outpath + "logfile.txt", "w")
	logfile.write("".join(VCF_preproc_log_in_short) + "\n" + "".join(VCF_preproc_log))
	logfile.close()
	print("Colliding variants (half of them removed): " +str(len(VCF_preproc_log_in_short) - 2))
	print("Writing log done: " + str(datetime.now() - starttime))


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
	Preprocessing_original_vcf_file(invcf, new1_vcf, new2_vcf, outpath)

# python3 navip.py --mode pre --invcf /grp/gf/Alle_temp_Ordner/datenaustausch_jan_sarah/merged_variants.g.vcf --outpath /prj/gf-arabseq/project_VariantAnnotation/navip_bugsearch/

#vcf_preprocessing("/grp/gf/Alle_temp_Ordner/datenaustausch_jan_sarah/merged_variants.g.vcf", "/prj/gf-arabseq/project_VariantAnnotation/navip_bugsearch/")



#if __name__ == '__main__':
#	vcf_preprocessing("/grp/gf/Alle_temp_Ordner/datenaustausch_jan_sarah/merged_variants.g.vcf",
#					  "/prj/gf-arabseq/project_VariantAnnotation/navip_bugsearch/")




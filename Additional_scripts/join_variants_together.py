
def join_variants_together(input_vcf_file:str, output_vcf_file:str):
	"""
	It is using the ID from ID_Handler to merge multiple variant lines. It is expecting, that the only different information
	in this lines are the navip tags. So everything else is ignored and will be replaced by the first line.
	:param input_vcf_file: VCF file with the NAVID tag.
	:param output_vcf_file: Merged VCF file, will be much smaller.
	:return:
	"""
	def join_them(all_of_them:dict)-> list:
		output_list = []
		for variant_id_list in all_of_them.values(): # dict with a list for every id and the values are these lists
			# in variant_id_list here we have a list of variants, with the same id
			# now all infos from the different transcripts will join under one entry(line)
			new_line = variant_id_list[0].split("\t")
			nav1_list = []
			nav2_list = []
			for other_variants in variant_id_list[1:]:
				more_split_lines = other_variants.split("\t")[7]
				for info in more_split_lines.split(";"):
					if info.startswith("NAV1="):
						nav1_list.append(info[5:])
					elif info.startswith("NAV2="):
						nav2_list.append(info[5:])
			create_new_info = []
			for info in new_line[7].split(";"):
				if info.startswith("NAV1="):
					create_new_info.append(info + "@".join(nav1_list))
				elif info.startswith("NAV2="):
					create_new_info.append(info + "@".join(nav2_list))
				else:
					create_new_info.append(info)
			new_line[7] = ";".join(create_new_info)
			output_list.append("\t".join(new_line))
		return output_list


	old_vcf = open(input_vcf_file, 'r')
	line = old_vcf.readline()
	current_id = -1
	ID_collection = {}

	new_vcf = open(output_vcf_file, 'w')
	chromflag = ""

	while line:
		if line.startswith('#'):
			new_vcf.write(line)
			line = old_vcf.readline()
			continue
		elif len(line) <= 4:
			new_vcf.write(line)
			line = old_vcf.readline()
			continue

		spline = line.split("\t")
		if chromflag == "":
			chromflag = spline[0]
			ID_collection[chromflag] = {}
		elif chromflag != spline[0] and spline[0] != "" and spline[0] != " ":
			#new chrom
			joined_variants_list = join_them(ID_collection[chromflag])
			joined_variants_list = sorted(joined_variants_list, key=lambda data_line: (int(data_line.split("\t")[1])))
			new_vcf.write("".join(joined_variants_list))
			joined_variants_list = []

			ID_collection[chromflag] = {}
			chromflag = spline[0]
			ID_collection[chromflag] = {}


		infoline = spline[7].split(";")
		for info in infoline:
			if info.startswith("NAVID="):
				try:
					ID_collection[chromflag][info[6:]].append(line)
				except KeyError:
					ID_collection[chromflag][info[6:]] = [line]

		line = old_vcf.readline()

	#last chrom after while, if not triggered
	if ID_collection[chromflag] != {}:
		joined_variants_list = join_them(ID_collection[chromflag])
		joined_variants_list = sorted(joined_variants_list, key=lambda data_line: (int(data_line.split("\t")[1])))
		new_vcf.write("".join(joined_variants_list))
		joined_variants_list = []

		ID_collection[chromflag] = {}




	old_vcf.close()
	new_vcf.close()
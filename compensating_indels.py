import Transcript

def find_all_cindels(vcf_file_link:str, mod_or_not: bool, outputfolder: str, max_bp_range: int):

	vcf_file = open(vcf_file_link, 'r')
	vcf_data = vcf_file.readlines()
	vcf_file.close()

	transcript_indels_dict = {}
	transcript_direction_dict = {}

	for line in vcf_data:
		if line.startswith('#'):
			continue
		spline = line.split('\t')
		infoline = spline[7]
		for info in infoline.split(";"):
			if info.startswith('NAV1'):
				# transcript	direction
				#NAV1=AT1G76520.2|FOR|SUB,Amino acid change|NONE|ATT/0|i|583|GTT/0|v|583;
				transcript_as_key = str(info.split("|")[0].split('=')[1])
				effect_annotation = str(info.split("|")[2])
				new_cds_position = int(info.split("|")[9])

				if Transcript.TranscriptEnum.FRAMESHIFT_2_DEL.value in effect_annotation:
					stuff = (new_cds_position,-2, str(spline[0]), int(spline[1]))
				elif Transcript.TranscriptEnum.FRAMESHIFT_2.value in effect_annotation:
					stuff = (new_cds_position, 2, str(spline[0]), int(spline[1]))
				elif Transcript.TranscriptEnum.FRAMESHIFT_1_DEL.value in effect_annotation:
					stuff = (new_cds_position, -1, str(spline[0]), int(spline[1]))
				elif Transcript.TranscriptEnum.FRAMESHIFT_1.value in effect_annotation:
					stuff = (new_cds_position, 1, str(spline[0]), int(spline[1]))
				elif Transcript.TranscriptEnum.FRAMESHIFT.value in effect_annotation:
					print('meh')
				else:continue
				transcript_direction_dict[transcript_as_key] = info.split("|")[1]
				if transcript_as_key in transcript_indels_dict:
					transcript_indels_dict[transcript_as_key].append(stuff)

				else:
					transcript_indels_dict[transcript_as_key] =  [stuff]
			else:
				continue


	possible_neutralizing_indel_transcripts = {}
	for transcripts in transcript_indels_dict.keys():
		t_entrys = transcript_indels_dict[transcripts]
		if len(t_entrys) > 1:
			shift =0
			sub_indel_list = []
			#testing rev-direction correction
			if transcript_direction_dict[transcripts] == "REV":
				t_entrys = t_entrys[::-1]
			for entry in t_entrys:
				shift += entry[1] # if this reaches 0 (or mod 3 == 0), the frameshift is (maybe, no stop detection here) compensated
				sub_indel_list.append(entry) # all entrys needed for one additional compensated fs (could be more than one)
				if mod_or_not:
					if shift % 3 == 0: # compensation
						for sub_indels in sub_indel_list: # all fs entrys
							if transcripts in possible_neutralizing_indel_transcripts:
								possible_neutralizing_indel_transcripts[transcripts].append(sub_indels)
							else:
								possible_neutralizing_indel_transcripts[transcripts] = [sub_indels]
						sub_indel_list = [] # search for next compensating fs
				if not mod_or_not:
					if shift == 0: # compensation
						for sub_indels in sub_indel_list: # all fs entrys
							if transcripts in possible_neutralizing_indel_transcripts:
								possible_neutralizing_indel_transcripts[transcripts].append(sub_indels)
							else:
								possible_neutralizing_indel_transcripts[transcripts] = [sub_indels]
						sub_indel_list = [] # search for next compensating fs


	#max_bp_range = 5000
	last_dict = {}
	table_output = []
	unique_table_output = []
	for maxbp_between_compensation_fs in range(1, max_bp_range+1):
		i_dict = {}
		unique_dict = {}
		for transcripts in possible_neutralizing_indel_transcripts.keys():
			#possible_neutralizing_indel_transcripts[transcripts] lloks like this:
			#<class 'list'>: [(435, 1, 'Chr1\t1463756'), (495, -1, 'Chr1\t1463815')]
			range_check = 0
			too_long = False # True, if length between indels is too large
			for i,entrys in enumerate(possible_neutralizing_indel_transcripts[transcripts][1:]):
				# i should be one less, then the current position of the list, because list starts at 1
				last_entry = possible_neutralizing_indel_transcripts[transcripts][i]
				range_check += abs(entrys[0] - last_entry[0])
				if range_check >= maxbp_between_compensation_fs:
					too_long = True
					break
			if too_long:
				continue
			unique_set = set()
			for entry in possible_neutralizing_indel_transcripts[transcripts]:
				# possible_neutralizing_indel_transcripts[transcripts] looks like this:
				# <class 'list'>: [(435, 1, 'Chr1\t1463756'), (495, -1, 'Chr1\t1463815')]
				# entry == (435, 1, 'Chr1', 1463756)
				unique_set.add(entry)
			unique_dict[transcripts.split('.')[0]] = list(unique_set)
			unique_set = set()
			if len(possible_neutralizing_indel_transcripts[transcripts]) in i_dict:
				i_dict[len(possible_neutralizing_indel_transcripts[transcripts])] += 1
			else:
				i_dict[len(possible_neutralizing_indel_transcripts[transcripts])] = 1
		if last_dict.items() == i_dict.items():
			continue
		else:
			last_dict = i_dict

		list_of_keys = i_dict.keys()
		list_of_keys = sorted(list_of_keys)
		number_of_zeros = 1
		table_output.append(str(maxbp_between_compensation_fs))
		for i in list_of_keys:
			#print("i:" + str(i) +"\tkeys:" + str(list_of_keys))
			number_of_zeros += 1
			if i != number_of_zeros:
				for zeros in range(0, abs(i-number_of_zeros)):
					#its for jumping from... maybe 5 introns to 7 introns without 6 introns in any transcript
					table_output.append("\t0")
					number_of_zeros += 1
			table_output.append("\t" + str(i_dict[i]))
		table_output.append("\n")

		count_indel_in_transcripts_unique_dict = {}
		for indellist in unique_dict.values():
			if len(indellist) in count_indel_in_transcripts_unique_dict:
				count_indel_in_transcripts_unique_dict[len(indellist)] += 1
			else:
				count_indel_in_transcripts_unique_dict[len(indellist)] = 1


		list_of_keys = count_indel_in_transcripts_unique_dict.keys()
		list_of_keys = sorted(list_of_keys)
		number_of_zeros = 1
		unique_table_output.append(str(maxbp_between_compensation_fs))
		for i in list_of_keys:
			#print("i:" + str(i) +"\tkeys:" + str(list_of_keys))
			number_of_zeros += 1
			if i != number_of_zeros:
				for zeros in range(0, abs(i-number_of_zeros)):
					# its for jumping from... maybe 5 introns to 7 introns without 6 introns in any transcript
					unique_table_output.append("\t0")
					number_of_zeros += 1
			unique_table_output.append("\t" + str(count_indel_in_transcripts_unique_dict[i]))
		unique_table_output.append("\n")






	max_key = max(last_dict.keys())
	key_output = []
	for key in range(2,max_key +1):
		key_output.append("\tX=" + str(key))

	formating_list1 = "".join(table_output)
	formating_list = formating_list1.split("\n") # every line should be here in a big list
	if formating_list[len(formating_list)-1] == "":
		formating_list.pop(len(formating_list)-1)
	maxlength = len(formating_list[len(formating_list)-1].split('\t'))# length of last entry
	new_output = []
	for line in formating_list:
		# line is a string for one line
		new_output.append(line)
		diff = abs(maxlength - len(line.split("\t"))) # if there are not enough zeros, there should be a difference
		if diff > 0:
			for i in range(0, diff):
				new_output.append("\t0")
		new_output.append("\n")



	mode_name = ""
	if mod_or_not:
		mode_name = "compInDels"
	else:
		mode_name = "compAA"

	table_outputname = "transcripts_isoform_" + mode_name +"_"+ str(max_bp_range) +'bpr' + '.txt'
	description = "#transcripts with X indels (cumultative) - " + str(mode_name) + ' - to max distance:' + str(max_bp_range) + "\n"
	description+= "#max_bp_between_compensation_fs" + "".join(key_output) + "\n"
	#print(description + "".join(table_output))

	table_output_file = open(outputfolder + table_outputname, 'w')
	table_output_file.write(description + "".join(new_output))
	table_output_file.close()

	formating_list1 = "".join(unique_table_output)
	formating_list = formating_list1.split("\n")  # every line should be here in a big list
	if formating_list[len(formating_list) - 1] == "":
		formating_list.pop(len(formating_list) - 1)
	maxlength = len(formating_list[len(formating_list) - 1].split('\t'))  # length of last entry
	new_output = []
	for line in formating_list:
		# line is a string for one line
		new_output.append(line)
		diff = abs(maxlength - len(line.split("\t")))  # if there are not enough zeros, there should be a difference
		if diff > 0:
			for i in range(0, diff):
				new_output.append("\t0")
		new_output.append("\n")

	unique_table_outputname = "transcripts_unique_" + mode_name +"_"+ str(max_bp_range) +'bpr' + '.txt'
	description2 = "#transcripts with X indels (cumultative) - unqiue dataset - " + str(mode_name) + ' - to max distance:' + str(max_bp_range) + "\n"
	description2 += "#max_bp_between_compensation_fs" + "".join(key_output) + "\n"
	#print(description2 + "".join(unique_table_output))

	table_output_file = open(outputfolder + unique_table_outputname, 'w')
	table_output_file.write(description2 + "".join(new_output))
	table_output_file.close()













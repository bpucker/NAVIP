import Transcript
import matplotlib.pyplot as plt
from matplotlib import rc

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
		isoform_dict = {}
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
			if len(possible_neutralizing_indel_transcripts[transcripts]) in isoform_dict:
				isoform_dict[len(possible_neutralizing_indel_transcripts[transcripts])] += 1
			else:
				isoform_dict[len(possible_neutralizing_indel_transcripts[transcripts])] = 1
		if last_dict.items() == isoform_dict.items():
			continue
		else:
			last_dict = isoform_dict

		list_of_keys = isoform_dict.keys()
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
			table_output.append("\t" + str(isoform_dict[i]))
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

def do_magic_plotting(data:list, outputfolder:str, orig_outputname:str, formats:str):
	# y-axis in bold
	rc('font', weight='bold')
	#'1\t30\t0\t0\t0'
	# Values of each group
	max_coloumns_so_bars = len(data[0].split("\t"))
	all_bars = {}

	for i in range(1,max_coloumns_so_bars):
		all_bars[i] = [] #initialization of bars

	x_axsis = []
	for dataline in data[0:30]:
		dataline = dataline.split("\t")
		x_axsis.append(int(dataline[0]))  # bpr
		for i in range(1, max_coloumns_so_bars):
			all_bars[i].append(int(dataline[i])) # involved transcripts with cindel_events
	###create all bars
	labels = []
	for i in range(1, max_coloumns_so_bars):
		plt.bar(x_axsis,all_bars[i],edgecolor='white', width=1 )
		labels.append(str(i+1) + " InDels")

	# for my data a bar would be the entire second, third [...] colomn, not row.
	#bars1 = [12, 28, 1, 8, 22]
	#bars2 = [28, 7, 16, 4, 10]
	#bars3 = [25, 3, 23, 25, 17]

	# Heights of bars1 + bars2
	#bars = np.add(bars1, bars2).tolist()

	# The position of the bars on the x-axis
	#r = [0, 1, 2, 3, 4]

	# Names of group and bar width
	#names = ['A', 'B', 'C', 'D', 'E']
	#barWidth = 1

	# Create brown bars
	#plt.bar(r, bars1, color='#7f6d5f', edgecolor='white', width=barWidth)
	# Create green bars (middle), on top of the firs ones
	#plt.bar(r, bars2, bottom=bars1, color='#557f2d', edgecolor='white', width=barWidth)
	# Create green bars (top)
	#plt.bar(r, bars3, bottom=bars, color='#2d7f5e', edgecolor='white', width=barWidth)

	# Custom X axis
	#plt.margins(x=0)
	plt.legend(labels)
	plt.xticks(x_axsis,x_axsis, fontweight='bold')
	plt.xscale('linear')
	plt.xlabel('distance between cInDels')
	plt.ylabel('number of transcripts with cInDel events')
	#plt.x
	#plt.xlabel("group")

	# Show graphic
	#plt.show()
	#"One of the file extensions supported by the active backend. Most backends support png, pdf, ps, eps and svg."
	probably_supported_formats = ["png","pdf","ps","eps","svg"]
	for picture_format in formats.split(','):
		if picture_format.lower() in probably_supported_formats:
			plt.savefig(outputfolder + orig_outputname.split(".")[0] + "." + picture_format.lower())
	plt.clf()

def find_all_cindels_v2(navip_vcf_file_link: str, mod_or_not: bool, outputfolder: str, formats:str):

		vcf_file = open(navip_vcf_file_link, 'r')
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
					# NAV1=AT1G76520.2|FOR|SUB,Amino acid change|NONE|ATT/0|i|583|GTT/0|v|583;
					transcript_as_key = str(info.split("|")[0].split('=')[1])
					effect_annotation = str(info.split("|")[2])
					new_cds_position = int(info.split("|")[9])

					if Transcript.TranscriptEnum.FRAMESHIFT_2_DEL.value in effect_annotation:
						stuff = (new_cds_position, -2, str(spline[0]), int(spline[1]))
					elif Transcript.TranscriptEnum.FRAMESHIFT_2.value in effect_annotation:
						stuff = (new_cds_position, 2, str(spline[0]), int(spline[1]))
					elif Transcript.TranscriptEnum.FRAMESHIFT_1_DEL.value in effect_annotation:
						stuff = (new_cds_position, -1, str(spline[0]), int(spline[1]))
					elif Transcript.TranscriptEnum.FRAMESHIFT_1.value in effect_annotation:
						stuff = (new_cds_position, 1, str(spline[0]), int(spline[1]))
					elif Transcript.TranscriptEnum.FRAMESHIFT.value in effect_annotation:
						print('There is somehow not the correct frameshift annotation in NAV1.')
					else:
						continue
					transcript_direction_dict[transcript_as_key] = info.split("|")[1]
					if transcript_as_key in transcript_indels_dict:
						transcript_indels_dict[transcript_as_key].append(stuff)
					else:
						transcript_indels_dict[transcript_as_key] = [stuff]

				else:
					continue

		possible_neutralizing_indel_transcripts_experimental = {}
		possible_neutralizing_indel_transcripts_experimental_unique = {}
		for transcripts in transcript_indels_dict.keys():
			t_entrys = transcript_indels_dict[transcripts]
			if len(t_entrys) > 1:
				shift = 0
				sub_indel_list = []
				# testing rev-direction correction
				if transcript_direction_dict[transcripts] == "REV":
					# because the list ist sorted after chr_pos in ascending order
					# for rev is the descending order needed
					t_entrys = t_entrys[::-1]
				for entry in t_entrys:
					shift += entry[1]  # if this reaches 0 (or mod 3 == 0), the frameshift is (maybe, no stop detection here) compensated
					sub_indel_list.append(entry)  # all entrys needed for one additional compensated fs (could be more than one)
					if mod_or_not:
						if shift % 3 == 0:  # compensation
							#sub_indel_list: <class 'list'>: [(16, 2, 'Chr1', 354830), (43, 1, 'Chr1', 354811)]
							max_dist = 0
							for i, entrys in enumerate(sub_indel_list[1:]):
								# i should be one less, then the current position of the list, because list starts at 1
								last_entry = sub_indel_list[i]
								max_dist = max((abs(entrys[0] - last_entry[0]), max_dist))
							# new_cindel : <class 'tuple'>: ('AT1G02020.1', 2, [(16, 2, 'Chr1', 354830), (43, 1, 'Chr1', 354811)], 27)
							new_cindel = (transcripts,len(sub_indel_list),sub_indel_list,max_dist)
							new_cindel_unique = (transcripts.split(".")[0],len(sub_indel_list),sub_indel_list,max_dist)

							if max_dist in possible_neutralizing_indel_transcripts_experimental:
								if len(sub_indel_list) in possible_neutralizing_indel_transcripts_experimental[max_dist]:
									possible_neutralizing_indel_transcripts_experimental[max_dist][len(sub_indel_list)].append(new_cindel)
									add_this_without_bug = str(new_cindel_unique[0]) +"\t" + str(new_cindel_unique[1]) +"\t" + str(new_cindel_unique[2]) +"\t" + str(new_cindel_unique[3])
									possible_neutralizing_indel_transcripts_experimental_unique[max_dist][len(sub_indel_list)].add(add_this_without_bug)
								else:
									possible_neutralizing_indel_transcripts_experimental[max_dist][len(sub_indel_list)] = [new_cindel]
									add_this_without_bug = str(new_cindel_unique[0]) + "\t" + str(new_cindel_unique[1]) + "\t" + str(new_cindel_unique[2]) + "\t" + str(new_cindel_unique[3])
									possible_neutralizing_indel_transcripts_experimental_unique[max_dist][len(sub_indel_list)] = set()
									possible_neutralizing_indel_transcripts_experimental_unique[max_dist][len(sub_indel_list)].add(add_this_without_bug)
							else:
								possible_neutralizing_indel_transcripts_experimental[max_dist] = {len(sub_indel_list):[new_cindel]}
								add_this_without_bug = str(new_cindel_unique[0]) + "\t" + str(new_cindel_unique[1]) + "\t" + str(new_cindel_unique[2]) + "\t" + str(new_cindel_unique[3])
								possible_neutralizing_indel_transcripts_experimental_unique[max_dist] = {len(sub_indel_list):set()}
								possible_neutralizing_indel_transcripts_experimental_unique[max_dist][len(sub_indel_list)].add(add_this_without_bug)
							sub_indel_list = [] # clear all entrys, because cindel combo is finished
					if not mod_or_not:
						if shift == 0:  # compensation
							max_dist = 0
							for i, entrys in enumerate(sub_indel_list[1:]):
								# i should be one less, then the current position of the list, because list starts at 1
								last_entry = sub_indel_list[i]
								max_dist = max((abs(entrys[0] - last_entry[0]), max_dist))

							new_cindel = (transcripts, len(sub_indel_list), sub_indel_list, max_dist)
							new_cindel_unique = (
							transcripts.split(".")[0], len(sub_indel_list), sub_indel_list, max_dist)

							if max_dist in possible_neutralizing_indel_transcripts_experimental:
								if len(sub_indel_list) in possible_neutralizing_indel_transcripts_experimental[
									max_dist]:
									possible_neutralizing_indel_transcripts_experimental[max_dist][
										len(sub_indel_list)].append(new_cindel)
									add_this_without_bug = str(new_cindel_unique[0]) + "\t" + str(
										new_cindel_unique[1]) + "\t" + str(new_cindel_unique[2]) + "\t" + str(
										new_cindel_unique[3])
									possible_neutralizing_indel_transcripts_experimental_unique[max_dist][
										len(sub_indel_list)].add(add_this_without_bug)
								else:
									possible_neutralizing_indel_transcripts_experimental[max_dist][
										len(sub_indel_list)] = [new_cindel]
									add_this_without_bug = str(new_cindel_unique[0]) + "\t" + str(
										new_cindel_unique[1]) + "\t" + str(new_cindel_unique[2]) + "\t" + str(
										new_cindel_unique[3])
									possible_neutralizing_indel_transcripts_experimental_unique[max_dist][
										len(sub_indel_list)] = set()
									possible_neutralizing_indel_transcripts_experimental_unique[max_dist][
										len(sub_indel_list)].add(add_this_without_bug)
							else:
								possible_neutralizing_indel_transcripts_experimental[max_dist] = {
									len(sub_indel_list): [new_cindel]}
								add_this_without_bug = str(new_cindel_unique[0]) + "\t" + str(
									new_cindel_unique[1]) + "\t" + str(new_cindel_unique[2]) + "\t" + str(
									new_cindel_unique[3])
								possible_neutralizing_indel_transcripts_experimental_unique[max_dist] = {
									len(sub_indel_list): set()}
								possible_neutralizing_indel_transcripts_experimental_unique[max_dist][
									len(sub_indel_list)].add(add_this_without_bug)
							sub_indel_list = []  # clear all entrys, because cindel combo is finished

		mode_name = ""
		if mod_or_not:
			mode_name = "cInDels"
		else:
			mode_name = "only_additive_cInDels"


		#################################
		##### experimental area     #####
		##### proceed with caution  #####
		#################################
		##### isoform management     ####
		normal_output = []
		normal_output_with_zeros = []
		involved_tid_list = []

		highest_cindel_combo = 0
		key_output = []
		for max_base_pairs in possible_neutralizing_indel_transcripts_experimental.keys():
			meh = [highest_cindel_combo]
			meh.extend(possible_neutralizing_indel_transcripts_experimental[max_base_pairs]) # returns void....
			highest_cindel_combo = max(meh)
		for key in range(2, highest_cindel_combo + 1):
			key_output.append("\tX=" + str(key))

		for max_base_pairs in sorted(list(possible_neutralizing_indel_transcripts_experimental.keys())):
			normal_output.append(str(max_base_pairs))

			i = 1
			for cindelsLength in sorted(list(possible_neutralizing_indel_transcripts_experimental[max_base_pairs])):
				i +=1

				while i != cindelsLength:
					normal_output.append("\t0")
					i +=1
					if i > cindelsLength+1:
						print("... ?") # should not happen, happened once >.<
				temp_list = []
				normal_output.append("\t" + str(len(possible_neutralizing_indel_transcripts_experimental[max_base_pairs][cindelsLength])))
				for cindel in possible_neutralizing_indel_transcripts_experimental[max_base_pairs][cindelsLength]:
					temp_list.append(str(cindel[0]))
				involved_tid_list.append((max_base_pairs, cindelsLength,temp_list))
			normal_output.append("\n")
		involved_tid_list = sorted(involved_tid_list)
		for lines in "".join(normal_output).split("\n"):
			if lines == "":
				continue
			dif = abs(highest_cindel_combo - len(lines.split("\t")))
			if dif > 0:
				for i in range(0, dif):
					lines += "\t0"
			normal_output_with_zeros.append(lines)
			#highest_cindel_combo

		table_outputname = "transcripts_isoform_" + mode_name + '.txt'
		description = "#transcripts with X indels - " + str(mode_name) + "\n"
		description += "#bp_between_compensation_fs" + "".join(key_output) + "\n"

		table_output_file = open(outputfolder + table_outputname, 'w')
		table_output_file.write(description + "\n".join(normal_output_with_zeros))
		table_output_file.close()

		do_magic_plotting(normal_output_with_zeros,outputfolder, table_outputname,formats )

		### write all tids and stuff
		involved_tid_list_output = []
		for mbp_clength_temp_entry in involved_tid_list:
			involved_tid_list_output.append(">" + str(mbp_clength_temp_entry[0]) + " ")
			involved_tid_list_output.append(str(mbp_clength_temp_entry[1]) + "\n")
			i = 0
			for tid in mbp_clength_temp_entry[2]:
				if i % 10 == 0:
					involved_tid_list_output.append(str(tid))
				else:
					involved_tid_list_output.append("," + str(tid))
				if i % 10 == 9:
					involved_tid_list_output.append("\n")
				i += 1
			if i % 10 != 0:
				involved_tid_list_output.append("\n")

		table_outputname = "transcripts_isoform_" + mode_name + '_TIDs' + '.txt'
		description = "#Header: > <bp> <quantity of involved InDels for one compensating InDel (cInDel) event>  \n"
		description += "#Data line: <tid|cInDel-event 1| cInDel-event2| ....> \n"
		description += "#cInDel-event:<Chr>,<Pos>,<CDS-Pos>,<frameshift-value>;<Chr>,<Pos>,<CDS-Pos>,<frameshift-value> [...]\n"
		description += "#frameshift-value: integer values: deletion: -1,-2 bases; insertion: 1,2 bases\n"

		table_output_file = open(outputfolder + table_outputname , 'w')
		table_output_file.write(description + "".join(involved_tid_list_output))
		table_output_file.close()


		#### unique management ####

		normal_output = []
		normal_output_with_zeros = []

		highest_cindel_combo = 0
		key_output = []
		for max_base_pairs in possible_neutralizing_indel_transcripts_experimental_unique.keys():
			meh = [highest_cindel_combo]
			meh.extend(possible_neutralizing_indel_transcripts_experimental_unique[max_base_pairs])  # returns void....
			highest_cindel_combo = max(meh)
		for key in range(2, highest_cindel_combo + 1):
			key_output.append("\tX=" + str(key))

		for max_base_pairs in sorted(list(possible_neutralizing_indel_transcripts_experimental_unique.keys())):
			normal_output.append(str(max_base_pairs))
			i = 1
			for cindelsLength in sorted(list(possible_neutralizing_indel_transcripts_experimental_unique[max_base_pairs])):
				i += 1

				while i != cindelsLength:
					normal_output.append("\t0")
					i += 1
					if i > cindelsLength + 1:
						print("... ?")  # should not happen, happened once >.<
				normal_output.append("\t" + str(
					len(possible_neutralizing_indel_transcripts_experimental_unique[max_base_pairs][cindelsLength])))
			normal_output.append("\n")
		for lines in "".join(normal_output).split("\n"):
			if lines == "":
				continue
			dif = abs(highest_cindel_combo - len(lines.split("\t")))
			if dif > 0:
				for i in range(0, dif):
					lines += "\t0"
			normal_output_with_zeros.append(lines)
		# highest_cindel_combo

		table_outputname = "transcripts_unique_" + mode_name + '.txt'
		description = "#transcripts with X indels - " + str(mode_name) + "\n"
		description += "#bp_between_compensation_fs" + "".join(key_output) + "\n"

		table_output_file = open(outputfolder + table_outputname, 'w')
		table_output_file.write(description + "\n".join(normal_output_with_zeros))
		table_output_file.close()
		do_magic_plotting(normal_output_with_zeros, outputfolder, table_outputname,formats)

		### write all tids and stuff
		involved_tid_list_unique = []
		involved_tid_list_unique_set = set()
		for mbp_clength_temp_entry in involved_tid_list:
			for tid in mbp_clength_temp_entry[2]:
				involved_tid_list_unique_set.add(str(mbp_clength_temp_entry[0]) + ' ' + str(mbp_clength_temp_entry[1]) + ' ' + str(tid.split(".")[0]))
		involved_tid_list_unique = list(involved_tid_list_unique_set)
		involved_tid_list_unique = sorted(involved_tid_list_unique, key= lambda x: (int(x.split(" ")[0]),int(x.split(" ")[1]))) # [000]:'1 2 AT1G22060'

		involved_tid_list_output = []
		old_header = "> " + str(involved_tid_list_unique[0].split(" ")[0]) + " " + str(involved_tid_list_unique[0].split(" ")[1] + "\n")
		involved_tid_list_output.append(old_header)
		i = 0
		for mbp_clength_temp_entry in involved_tid_list_unique:
			tid = mbp_clength_temp_entry.split(" ")[2]
			new_header = "> " + str(mbp_clength_temp_entry.split(" ")[0]) + " " + str(mbp_clength_temp_entry.split(" ")[1] + "\n")
			if new_header != old_header: # new header == new tids
				if i % 10 != 0:
					involved_tid_list_output.append("\n")
				old_header = new_header
				involved_tid_list_output.append(new_header)
				i = 0
			if i % 10 == 0:
				involved_tid_list_output.append(str(tid))
			else:
				involved_tid_list_output.append("," + str(tid))
			if i % 10 == 9:
				involved_tid_list_output.append("\n")
			i += 1

		description = "#Header: > <bp> <quantity of involved InDels for one compensating InDel (cInDel) event>  \n"
		description += "#Data line: <tid|cInDel-event 1| cInDel-event2| ....> \n"
		description += "#cInDel-event:<Chr>,<Pos>,<CDS-Pos>,<frameshift-value>;<Chr>,<Pos>,<CDS-Pos>,<frameshift-value> [...]\n"
		description += "#raster-change-value: integer values: deletion: -1,-2 bases; insertion: 1,2 bases\n"

		table_outputname = "transcripts_unique_" + mode_name + '_TIDs' + '.txt'
		table_output_file = open(outputfolder + table_outputname , 'w')
		table_output_file.write(description + "".join(involved_tid_list_output))
		table_output_file.close()


		detailed_output_stuff_dict = {}
		output_sorted_by_tid_dict = {}
		bpr_quantity_pairs = [] # can be sorted and this way i have the needed entrys for everything
		for bpr in possible_neutralizing_indel_transcripts_experimental:
			for quantity_of_involved_indels in possible_neutralizing_indel_transcripts_experimental[bpr]:
				bpr_quantity_pairs.append((bpr,quantity_of_involved_indels))
				for cindel_event in possible_neutralizing_indel_transcripts_experimental[bpr][quantity_of_involved_indels]:
					#<class 'tuple'>: ('AT1G06220.1', 2, [(350, 2, 'Chr1', 1900873), (353, 1, 'Chr1', 1900874)], 3)
					"""
					### planned output ###
					#Header: > <bp> <quantity of involved InDels for one compensating InDel (cInDel) event>
					#Data line: <tid|cInDel-event 1| cInDel-event2| ....>
					#cInDel-event:<Chr>,<Pos>,<CDS-Pos>,<frameshift-value>;<Chr>,<Pos>,<CDS-Pos>,<frameshift-value> [...]
					> 3 2
					AT1G27170.1|Chr1,9436258,1442,-1;Chr1,9436262,1445,-2|Chr1,9436894,1994,2;Chr1,9436895,1997,-2
					AT1G27170.2|Chr1,9436258,1442,-1;Chr1,9436262,1445,-2|Chr1,9436894,1994,2;Chr1,9436895,1997,-2 
					"""

					tid = cindel_event[0]
					cindel_details = ""
					for data_tuple in cindel_event[2]:
						# <class 'tuple'>: ('AT1G06220.1', 2, [(350, 2, 'Chr1', 1900873), (353, 1, 'Chr1', 1900874)], 3)
						cindel_details += str(data_tuple[2]) + "," + str(data_tuple[3]) + "," + str(data_tuple[0]) + ","+ str(data_tuple[1]) + ";"
					#cindel_details += "|"
					if tid in output_sorted_by_tid_dict.keys():
						output_sorted_by_tid_dict[tid].append(cindel_details[0:len(cindel_details)-1])
					else:
						output_sorted_by_tid_dict[tid] = [cindel_details[0:len(cindel_details)-1]]

					if (bpr,quantity_of_involved_indels) in detailed_output_stuff_dict.keys():
						if tid in detailed_output_stuff_dict[(bpr,quantity_of_involved_indels)].keys():
							detailed_output_stuff_dict[(bpr, quantity_of_involved_indels)][tid].append(cindel_details)
						else:
							detailed_output_stuff_dict[(bpr, quantity_of_involved_indels)][tid] = [cindel_details]
					else:
						detailed_output_stuff_dict[(bpr, quantity_of_involved_indels)] = {tid: [cindel_details]}
		bpr_quantity_pairs = sorted(bpr_quantity_pairs)

		involved_tid_list_detailed_output = []
		for bpr,quantity in bpr_quantity_pairs:
			#print(bpr)
			#print(quantity)
			involved_tid_list_detailed_output.append(">" + str(bpr) +" " + str(quantity) + "\n")
			for tid in detailed_output_stuff_dict[(bpr,quantity)]:
				involved_tid_list_detailed_output.append(str(tid) + "|")
				for cindel_details_list in detailed_output_stuff_dict[(bpr,quantity)][tid]:
					cindel_events_list = cindel_details_list.split("|")
					i = len(detailed_output_stuff_dict[(bpr,quantity)][tid])
					for cindel_event in cindel_events_list:
						#cindel_event should be a string
						i -= 1
						involved_tid_list_detailed_output.append(cindel_event[0:len(cindel_event)-1])# -2 because of ignoring the last ";"
						if i != 0:
							involved_tid_list_detailed_output.append("|")
				involved_tid_list_detailed_output.append("\n")

				"""
				### planned output ###
				#Header: > <bp> <quantity of involved InDels for one compensating InDel (cInDel) event>
				#Data line: <tid|cInDel-event 1| cInDel-event2| ....>
				#cInDel-event:<Chr>,<Pos>,<CDS-Pos>,<frameshift-value>;<Chr>,<Pos>,<CDS-Pos>,<frameshift-value> [...]
				> 3 2
				AT1G27170.1|Chr1,9436258,1442,-1;Chr1,9436262,1445,-2|Chr1,9436894,1994,2;Chr1,9436895,1997,-2
				AT1G27170.2|Chr1,9436258,1442,-1;Chr1,9436262,1445,-2|Chr1,9436894,1994,2;Chr1,9436895,1997,-2 
				"""



		description = "#Header: > <bp> <quantity of involved InDels for one compensating InDel (cInDel) event>  \n"
		description += "#Data line: <tid|cInDel-event 1| cInDel-event2| ....> \n"
		description += "#cInDel-event:<Chr>,<Pos>,<CDS-Pos>,<frameshift-value>;<Chr>,<Pos>,<CDS-Pos>,<frameshift-value> [...]\n"
		description += "#raster-change-value: integer values: deletion: -1,-2 bases; insertion: 1,2 bases\n"

		table_outputname = "transcripts_unique_" + mode_name + '_TIDs_detailed' + '.txt'
		table_output_file = open(outputfolder + table_outputname, 'w')
		table_output_file.write(description + "".join(involved_tid_list_detailed_output))
		table_output_file.close()

		output_sorted_by_tid = []
		sort_this_tids = []
		unique_set = set()

		for tid in output_sorted_by_tid_dict.keys():
			sort_this_tids.append(tid)
			unique_set.add(tid.split(".")[0])
		sort_this_tids = sorted(sort_this_tids)
		counti_all = 0
		counti_unique_list = []
		for tid in sort_this_tids:
			entry_i = len(output_sorted_by_tid_dict[tid])
			if entry_i > 1:
				counti_all += 1
				asdasdasdasd = tid.split(".")[0]
				if asdasdasdasd not in counti_unique_list:
					counti_unique_list.append(asdasdasdasd)

			output_sorted_by_tid.append(str(tid) + "," + str(len(output_sorted_by_tid_dict[tid])) + "|")
			for value in sorted(output_sorted_by_tid_dict[tid], key= lambda entry_in_here: int(entry_in_here.split(",")[2])):
				#value = sorted(value, key= lambda entry: int(entry.split(",")[2]))
				entry_i -= 1
				indel_list = value.split(";")
				max_bpr = 0
				old_indel = indel_list[0]
				for indel in indel_list[1:]:
					max_bpr = max([abs(int(old_indel.split(",")[2]) - int(indel.split(",")[2])), max_bpr])
					old_indel = indel
				output_sorted_by_tid.append(str(max_bpr) + "," + str(len(indel_list)) +",")
				i = len(indel_list)
				for stuff in indel_list:
					i -=1
					#output_sorted_by_tid.append(str(stuff[0]) + str(stuff[1]) + str(stuff[2]) +str(stuff[3]) )
					output_sorted_by_tid.append(stuff)
					if i != 0:
						output_sorted_by_tid.append(";")
				if entry_i != 0:
					output_sorted_by_tid.append("|")
			output_sorted_by_tid.append("\n")





		description = ""
		description += "#Data line: <tid>,<quantity cindel-events>|<max_bpr>,<InDel_quantity>,cInDel-event 1|<bpr><InDel_quantity>,cInDel-event2| ....> \n"
		description += "#cInDel-event:<Chr>,<Pos>,<CDS-Pos>,<frameshift-value>;<Chr>,<Pos>,<CDS-Pos>,<frameshift-value> [...]\n"

		table_outputname = "transcripts_unique_" + mode_name + '_TIDs_detailed_sorted_by_TID' + '.txt'
		table_output_file = open(outputfolder + table_outputname, 'w')
		table_output_file.write(description + "".join(output_sorted_by_tid))
		table_output_file.close()

		print("### Overview ###")
		print("Number of transcripts (+isoforms) with at least one compensation InDel (cInDel) event: " + str(len(sort_this_tids)))
		print("Number of transcripts (unique) with at least one compensation InDel (cInDel) event: " +str(len(unique_set)))
		print("Number of transcripts (+isoforms) with more then one compensation InDel (cInDel) event: " + str(counti_all))
		print("Number of transcripts (unique) with more then one compensation InDel (cInDel) event: " +  str(len(counti_unique_list)))











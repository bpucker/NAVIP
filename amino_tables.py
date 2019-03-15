import re

def different_AA(entry: str) -> (str, str):
	"""
	:param entry: Entry betweeen pipes: |p.Thr56Ser|
	# p.Thr56Ser	p.Thr56Thr
	:return: (Amino_from , Amino_to)
	"""
	entry = entry.replace("p.", "")
	# spentry = re.split('[a-z]+(?:[0-9]])',entry)
	# spentry = re.split('[A-Z]{1}{a-z}{2}\d+[A-Z]{1}{a-z}{2}',entry)
	spentry = re.split('\d+', entry)
	# print(entry)
	# print(spentry)
	"""
		Pro119Asp
		['Pro', '', '', 'Asp']
	"""
	# spentry = re.split('\D', entry)
	# print(entry)
	# print(spentry)
	"""
		Pro119Asp
		['', '', '', '119', '', '', '']
	"""
	return spentry


aminodict = {}
aminodict["Ala"] = "A"
aminodict["A"] = "Ala"
aminodict["Cys"] = "C"
aminodict["C"] = "Cys"
aminodict["Asp"] = "D"
aminodict["D"] = "Asp"
aminodict["Glu"] = "E"
aminodict["E"] = "Glu"
aminodict["Phe"] = "F"
aminodict["F"] = "Phe"
aminodict["Gly"] = "G"
aminodict["G"] = "Gly"
aminodict["His"] = "H"
aminodict["H"] = "His"
aminodict["Ile"] = "I"
aminodict["I"] = "Ile"
aminodict["Lys"] = "K"
aminodict["K"] = "Lys"
aminodict["Leu"] = "L"
aminodict["L"]  = "Leu"
aminodict["Met"] = "M"
aminodict["M"] = "Met"
aminodict["Asn"] = "N"
aminodict["N"] = "Asn"
aminodict["Pro"] = "P"
aminodict["P"] = "Pro"
aminodict["Gln"] = "Q"
aminodict["Q"] = "Gln"
aminodict["Arg"] = "R"
aminodict["R"] = "Arg"
aminodict["Ser"] = "S"
aminodict["S"] = "Ser"
aminodict["Thr"] = "T"
aminodict["T"] = "Thr"
aminodict["Val"] = "V"
aminodict["V"] = "Val"
aminodict["Trp"] = "W"
aminodict["W"] = "Trp"
aminodict["Tyr"] = "Y"
aminodict["Y"] = "Tyr"
aminodict["*"] = "stop_gained"
aminodict["stop_gained"] = "*"
aminodict["X"] = "Xaa"
aminodict["Xaa"] = 'X'

#vcf_file_link = "/homes/janbaas/NAVIP_prj/members/Snpeff_VS_NAVIP/2019-02-13/better_format/All_VCF.vcf"
vcf_file_link = "/homes/janbaas/NAVIP_prj/members/Snpeff_VS_NAVIP/2019-02-13/neu-validiert/first/All_VCF.vcf"
vcf_file = open(vcf_file_link, 'r')
vcf_data = vcf_file.readlines()
vcf_file.close()

some_dict = {}
for line in vcf_data:
	if line.startswith('#'):
		continue
	spline = line.split("\t")
	if (spline[0],int(spline[1])) in some_dict.keys():
		some_dict[(spline[0], int(spline[1]))] += ";" + spline[7].split(";")[0] + ";" + spline[7].split(";")[1] #+ nav1, +nav 2
		#0		1		2	3	4	5		6		7
		#Chr1	1219320	.	T	C	4083.77	PASS	NAV2=C|SUB|LOW|AT1G04490|AT1G04490|transcript|AT1G04490.2|protein_coding||c.6A>G|p.Glu2Glu|/|6/1206|2/402||;NAV1=AT1G04490.2|REV|SUB|NONE|TTC/2|e|6|CTC/2|e|
	else:
		some_dict[(spline[0], int(spline[1]))] = line

only_my_keys = []
somemore_transcripts_dict = {}
for line in some_dict.values():
	#NAV2=C|SUB|LOW|AT1G04490|AT1G04490|transcript|AT1G04490.2|protein_coding||c.6A>G|p.Glu2Glu|/|6/1206|2/402||
	#NAV1=AT1G04490.2|REV|SUB|NONE|TTC/2|e|6|CTC/2|e|
	#ANN=C|missense_variant|MODERATE|AT1G04500|AT1G04500|transcript|AT1G04500.1|protein_coding|3/5|c.775C>G|p.Gln259Glu|906/1616|775/1161|259/386||
	##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type |
	# Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	# CHR,POS,TRANSCRIPT,NAV1
	# CHR,POS,TRANSCRIPT,NAV2
	# CHR,POS,TRANSCRIPT,ANN
	spline = line.split("\t")
	for infolineentry in spline[7].replace("\n","").split(";"):
		if infolineentry.startswith("NAV2="):
			somemore_transcripts_dict[(spline[0], int(spline[1]), infolineentry.split("|")[6] , 'NAV2')] = infolineentry
			only_my_keys.append((spline[0], int(spline[1]), infolineentry.split("|")[6])) # only for keys
		elif infolineentry.startswith("NAV1="):
			somemore_transcripts_dict[(spline[0], int(spline[1]), infolineentry.split("|")[0].split("=")[1], 'NAV1')] = infolineentry
		elif infolineentry.startswith("ANN="):
			for spANNline in infolineentry.split(","):
				somemore_transcripts_dict[(spline[0], int(spline[1]), spANNline.split("|")[6], 'ANN')] = spANNline
outsorted = 0
identical = 0
not_found = 0
navip_dict_every_non_indel_prot_substitution = {}
snpeff_dict_every_non_indel_prot_substitution = {}

stop_differences = {}

navip_dict_only_dna_subs = {}
snpeff_dict_only_dna_subs = {}

navip_dict_only_dna_subs_and_mvc = {} #mvc == multiple variant codon
snpeff_dict_only_dna_subs_and_mvc = {}


for from_aa_init in aminodict:
	for to_aa_init in aminodict:
		if from_aa_init == '*':
			from_aa_init = 'Ter'
		if to_aa_init == '*':
			to_aa_init = 'Ter'

		navip_dict_every_non_indel_prot_substitution[from_aa_init, to_aa_init] = 0
		snpeff_dict_every_non_indel_prot_substitution[from_aa_init, to_aa_init] = 0
		navip_dict_only_dna_subs[from_aa_init, to_aa_init] = 0
		snpeff_dict_only_dna_subs[from_aa_init, to_aa_init] = 0
		navip_dict_only_dna_subs_and_mvc[from_aa_init, to_aa_init] = 0.0
		snpeff_dict_only_dna_subs_and_mvc[from_aa_init, to_aa_init] = 0.0
		stop_differences[from_aa_init] = 0.0

only_my_keys_sorted_and_sub = []
only_my_keys_sorted_and_sub_dict = {}

test1 = 0
test2 = 0

for keys in only_my_keys:
	chr =keys[0]
	pos =keys[1]
	tra =keys[2]
	try:
		nav1 = somemore_transcripts_dict[(chr, int(pos), tra, "NAV1")]
		nav2 = somemore_transcripts_dict[(chr, int(pos), tra, "NAV2")]
		speff = somemore_transcripts_dict[(chr, int(pos), tra, "ANN")]
	except KeyError:
		not_found +=1
		continue
	"""
	NAV1=AT1G01030.1|REV|SUB|NONE|AAC/2|v|357|CAC/2|v|357
	NAV2=C|SUB|LOW|AT1G01030|AT1G01030|transcript|AT1G01030.1|protein_coding||c.357T>G|p.Val119Val|/|357/1077|119/359||
	ANN=C|synonymous_variant|LOW|AT1G01030|AT1G01030|transcript|AT1G01030.1|protein_coding|2/2|c.357T>G|p.Val119Val|970/1905|357/1077|119/358||
	"""

	nav2_aa_dup = different_AA(nav2.split("|")[10].replace("p.",""))
	snpeff_aa_dup = different_AA(speff.split("|")[10].replace("p.",""))

	if nav2_aa_dup == ['']:
		outsorted +=1
		continue
	nav2_aa_dup[0] = nav2_aa_dup[0].replace("*", "Ter")
	nav2_aa_dup[1] = nav2_aa_dup[1].replace("*", "Ter")

	if snpeff_aa_dup == [''] or snpeff_aa_dup == ['',''] or snpeff_aa_dup[0] == '' or snpeff_aa_dup[1] == '' :
		outsorted +=1
		continue
	snpeff_aa_dup[0] = snpeff_aa_dup[0].replace("*", "Ter")
	snpeff_aa_dup[1] = snpeff_aa_dup[1].replace("*", "Ter")

	#if nav2_aa_dup == snpeff_aa_dup:
	#	identical += 1
	#	continue
	if "del" == nav2_aa_dup[1] \
			or "ins" == nav2_aa_dup[1] \
			or 'dup' == nav2_aa_dup[1] \
			or nav2_aa_dup[1] == 'fs' \
			or nav2_aa_dup[1] == 'Ter?' \
			or len(nav2_aa_dup[0]) > 3 \
			or len(nav2_aa_dup[1]) > 3:
		outsorted +=1
		continue
	elif len (snpeff_aa_dup[0]) > 3 \
			or len (snpeff_aa_dup[1]) > 3 \
			or snpeff_aa_dup[1] == 'del' \
			or snpeff_aa_dup[1] == 'dup' \
			or snpeff_aa_dup[1] == 'fs' \
			or snpeff_aa_dup[1] == 'Ter?' \
			or snpeff_aa_dup[1] == 'ins' \
			or '?' in snpeff_aa_dup[1]:
		outsorted +=1
		continue
	#print(speff)
	navip_dict_every_non_indel_prot_substitution[nav2_aa_dup[0], nav2_aa_dup[1]] += 1  # (1/len(entry))
	snpeff_dict_every_non_indel_prot_substitution[snpeff_aa_dup[0], snpeff_aa_dup[1]] += 1  # (1/len(entry))
	test1 += 1
	if 'SUB' in nav1:
		test2 += 1
		navip_dict_only_dna_subs[nav2_aa_dup[0], nav2_aa_dup[1]] += 1
		snpeff_dict_only_dna_subs[snpeff_aa_dup[0], snpeff_aa_dup[1]] += 1

		only_my_keys_sorted_and_sub.append((chr, pos, tra))
		only_my_keys_sorted_and_sub_dict[(chr, int(pos), tra, 'NAV1')] = nav1
		only_my_keys_sorted_and_sub_dict[(chr, int(pos), tra, 'NAV2')] = nav2
		only_my_keys_sorted_and_sub_dict[(chr, int(pos), tra, 'ANN')] = speff
list_of_stop_codons = []
more_stop_codon_action = []
for keys in only_my_keys_sorted_and_sub:
	"""
	NAV1=AT1G01030.1|REV|SUB|NONE|AAC/2|v|357|CAC/2|v|357
	NAV2=C|SUB|LOW|AT1G01030|AT1G01030|transcript|AT1G01030.1|protein_coding||c.357T>G|p.Val119Val|/|357/1077|119/359||
	ANN=C|synonymous_variant|LOW|AT1G01030|AT1G01030|transcript|AT1G01030.1|protein_coding|2/2|c.357T>G|p.Val119Val|970/1905|357/1077|119/358||
	"""

	chr =keys[0]
	pos =int(keys[1])
	tra =keys[2]

	nav1 = only_my_keys_sorted_and_sub_dict[(chr, pos, tra, "NAV1")]
	nav2 = only_my_keys_sorted_and_sub_dict[(chr, pos, tra, "NAV2")]
	speff = only_my_keys_sorted_and_sub_dict[(chr, pos, tra, "ANN")]

	aa_pos = nav2.split("|")[13].split("/")[0]

	nav2_aa_dup = different_AA(nav2.split("|")[10].replace("p.", ""))
	snpeff_aa_dup = different_AA(speff.split("|")[10].replace("p.", ""))


	if nav2_aa_dup == ['']:
		outsorted += 1
		continue
	nav2_aa_dup[0] = nav2_aa_dup[0].replace("*", "Ter")
	nav2_aa_dup[1] = nav2_aa_dup[1].replace("*", "Ter")

	if snpeff_aa_dup == ['']:
		outsorted += 1
		continue
	snpeff_aa_dup[0] = snpeff_aa_dup[0].replace("*", "Ter")
	snpeff_aa_dup[1] = snpeff_aa_dup[1].replace("*", "Ter")


	nav2_compare_list = []
	#try:
	#	nav2_compare_list.append(only_my_keys_sorted_and_sub_dict[(chr, pos +1 , tra, "NAV2")])
	#except KeyError:
	#	nav2_compare_list = []
	try:
		nav2_compare_list.append(only_my_keys_sorted_and_sub_dict[(chr, pos + 2, tra, "NAV2")])
		nav1_2 = only_my_keys_sorted_and_sub_dict[(chr, pos+2, tra, "NAV1")]

		nav2_2 = nav2_compare_list[0]
		speff_2 = only_my_keys_sorted_and_sub_dict[(chr, pos + 2, tra, "ANN")]

		nav2_aa_dup2 = different_AA(nav2_2.split("|")[10].replace("p.", ""))
		snpeff_aa_dup2 = different_AA(speff_2.split("|")[10].replace("p.", ""))

		if nav2_aa_dup2 == ['']:
			outsorted += 1
			continue
		nav2_aa_dup2[0] = nav2_aa_dup2[0].replace("*", "Ter")
		nav2_aa_dup2[1] = nav2_aa_dup2[1].replace("*", "Ter")

		if snpeff_aa_dup2 == ['']:
			outsorted += 1
			continue
		snpeff_aa_dup2[0] = snpeff_aa_dup2[0].replace("*", "Ter")
		snpeff_aa_dup2[1] = snpeff_aa_dup2[1].replace("*", "Ter")
	except KeyError:
		if nav2_compare_list == []:
			continue

	#navip_dict_only_dna_subs_and_mvc = {}  # mvc == multiple variant codon
	#snpeff_dict_only_dna_subs_and_mvc = {}

	if len(nav2_compare_list) == 1:
		if aa_pos == nav2_compare_list[0].split("|")[13].split("/")[0]:
			if nav2.split("|")[10] != nav2_compare_list[0].split("|")[10]:
				#print(str(chr) +"\t"+  str(pos) +"\t"+ str(tra))
				#print("impossible1")
				continue
			else:
				if 'Ter' in snpeff_aa_dup[0] or 'Ter' in snpeff_aa_dup[1] or 'Ter' in snpeff_aa_dup2[0] or 'Ter' in snpeff_aa_dup2[1]:
					if 'Ter' in nav2_aa_dup[0] or 'Ter' in nav2_aa_dup[1] or 'Ter' in nav2_aa_dup2[0] or 'Ter' in nav2_aa_dup2[1]:
						list_of_stop_codons.append((chr, pos, tra, 'navip+snpeff'))
						more_stop_codon_action.append("\t".join([str(chr), str(pos), str(tra),'navip', str(nav2_aa_dup[0]), str(nav2_aa_dup[1]), str(nav1.split('|')[4]), str(nav1.split('|')[7]) ]))
						more_stop_codon_action.append("\t".join([str(chr), str(pos+2), str(tra),'navip', str(nav2_aa_dup2[0]), str(nav2_aa_dup2[1]), str(nav1_2.split('|')[4]), str(nav1_2.split('|')[7]) ]))
						more_stop_codon_action.append("\t".join([str(chr), str(pos), str(tra), 'snpeff', str(snpeff_aa_dup[0]), str(snpeff_aa_dup[1])]))
						more_stop_codon_action.append("\t".join([str(chr), str(pos+2), str(tra), 'snpeff', str(snpeff_aa_dup2[0]), str(snpeff_aa_dup2[1])]))
						stop_differences[nav2_aa_dup[1]] += 1
					else:
						list_of_stop_codons.append((chr, pos, tra, 'snpeff'))
						more_stop_codon_action.append("\t".join([str(chr), str(pos), str(tra),'navip', str(nav2_aa_dup[0]), str(nav2_aa_dup[1]), str(nav1.split('|')[4]), str(nav1.split('|')[7]) ]))
						more_stop_codon_action.append("\t".join([str(chr), str(pos+2), str(tra),'navip', str(nav2_aa_dup2[0]), str(nav2_aa_dup2[1]), str(nav1_2.split('|')[4]), str(nav1_2.split('|')[7]) ]))
						more_stop_codon_action.append("\t".join([str(chr), str(pos), str(tra), 'snpeff', str(snpeff_aa_dup[0]), str(snpeff_aa_dup[1])]))
						more_stop_codon_action.append("\t".join([str(chr), str(pos+2), str(tra), 'snpeff', str(snpeff_aa_dup2[0]), str(snpeff_aa_dup2[1])]))
						stop_differences[nav2_aa_dup[1]] += 1
				elif 'Ter' in nav2_aa_dup[0] or 'Ter' in nav2_aa_dup[1]:
					list_of_stop_codons.append((chr, pos, tra, 'navip'))
					#more_stop_codon_action.append("\t".join([str(chr), str(pos), str(tra),'navip', str(nav2_aa_dup[0]), str(nav2_aa_dup[1]), str(nav1.split('|')[4]), str(nav1.split('|')[7]) ]))
					#more_stop_codon_action.append("\t".join([str(chr), str(pos+2), str(tra),'navip', str(nav2_aa_dup2[0]), str(nav2_aa_dup2[1]), str(nav1_2.split('|')[4]), str(nav1_2.split('|')[7]) ]))
					#more_stop_codon_action.append("\t".join([str(chr), str(pos), str(tra), 'snpeff', str(snpeff_aa_dup[0]), str(snpeff_aa_dup[1])]))
					#more_stop_codon_action.append("\t".join([str(chr), str(pos+2), str(tra), 'snpeff', str(snpeff_aa_dup2[0]), str(snpeff_aa_dup2[1])]))
				if nav2_aa_dup[0] != nav2_aa_dup2[0] or nav2_aa_dup[1] != nav2_aa_dup2[1]:
					print(str(tra)+"\t"+str(nav2_aa_dup[0])+"\t"+str(nav2_aa_dup[1])+"\t"+str(snpeff_aa_dup[0])+"\t"+str(snpeff_aa_dup[1]))
					print(str(tra)+"\t" + str(nav2_aa_dup2[0])+"\t"+str( nav2_aa_dup2[1])+"\t"+str(snpeff_aa_dup2[0])+"\t"+str( snpeff_aa_dup2[1]))
				navip_dict_only_dna_subs_and_mvc[nav2_aa_dup[0],nav2_aa_dup[1]] += 0.5
				snpeff_dict_only_dna_subs_and_mvc[snpeff_aa_dup[0],snpeff_aa_dup[1]] += 0.5

				navip_dict_only_dna_subs_and_mvc[nav2_aa_dup2[0], nav2_aa_dup2[1]] += 0.5
				snpeff_dict_only_dna_subs_and_mvc[snpeff_aa_dup2[0], snpeff_aa_dup2[1]] += 0.5

		else:
			continue
	if len(nav2_compare_list) == 2:
		print('hmmmh')
		if aa_pos == nav2_compare_list[0].split("|")[13].split("/")[0] :
			if nav2.split("|")[10] != nav2_compare_list[0].split("|")[10]:
				asd = nav2.split("|")[10]
				asdasdasd = nav2_compare_list[0].split("|")[10]
				print("AA_switching because of frameshift")
			else:
				if 'Ter' in snpeff_aa_dup[0] or 'Ter' in snpeff_aa_dup[1]:
					if 'Ter' in nav2_aa_dup[0] or 'Ter' in nav2_aa_dup[1]:
						list_of_stop_codons.append((chr, pos, tra, 'navip+snpeff'))
					else:
						list_of_stop_codons.append((chr, pos, tra, 'snpeff'))
				elif 'Ter' in nav2_aa_dup[0] or 'Ter' in nav2_aa_dup[1]:
					list_of_stop_codons.append((chr, pos, tra, 'navip'))
					print("aaaaaah")
				navip_dict_only_dna_subs_and_mvc[nav2_aa_dup[0], nav2_aa_dup[1]] += 0.5
				snpeff_dict_only_dna_subs_and_mvc[snpeff_aa_dup[0], snpeff_aa_dup[1]] += 0.5

				navip_dict_only_dna_subs_and_mvc[nav2_aa_dup2[0], nav2_aa_dup2[1]] += 0.5
				snpeff_dict_only_dna_subs_and_mvc[snpeff_aa_dup2[0], snpeff_aa_dup2[1]] += 0.5
		elif aa_pos == nav2_compare_list[1].split("|")[13].split("/")[0]:
			if nav2.split("|")[10] != nav2_compare_list[1].split("|")[10]:
				print("impossible3")
			else:
				if 'Ter' in snpeff_aa_dup[0] or 'Ter' in snpeff_aa_dup[1]:
					if 'Ter' in nav2_aa_dup[0] or 'Ter' in nav2_aa_dup[1]:
						list_of_stop_codons.append((chr, pos, tra, 'navip+snpeff'))
					else:
						list_of_stop_codons.append((chr, pos, tra, 'snpeff'))
				elif 'Ter' in nav2_aa_dup[0] or 'Ter' in nav2_aa_dup[1]:
					list_of_stop_codons.append((chr, pos, tra, 'navip'))
					print("bbbbb")
				navip_dict_only_dna_subs_and_mvc[nav2_aa_dup[0], nav2_aa_dup[1]] += 0.5
				snpeff_dict_only_dna_subs_and_mvc[snpeff_aa_dup[0], snpeff_aa_dup[1]] += 0.5

				navip_dict_only_dna_subs_and_mvc[nav2_aa_dup2[0], nav2_aa_dup2[1]] += 1
				snpeff_dict_only_dna_subs_and_mvc[snpeff_aa_dup2[0], snpeff_aa_dup2[1]] += 0.5
		else:
			continue
	if len(nav2_compare_list) > 2:
		print("dafuk?")



AA = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'Ter']

output = ".\t" + "\t".join(AA) + "\n"
output2 = ".\t" + "\t".join(AA) + "\n"

for from_aa in AA:
	output += from_aa
	output2 +=from_aa
	for to_aa in AA:
		#for number in numbers:
		#	print(f'{number:9.4f}')
		#output += "\t" + str("{:10.2f}".format(navip_dict[my_aa_dup, snpeff_aa_dup]))
		#output2+= "\t" + str("{:10.2f}".format(snpeff_dict[my_aa_dup, snpeff_aa_dup]))
		output += "\t" + str(navip_dict_every_non_indel_prot_substitution[from_aa, to_aa])
		output2 += "\t" + str(snpeff_dict_every_non_indel_prot_substitution[from_aa, to_aa])
	output += "\n"
	output2 += "\n"
print("navip_dict_every_non_indel_prot_substitution\n")
print(output)
print("snpeff_dict_every_non_indel_prot_substitution\n")
print(output2)

""" identical to first output
output = ".\t" + "\t".join(AA) + "\n"
output2 = ".\t" + "\t".join(AA) + "\n"

for from_aa in AA:
	output += from_aa
	output2 +=from_aa
	for to_aa in AA:
		#for number in numbers:
		#	print(f'{number:9.4f}')
		#output += "\t" + str("{:10.2f}".format(navip_dict[my_aa_dup, snpeff_aa_dup]))
		#output2+= "\t" + str("{:10.2f}".format(snpeff_dict[my_aa_dup, snpeff_aa_dup]))
		output += "\t" + str(navip_dict_only_dna_subs[from_aa, to_aa])
		output2 += "\t" + str(navip_dict_only_dna_subs[from_aa, to_aa])
	output += "\n"
	output2 += "\n"
print("navip_dict_only_dna_subs\n")
print(output)
print("snpeff_dict_only_dna_subs\n")
print(output2)
"""



output = ".\t" + "\t".join(AA) + "\n"
output2 = ".\t" + "\t".join(AA) + "\n"
output3 = ".\t" + "\t".join(AA) + "\n"

for from_aa in AA:
	output += from_aa
	output2 +=from_aa
	output3 += "\t" + str("{:10.2f}".format(stop_differences[from_aa]))
	for to_aa in AA:
		#for number in numbers:
		#	print(f'{number:9.4f}')
		output += "\t" + str("{:10.2f}".format(navip_dict_only_dna_subs_and_mvc[from_aa, to_aa]))
		output2+= "\t" + str("{:10.2f}".format(snpeff_dict_only_dna_subs_and_mvc[from_aa, to_aa]))
		#output += "\t" + str(navip_dict_only_dna_subs_and_mvc[from_aa, to_aa])
		#output2 += "\t" + str(snpeff_dict_only_dna_subs_and_mvc[from_aa, to_aa])

	output += "\n"
	output2 += "\n"
print("navip_dict_only_dna_subs_and_mvc, always counted only one of the multiple snips per codon\n")
print(output)
print("snpeff_dict_only_dna_subs_and_mvc, only first snip used, but all are different\n")
print(output2)



print(list_of_stop_codons)
print(len(list_of_stop_codons))
"""
smalerdict = {}
for a in more_stop_codon_action:
	if a in smalerdict:
		smalerdict[a] += 1
	else:
		smalerdict[a] = 1
for a in smalerdict.keys():
	print(str(smalerdict[a]) +" "+ str(a))
"""
print("\n".join(more_stop_codon_action))
print(output3)
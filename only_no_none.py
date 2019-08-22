
def only_no_none(file_path:str, file_out:str):
	all_vcf_file = open(file_path, 'r')
	all_lines = all_vcf_file.readlines()
	all_vcf_file.close()

	all_no_none_lines = []

	for line in all_lines:
		if line.startswith("#"):
			continue
		elif len(line.split("\t")) <= 4:
			continue
		spline = line.split("\t")
		infoline = spline[7]
		for info in infoline.split(";"):
			if info.startswith("NAV1="):
				shared_keys = info.split("|")[3]
				if shared_keys == "NONE":
					continue
				else:
					all_no_none_lines.append(line)

	##Info=<ID=NAV1, Type=String,Number=.,Values=[TranscriptID|Strand_Direction|Variant_Classification1,Variant_Classification2,...|Shared_EffKey(s)|REF_Codon(s)/Variant_Position_in_Codon|REF_AA|old_CDS_Position|ALT_Codon(s)/Variant_Position_in_Codon|ALT_AA|new_CDS_Position]>
	only_no_none_file = open(file_out, 'w')
	only_no_none_file.write("".join(all_no_none_lines))
	only_no_none_file.close()


def add_stable_id_tag(vcf_file:str,output_vcf_file):
	"""
	This def will create a <NAVID=ID> tag inside the infoline in the vcf-file.
	The ID will be the position of the variant. This is needed for a few special cases, in which the variant-position
	will change after the pre-processing.
	:param vcf_file: incomming VCF-File.
	:param output_vcf_file: Outgoing vcf-file.
	:return: nothing
	"""

	all_vcf_file = open(vcf_file, 'r')
	output_vcf = open(output_vcf_file, 'w')
	line = all_vcf_file.readline()
	while line:
		if line.startswith("#"):
			output_vcf.write(line)
			line = all_vcf_file.readline()
			continue
		elif len(line.split("\t")) <= 4:
			output_vcf.write(line)
			line = all_vcf_file.readline()
			continue
		spline = line.split("\t")
		#infoline = spline[7]
		newinfo = "NAVID=" + str(spline[1])
		spline[7] = newinfo + ";" + spline[7]
		output_vcf.write("\t".join(spline))
		line = all_vcf_file.readline()
	all_vcf_file.close()
	output_vcf.close()

def remove_stable_id_tag(vcf_file:str,output_vcf_file):
	"""
	This will def will simply remove the <NAVID=ID> Tag. Useful, because the Tag is a solution for a minor problem.
	:param vcf_file: incomming VCF-File.
	:param output_vcf_file: Outgoing vcf-file.
	:return: nothing
	"""

	all_vcf_file = open(vcf_file, 'r')
	output_vcf = open(output_vcf_file, 'w')
	line = all_vcf_file.readline()
	while line:
		if line.startswith("#"):
			output_vcf.write(line)
			line = all_vcf_file.readline()
			continue
		elif len(line.split("\t")) <= 4:
			output_vcf.write(line)
			line = all_vcf_file.readline()
			continue
		spline = line.split("\t")
		infoline = spline[7]
		newinfo = []
		for info in infoline.split(";"):
			if not info.startswith("NAVID="):
				newinfo.append( info)

		spline[7] = ";".join(newinfo)
		output_vcf.write("\t".join(spline))
		line = all_vcf_file.readline()
	all_vcf_file.close()
	output_vcf.close()
__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import VCF_Variant

class VCF_Handler:
	"""
	This class will gather all information from a given VCF File together and will
	store every data line into a VCF_Variant object.
	"""

	# JustXChromosomes = 0 -> all
	def __init__ (self, VCF_data_path: str, justx_chromosomes: int):
		"""
		Will read an entire VCF file, handle and store all its data.
		:param VCF_data_path: Path to (including) the file.
		:param justx_chromosomes: Zero for all chromosomes, Any other number for the first X chromosomes.
		"""
		#Initialization
		self.VCF_list_chromosomes = [[[]]]
		self.VCF_variant_list = [[]] # => Genome: [ Chromosome: [Variant(s)] ]
		self.dict_chr_names = {}

		count_chr = 0

		with open(VCF_data_path, "r") as DataFile :
			id = 0
			lines = DataFile.readline()
			while lines.startswith("#"):
				lines = DataFile.readline()
			chromosome_flag =  lines.split('\t' , 1)[0]
			current_list = [[]]
			current_variant_list = []
			count_variant = 0

			# for easy name-id relation
			self.dict_chr_names[id] = chromosome_flag
			self.dict_chr_names[chromosome_flag] = id


			while lines:

				while lines.startswith("#"):
					lines = DataFile.readline()
					continue

				if lines.split("\t")[0] == chromosome_flag:
					now = lines.split('\t')
					if len(now) >8:
						infoline = "\t".join(now[7:])
					elif len(now) == 8:
						infoline = now[7]
					current_list.append(now)
					variant = VCF_Variant.Variant(
											  now[0],		 #Chr
											  int(now[1]),	 #Pos
											  now[2],		 #ID
											  count_variant, #useful ID (not VCF_File)
											  now[3],		 #Ref
											  now[4],		 #Alternate
											  now[5],		 #Qual
											  now[6],		 #Filter
										 	  infoline		 #Info
											  )

					current_variant_list.append(variant)

					count_variant += 1

					lines = DataFile.readline()
				elif lines.startswith("###"):
					continue # because info-field
				else :
					count_chr += 1

					chromosome_flag =  lines.split('\t' , 1)[0]
					if not current_list[0]:
						current_list.pop(0) # To Remove the [[]]-Entry
					self.VCF_list_chromosomes.append(current_list)
					current_list = [[]] # That's current_list.clear in Python 2.7

					if not current_variant_list[0]:
						current_variant_list.pop(0)
					self.VCF_variant_list.append(current_variant_list)
					current_variant_list = []


					if count_chr >= justx_chromosomes != 0:
						break

					# for easy name-id relation
					id +=1
					self.dict_chr_names[id] = chromosome_flag
					self.dict_chr_names[chromosome_flag] = id

			if justx_chromosomes == 0:
				self.VCF_list_chromosomes.append(current_list)
				current_list = [[]]
				self.VCF_variant_list.append(current_variant_list)



			if not current_list[0]:
				current_list.pop(0) # To Remove the [[]]-Entry
				self.VCF_list_chromosomes.append(current_list)


			if self.VCF_list_chromosomes[0] == [[]]:
				self.VCF_list_chromosomes.pop(0) # To Remove the [[]]-Entry

			if not self.VCF_variant_list[0]:
				self.VCF_variant_list.pop(0)


			DataFile.close()

	def get_chromosome_names(self):
		"""
		Returns a list of all chromosome names.
		:return: List of all chromosome names.
		"""
		count_names = len(self.dict_chr_names) / 2
		i = 0
		name_list = []
		while count_names != i:
			name_list.append(self.dict_chr_names[i])
			i += 1
		if not name_list[0]:
			name_list.pop(0)

		return name_list

	def get_chr_VCF_variant_list(self, chr_name: str):
		"""
		Returns a list of all variants inside the specified chromosome.
		:param chr_name: Chromosome name.
		:return: List of all variants inside the specified chromosome.
		"""
		try:
			return self.VCF_variant_list[self.dict_chr_names[chr_name]]
		except KeyError:
			return []

	def free_RAM(self, chr_name:str):

		self.VCF_variant_list[self.dict_chr_names[chr_name]] = []


class NAVIP_VCF_File_Manager:
	def __init__(self, navip_file_path: str):
		vcf_handler = VCF_Handler(navip_file_path, 0)
		self.chr_list = vcf_handler.get_chromosome_names()
		self.variants_in_list_in_dict = {}
		for name in self.chr_list:
			if name not in self.variants_in_list_in_dict.keys():
				self.variants_in_list_in_dict[name] = []
			for v_list in vcf_handler.get_chr_VCF_variant_list(name):
				for variant in v_list:
					self.variants_in_list_in_dict[name].add(VCF_Variant.Variant_NAVIP(variant))

	def merge_with_another_navip_file_manager(self,navip_file_manager):
		# needs a check, if the new file manager has chr entries, which the current does not have.
		# they can easily be added
		merged_dict = {}
		for name in self.chr_list:
			if name not in navip_file_manager.chr_list:
				continue

			v_in_list_in_dict = {}
			for variant in self.variants_in_list_in_dict[name]:
				try:
					v_in_list_in_dict[variant.Chromosome, variant.Position].add(variant)
				except KeyError:
					v_in_list_in_dict[variant.Chromosome, variant.Position] = [variant]
			for variant in self.chr_list[name]:
				try:
					v_in_list_in_dict[variant.Chromosome, variant.Position].add(variant)
				except KeyError:
					v_in_list_in_dict[variant.Chromosome, variant.Position] = [variant]
			for position in v_in_list_in_dict.keys():
				v_dict = {}
				for variant in v_in_list_in_dict[position]:
					v_dict[variant.Info] = variant #overwrite, if present == identical
				v_in_list_in_dict[position] = list(v_dict)
			merged_dict[name] = v_in_list_in_dict
			# so it's here a dict with entries per chromosome.
			# in these entries are dicts with entries for every position
			# in these entries are lists with all unique, merged variants

		for name in merged_dict.keys(): # name keys
			for variant_in_list_in_dict_per_position in merged_dict[name]: #dict with list-entry for positions
				for position in variant_in_list_in_dict_per_position.keys(): #position keys
					print ("experimental")

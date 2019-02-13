__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import VCF_Variant

class VCF_HANDLER:
	"""
	This class will gather all information from a given VCF File together and will
	store every data line into a VCF_Variant object.
	"""

	# JustXChromosomes = 0 -> all
	def __init__ (self, VCF_DATA_PATH: str,JustXChromosomes: int):
		"""
		Will read an entire VCF file, handle and store all its data.
		:param VCF_DATA_PATH: Path to (including) the file.
		:param JustXChromosomes: Zero for all chromosomes, Any other number for the first X chromosomes.
		"""
		#Initialization
		self.VCF_ListChromosomes = [[[]]]
		self.VCF_Variant_List = [[]] # => Genom: [ Chromosome: [Variant(s)] ]
		self.dictChrNames = {}

		countChr = 0

		with open(VCF_DATA_PATH, "r") as DataFile :
			id = 0;
			lines = DataFile.readline()
			while (lines.startswith("#")):
				lines = DataFile.readline()
			ChromosomFlag =  lines.split('\t' , 1)[0]
			CurrentList = [[]]
			CurrentVariantList = []
			CountVariant = 0

			# for easy name-id relation
			self.dictChrNames[id] = ChromosomFlag
			self.dictChrNames[ChromosomFlag] = id


			while (lines):

				while (lines.startswith("#")):
					lines = DataFile.readline()
					continue

				if (lines.split("\t")[0] == ChromosomFlag) :
					Now = lines.split('\t')
					if len(Now) >8:
						infoline = "\t".join(Now[7:])
					elif len(Now) == 8:
						infoline = Now[7]
					CurrentList.append(Now)
					varianti = VCF_Variant.Variant(
											  Now[0],		# Chr
											  int(Now[1]), # Pos
											  Now[2],		#ID
											  CountVariant, #usefull ID (not VCF_File)
											  Now[3],		#Ref
											  Now[4],		#Alternate
											  Now[5],		#Qual
											  Now[6],		#Filter
										 	  infoline		#Info
											  )

					CurrentVariantList.append(varianti)

					CountVariant += 1

					lines = DataFile.readline()
				elif (lines.startswith("###")):
					continue # because info-field
				else :
					countChr = +1

					ChromosomFlag =  lines.split('\t' , 1)[0]
					if (CurrentList[0] == []):
						CurrentList.pop(0) # To Remove the [[]]-Entry
					self.VCF_ListChromosomes.append(CurrentList)
					CurrentList = [[]] # Thats CurrentList.clear in Python 2.7

					if (CurrentVariantList[0] == []):
						CurrentVariantList.pop(0)
					self.VCF_Variant_List.append(CurrentVariantList)
					CurrentVariantList = []


					if (countChr >= JustXChromosomes and  JustXChromosomes != 0 ) :
						break

					# for easy name-id relation
					id +=1
					self.dictChrNames[id] = ChromosomFlag
					self.dictChrNames[ChromosomFlag] = id

			if JustXChromosomes == 0:
				self.VCF_ListChromosomes.append(CurrentList)
				CurrentList = [[]]
				self.VCF_Variant_List.append(CurrentVariantList)
				CurrentVariantList = []



			if (CurrentList[0] == []):
				CurrentList.pop(0) # To Remove the [[]]-Entry
				self.VCF_ListChromosomes.append(CurrentList)


			if (self.VCF_ListChromosomes[0] == [[]] ):
				self.VCF_ListChromosomes.pop(0) # To Remove the [[]]-Entry

			if self.VCF_Variant_List[0] == []:
				self.VCF_Variant_List.pop(0)


			DataFile.close()

	def GetChromosomeNames(self):
		"""
		Returns a list of all chromosome names.
		:return: List of all chromosome names.
		"""
		countNames = len(self.dictChrNames) / 2
		i = 0
		NameList = []
		while countNames != i:
			NameList.append(self.dictChrNames[i])
			i += 1
		if NameList[0] == []:
			NameList.pop(0)

		return NameList

	def GetChr_VCF_Variant_List(self, ChrName: str):
		"""
		Returns a list of all variants inside the specified chromosome.
		:param ChrName: Chromosome name.
		:return: List of all variants inside the specified chromosome.
		"""
		try:
			return self.VCF_Variant_List[self.dictChrNames[ChrName]]
		except KeyError as e:
			return []


class NAVIP_VCF_File_Manager():
	def __init__(self, navip_file_path: str):
		vcf_handler = VCF_HANDLER(navip_file_path,0)
		self.chr_list = vcf_handler.GetChromosomeNames()
		self.variants_in_list_in_dict = {}
		for name in self.chr_list:
			if name not in self.variants_in_list_in_dict.keys():
				self.variants_in_list_in_dict[name] = []
			for v_list in vcf_handler.GetChr_VCF_Variant_List(name):
				for variant in v_list:
					self.variants_in_list_in_dict[name].add(VCF_Variant.Variant_NAVIP(variant))

	def merge_with_another_navip_file_manager(self,navip_file_manager):
		# needs a check, if the new file manager has chr entrys, which the current does not have.
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
			# so its here a dict with entrys per chromosome.
			# in this entrys are dicts with entrys for every position
			# in this entrys are lists with all unique, merged variants

		final_new_variants_in_list_in_dict = {}
		for name in merged_dict.keys(): # name keys
			for variant_in_list_in_dict_per_position in merged_dict[name]: #dict with list-entry for positions
				for postion in variant_in_list_in_dict_per_position.keys(): #position keys
					print ("experimental")



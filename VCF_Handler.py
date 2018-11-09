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
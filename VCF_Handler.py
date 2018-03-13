import VCF_Variant
from enum import Enum

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
		self.VCF_Dict_Variant_List = [{}] # For each Chromosome one dictionary will all its variants. the dcits in a list
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
			CurrentVariantDict = {}
			CountVariant = 0

			# for easy name-id relation
			self.dictChrNames[id] = ChromosomFlag
			self.dictChrNames[ChromosomFlag] = id


			while (lines):

				if (lines.split("\t")[0] == ChromosomFlag) :
					Now = lines.split('\t')
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
											  Now[7]		#Info
											  )

					CurrentVariantList.append(varianti)
					CurrentVariantDict[CountVariant] = varianti
					CurrentVariantDict[varianti.Chromosome + "." + str(varianti.Position)] = varianti

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

					self.VCF_Dict_Variant_List.append(CurrentVariantDict)
					CurrentVariantDict = {}

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
				self.VCF_Dict_Variant_List.append(CurrentVariantDict)
				CurrentVariantDict = {}



			if (CurrentList[0] == []):
				CurrentList.pop(0) # To Remove the [[]]-Entry
				self.VCF_ListChromosomes.append(CurrentList)


			if (self.VCF_ListChromosomes[0] == [[]] ):
				self.VCF_ListChromosomes.pop(0) # To Remove the [[]]-Entry

			if self.VCF_Variant_List[0] == []:
				self.VCF_Variant_List.pop(0)

			if self.VCF_Dict_Variant_List[0] == {}:
				self.VCF_Dict_Variant_List.pop(0)

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

	def GetChrDict_Variants(self, ChrName: str):
		"""
		Returns a dictionary of variants from the specified chromosome.
		:param ChrName: Chromosome name.
		:return: A dictionary of variants from the specified chromosome.
		"""
		if (self.VCF_Dict_Variant_List[0] == {}) :
			self.VCF_Dict_Variant_List.pop(0)
		return self.VCF_Dict_Variant_List[self.dictChrNames[ChrName]]

	def GetOneVariant (self, ChrName: str, VariantID: int) -> VCF_Variant:
		"""
		Returns one variant from one chromosome.
		:param ChrName: Chromosome name.
		:param VariantID: Usefull ID from the variant.
		:return: The VCF_Variant object.
		"""
		if (self.VCF_Dict_Variant_List[0] == {}) :
			self.VCF_Dict_Variant_List.pop(0)
		return self.VCF_Dict_Variant_List[self.dictChrNames[ChrName]][VariantID]
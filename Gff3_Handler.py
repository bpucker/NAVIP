__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import Transcript
import sys

class GFF3_Handler:
	"""
	This class reads, contains and manages all used data from the gff3 file.
	Furthermore it creates and contains all transcript objects.
	"""

	def __init__ (self, GFF3_DATA_PATH: str, JustXChromosomes: int) :
		"""
		Reads the entire Gff3 data or the first x chromosomes.
		Searches for "gene", "UTR", and "RNA" to create transcript objects.
		Works for the gff3 augustus format, too. In this case Parend-IDs will be used, too.
		:param GFF3_DATA_PATH: Path and filename of the gff3 file.
		:param JustXChromosomes: The number of the first X chromosomes to read. 0 for all chromosomes available.
		"""
		#Initialization
		self.InfoField = "" # gene info field
		self.ListChromosomes =  [[[]]] # same as dictListOfTranscripts, but as a list
		self.ListOfTranscripts = [[]] # temporary list for all transcripts of a chromosome
		self.dictListOfTranscripts = [{}] # List for chromosomes, dict for normal entrys
		self.dictChrNames = {}
		self.ListOfSpecialLNCTranscripts = [[]] #long non coding transcripts
		self.transIndex = 0 #id for every transcript
		self.TID_DICT = {} # acces to transcripts over the TID

		with open(GFF3_DATA_PATH, "r") as DataFile :
			# Version of GFF
			lines = DataFile.readline() # First Data-Line
			while lines.startswith("#"):
				if self.InfoField == "":
					self.InfoField = DataFile.readline()
				lines = DataFile.readline()
			ChromosomFlag =  lines.split('\t' , 1)[0]
			countChromosomes = 0
			self.dictChrNames[countChromosomes] = ChromosomFlag
			self.dictChrNames[ChromosomFlag] = countChromosomes
			CurrentList = [[]]
			CurrentTransList = []
			countlines = 0
			self.transIndex = 0
			littleDict = {}
			CurrentSpeciallncTransList = []
			Gene_Info_String = ""
			Gene_Start_End_Position = (0,0)
			UTR_Description = []

			while (lines):
				while lines.startswith("#"):
					lines = DataFile.readline()
					continue


				if (lines.split("\t")[0] == ChromosomFlag) :
					CurrentLine = lines.split('\t')
					if "gene" in CurrentLine[2]:
						Gene_Info_String = CurrentLine[8]
						Gene_Start_End_Position = (int(CurrentLine[3]), int(CurrentLine[4]))
						UTR_Description = []
					###
					# sorting out after locus_type
					###
					Key_Error_Unknown_type = True
					if "locus_type=novel_transcribed_region" in Gene_Info_String \
							or "locus_type=pseudogene" in Gene_Info_String:
						Key_Error_Unknown_type = False

					if "UTR" in CurrentLine[2] \
							or "transcription_end_site" in CurrentLine[2] \
							or "transcription_start_site" in CurrentLine[2]:
						UTR_Description.append("".join(CurrentLine[2:]))

					if "gene" not in CurrentLine[2]:
						parent_id_test = CurrentLine[8].split(";")
						for testForId in parent_id_test:

							if "Parent" not in testForId:
								continue
							if "," in testForId and "UTR" not in CurrentLine[2] and "exon" not in CurrentLine[2]:
								print("(Maybe)Can't handle entry(critical with CDS-entry): " + str(CurrentLine[8]))
								break
							elif "UTR" in CurrentLine[2]:
								#usefull information, if the transcript is somehow broken
								#Parent=AT1G67300.5,AT1G67300.7;
								#0123456
								IDS = testForId[7:]
								IDS = IDS.split(",")
								for ID in IDS:
									ID = ID.replace("\n", "")
									try:
										self.TID_DICT[ID].AddUTR_Description(lines)
									except KeyError:
										if Key_Error_Unknown_type:
											print("No transcript entry: " + str(ID) + "\t" + CurrentLine[2] + CurrentLine[8])
										pass

							elif "exon" in CurrentLine[2]:
								IDS = testForId[7:]
								IDS = IDS.split(",")
								for ID in IDS:
									ID = ID.replace("\n","")
									try:
										self.TID_DICT[ID].AddEXON_Descriptin(lines)
									except KeyError:
										if Key_Error_Unknown_type:
											print("No transcript entry: " +str(ID)+ "\t" + CurrentLine[2] + CurrentLine[8])
										pass

					if "RNA" in CurrentLine[2] or "transcript" == str(CurrentLine[2]):
						TID = CurrentLine[8][3:CurrentLine[8].find(";")]
						#def __init__ (self, IndexKey: int, TID: int, StartOfRNA: int, EndOfRNA: int, Forward: bool):

						if (CurrentLine[6] == "-"):
							ForwardOrReverse = Transcript.TranscriptEnum.REVERSE
						elif (CurrentLine[6] == "+"):
							ForwardOrReverse = Transcript.TranscriptEnum.FORWARD
						else:
							ForwardOrReverse = Transcript.TranscriptEnum.UNKNOWN_STRAND_DIRECTION
						trans = Transcript.Transcript(int(self.transIndex), str(TID), int(CurrentLine[3]), int(CurrentLine[4]), ForwardOrReverse)
						trans.SetGene_Info_String(Gene_Info_String)
						trans.SetGene_Start_Position(Gene_Start_End_Position[0])
						trans.SetGene_End_Position(Gene_Start_End_Position[1])
						self.transIndex +=1
						CurrentTransList.append(trans)
						littleDict[(self.transIndex-1)] = trans
						self.TID_DICT[trans.TID] = trans


						if ("lnc" in CurrentLine[2]): #long non coding
							CurrentSpeciallncTransList.append(trans)

					elif ("CDS" in CurrentLine[2]):
						parent_id_test = CurrentLine[8].split(";")
						for testForId in parent_id_test:

							if "Parent" not in testForId:
								continue

							IDS = testForId[7:]
							IDS = IDS.split(",")
							for ID in IDS:
								ID = ID.replace("\n", "")
								try:
									self.TID_DICT[ID].addCDS(int(CurrentLine[3]), int(CurrentLine[4]), str(CurrentLine[7]))
								except KeyError:
									if Key_Error_Unknown_type:
										print("No transcript entry: " + str(ID) + "\t" + CurrentLine[2] + CurrentLine[8])
									pass

						#trans.addCDS(int(CurrentLine[3]), int(CurrentLine[4]), int(CurrentLine[7]))


					CurrentList.append(CurrentLine)
					lines = DataFile.readline()
					countlines += 1
				elif lines.startswith("###") or lines.startswith("#") :
					lines = DataFile.readline()
				else :
					######
					## For next chromosome data
					######
					if (self.dictListOfTranscripts[0] == {}):
						self.dictListOfTranscripts.pop(0)
					self.dictListOfTranscripts.append(littleDict)
					littleDict = {}

					ChromosomFlag =  lines.split('\t' , 1)[0]
					countChromosomes += 1


					######
					## For the next transcripts (in next chromosome)
					######

					if (CurrentTransList[0] ==  []):
						CurrentTransList.pop(0)# To Remove the []-Entry
					if (CurrentSpeciallncTransList == [] and len(CurrentSpeciallncTransList) > 0):
						CurrentSpeciallncTransList.pop(0)

					CurrentTransList = sorted(CurrentTransList, key = lambda sTranscript: sTranscript.StartOfRNA) #todo test me
					self.ListOfTranscripts.append(CurrentTransList)
					CurrentTransList = []
					if (self.ListOfSpecialLNCTranscripts[0] == []):
						self.ListOfSpecialLNCTranscripts.pop(0)
					self.ListOfSpecialLNCTranscripts.append(CurrentSpeciallncTransList)


					if(countChromosomes >= JustXChromosomes and (JustXChromosomes != 0)):
						break

					self.dictChrNames[countChromosomes] = ChromosomFlag
					self.dictChrNames[ChromosomFlag] = countChromosomes
					if(CurrentList[0] == []):
						CurrentList.pop(0) # To ReNdChr2.g10821.t1move the [[]]-Entry
					self.ListChromosomes.append(CurrentList)
					CurrentList = [[]] # Thats CurrentList.clear in Python 2.7
					#CurrentTransList = []


			# after while
			# for the last entry

			######
			## For next chromosome data
			######
			if littleDict != {}:
				if (self.dictListOfTranscripts[0] == {}):
					self.dictListOfTranscripts.pop(0)
				self.dictListOfTranscripts.append(littleDict)



			if(CurrentList[0] == []):
				CurrentList.pop(0) # To Remove the [[]]-Entry
				self.ListChromosomes.append(CurrentList)

			if (self.ListChromosomes[0] == [[]] ):
				self.ListChromosomes.pop(0) # To Remove the [[]]-Entry

			#if (CurrentTransList[0] ==  []):
			#	CurrentTransList.pop(0)# To Remove the []-Entry

			if (self.ListOfTranscripts[0] == []):
				self.ListOfTranscripts.pop(0)# To Remove the [[]]-Entry
			self.ListOfTranscripts.append(CurrentTransList)


		DataFile.close()

	def GetNextTranscriptIndex (self):
		"""
		For creating new transcripts (need for a new ID).
		:return: Next integer ID.
		"""
		self.transIndex +=1
		return self.transIndex

	def GetChromosomeNames (self)->list:
		"""
		Returns all names of the chromosomes inside a list.
		:return: List of all chromosome names.
		"""
		countNames = len(self.dictChrNames)/2
		i = 0
		NameList = []
		while countNames != i:
			NameList.append(self.dictChrNames[i])
			i +=1
		if NameList[0] == []:
			NameList.pop(0)

		return NameList

	def GetChrTranscripts(self, ChrName: str):
		"""
		Returns a dictionary with all transcripts from the choosen chromosome.
		:param ChrName: Name of the chromosome.
		:return: Dictionary with all transcripts inside this chromosome.
		"""
		return self.dictListOfTranscripts[self.dictChrNames[ChrName]]


class GFF3_Handler_V2:
	"""
	Working on a better version, using only the ids.
	No idea, why i didn't do this before.
	"""
	def create_Transcripts(self,Chromosomflag:str):
		TransList = []
		if Chromosomflag not in self.dict_Chr_TranscriptList.keys():
			for entry in self.gff3_ID_Dict.values():
				spline_mRNA = entry[0]
				CDSlist = entry[1:]
				CDSlist = sorted(CDSlist, key=lambda CDSlist: CDSlist[0]) # sorting after startposition

				TID = spline_mRNA[8][3:spline_mRNA[8].find(";")]

				if (spline_mRNA[6] == "-"):
					ForwardOrReverse = Transcript.TranscriptEnum.REVERSE
				elif (spline_mRNA[6] == "+"):
					ForwardOrReverse = Transcript.TranscriptEnum.FORWARD
				else:
					ForwardOrReverse = Transcript.TranscriptEnum.UNKNOWN_STRAND_DIRECTION

				# def __init__ (self, IndexKey: int, TID: int, StartOfRNA: int, EndOfRNA: int, Forward: bool):
				trans = Transcript.Transcript(int(self.transIndex), str(TID), int(spline_mRNA[3]), int(spline_mRNA[4]),
											  ForwardOrReverse)
				self.transIndex += 1

				for CDS in CDSlist:
					trans.addCDS(CDS[0],CDS[1],CDS[2])

				TransList.append(trans)
			TransList = sorted(TransList, key=lambda TransList: TransList.StartOfRNA) # sorting after startposition
			self.dict_Chr_TranscriptList[Chromosomflag] = TransList
		else:
			print(str(Chromosomflag) + " liked to be two times inside the gff3 handler. " +
									"Maybe the gff3 file is not completely sorted after seqid (or some bug occurred).")

	def __init__(self, GFF3_DATA_PATH: str):
		"""
		Reads the gff3 file and creates the transcript objects.

		:param GFF3_DATA_PATH: Path to file, including file.
		"""

		self.gff3_ID_Dict = {}
		self.dict_Chr_Names = {}
		self.dict_Chr_TranscriptList = {}
		self.transIndex = 0



		gff3_file = open(GFF3_DATA_PATH)
		line = gff3_file.readline()
		chromosom_Flag = ""
		current_Chr = 0


		while line:
			#for the very first information

			if line.startswith("#"):
				line = gff3_file.readline()
				continue

			spline = line.split("\t")
			chromosom_Flag = spline[0]
			self.dict_Chr_Names[current_Chr] = chromosom_Flag
			self.dict_Chr_Names[chromosom_Flag] = current_Chr
		while line:
			if line.startswith("#"):
				line = gff3_file.readline()
				continue
			spline = line.split("\t")
			#0:	seqid(chr)
			#1:	source
			#2: type
			#3: start(int)
			#4:	end(int)
			#5:	score
			#6: strand(direction, +,-)
			#7: phase(+0,1,2 to next codon)
			#8:	attributes (warning, not always all attributes here, not even the first few)
			#8.0: ID
			#8.1: Name
			#8.2: Alias
			#8.3: Parent
			#8.4-...: Other stuff
			if spline[0] == chromosom_Flag: #same chromosome
				if "mRNA" in spline[2]:
					id = spline[8].split(";")[0] #ID=....;
					id = id.replace("\n","") # todo try id, if there is no id --> old gff version
					self.gff3_ID_Dict[id[3:]] = [spline] #ID=...
				elif "CDS" in spline[2]:
					try:
						pids = spline[8].split("Parent=")[1] #Parent=....;
						pids = pids.split(";")[0] #ignore all attributes after parents
						pids = pids.replace("\n","")
						for pid in pids.split(","): #if there are multiple parents
							self.gff3_ID_Dict[pid].append([int(spline[3]),int(spline[4]),spline[7]])#[Start,End,phase]
					except KeyError as e:
						print(e)
			else: # new chromosome
				if len(spline) > 1 :# not empty line
					self.create_Transcripts(chromosom_Flag)
					chromosom_Flag = spline[0]
					current_Chr += 1
					self.dict_Chr_Names[current_Chr] = chromosom_Flag
					self.dict_Chr_Names[chromosom_Flag] = current_Chr


				else: #empty line
					line= gff3_file.readline()
					continue
		self.create_Transcripts(chromosom_Flag)


		gff3_file.close()


	def GetNextTranscriptIndex (self):
		"""
		For creating new transcripts (need for a new ID).
		:return: Next integer ID.
		"""
		self.transIndex +=1
		return self.transIndex

	def GetChromosomeNames (self)->list:
		"""
		Returns all names of the chromosomes inside a list.
		:return: List of all chromosome names.
		"""
		countNames = len(self.dict_Chr_Names)/2
		i = 0
		NameList = []
		while countNames != i:
			NameList.append(self.dict_Chr_Names[i])
			i +=1
		if NameList[0] == []:
			NameList.pop(0)

		return NameList

	def GetChrListTranscripts(self, ChrName: str):
		"""
		Returns a list with all transcripts from the choosen chromosome.
		:param ChrName: Name of the chromosome.
		:return: List with all transcripts inside this chromosome.
		"""
		return self.dict_Chr_TranscriptList[self.dict_Chr_Names[ChrName]]


class GFF3_Handler_V3:


	def __init__(self, GFF3_DATA_PATH: str):
		self.gff3_ID_Dict = {}
		self.dict_gff_for_parents = {}
		self.dict_Chr_Names = {}
		self.dict_Chr_dict_Transcript = {}
		self.List_Of_Transcripts = {}
		self.nr_chroms = -1
		self.transIndex = -1

		self.readGFF3(GFF3_DATA_PATH)
		self.createTranscripts()
		self.sortTranscripts()

	def readGFF3(self,GFF3_DATA_PATH: str):
		gff3 = open(GFF3_DATA_PATH, 'r')
		lines = gff3.readlines()
		gff3.close()

		for line in lines:
			if line.startswith('#') or len(line) <=2:
				continue
			spline = line.split('\t')
			seqid = spline[0]
			if seqid not in self.dict_Chr_Names.keys():
				self.nr_chroms += 1
				self.dict_Chr_Names[seqid] = self.nr_chroms
				self.dict_Chr_Names[self.nr_chroms] = seqid
				self.List_Of_Transcripts[seqid] = []

			atts = spline[8].split(';')
			gff3_id = ""
			parent = ""
			for att in atts:
				if 'ID=' in att:
					gff3_id = att[3:].replace("\n","")
				elif 'Parent=' in att:
					parent = att[7:].replace("\n","")
					if ',' in parent:
						parent = parent.split(',')
					else:
						parent = [parent]
			if gff3_id == "":
				print(line)
				print("Item(s) without ID - is it really gff version 3?(or empty lines....)")
				sys.exit()
			self.gff3_ID_Dict[seqid,gff3_id] = spline

			if len(parent) == 0:
				#when its a gene entry
				continue
			for pp in parent:
				try:
					self.dict_gff_for_parents[seqid, pp].append(spline)
				except KeyError:
					self.dict_gff_for_parents[seqid, pp] = [spline]

	def createTranscripts(self):

		for key, spline in self.gff3_ID_Dict.items():

			if spline[2] != "gene":
				continue

			gene_seqid = spline[0]
			gene_gff3_ID = key[1] # because of (seqid,gff3_id)
			gene_info_string = spline[8]

			if gene_seqid not in self.dict_Chr_dict_Transcript:
				self.dict_Chr_dict_Transcript[gene_seqid] = {}


			for spline_child in self.dict_gff_for_parents[gene_seqid,gene_gff3_ID]:
				gfftype = spline_child[2]
				if 'mRNA' not in gfftype:
					continue

				start = int(spline_child[3])
				end = int(spline_child[4])
				if spline_child[6] == '+':
					strand = Transcript.TranscriptEnum.FORWARD
				elif spline_child[6] == '-':
					strand = Transcript.TranscriptEnum.REVERSE
				else:
					strand = Transcript.TranscriptEnum.UNKNOWN_STRAND_DIRECTION
				phase = spline_child[7]

				atts = spline_child[8].split(';')
				TID = ""
				for att in atts:
					if 'ID=' in att:
						TID = att[3:]
						break
				if TID == "":
					print("No ID:\n" + spline_child)
					exit()

				# IndexKey: int, TID: str, StartOfRNA: int, EndOfRNA: int, ForwardDirection: TranscriptEnum.REVERSE):
				transcript = Transcript.Transcript(self.GetNextTranscriptIndex(),
										TID,
										start,
										end,
										strand)
				cdslist = []
				for rna_child_spline in self.dict_gff_for_parents[gene_seqid,TID]:
					if "CDS" in rna_child_spline[2]:
						cdslist.append((int(rna_child_spline[3]),int(rna_child_spline[4]),rna_child_spline[7]))
					elif "exon" in rna_child_spline[2]:
						transcript.AddEXON_Descriptin("\t".join(rna_child_spline))
					elif "utr" in rna_child_spline[2]:
						transcript.AddUTR_Description("\t".join(rna_child_spline))

				cdslist = sorted(cdslist)
				for cds in cdslist:
					transcript.addCDS(cds[0],cds[1],cds[2])
				transcript.SetGene_Info_String(gene_info_string)
				self.dict_Chr_dict_Transcript[gene_seqid][transcript.IndexKey] = transcript
				self.List_Of_Transcripts[gene_seqid].append(transcript)

	def GetNextTranscriptIndex (self):
		"""
		For creating new transcripts (need for a new ID).
		:return: Next integer ID.
		"""
		self.transIndex +=1
		return self.transIndex

	def GetChromosomeNames (self)->list:
		"""
		Returns all names of the chromosomes inside a list.
		:return: List of all chromosome names.
		"""
		countNames = len(self.dict_Chr_Names)/2
		i = 0
		NameList = []
		while countNames != i:
			NameList.append(self.dict_Chr_Names[i])
			i +=1
		if NameList[0] == []:
			NameList.pop(0)

		return NameList

	def GetChromosomeID (self, ChrName:str)->int:
		return self.dict_Chr_Names[ChrName]

	def GetChrTranscriptsDict(self, ChrName: str):
		"""
		Returns a dictionary with all transcripts from the choosen chromosome.
		:param ChrName: Name of the chromosome.
		:return: Dictionary with all transcripts inside this chromosome.
		"""
		return self.dict_Chr_dict_Transcript[ChrName]

	def GetChrTranscriptList(self,ChrName:str):
		return self.List_Of_Transcripts[ChrName]


	def sortTranscripts(self):
		for chrName in self.List_Of_Transcripts.keys():
			self.List_Of_Transcripts[chrName] = sorted(self.List_Of_Transcripts[chrName], key=lambda sTranscript: sTranscript.StartOfRNA)

	def updateTransripts(self,transcriptList: list , ChrName: str):
		#i = self.dict_Chr_Names[ChrName]
		#self.List_Of_Transcripts[i] = transcriptList
		self.List_Of_Transcripts[ChrName] = sorted(transcriptList, key=lambda sTranscript: sTranscript.StartOfRNA)

	def AddNewTranscriptToDict(self, ChrName: str, transcript: Transcript):
		self.dict_Chr_dict_Transcript[ChrName][transcript.TID] = transcript

# 0:	seqid(chr)
# 1:	source
# 2: type
# 3: start(int)
# 4:	end(int)
# 5:	score
# 6: strand(direction, +,-)
# 7: phase(+0,1,2 to next codon)
# 8:	attributes (warning, not always all attributes here, not even the first few)
# 8.0: ID
# 8.1: Name
# 8.2: Alias
# 8.3: Parent
# 8.4-...: Other stuff




















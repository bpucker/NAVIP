import Transcript

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

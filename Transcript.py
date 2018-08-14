__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


from enum import Enum, unique
from VCF_Variant import Variant,VariantEnum
from LogOrganizer import LogEnums, LogOrganizer



@unique
class TranscriptEnum (Enum):
	"""
	This class supports enums for transcripts and the classification of variants inside these transcripts.
	"""
	START_NOT_IN_CDS = "Error - Startposition is not in the CDS"
	END_NOT_IN_CDS = "Error - Endposition is not in the CDS"
	POSITION_NOT_IN_CDS = -1
	FORWARD = True
	REVERSE = False
	UNKNOWN_STRAND_DIRECTION = "No Direction for this transcript available" #will hopefully never happen

	# Stuff for variants inside the transcript
	INSERTION = "INS"
	DELETION = "DEL"
	SUBSTITUTION = "SUB"
	FRAMESHIFT = "Frameshift"
	FRAMESHIFT_1 = "Frameshift+1"
	FRAMESHIFT_2 = "Frameshift+2"
	FRAMESHIFT_1_DEL = "Frameshift-1"
	FRAMESHIFT_2_DEL = "Frameshift-2"
	STOP_GAINED = "Stop gained"
	STOP_LOST = "Stop lost"
	STOP_CHANGED = "Stop changed"
	AA_CHANGE = "Amino acid change"
	START_LOST = "Start lost"

	# when the transcript structure is damaged
	CDS_START_INVOLVED = "Hits before Start into CDS"
	CDS_STOP_INVOLVED = "Hits CDS Stop and after"
	CDS_INTRON_TO_EXON = "Hits from intron to exon"
	CDS_EXON_TO_INTRON = "Hits from exon to intron"


class Variant_Information_Storage:
	"""
	A storage class for in-transcript use only. It contains all necessary variant information.
	"""

	def __init__(self, ChrPosition: int, Ref: str, Alt: str, Unchanged_CDS_Position: int, ID: str, qual:str, filter:str, info:str):
		"""
		Initialization of all necessary information.
		:param ChrPosition:  Position of the variant.
		:param Ref: VCF REF entry.
		:param Alt: VCF ALT entry.
		:param Unchanged_CDS_Position: Position inside the original CDS.
		:param ID: VCF ID entry.
		:param qual: VCF QUAL entry.
		:param filter: VCF FILTER entry.
		:param info: VCF INFO entry.
		"""
		self.ChrPosition = ChrPosition
		self.Unchanged_CDS_Position = Unchanged_CDS_Position #without variant effects
		self.Changed_CDS_Position = 0 #with variant effects
		self.Ref = Ref
		self.Alt = Alt
		self.ID = ID
		self.Qual = qual
		self.Filter = filter
		self.OLD_Info = info

		self.ReverseRef = For_Type_Safety_and_statics.ReverseSeq(Ref)
		self.ReverseAlt = For_Type_Safety_and_statics.ReverseSeq(Alt)
		if self.Unchanged_CDS_Position != TranscriptEnum.POSITION_NOT_IN_CDS:
			self.OrigRaster = For_Type_Safety_and_statics.calculateRaster(self.Unchanged_CDS_Position)
		else:
			self.OrigRaster = self.Unchanged_CDS_Position
		self.Classification = []
		self.StartOfOwnEffect = 0
		self.EndOfOwnEffect = 0
		self.OrigTriplets = ""
		self.OrigRevTriplets = ""
		self.OrigAmino = ""
		# after variant effects
		self.Changed_Raster = 0
		self.ChangedTriplets = ""
		self.ChangedRevTriplets = ""
		self.NewAmino = ""

		self.SharedEffectsWith = []

	def SetChanged_CDS_Position(self, Changed_CDS_Position: int):
		"""
		Sets the value for the changed cds position and calculates the changed raster, too.
		But only, if possible, this position is maybe outside the cds and will note it.
		:param Changed_CDS_Position: New CDS position after all variant effects.
		:return: Nothing.
		"""
		self.Changed_CDS_Position = Changed_CDS_Position
		if type(Changed_CDS_Position) == int:
			self.Changed_Raster = For_Type_Safety_and_statics.calculateRaster(Changed_CDS_Position)
		else:
			self.Changed_Raster = Changed_CDS_Position






class Transcript:
	"""
	Contains all information about a single transcript, including all variants in it's range.
	Contains all necessary functions to calculate classification of the variants inside a transcript.
	"""


	def __init__ (self, IndexKey: int, TID: str, StartOfRNA: int, EndOfRNA: int, ForwardDirection: TranscriptEnum.REVERSE):
		"""
		Initialization of the transcript-class.
		:param IndexKey: Unique number, for identifying this transcript: 0...n.
		:param TID: Not unique(tri-allele)! Transcripts identification string from the GFF3-File: For example: "AT5G40340.1"
		:param StartOfRNA: Start position inside the genome/dna-data. Not String-Position.
		:param EndOfRNA: End position inside genome/dna-data. Not String-Position.
		:param ForwardDirection: TranscriptEnum.FORWARD or TranscriptEnum.REVERSE for transcript/strand orientation.
		"""

		self.IndexKey = IndexKey
		self.TID = TID
		self.StartOfRNA = int(StartOfRNA)
		self.EndOfRNA = int(EndOfRNA)
		self.ListofCDS = [] # [(StartPos, EndPos, Raster)]
		self.ForwardDirection = ForwardDirection
		self.Complete_CDS = ""
		self.Rev_CDS = ""
		self.ListofVariants = [] # [VariantIDs]
		self.LastCDSPosition = 0
		self.Gene_Info_String = ""
		self.UTR_Description = []
		self.EXON_Descriptin = []
		self.Gene_Start_Position = 0
		self.Gene_End_Position = 0



		# for integrated variants action:
		self.IntegratedVariantObjects_CDS_Hits = []
		self.IntegratedVariantObjects_NotCDS = []
		self.IV_Changed_DNA_CDS_Seq = ""
		self.IV_Ready = False
		self.IV_NewPositionList = [0]
		self.IV_OriginalTranslation = ""
		self.IV_ChangedTranslation = ""
		self.IV_Count_Stops_in_New_AA = -1

		# state's of transcript
		self.TID_locked = False  # for changing TID because of more than one multi_allel_variants
		self.Transcript_CDS_damaged = False # exon or intron damaged
		self.CDS_Exist = False
		self.Rev_CDS_Exist = False
		self.MultiAllelVariants = False
		self.Lost_Stop = False
		self.Found_New_Stop = False

	def SetGene_Start_Position(self,Gene_Start_Position:int):
		"""
		Set the start position of the gene.
		:param Gene_Start_Position: Position in the chromosome.
		:return: Nothing.
		"""
		self.Gene_Start_Position = Gene_Start_Position
	def SetGene_End_Position(self,Gene_End_Position:int):
		"""
		Set the end position of the gene.
		:param Gene_End_Position: Position in the chromosome.
		:return: Nothing.
		"""
		self.Gene_End_Position = Gene_End_Position
	def SetGene_Info_String(self,Gene_Info_String:str):
		"""
		Normaly there exist an information string about the classification of the certain gene in the gff3 file.
		Especially hypothetical genes can make some problems, so the information can be stored for explanation,
		if something went wrong.
		:param Gene_Info_String: String description of the gene, if the gff3 file holds this information.
		:return: Nothing.
		"""
		self.Gene_Info_String = Gene_Info_String
	def AddUTR_Description(self,UTR_Description:list ):
		"""
		Sometimes the gff3 file contains information about the UTR in some genes.
		It can be usefull, if the transcript is somehow damaged from variants.
		:param UTR_Description: String from the gff3 file with the UTR description.
		:return: Nothing.
		"""
		self.UTR_Description.append(UTR_Description)
	def AddEXON_Descriptin(self,EXON_Descriptin:list):
		"""
		Sometimes the gff3 file contains information about the exon in some genes.
		It can be usefull, if the transcript is somehow damaged from variants.
		:param EXON_Descriptin: String from the gff3 file with the exon description.
		:return: Nothing
		"""
		self.EXON_Descriptin.append(EXON_Descriptin)

	def Create_IV_OriginalTranslation(self, genetic_code:dict) -> bool:
		"""
		Create the AminoAcid sequence from the original dna source without any variant effects.
		:param genetic_code: A dictionary which translates lowercase triplet RNA to AA. Example: genetic_code = {'agg':'r'}
		:return: Bool value if the transcription and translation of the DNA was successful.
		"""

		if self.ForwardDirection == TranscriptEnum.FORWARD and self.Complete_CDS != "":
			self.IV_OriginalTranslation = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(self.Complete_CDS), genetic_code)
			return True
		elif self.ForwardDirection == TranscriptEnum.REVERSE and self.Rev_CDS != "":
			self.IV_OriginalTranslation = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(self.Rev_CDS), genetic_code)
			return True
		else:
			if self.Complete_CDS != "":
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_CDS_CREATION_LOG, "No CDS: " + str(self.TID) + "\n")
				#print("No CDS: " + str(self.TID))
				return False
			elif self.Rev_CDS != "":
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_CDS_CREATION_LOG, "No Rev CDS: " + str(self.TID) + "\n")
				#print("No Rev CDS: " + str(self.TID))
				return False
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_CDS_CREATION_LOG, "No direction?: " + str(self.TID) + "\n")
				#print ("No direction?: " + str(self.TID))
				return False

	def Create_IV_ChangedTranslation (self, genetic_code: dict) -> bool:

		"""
		Create the AminoAcid sequence with the original dna source AND with variant effects inside the CDS.
		:param genetic_code: A dictionary which translates lowercase triplet RNA to AA. Example: genetic_code = {'agg':'r'}
		:return :Bool value if the transcription and translation of the DNA was successful.
		"""


		if self.ForwardDirection == TranscriptEnum.FORWARD and self.IV_Changed_DNA_CDS_Seq != "":
			self.IV_ChangedTranslation = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(self.IV_Changed_DNA_CDS_Seq), genetic_code)
			return True
		elif self.ForwardDirection == TranscriptEnum.REVERSE and self.IV_Changed_DNA_CDS_Seq != "":
			self.IV_ChangedTranslation = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(self.IV_Changed_DNA_CDS_Seq), genetic_code)
			return True
		else:
			if self.IV_Changed_DNA_CDS_Seq == "" and self.Complete_CDS == "":
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No CDS: " + str(self.TID) + "\n" )
				#print("No CDS: " + str(self.TID))
				return False
			elif self.IV_Changed_DNA_CDS_Seq == "":
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No Changed_DNA_CDS_Seq, but CDS without changes exist: " + str(self.TID) + "\n")
				#For_Type_Safety_and_statics.log.append("No Changed_DNA_CDS_Seq, but CDS without changes exist: " + str(self.TID) + "\n")
				return False
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No direction?: " + str(self.TID)+ "\n")
				#print("No direction?: " + str(self.TID))
				return False

	def Create_IV_Changed_DNA_CDS_Seq (self, genetic_code: dict, IntegratedVariantObjects: list, stopcodon: str):
		"""
		First Step:
			Simple classification for every Variant_Information_Storage object in SUB, DEl and INS.
			While doing this classification, the new CDS position for every variant will be calculated.
			In this step the new DNA sequence of the transcript will be created
		Second Step:
			 For every variant the IV_Local_Classification will be used to create the exact changes inside the CDS.
			 For every variant the IV_Local_Effect_Length will be used to calculate the length of the effects.
		Third Step.
			For all variants SetKeysToCombinedEffects will be used to add the information about side effects.
			-> After this step every variant contains the information with which other variant it shares an effect.
		:param genetic_code: The genetic code as a dictionary.
		:param IntegratedVariantObjects: A list from all variants in Variant_Information_Storage format, which all hits the CDS.
		:param stopcodon: Which letter stands for the stop codon.
		:return: Nothing.
		"""
		if self.ForwardDirection == TranscriptEnum.FORWARD and not self.IV_Ready:
			self.IV_Changed_DNA_CDS_Seq = self.Complete_CDS
			self.IV_Ready = True
		elif self.ForwardDirection == TranscriptEnum.REVERSE and not self.IV_Ready:
			self.IV_Changed_DNA_CDS_Seq = self.Rev_CDS
			self.IV_Ready = True


		IntegratedVariantObjects = sorted(IntegratedVariantObjects, key = lambda variant_information: variant_information.Unchanged_CDS_Position)
		# list needs to be ordered after position, so the new cds position can be calculated.
		# because every variants cds position can be changed with previous indels.
		for vinfo in IntegratedVariantObjects:

			if self.ForwardDirection == TranscriptEnum.FORWARD:
				alt = vinfo.Alt
				ref = vinfo.Ref
			elif self.ForwardDirection == TranscriptEnum.REVERSE:
				alt = vinfo.ReverseAlt
				ref = vinfo.ReverseRef
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Transcript without Direction: " + str(self.TID) + "\n")
				#print("Transcript without Direction: " + str(self.TID))
				continue
			cds_position = vinfo.Unchanged_CDS_Position
			current_additional_position = self.IV_NewPositionList[len(self.IV_NewPositionList)-1]
			vinfo.SetChanged_CDS_Position(cds_position + current_additional_position)
			if len(ref) == 1 and len(alt) == 1:
				firstkoord = cds_position + current_additional_position -1
				###
				# important, because seq[0:0] does not work.
				# but there can be the case, that the first position will be substituted
				###
				if firstkoord !=  0:
					first = self.IV_Changed_DNA_CDS_Seq[0:firstkoord]
				else:
					first = ""

				test = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position - 1]
				if test != ref:
					LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
																"Error 1: SUB not identical with Ref: " + str(
																	self.TID) + " " + str(vinfo.ChrPosition) + "\n")
					#print("Error 1: SUB not identical with Ref: " + str(self.TID) + " " + str(vinfo.ChrPosition))
				substitution = alt
				second = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position:]
				self.IV_Changed_DNA_CDS_Seq = first + substitution + second
				vinfo.Classification.append(TranscriptEnum.SUBSTITUTION)

				self.IV_NewPositionList.append(current_additional_position)
			elif len(ref) > 1 and len(alt) == 1: # del
				dellength = (len(ref)-1)

				if self.ForwardDirection == TranscriptEnum.REVERSE:
					if cds_position + current_additional_position + dellength <= 3 \
							and cds_position + current_additional_position + dellength  > 0:
						first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position -1 ]
						LogOrganizer.addToLog(LogEnums.TRANSCRIPT_ADDITIONAL_INFO_LOG, "Note 1: DEL removed start: " + str(self.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
						#print("Note 1: DEL removed start: " + str(self.TID) + "\t" + str(vinfo.ChrPosition))
						vinfo.Classification.append(TranscriptEnum.START_LOST)
					elif cds_position + current_additional_position + dellength <= 0:
						LogOrganizer.addToLog(LogEnums.TRANSCRIPT_ADDITIONAL_INFO_LOG,"Note 1.1: DEL destroyed start exon: " + str(self.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
						#print("Note 1.1: DEL destroyed start exon: " + str(self.TID) + "\t" + str(vinfo.ChrPosition))
						vinfo.Classification.append(TranscriptEnum.CDS_START_INVOLVED)
						first = ""
					else:
						first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position -1 ]

					second = self.IV_Changed_DNA_CDS_Seq[cds_position + dellength + current_additional_position -1:]
					if second[0] != ref[len(ref) - 1]:
						LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,"Error 2: DEL not identical with Ref: " + str(self.TID) + str(vinfo.ChrPosition) + "\n" )
						#print("Error 2: DEL not identical with Ref: " + str(self.TID) + str(vinfo.ChrPosition))

				else:
					first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position]

					second = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position + dellength:]
					if first != "":
						if first[len(first) - 1] != ref[0]:
							LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Error 3: DEL not identical with Ref: " + str(self.TID) + str(vinfo.ChrPosition) + "\n")
							#print("Error 3: DEL not identical with Ref: " + str(self.TID) + str(vinfo.ChrPosition))

				self.IV_Changed_DNA_CDS_Seq = first + second
				self.IV_NewPositionList.append(current_additional_position - (-1 + len(ref)))
				vinfo.Classification.append(TranscriptEnum.DELETION)
			elif len(ref) == 1 and len(alt) > 1: #insert
				test = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position - 1]
				first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position - 1]
				insert = alt
				second = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position:]
				self.IV_Changed_DNA_CDS_Seq = first + insert + second
				vinfo.Classification.append(TranscriptEnum.INSERTION)

				self.IV_NewPositionList.append(current_additional_position + (-1 + len(alt)))
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Create_IV_IV_Changed_DNA_CDS_Seq - error?"+ str(self.TID) + str(vinfo.ChrPosition) + "\n")
				#print("Create_IV_IV_Changed_DNA_CDS_Seq - error?"+ str(self.TID) + str(vinfo.ChrPosition) + "\n")
				self.IV_NewPositionList.append(current_additional_position)
				continue
		###
		# after every mutation is inside the new cds:
		# start with classification, its important, that every entry is inside, because of sideeffects
		#
		###
		for vinfo in IntegratedVariantObjects:
			self.IV_Local_Classification(vinfo, genetic_code, stopcodon)
			self.IV_Local_Effect_Length(vinfo)

		#IntegratedVariantObjects = self.updateIntegratedVariantObjectsList(IntegratedVariantObjects)
		self.SetKeysToCombinedEffects(IntegratedVariantObjects)

	def SetKeysToCombinedEffects(self, IntegratedVariantObjects: list):
		"""
		Every variant, which has an effect to another variant will be calculated inside this function.
		Example: A variant which will cause a frameshift get the keys from every variant after this frameshift.
				 Every variant inside this frameshift will get the key from the one it caused it.
		:param IntegratedVariantObjects: List of all classified variants. They all needs to have the effect length already.
		:return: Nothing.
		"""

		if len(IntegratedVariantObjects) > 0:

			vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(IntegratedVariantObjects[0])
			listWithVariants = [vinfo]

			#intForFrameshifts = 0
			Marked = False

			in_the_end = False
			for i in range(1, len(IntegratedVariantObjects)):
				current_vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(IntegratedVariantObjects[i])
				listed_for_deletion = []
				for old_vinfo in listWithVariants:
					if old_vinfo.StartOfOwnEffect == current_vinfo.StartOfOwnEffect:
						Marked = True
					elif old_vinfo.EndOfOwnEffect == VariantEnum.NO_EndOfOwnEffect:
						if TranscriptEnum.STOP_GAINED in old_vinfo.Classification:
							#if TranscriptEnum.STOP_GAINED in current_vinfo.Classification:
								if TranscriptEnum.DELETION in old_vinfo.Classification:
									endeffect_without_stop = old_vinfo.Changed_CDS_Position + (2 - old_vinfo.Changed_Raster)
								elif TranscriptEnum.INSERTION in old_vinfo.Classification:
									endeffect_without_stop = old_vinfo.Changed_CDS_Position + (len(old_vinfo.Alt) - 1) + (
									2 - old_vinfo.Changed_Raster)
								elif TranscriptEnum.SUBSTITUTION in old_vinfo.Classification \
										or TranscriptEnum.AA_CHANGE in old_vinfo.Classification:
									endeffect_without_stop = old_vinfo.Changed_CDS_Position + (2 - old_vinfo.Changed_Raster)
								else:
									LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Error: SetKeysToCombinedEffects:" + str(self.TID)+
										  "\t" + str(old_vinfo.ChrPosition) +
										  "\t" + str(current_vinfo.ChrPosition) + "\n")
									#print("Error: SetKeysToCombinedEffects:" + str(self.TID)+
									#	  "\t" + str(old_vinfo.ChrPosition) +
									#	  "\t" + str(current_vinfo.ChrPosition) + "\n")
									break
								#if endeffect_without_stop >= current_vinfo.StartOfOwnEffect:
								Marked = True
							#else:
							#	Marked = False
							#	in_the_end = True # after stop nothing matters
							#	break
						else:
							Marked = True
					elif old_vinfo.EndOfOwnEffect >= current_vinfo.StartOfOwnEffect:
						Marked = True
					else:
						listed_for_deletion.append(old_vinfo)
						#listWithVariants.remove(old_vinfo)
					if Marked:
						old_vinfo.SharedEffectsWith.append(current_vinfo)
						current_vinfo.SharedEffectsWith.append(old_vinfo)
						Marked = False
				if len(listed_for_deletion) > 0:
					for deleteit in listed_for_deletion:
						if deleteit == []:
							pass
						listWithVariants.remove(deleteit)
				if in_the_end: # stop == nothing else matters
					# Look 14 lines above, to reactivate this
					if len(current_vinfo.SharedEffectsWith) > 0: # if there are more effects before the stop, they have to be removed
						for vinfo_in_shared_effects in current_vinfo.SharedEffectsWith:
							vinfo_in_shared_effects.SharedEffectsWith.remove(current_vinfo)
							current_vinfo.SharedEffectsWith.remove(vinfo_in_shared_effects)

					break
				listWithVariants.append(current_vinfo)

	def RasterZero (self, dna:str, vinfo):
		return dna[vinfo.Unchanged_CDS_Position-1:vinfo.Unchanged_CDS_Position+2]
	def RasterOne (self, dna:str, vinfo):
		return dna[vinfo.Unchanged_CDS_Position-2:vinfo.Unchanged_CDS_Position+1]
	def RasterTwo (self, dna:str, vinfo):
		return dna[vinfo.Unchanged_CDS_Position-3:vinfo.Unchanged_CDS_Position]
	def RasterZeroChanged (self, dna:str, vinfo):
		return dna[vinfo.Changed_CDS_Position-1:vinfo.Changed_CDS_Position+2]
	def RasterOneChanged (self, dna:str, vinfo):
		return dna[vinfo.Changed_CDS_Position-2:vinfo.Changed_CDS_Position+1]
	def RasterTwoChanged (self, dna:str, vinfo):
		return dna[vinfo.Changed_CDS_Position-3:vinfo.Changed_CDS_Position]

	def IV_Local_Classification(self,vinfo: Variant_Information_Storage, genetic_code: dict, stopcodon: str):
		"""
		Adds (and calculates) following information to the Variant_Information_Storage object:
		OrigTriplets -> original codons for the position of the variant.
		ChangedTriplets -> codons for the NEW position in the NEW transcript.
		Classifications -> FRAMESHIFT_1; FRAMESHIFT_2; FRAMESHIFT_1_DEL; FRAMESHIFT_2_DEL;
						-> STOP_LOST; STOP_GAINED; AA_CHANGE

		:param vinfo: A single Variant_Information_Storage object.
		:param genetic_code: Dictionary of the genetic code.
		:param stopcodon: Character of the stopcodon.
		:return: Nothing.
		"""

		chrPos = vinfo.ChrPosition
		new_Raster = vinfo.Changed_Raster
		OrigRaster = vinfo.OrigRaster

		cdsPosition = vinfo.Unchanged_CDS_Position
		cdsPositionChanged = vinfo.Changed_CDS_Position

		if self.ForwardDirection == TranscriptEnum.FORWARD:
			Forward = True
			current_CDS = self.Complete_CDS
		elif self.ForwardDirection == TranscriptEnum.REVERSE:
			Forward = False
			current_CDS = self.Rev_CDS
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Error: Direction: " + str(self.TID) +" " + str(chrPos) + "\n")
			#print ("Error: Direction: " + str(self.TID) +" " + str(chrPos) + "\n")
			return False

		if TranscriptEnum.SUBSTITUTION in vinfo.Classification:



			if OrigRaster == 0:
				OrigTriplets = self.RasterZero(current_CDS, vinfo)
			elif OrigRaster == 1:
				OrigTriplets = self.RasterOne(current_CDS, vinfo)
			elif OrigRaster == 2:
				OrigTriplets = self.RasterTwo(current_CDS, vinfo)
			else:
				#print("(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(OrigRaster) + "\n")
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(OrigRaster) + "\n")
				return False

			if new_Raster == 0:
				ChangedTriplets = self.RasterZeroChanged(self.IV_Changed_DNA_CDS_Seq, vinfo)
			elif new_Raster == 1:
				ChangedTriplets = self.RasterOneChanged(self.IV_Changed_DNA_CDS_Seq, vinfo)
			elif new_Raster == 2:
				ChangedTriplets = self.RasterTwoChanged(self.IV_Changed_DNA_CDS_Seq, vinfo)
			else:
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  new_Raster) + "\n")
				return False

			if Forward:
				vinfo.OrigTriplets = OrigTriplets
				vinfo.OrigAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(OrigTriplets), genetic_code)
			else:
				vinfo.OrigTriplets = For_Type_Safety_and_statics.ReverseSeq(OrigTriplets)
				vinfo.OrigAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(OrigTriplets), genetic_code)
				vinfo.OrigAmino = vinfo.OrigAmino[::-1]

			if Forward:
				vinfo.ChangedTriplets = ChangedTriplets
				vinfo.NewAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(ChangedTriplets), genetic_code)
			else:
				vinfo.ChangedTriplets = For_Type_Safety_and_statics.ReverseSeq(ChangedTriplets)
				vinfo.NewAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(ChangedTriplets), genetic_code)
				vinfo.NewAmino = vinfo.NewAmino[::-1]

		elif TranscriptEnum.INSERTION in vinfo.Classification:
			###
			#
			# Rasterchange = (len(Alt) -1) % 3
			# Raster	Rasterchange		(Normal without change)	Before	Next
			# 0				0				2						0		2
			# 0				1				2						0		1
			# 0				2				2						0		0
			#
			# 1				0				1						1		1
			# 1				1				1						1		0
			# 1				2				1						1		2
			#
			# 2				0				0						2		0
			# 2				1				0						2		2
			# 2				2				0						2		1
			###


			RasterChange = (len(vinfo.Alt) -1) % 3

			if RasterChange == 1:
				vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_1)
			elif RasterChange == 2:
				vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_2)

			if OrigRaster == 0 :
				#before = ""
				before = current_CDS[cdsPosition-1:cdsPosition+2]
				next = ""
			elif OrigRaster == 1:
				#before = current_CDS[cdsPosition-2:cdsPosition-1]
				before = current_CDS[cdsPosition-2:cdsPosition + 1]
				next = ""
			elif OrigRaster == 2:
				#before = current_CDS[cdsPosition-3:cdsPosition-1]
				before = current_CDS[cdsPosition-3:cdsPosition]
				next = ""
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  OrigRaster) + "\n")
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				return False


			if new_Raster == 0:
				before2 = ""
			elif new_Raster == 1:
				before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged-2:cdsPositionChanged-1]
			elif new_Raster == 2:
				before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged-3:cdsPositionChanged-1]
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  new_Raster) + "\n")
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				return False

			if new_Raster == 0 and RasterChange == 0:
				next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged  + len(vinfo.Alt) - 1: cdsPositionChanged + 2 + len(vinfo.Alt) - 1]
			elif new_Raster == 0 and RasterChange == 1:
				next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged  + len(vinfo.Alt) - 1: cdsPositionChanged + 1 + len(vinfo.Alt) - 1]
			elif new_Raster == 0 and RasterChange == 2:
				next2 = ""
			elif new_Raster == 1 and RasterChange == 0:
				next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged  + len(vinfo.Alt) - 1: cdsPositionChanged + 1 + len(vinfo.Alt) - 1]
			elif new_Raster == 1 and RasterChange == 1:
				next2 = ""
			elif new_Raster == 1 and RasterChange == 2:
				next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged  + len(vinfo.Alt) - 1: cdsPositionChanged + 2 + len(vinfo.Alt) - 1]
			elif new_Raster == 2 and RasterChange == 0:
				next2 = ""
			elif new_Raster == 2 and RasterChange == 1:
				next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged  + len(vinfo.Alt) - 1: cdsPositionChanged + 2 + len(vinfo.Alt) - 1]
			elif new_Raster == 2 and RasterChange == 2:
				next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged  + len(vinfo.Alt) - 1: cdsPositionChanged + 1 + len(vinfo.Alt) - 1]
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  new_Raster) + " " + str(RasterChange) + "\n")
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				return False

			if Forward:
				vinfo.OrigTriplets = before + next
				vinfo.OrigAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(vinfo.OrigTriplets), genetic_code)
			else:
				vinfo.OrigTriplets = For_Type_Safety_and_statics.ReverseSeq(before + next)
				vinfo.OrigAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(before + next), genetic_code)
				vinfo.OrigAmino = vinfo.OrigAmino[::-1]

			if Forward:
				vinfo.ChangedTriplets = before2 + vinfo.Alt + next2
				vinfo.NewAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(vinfo.ChangedTriplets), genetic_code)
			else:
				vinfo.ChangedTriplets = For_Type_Safety_and_statics.ReverseSeq(before2 + For_Type_Safety_and_statics.ReverseSeq(vinfo.Alt) + next2)
				vinfo.NewAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(before2 + For_Type_Safety_and_statics.ReverseSeq(vinfo.Alt) + next2), genetic_code)
				vinfo.NewAmino = vinfo.NewAmino[::-1]

		elif TranscriptEnum.DELETION in vinfo.Classification:

			RasterChange = (len(vinfo.Ref) - 1) % 3

			if RasterChange == 1:
				vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_1_DEL)
			elif RasterChange == 2:
				vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_2_DEL)

			if OrigRaster == 0:
				#before = current_CDS[cdsPosition - 1:cdsPosition] # del position
				#next = current_CDS[cdsPosition + len(vinfo.Ref) -1 :cdsPosition + len(vinfo.Ref) -1  +2]
				before = ""
				next = ""
			elif OrigRaster == 1:
				#before = current_CDS[cdsPosition - 2:cdsPosition] #del position, too
				#next =  current_CDS[cdsPosition +len(vinfo.Ref) -1 :cdsPosition +len (vinfo.Ref)-1 +1]
				before = current_CDS[cdsPosition - 2:cdsPosition-1]
				next = ""
			elif OrigRaster == 2:
				before = current_CDS[cdsPosition - 3:cdsPosition-1]
				next = ""
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  OrigRaster)  + "\n")
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				return False

			if OrigRaster == 0 and RasterChange == 0:
				next = current_CDS[cdsPosition + len(vinfo.Ref) - 1: cdsPosition + len(vinfo.Ref) - 1 + 2]
			elif OrigRaster == 0 and RasterChange == 1:
				next = current_CDS[cdsPosition  + len(vinfo.Ref) - 1: cdsPosition + len(vinfo.Ref) - 1 + 1]
			elif OrigRaster == 0 and RasterChange == 2:
				next = ""
			elif OrigRaster == 1 and RasterChange == 0:
				next = current_CDS[cdsPosition + len(vinfo.Ref) - 1: cdsPosition + len(vinfo.Ref) - 1 + 1]
			elif OrigRaster == 1 and RasterChange == 1:
				next = ""
			elif OrigRaster == 1 and RasterChange == 2:
				next = current_CDS[cdsPosition + len(vinfo.Ref) - 1: cdsPosition + len(vinfo.Ref) - 1 + 2]
			elif OrigRaster == 2 and RasterChange == 0:
				next = ""
			elif OrigRaster == 2 and RasterChange == 1:
				next = current_CDS[cdsPosition + len(vinfo.Ref) - 1: cdsPosition + len(vinfo.Ref) - 1 + 2]
			elif OrigRaster == 2 and RasterChange == 2:
				next = current_CDS[cdsPosition + len(vinfo.Ref) - 1: cdsPosition + len(vinfo.Ref) - 1 + 1]
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  OrigRaster) + " " + RasterChange + "\n")
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				return False


			if new_Raster == 0:
				#before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged - 1: cdsPositionChanged + 2]
				# next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged + len(vinfo.Ref) -1: cdsPositionChanged +len(vinfo.Ref) -1 +2]
				before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged - 1: cdsPositionChanged + 2]
				next2 = ""
			elif new_Raster == 1:
				#before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged-2 : cdsPositionChanged + 1]
				#next2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged + len(vinfo.Ref) -1 : cdsPositionChanged + len(vinfo.Ref)-1 +1]
				before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged-2 : cdsPositionChanged + 1]
				next2 = ""
			elif new_Raster == 2:
				#before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged-3: cdsPositionChanged]
				before2 = self.IV_Changed_DNA_CDS_Seq[cdsPositionChanged-3: cdsPositionChanged]
				next2 = ""
			else :
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
									  "(Impossible) Raster Error: " + str(self.TID) + " " + str(chrPos) + " " + str(
										  new_Raster)  + "\n")
				#print("Error: " + str(self.TID) + " " + str(chrPos))
				return False

			if Forward:
				vinfo.OrigTriplets = before + vinfo.Ref +next
				vinfo.OrigAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(vinfo.OrigTriplets), genetic_code)
			else:
				vinfo.OrigTriplets = For_Type_Safety_and_statics.ReverseSeq(before + For_Type_Safety_and_statics.ReverseSeq(vinfo.Ref) + next)
				vinfo.OrigAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(before + For_Type_Safety_and_statics.ReverseSeq(vinfo.Ref) + next), genetic_code)
				vinfo.OrigAmino = vinfo.OrigAmino[::-1]

			if Forward:
				vinfo.ChangedTriplets = before2 + next2
				vinfo.NewAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(vinfo.ChangedTriplets), genetic_code)
			else:
				vinfo.ChangedTriplets = For_Type_Safety_and_statics.ReverseSeq(before2 + next2)
				vinfo.NewAmino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(before2 + next2), genetic_code)
				vinfo.NewAmino = vinfo.NewAmino[::-1]

		origAmino = vinfo.OrigAmino
		newAmino = vinfo.NewAmino
		vinfo.OrigRevTriplets = For_Type_Safety_and_statics.ReverseSeq(vinfo.OrigTriplets)

		if stopcodon in origAmino and not stopcodon in newAmino:
			vinfo.Classification.append(TranscriptEnum.STOP_LOST)
			self.Lost_Stop = True
			self.Calculate_Last_CDS_Position()
		elif stopcodon in newAmino and not stopcodon in origAmino:
			vinfo.Classification.append(TranscriptEnum.STOP_GAINED)
		elif stopcodon in newAmino and stopcodon in origAmino:
			if len(newAmino) == len(origAmino):
				vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
			elif TranscriptEnum.SUBSTITUTION in vinfo.Classification:
				vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
			elif TranscriptEnum.DELETION in vinfo.Classification:
				vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
			elif TranscriptEnum.INSERTION in vinfo.Classification:
				if vinfo.NewAmino[0] == stopcodon:
					vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
				else:
					vinfo.Classification.append(TranscriptEnum.STOP_LOST)
					vinfo.Classification.append(TranscriptEnum.STOP_GAINED)
		elif origAmino != newAmino:
			vinfo.Classification.append(TranscriptEnum.AA_CHANGE)

	def IV_Local_Effect_Length(self, vinfo: Variant_Information_Storage):
		"""
		Calculates the effect length and add it to a single Variant_Information_Storage (object).
		The starteffect is always a position, because every variant in here is inside the transcript.
		The endeffect can be a position and the NO_EndOfOwnEffect enum value.
		Frameshifts or other high impact effects (Stop) will effect every following variant.
		:param vinfo: A single Variant_Information_Storage object, already containing its classifications.
		:return: False, if anything went wrong.
		"""
		starteffect = vinfo.Changed_CDS_Position - vinfo.Changed_Raster
		classification = vinfo.Classification

		if TranscriptEnum.STOP_GAINED in classification or TranscriptEnum.STOP_LOST in classification or TranscriptEnum.FRAMESHIFT in classification or TranscriptEnum.FRAMESHIFT_1 in classification or TranscriptEnum.FRAMESHIFT_2 in classification:
			if TranscriptEnum.FRAMESHIFT in classification:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Still FRAMESHIFT in use." + "\n" )
			endeffect = VariantEnum.NO_EndOfOwnEffect
		elif TranscriptEnum.DELETION in classification:
			endeffect = vinfo.Changed_CDS_Position + (2-vinfo.Changed_Raster)
		elif TranscriptEnum.INSERTION in classification:
			endeffect = vinfo.Changed_CDS_Position + (len(vinfo.Alt)-1) + (2-vinfo.Changed_Raster)
		elif TranscriptEnum.SUBSTITUTION in classification or TranscriptEnum.AA_CHANGE in classification:
			endeffect = vinfo.Changed_CDS_Position + (2-vinfo.Changed_Raster)
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Error in IV_Local_Effect_Length: " + str(self.TID) + " " + str(vinfo.ChrPosition) + "\n")
			#print("Error in IV_Local_Effect_Length: " + str(self.TID) + " " + str(vinfo.ChrPosition))
			return False


		vinfo.StartOfOwnEffect = starteffect
		vinfo.EndOfOwnEffect = endeffect

	def Calculate_Last_CDS_Position(self):
		"""
		The last position of the CDS in the chromosome may be needed to extend the cds.
		:return: The last position of the CDS in the chromosome. In reverse direction it
				is of course the position on the left side.
		"""
		if self.LastCDSPosition != 0 or len(self.ListofCDS) == 0:
			return False
		if self.ForwardDirection == TranscriptEnum.FORWARD:
			self.LastCDSPosition = self.ListofCDS[len(self.ListofCDS)-1][1]
		elif self.ForwardDirection == TranscriptEnum.REVERSE:
			self.LastCDSPosition = self.ListofCDS[0][0]
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No Direction: " + str(self.TID) + "\n")
			#print ("No Direction: " + str(self.TID))

	def update_Last_CDS_Position(self):
		"""
		Updates the CDS Position inside the ListofCDS list.
		:return: Nothing.
		"""
		if self.ForwardDirection == TranscriptEnum.FORWARD:
			self.ListofCDS[len(self.ListofCDS)-1][1] = self.LastCDSPosition
		elif self.ForwardDirection == TranscriptEnum.REVERSE:
			self.ListofCDS[0][0] = self.LastCDSPosition
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No Direction: " + str(self.TID) + "\n")

	def Find_New_Stop(self, nextDNA:str , genetic_code:dict , stopcodon: str):
		"""
		Uses the IV_Changed_DNA_CDS_Seq and add the nextDNA, translates it and search for a new stopcodon.
		:param nextDNA: DNA, which will be added to the existing transcript.
		:param genetic_code: Dictionary, used for translation.
		:param stopcodon: Char, which will be used as stop.
		:return: False, if there is no new stop.
		"""
		currentRaster = For_Type_Safety_and_statics.calculateRaster(len(self.IV_Changed_DNA_CDS_Seq))
		###
		# len -1 % 3 --> last position is inside raster x
		###
		# 0 -> last position + new
		# 1 -> last -1 + last + new
		# 2 -> last -2 + last -1 + last + new -----> only new
		###
		# new stop position
		# dna to this code -< append changed_cds
		###

		if currentRaster == 0:
			dnaToTest = self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq)-1] + nextDNA
			oldDNA = 1
		elif currentRaster == 1:
			dnaToTest = self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq)-2] + self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq)-1] + nextDNA
			oldDNA = 2
		elif currentRaster == 2:
			dnaToTest = nextDNA
			oldDNA = 0
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Error in find_New_Stop: " + str(self.TID) + " " + str(currentRaster) + "\n" )
			#print ("Error in find_New_Stop: " + str(self.TID) + " " + str(currentRaster) + "\n" )
			return False

		new_Amino = For_Type_Safety_and_statics.Translation(For_Type_Safety_and_statics.Transcription(dnaToTest), genetic_code)
		position_in_string = new_Amino.find(stopcodon)

		if position_in_string != -1:
			###
			# dna_position_in_string = 3*Amino_Acid_Position_in_String
			# its alsways the start position -> Position 5 in amino == Stop is in dna_position 15,16,17
			###

			self.Lost_Stop = False
			self.Found_New_Stop = True # maybe i can use this later
			if self.ForwardDirection == TranscriptEnum.FORWARD:
				self.IV_Changed_DNA_CDS_Seq += nextDNA[0:(1 + position_in_string) * 3 - oldDNA]
				self.LastCDSPosition += (1+position_in_string)*3 -oldDNA
				self.Complete_CDS += nextDNA[0:(1+position_in_string)*3 - oldDNA]
				self.Rev_CDS += For_Type_Safety_and_statics.ReverseSeq(nextDNA[0:(1 + position_in_string) * 3 - oldDNA])
			else:
				revdna = For_Type_Safety_and_statics.ReverseSeq(nextDNA[0:(1 + position_in_string) * 3 - oldDNA])
				self.IV_Changed_DNA_CDS_Seq += revdna
				self.Rev_CDS += revdna
				self.Complete_CDS += nextDNA[0:(1 + position_in_string) * 3 - oldDNA]
				self.LastCDSPosition -= (1 + position_in_string) * 3 - oldDNA
			self.IV_Check_For_New_Variants(genetic_code, stopcodon)
		else:

			self.Rev_CDS += For_Type_Safety_and_statics.ReverseSeq(nextDNA)
			self.Complete_CDS += nextDNA
			#self.IV_Check_For_New_Variants(genetic_code, stopcodon)
			if self.ForwardDirection == TranscriptEnum.FORWARD:
				self.LastCDSPosition += len(nextDNA)
				self.IV_Changed_DNA_CDS_Seq += nextDNA
			else:
				self.LastCDSPosition -= len(nextDNA)
				self.IV_Changed_DNA_CDS_Seq += For_Type_Safety_and_statics.ReverseSeq(nextDNA)
			return False

	def resetTranscript(self):
		self.IV_Changed_DNA_CDS_Seq = ""
		self.IV_Ready = False
		self.IV_NewPositionList = [0]
		self.IV_OriginalTranslation = ""
		self.IV_ChangedTranslation = ""
		self.IV_Count_Stops_in_New_AA = -1

	def IV_Check_For_New_Variants(self, genetic_code: dict, stopcodon:str):
		"""
		Prototype function: It is more a reminder, that it is possible for longer transcripts, to have more variants.
		Not completed.
		:param genetic_code: Dictionary for the genetic code.
		:param stopcodon: Character for the stopcodon.
		:return: Nothing.
		"""
		###
		# need a new function for finding cds position in the new, longer cds.
		# current cds function needs old cds from gff3 file - maybe this can be modified
		# or the last entry can be used for this
		#
		###
		self.update_Last_CDS_Position()
		for vinfo in self.IntegratedVariantObjects_NotCDS:
			vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
			if self.ForwardDirection == TranscriptEnum.FORWARD:
				if self.SearchPositionInCDS(vinfo.ChrPosition) == TranscriptEnum.POSITION_NOT_IN_CDS:
					continue
			elif self.ForwardDirection == TranscriptEnum.REVERSE:
				if self.SearchPositionInCDSReverse(vinfo.ChrPosition) == TranscriptEnum.POSITION_NOT_IN_CDS:
					continue
			variant = Variant("not needed here",
							  vinfo.ChrPosition,
							  vinfo.ID,
							  -1,
							  vinfo.Ref,
							  vinfo.Alt,
							  vinfo.Qual,
							  vinfo.Filter,
							  vinfo.OLD_Info)
			self.IntegratedVariantObjects_NotCDS.remove(vinfo)
			self.Add_Variant_Information(variant)
		self.resetTranscript()
		self.Create_IV_Changed_DNA_CDS_Seq(genetic_code,self.IntegratedVariantObjects_CDS_Hits, stopcodon )
		"""
		chromosome: str,
		  position: int,
		  ID: int, #from vcf-file
		  usefullID: int, #not from vcf-file (count variants)
		  reference: str,
		  alternate: str,
		  qual: str,
		  filter_: str,
		  info: str) :
		"""

	def IV_CDS_Damaging_Variant_Handler(self, vinfo: Variant_Information_Storage):
		"""
		Variants, which are classified as deletions do have two different positions:
			Position of the variant and endposition of the deletion.
		If one of them is inside the CDS and one outside the transcript will be marked as damaged.
		The classification of the 'damage' will be made inside this function and added to the
		incomming Variant_Information_Storage object.
		:param vinfo: Variant_Information_Storage object with a transcript 'damaging' deletion
		:return: Nothing.
		"""
		if self.ForwardDirection == TranscriptEnum.FORWARD:
			normal_cds = vinfo.Unchanged_CDS_Position
			###
			# 4 additional cases DELETION (Forward):
			# Del-length + chrPos == Inside CDS :;: BUT ChrPos alone == Outside CDS
			#	A1-> before transcript starts, into transcript -> Start lost
			#	A2-> starts inside intron, into exon -> Spliceside destroyed
			# Del-length + chrPos == Outside CDS :;: BUT ChrPos alone == Inside CDS
			#	B1-> Start: Last-Exon; End: After-Last Exon -> Stop destroyed
			#	B2-> Start: In one Exon; End: In one Intron -> Spliceside destroyed
			###
			# While bughunting, new cases appeared:
			#  -> before transcript starts into intron -> Start + a lot more lost
			#  -> starts inside intron, into next intron -> exon or more lost + Spliceside(s) destroyed
			#  -> starts inside intron, ends after stop -> Stop destroyed + more
			#  -> starts inside exon, ends in next exon -> intron lost + Spliceside(s) destroyed

			second_position = vinfo.ChrPosition + (len(vinfo.Ref) -1)
			cds_2 = self.SearchPositionInCDS(second_position)

			firstCDS = self.ListofCDS[0]
			lastCDS = self.ListofCDS[len(self.ListofCDS)-1]

			# case A
			if normal_cds == TranscriptEnum.POSITION_NOT_IN_CDS and type(cds_2) == int:

				# case A1
				if firstCDS[0]  <= second_position <= firstCDS[1]:
					vinfo.Classification.append(TranscriptEnum.CDS_START_INVOLVED)
				#case A2
				else:
					vinfo.Classification.append(TranscriptEnum.CDS_INTRON_TO_EXON)
			#case B
			elif type(normal_cds) == int and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:

				#case B1
				if lastCDS[0] <= vinfo.ChrPosition <= lastCDS[1]:
					vinfo.Classification.append(TranscriptEnum.CDS_STOP_INVOLVED)
				#case B2
				else:
					vinfo.Classification.append(TranscriptEnum.CDS_EXON_TO_INTRON)
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,"Error - IV_CDS_Damaging_Variant_Handler?,Forward: " + str(self.TID) + "\t" + str(vinfo.ChrPosition) + "\n" )
				#print ("Error - IV_CDS_Damaging_Variant_Handler?,Forward: " + str(self.TID) + "\t" + str(vinfo.ChrPosition) + "\n")


		else:
			normal_cds = self.SearchPositionInCDSReverse(vinfo.ChrPosition)
			###
			# 4 additional cases DELETION (Reverse):
			# Del-Length + ChrPos == Inside CDS :;: BUT ChrPos alone == Outside CDS
			# 	A1 -> before transcript starts, into transcript -> Stop lost
			#   A2 -> starts inside intron, into exon -> Spliceside destroyed
			# Del-length + chrPos == Outside CDS :;: BUT ChrPos alone == Inside CDS
			#	B1 -> Start: Last-Exon(Forward Direction); End: After-Last Exon -> Start destroyed
			#	B2 -> Start: In one Exon; End: In one Intron -> Spliceside destroyed
			###
			second_position = vinfo.ChrPosition + (len(vinfo.Ref) - 1)
			cds_2 = self.SearchPositionInCDSReverse(second_position)

			firstCDS = self.ListofCDS[0] # stop in here for reverse transcripts
			lastCDS = self.ListofCDS[len(self.ListofCDS)-1] # start in here for reverse transcripts

			#case A
			if normal_cds == TranscriptEnum.POSITION_NOT_IN_CDS and type(cds_2) == int:
				#case A1
				if firstCDS[0] <= second_position <= firstCDS[1]:
					vinfo.Classification.append(TranscriptEnum.CDS_STOP_INVOLVED)
				#case A2
				else:
					vinfo.Classification.append(TranscriptEnum.CDS_INTRON_TO_EXON)
			#case B
			elif type(normal_cds) == int and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:
				#case B1
				if lastCDS[0] <= vinfo.ChrPosition <= lastCDS[1]:
					vinfo.Classification.append(TranscriptEnum.CDS_START_INVOLVED)
				#case B2
				else:
					vinfo.Classification.append(TranscriptEnum.CDS_EXON_TO_INTRON)
			else:
				LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Error - IV_CDS_Damaging_Variant_Handler?,Reverse: " + str(self.TID) + "\t" + str(vinfo.ChrPosition) +"\n")
				#print ("Error - IV_CDS_Damaging_Variant_Handler?,Reverse: " + str(self.TID) + "\t" + str(vinfo.ChrPosition) +"\n")

	def Add_Variant_Information (self, variant: Variant):
		"""
		The raw data object VCF_Variant contains all information from the original vcf file.
		But one VCF variant object can have multiple classifications and modifications in different transcripts.
		The solution of this was: A new Variant_Information_Storage object for every variant which is inside a transcript.
		--> VCF_Variant goes in
		 	-> deletion will be checked, if they delete splice sides
		 		-> IV_CDS_Damaging_Variant_Handler -> classification
			-> store and save data in a Variant_Information_Storage object for this transcript.
				-> list for variants, inside cds
				-> list for variants, outside the cds
		:param variant: VCF_Variant class object, which starting position is between the rna start <-> end (+ range)
		:return: Variant_Information_Storage.
		"""


		if self.ForwardDirection == TranscriptEnum.FORWARD:
			CDS_Position = self.SearchPositionInCDS(variant.Position)
			newEntry = Variant_Information_Storage(variant.Position,
												   variant.Reference,
												   variant.Alternate,
												   CDS_Position,
												   variant.ID,
												   variant.Qual,
												   variant.Filter,
												   variant.Info
												   )
			if len(variant.Reference) > 1:  # DEL
				cds_2 = self.SearchPositionInCDS(variant.Position + len(variant.Reference) - 1)

		else: # self.ForwardDirection == TranscriptEnum.REVERSE:
			CDS_Position = self.SearchPositionInCDSReverse(variant.Position)
			if len(variant.Reference) > 1:  # DEL
				CDS_Position = self.SearchPositionInCDSReverse(variant.Position + (len(variant.Reference) - 1))
				cds_2 = self.SearchPositionInCDSReverse(variant.Position)


			newEntry = Variant_Information_Storage(variant.Position,
												   variant.Reference,
												   variant.Alternate,
												   CDS_Position,
												   variant.ID,
												   variant.Qual,
												   variant.Filter,
												   variant.Info
												   )
		if len(variant.Reference) > 1:  # DEL

			if type(CDS_Position) == int and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:
				self.Transcript_CDS_damaged = True
				self.IntegratedVariantObjects_CDS_Hits.append(newEntry)
				self.IV_CDS_Damaging_Variant_Handler(newEntry)
			elif CDS_Position == TranscriptEnum.POSITION_NOT_IN_CDS and type(cds_2) == int:
				self.Transcript_CDS_damaged = True
				self.IntegratedVariantObjects_CDS_Hits.append(newEntry)
				self.IV_CDS_Damaging_Variant_Handler(newEntry)

		if CDS_Position == TranscriptEnum.POSITION_NOT_IN_CDS:
			self.IntegratedVariantObjects_NotCDS.append(newEntry)
		elif type(CDS_Position) == int :
			self.IntegratedVariantObjects_CDS_Hits.append(newEntry)
			return newEntry
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Possible bug in Add_Variant_Information: " + str(self.TID) +"\t"+ str(variant.Position) + "\n" )
			#print("Possible bug in Add_Variant_Information: " + str(self.TID) +"\t"+ str(variant.Position) + "\n")

	def Remove_Mult_Allel_Entry_In_All_Variant_Information(self, zero_or_one: int):
		"""
		When no vcf preprocessing is used there may exist triallel variants.
		Here they will be split into the first and second entry.
		One transcript will only get all first entrys, the other all second entrys.
		Will add an "A" or a "B" to the TID.
		:param zero_or_one: Decides if this transcript will get the first or the second entrys.
		:return: Nothing.
		"""
		for multiAllelVariant in self.IntegratedVariantObjects_CDS_Hits:
			if "," in multiAllelVariant.Alt:

				self.MultiAllelVariants = True
				multiAllelVariant.Alt = multiAllelVariant.Alt.split(",")[zero_or_one]

				if zero_or_one == 0 and not self.TID_locked:
					self.TID = str(self.TID) + "B"
					self.TID_locked = True
				elif zero_or_one == 1 and not self.TID_locked:
					self.TID = str(self.TID) + "A"
					self.TID_locked = True
				else:
					LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Maybe bug in :Remove_Mult_Allel_Entry_In_All_Variant_Information with :" +self.TID + "\n")
					#print("Maybe bug in :Remove_Mult_Allel_Entry_In_All_Variant_Information with :" +self.TID + "\n")

	def addCDS(self, StartPosition: int, EndPosition: int, Raster: str):
		"""
		Adds the CDS from the gff3 file to this transcript, one entry at the time.
		Every entry stands for an exon.
		The raster argument will not be used.
		:param StartPosition: Startposition of the exon in the chromosome.
		:param EndPosition: Endposition of the exon in the chromosome.
		:param Raster: The gff3 files often contains a raster parameter.
		:return: Nothing.
		"""
		self.ListofCDS.append(
			[int(StartPosition), int(EndPosition), str(Raster)]
		)
		if (self.ListofCDS[0] == []):
			self.ListofCDS.pop(0)

	def CompleteTheCDS (self, SeqOfTranscript: str) :
		"""
		After the CDS list is complete: This function will calculate the DNA CDS from a given sequence.
		The sequence, which will given to this function has to start at the startposition of the transcript.
		(StartOfRNA argument)
		:param SeqOfTranscript: String DNA sequence from the chromosome.
		:return: Nothing.
		"""
		completetheCDS = []
		for CDS in self.ListofCDS:
			completetheCDS.append( # short list of all CDS seqs
				SeqOfTranscript[((CDS[0] - self.StartOfRNA)):
				((abs(CDS[0]-self.StartOfRNA) + abs(CDS[1]-CDS[0])))+1])
		self.Complete_CDS = "".join(completetheCDS)
		self.Change_CDS_Existence(True)
		if (self.ListofCDS == []):
			self.Change_CDS_Existence(False)
		if (self.Complete_CDS == ""):
			self.Change_CDS_Existence(False)
		self.Calculate_Last_CDS_Position()

	def ReverseTheCDS(self) -> bool:
		"""
		This will reverse the existing Complete_CDS for reverse direction transcripts.
		:return: False, if there is no Complete_CDS. True  if it succeded.
		"""
		if self.Rev_CDS_Exist:
			return True
		if(self.Complete_CDS == ""):
			return False
		else:
			# A - C and G - T
			# A To c, C to a, G to t, T to g, then uppercase
			self.Rev_CDS = self.Complete_CDS.replace("A", "t")
			self.Rev_CDS = self.Rev_CDS.replace("C", "g")
			self.Rev_CDS = self.Rev_CDS.replace("G", "c")
			self.Rev_CDS = self.Rev_CDS.replace("T", "a")
			self.Rev_CDS = self.Rev_CDS.upper()
			self.Rev_CDS = self.Rev_CDS[::-1] #reverse #
			self.Rev_CDS_Exist = True

			return True

	def SearchPositionInCDS(self, PositionInChr: int) -> int:
		"""
		Calculates the position inside the CDS.
		Only in forward direction.
		:param PositionInChr: Position inside the chromosome.
		:return: The CDS Position or POSITION_NOT_IN_CDS enum.
		"""
		PositionInCDS = TranscriptEnum.POSITION_NOT_IN_CDS
		CDS_length = 0
		for CDS in self.ListofCDS:
			if (CDS[0] <= PositionInChr <= CDS[1]):
				PositionInCDS = CDS_length + abs(PositionInChr - CDS[0])  +1
				break #break, when you found it
			else:
				#+1, because... from position 5 to 10 are not 5 positions, it's 6: P:[5,6,7,8,9,10]
				CDS_length += abs(CDS[1] - CDS[0]) +1#EndPosition - StartPosition = length


		return PositionInCDS

	def SearchPositionPairAroundTranscript(self, StartpositionInChr:int , EndpositionInChr:int):
		"""
		CDS are numbered from 0 up. The output of this funtcion is a list with two values from -0.5 to the maximal CDS number +0.5 to
		prescind where the Positions are.
		:param StartpositionInChr: The position of the variant.
		:param EndpositionInChr: The position of the variant + length of deletion.
		:return: List of two floats.
		"""
		result = [float(0.0),float(0.0)]
		for i,CDS in enumerate(self.ListofCDS):
			if StartpositionInChr < CDS[0]:
				result[0] = float(i) - float(0.5)
			elif CDS[0] <= StartpositionInChr <= CDS[1]:
				result[0] = float(i)
			elif StartpositionInChr > CDS[1] and i -1 == len(self.ListofCDS):
				result[0] = float(i) + float(0.5)

			if EndpositionInChr < CDS[0]:
				result[1] = float(i) - float(0.5)
			elif CDS[0] <= EndpositionInChr <= CDS[1]:
				result[1] = float(i)
			elif EndpositionInChr > CDS[1] and i - 1 == len(self.ListofCDS):
				result[1] = float(i) + float(0.5)
		return result

	def SearchPositionInCDSReverse(self, PositionInChr: int) -> int:
		"""
		Calculates the position inside the CDS.
		Only in reverse direction.
		:param PositionInChr: Position inside the chromosome.
		:return: The CDS Position or POSITION_NOT_IN_CDS enum.
		"""
		PositionInCDS = TranscriptEnum.POSITION_NOT_IN_CDS
		CDS_length = 0
		for CDS in self.ListofCDS[::-1]:
			#CDS[0] == Start of CDS
			#CDS[1] == End of CDS
			#CDS is always sorted in direction of forward strand DNA ... for reverse strand DNA we have to reverse the positions, too
			# here in reverse is  CDS[1] the starting Position, not the End, but still have the highest position, because its sorted forward
			if CDS[0] <= PositionInChr <= CDS[1]:
				# End |CDS[0] - PositionInChr| + current length inside CDS
				PositionInCDS = CDS_length + abs(PositionInChr - CDS[1]) +1
				break
			else :
				CDS_length += abs(CDS[1] - CDS[0]) +1
		return PositionInCDS

	def SeqInCDS(self, StartPosChr: int, EndPosChr: int) -> str:
		"""
		Calculates a sequence between the start and end position inside the cds.
		Works only in forward strand direction.
		:param StartPosChr: Position inside the chromosome.
		:param EndPosChr: Position inside the chromosome.
		:return: The CDS position or START_NOT_IN_CDS/END_NOT_IN_CDS enums.
		"""
		StartInCDS = self.SearchPositionInCDS(StartPosChr)
		EndInCDS = self.SearchPositionInCDS(EndPosChr)

		if (StartInCDS == TranscriptEnum.POSITION_NOT_IN_CDS):
			return TranscriptEnum.START_NOT_IN_CDS
		elif (EndInCDS == TranscriptEnum.POSITION_NOT_IN_CDS):
			return TranscriptEnum.END_NOT_IN_CDS
		else:
			if (StartPosChr == EndPosChr):
				return self.Complete_CDS[StartInCDS -1]
			else :
				return self.Complete_CDS[(StartInCDS-1) : StartInCDS + abs(EndInCDS-StartInCDS)]

	def SeqInRevCDS(self, StartPosChr: int, EndPosChr: int) -> str:
		"""
		Calculates a sequence between the start and end position inside the cds.
		Works only in reverse strand direction.
		:param StartPosChr: Position inside the chromosome.
		:param EndPosChr: Position inside the chromosome.
		:return: The CDS position or START_NOT_IN_CDS/END_NOT_IN_CDS enums.
		"""
		StartInCDS = self.SearchPositionInCDSReverse(StartPosChr)
		EndInCDS = self.SearchPositionInCDSReverse(EndPosChr)

		if(StartInCDS == -1):
			return TranscriptEnum.START_NOT_IN_CDS

		elif (EndInCDS == -1):
			return TranscriptEnum.END_NOT_IN_CDS

		else:
			if (StartPosChr == EndPosChr):
				return self.Rev_CDS[StartInCDS -1]
			else :
				return self.Rev_CDS[StartInCDS-1 : StartInCDS + abs(EndInCDS-StartInCDS)]

	def SeqInStrandDirection(self, StartPosChr: int, EndPosChr: int) -> str:
		"""
		Calculates a sequence between the start and end position inside the cds.
		Works in both strand direction, using the other two direction functions.
		:param StartPosChr: Position inside the chromosome.
		:param EndPosChr: Position inside the chromosome.
		:return: The CDS position or START_NOT_IN_CDS/END_NOT_IN_CDS enums.
		"""
		if self.ForwardDirection == TranscriptEnum.FORWARD:
			return self.SeqInCDS(StartPosChr, EndPosChr)
		elif self.ForwardDirection == TranscriptEnum.REVERSE:
			return self.SeqInRevCDS(StartPosChr, EndPosChr)
		else:
			LogOrganizer.addToLog(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,"Error: No Strand Direction: " + str(self.TID) + "\n" )
			#print("Error: No Strand Direction: " + str(self.TID) + "\n")
			return False

	def Change_CDS_Existence (self, Does_it_Exist: bool):
		"""
		Setter method for the bool value CDS_Exist.
		:param Does_it_Exist:
		:return:
		"""
		self.CDS_Exist = Does_it_Exist

	def SeqInCDS_OverCDS_Position(self, StartPosCDS: int, EndPosCDS: int)-> str:
		"""
		Calculates the sequence between two CDS positions.
		Works only in forward strand direction.
		:param StartPosCDS: Position inside the CDS.
		:param EndPosCDS: Position inside the CDS.
		:return: Sequence between the two Positions. Or one single base, if start and end position are identical.
		"""
		if StartPosCDS == EndPosCDS:
			return self.Complete_CDS[StartPosCDS-1]
		else:
			return self.Complete_CDS[
				   StartPosCDS -1 :
				   StartPosCDS + abs(StartPosCDS - EndPosCDS)]
	def SeqInRevCDS_OverCDS_Position(self, StartPosCDS: int, EndPosCDS: int):
		if StartPosCDS == EndPosCDS:
			return self.Rev_CDS[StartPosCDS-1]
		else:
			return self.Rev_CDS[
				   StartPosCDS -1 :
				   StartPosCDS + abs(StartPosCDS - EndPosCDS)]
	def SeqInStrandDirectionCDS_OverCDS_Position(self,StartPosCDS:int, EndPosCDS:int):
		"""
		Calculates the sequence between two CDS positions.
		Works only in reverse strand direction.
		:param StartPosCDS: Position inside the CDS.
		:param EndPosCDS: Position inside the CDS.
		:return: Sequence between the two Positions. Or one single base, if start and end position are identical.
		"""
		if self.ForwardDirection == TranscriptEnum.FORWARD:
			return self.SeqInCDS_OverCDS_Position(StartPosCDS, EndPosCDS)
		elif self.ForwardDirection == TranscriptEnum.REVERSE:
			return self.SeqInRevCDS_OverCDS_Position(StartPosCDS, EndPosCDS)

	def SeqInIV_Changed_DNA_CDS_Seq(self, StartPosNewCDS:int, EndPosNewCDS:int) -> str:
		"""
		Calculates the sequence between two CDS positions.
		Works in both strand direction, because the IV_Changed_DNA_CDS_Seq
		will only calculated for the current direction.
		:param StartPosCDS: Position inside the CDS.
		:param EndPosCDS: Position inside the CDS.
		:return: Sequence between the two Positions. Or one single base, if start and end position are identical.
		"""
		if StartPosNewCDS == EndPosNewCDS:
			return self.IV_Changed_DNA_CDS_Seq[StartPosNewCDS-1]
		else:
			return self.IV_Changed_DNA_CDS_Seq[
				   StartPosNewCDS -1 :
				   StartPosNewCDS + abs(StartPosNewCDS - EndPosNewCDS)]

class For_Type_Safety_and_statics:
	"""
	This class exist for supporting while programming.
	It is more comfortable to use python, when the IDE knows, which type the current variables are.
	Furthermore it has a variety of methods which are used often.
	"""

	@staticmethod
	def Variant_Information_Storage_Type_Safety(vinfo: Variant_Information_Storage) -> Variant_Information_Storage:
		"""
		Only here for a trick: Get a Variant_Information_Storage type in, get a Variant_Information_Storage
		out and python and IDE will "know" its type.
		(So you can use all programming functions like usage, refactor and auto extension.)
		:param vinfo: will be returned without change
		:return: Variant_Information_Storage will be returned without change.
		"""
		return vinfo

	@staticmethod
	def Transcript_Type_Safety(transcript: Transcript)->Transcript:
		"""
		Only here for a trick: Get a Transcript type in, get a Transcript
		out and python and IDE will "know" its type.
		(So you can use all programming functions like usage, refactor and auto extension.)
		:param transcript: will be returned without change
		:return: Transcript will be returned without change
		"""
		return transcript

	@staticmethod
	def ReverseSeq (dna: str):
		"""
		It will build a biologically reversed dna strand.
		:param dna: string
		:return: biologically reversed dna string
		"""
		revdna = dna.upper()
		revdna = revdna.replace("A","t")
		revdna = revdna.replace("T","a")
		revdna = revdna.replace("C","g")
		revdna = revdna.replace("G","c")
		revdna = revdna.upper()
		revdna = revdna[::-1]
		return revdna

	@staticmethod
	def calculateRaster (VariantPositionInsideCDS: int) -> int:
		"""
		A simple calculation of the raster with the currently position inside the CDS.
		Zero == Start of the triplet
		One == In the middle of the triplet
		Two == Last Position of the triplet
		:param VariantPositionInsideCDS: self-explanatory.
		:return: 0,1 or 2 for current position inside triplet.
		"""
		return (VariantPositionInsideCDS-1 )% 3

	@staticmethod
	def Translation(rna:str, genetic_code:dict):
		"""
		Translating the RNA to AA with the given genetic code.
		Will not translate incomplete triplets, like the last 2 bp left.
		"x" for every match, which is not inside the genetic code.
		:param rna: string, lower or upper case
		:param genetic_code: dictionary, all in lowercase
		:return: amino acid string, lower case
		"""
		rna = rna.lower()
		amino = ""
		rnalen = len(rna)
		position = 0

		while ((position + 3) <= rnalen):
			try:
				amino += genetic_code[rna[position:(position + 3)]]
			except:
				amino += 'x'
			position += 3
		return amino

	@staticmethod
	def Transcription(dna:str):
		"""
		Simple replacement from 't' to 'u'.
		:param dna: string, lower or upper case
		:return: string in lowercase
		"""
		return (str(dna).lower()).replace("t", "u").upper()



#1280
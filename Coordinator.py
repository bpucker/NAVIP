__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import copy #to copy objects
from GenomeHandler import GenomeHandler, Fasta_Enum, SequenceHandlingError
from VCF_Handler import VCF_HANDLER
from VCF_Variant import Variant, VariantEnum
from Transcript import Transcript, TranscriptEnum, For_Type_Safety_and_statics
from Gff3_Handler import GFF3_Handler_V3
from datetime import datetime
from LogOrganizer import LogOrganizer, LogEnums
from snpeff_hgvs_converter import snpeff_hgvs_converter
import os
import resource
import gc

def navip_main_coordinator(invcf, ingff, infasta, outpath):
	"""
	Coordinates all functions and classes of the main module.
	Also contains a lot of functions, but naming this module "main" or "navip" isn't
	an option, because it would be confusing with the different modules.
	:param invcf: Path and name to a vcf file. Preferred the preprocessed file.
	:param ingff: Path and name to a gff3 file.
	:param infasta: Path and name to a fasta file.
	:param outpath: Path to and including the output folder.
	:return: Nothing.
	"""

	def Write_All_VCF_File(data_path: str,
						   gff3: GFF3_Handler_V3,
						   Old_Info: bool,
						   Ref_Codons: bool,
						   Ref_AA: bool,
						   Alt_Codons: bool,
						   Alt_AA: bool,
						   extend_file:bool,
						   chrName:str):
		"""
		In this function the VCF-File with the additional neighbourhood awareness information will be written.
		Furthermore it is possible (and standard) to write a few more Information into the NAVIP info:
		TranscriptID, strand direction, classifications, keys with shared effect, ref and alt codon
		with the variant position inside this codon. This information will be written into every new VCF file, too.
		Transcripts without variants inside the CDS, without CDS and those, which are (possible) damaged within the
		splicing structure will be skipped.

		:param data_path: Path into the existing output folder.
		:param gff3: The gff3 data structure with all transcripts and their variants/sequences.
		:param Old_Info: The old VCF info column will be placed after the NAVIP info, if this is true.
		:param Ref_Codons: The REF DNA codons will be written, if true.
		:param Ref_AA: The REF AA will be written, if true.
		:param Alt_Codons: The ALT codons will be written, if true.
		:param Alt_AA: The ALT AA will be written, if true.
		:return: Nothing.
		"""
		# chrom	Pos	ID	Ref	Alt Info
		if not extend_file:
			vcf_file = open(data_path + "All_VCF" + ".vcf", "w")
			vcf_file.write("##Chrom\tPos\tID\tRef\tAlt\tQual\tFilter\tInfo\n")
			vcf_file.write("##NAVIP: All Data Output\n")
			vcf_file.write("##Please note, that the Variant_Position_in_Codon is read from left to right in forward "
						   "and rigth to left in reverse strand direction.\n")
			vcf_file.write("##If there are no shared effect keys, the value is:\"NONE\".\n")
			vcf_file.write("##Info=<ID=NAV1, Type=String,Number=.,Values=[TranscriptID|"
						   "Strand_Direction|"
						   "Variant_Classification1,Variant_Classification2,...|"
						   "Shared_EffKey(s)|"
						   "REF_Codon(s)/Variant_Position_in_Codon|"
						   "REF_AA|"
						   "old_CDS_Position|"
						   "ALT_Codon(s)/Variant_Position_in_Codon|"
						   "ALT_AA|"
						   "new_CDS_Position]>\n")

			vcf_file.write("##Info=<ID=NAV2, Type=String,Number=.,Values=[ALT|ANNOTATIONS|ANNOTATION_IMPACT|GENE_NAME|"
						   "GENE_ID|FEATURE_TYPE|FEATURE_ID|TRANSCRIPT_BIOTYPE|RANK|HGVS_C|HGVS_P|cDNA_pos/cDNA_length|"
							"CDS_pos/CDS_DNA_length|AA_pos/AA_length|distance|errors_warnings]>\n")
		else:
			vcf_file = open(data_path + "All_VCF" + ".vcf", "a")
		data_to_write = []
		infoline_parser = snpeff_hgvs_converter
		#for name in gff3.GetChromosomeNames():
		name = chrName
		for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
			transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)

			if not transcriptHier.CDS_Exist or transcriptHier.Transcript_CDS_damaged:
				#no CDS or transcript exons are damaged
				continue
			if transcriptHier.IntegratedVariantObjects_CDS_Hits == []:
				#no variants inside the cds
				continue

			if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
				direction = "FOR"
			elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
				direction = "REV"
			else:
				#will never happen
				LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Error:NOT Forward, NOT Reverse:" + str(transcriptHier.TID) + "\n")
				direction = "Error:NOT Forward, NOT Reverse."

			for vinfo in transcriptHier.IntegratedVariantObjects_CDS_Hits:
				vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
				if vinfo.ChrPosition == 5921700 and transcriptHier.TID == "Ma09_t08910.1":
					print("bugsearch")
				if vinfo.ChrPosition == 3917462 and transcriptHier.TID == "Ma00_t01100.1":
					print("bugsearch")
				try:
					snpeff_string = infoline_parser.convert_main(infoline_parser,transcriptHier,vinfo)
				except:
					snpeff_string = ""
					LogOrganizer.addToLog(LogEnums.COORDINATOR_VARIANT_LOG,
										  str(vinfo.ChrPosition) + "\t" + str(currentTranscript.TID))
					print('Critical error with variant: ' + str(vinfo.ChrPosition) + '\t' + str(transcriptHier.TID))
				classificationstring = ""
				class_list = vinfo.Classification
				class_list_length = len(class_list)
				for i in range(0, class_list_length):
					classificationstring += class_list[i].value
					if class_list[i] == TranscriptEnum.STOP_CAUSED_IN:
						classificationstring += str(vinfo.STOP_CAUSED_IN)
					if i != (class_list_length - 1):  # no tab after last entry
						classificationstring += ","
					elif i == (class_list_length - 1):
						classificationstring += "|"
					else:
						LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Bug Write_All_VCF_File: classification list length:" + str(transcriptHier.TID) + "\n")

				shared_effect_list = ""
				if len(vinfo.SharedEffectsWith) > 0:
					shared_effect_list = ""  # "\tShared Effects with:"
					shared_effect_list_len = len(vinfo.SharedEffectsWith)
					#creates the string with all keys
					for i in range(0, shared_effect_list_len):
						vinfo_effects = vinfo.SharedEffectsWith[i]
						vinfo_effects = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo_effects)
						shared_effect_list += str(vinfo_effects.ChrPosition)
						if i != (shared_effect_list_len - 1):
							shared_effect_list += ","
						#seperate entrys
						elif i == shared_effect_list_len - 1:
							shared_effect_list += "|"
						#after last entry
						else:
							LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Bug Write_All_VCF_File: shared_effect_list:" + str(transcriptHier.TID) + "\n")

				if len(vinfo.OrigTriplets) % 3 != 0:
					bug = "OrigTriplets: " + str(vinfo.OrigTriplets) \
						  + "\tFrom: " \
						  + transcriptHier.TID \
						  + ": " \
						  + str(vinfo.ChrPosition) \
						  + " " \
						  + classificationstring
					LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,bug + "\n")
				elif len(vinfo.ChangedTriplets) % 3 != 0:
					bug = "ChangedTriplets: " + str(vinfo.ChangedTriplets) \
						  + "\tFrom: " \
						  + transcriptHier.TID \
						  + ": " \
						  + str(vinfo.ChrPosition) \
						  + " " \
						  + classificationstring
					LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,bug + "\n")
				# Info:TranscriptID|"
				#   "Strand_Direction|"
				#   "Variant_Classification1,Variant_Classification2,...|"
				#   "Shared_EffKey(s)|"
				#   "REF_Codon(s);Variant_Position_in_Codon|"
				#   "REF_AA|"
				#   "ALT_Codon(s);Variant_Position_in_Codon|"
				#   "ALT_AA|"
				#   "NAVIP_END|"
				#   "<old info field>"
				NAVIP_INFO = ""
				NAVIP_INFO_LIST = ["NAV1=" + str(transcriptHier.TID) + "|", str(direction) + "|", classificationstring]
				if shared_effect_list == "":
					NAVIP_INFO_LIST.append("NONE|")
				else:
					NAVIP_INFO_LIST.append(str(shared_effect_list))
				if Ref_Codons:
					NAVIP_INFO_LIST.append(str(vinfo.OrigTriplets) + "/" + str(vinfo.OrigRaster) + "|")
				if Ref_AA:
					NAVIP_INFO_LIST.append(str(vinfo.OrigAmino) + "|")
				NAVIP_INFO_LIST.append(str(vinfo.Unchanged_CDS_Position) + "|")
				if Alt_Codons:
					NAVIP_INFO_LIST.append(str(vinfo.ChangedTriplets) + "/" + str(vinfo.Changed_Raster) + "|")
				if Alt_AA:
					NAVIP_INFO_LIST.append(str(vinfo.NewAmino) + "|")
				NAVIP_INFO_LIST.append(str(vinfo.Changed_CDS_Position) + ";")
				if Old_Info:
					if "\n" in vinfo.OLD_Info:
						NAVIP_INFO_LIST.append(str(vinfo.OLD_Info))
					else:
						NAVIP_INFO_LIST.append(str(vinfo.OLD_Info) + "\n")

				NAVIP_INFO = "".join(NAVIP_INFO_LIST)
				NAVIP_INFO_LIST = [] # clear

				data_to_write.append(
					str(name) + "\t"
					+ str(vinfo.ChrPosition) + "\t"
					+ str(vinfo.ID) + "\t"
					+ str(vinfo.Ref) + "\t"
					+ str(vinfo.Alt) + "\t"
					+ str(vinfo.Qual) + "\t"
					+ str(vinfo.Filter) + "\t"
					+ snpeff_string
					+ str(NAVIP_INFO))
		data_to_write = sorted(data_to_write,
							   key=lambda data_line: (int(data_line.split("\t")[1]), str(data_line.split("\t")[7])))
		vcf_file.write("".join(data_to_write))
		data_to_write = []
		vcf_file.close()

	def Create_Sequences(gff3: GFF3_Handler_V3, Orig_AA: bool, New_AA: bool, genetic_code: dict, chrName:str):
		"""
		This function calculates the old and new AA CDS to be sure all transcripts are having the AA sequences.
		There are a few possibilities to this point, which can lead to skip a few or all transcripts AA sequences.
		Transcripts without CDS are excluded. Transcripts with a damaged CDS do not get the new AA sequence.
		:param gff3: The gff3 data structure with all transcripts and their variants/sequences.
		:param Orig_AA: Will create the original AA sequence, if true.
		:param New_AA: Will create the new AA sequence, if true.
		:param genetic_code: The dictionary with the genetic code.
		:return: Nothing.
		"""
		#for chrName in gff3.GetChromosomeNames():
		#print(chrName)
		for currentTranscript in gff3.GetChrTranscriptsDict(chrName).values():
			try:
				currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)

				if not currentTranscript.CDS_Exist:
					continue
				if Orig_AA:
					currentTranscript.Create_IV_OriginalTranslation(genetic_code)

				if currentTranscript.Transcript_CDS_damaged:
					continue
				if New_AA:
					currentTranscript.Create_IV_ChangedTranslation(genetic_code)
			except:
				LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG,"\t" + str(currentTranscript.TID))

	def Complete_Check(gff3: GFF3_Handler_V3, ghandler: GenomeHandler, chrName:str):
		"""
		Its always a good idea to test the data in the end.
		Will print warnings, if there is incorrect data.
		:param gff3: The gff3 data structure with all transcripts and their variants/sequences.
		:param ghandler: The reference genome sequences for every chromosome is inside.
		:return: Nothing.
		"""
		#for chrName in Shared_Chromosomes_FA_GFF3:
		Chr_Transcript_List = gff3.GetChrTranscriptsDict(chrName).values()
		for currentTranscript in Chr_Transcript_List:
			currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)
			if not currentTranscript.CDS_Exist or currentTranscript.Transcript_CDS_damaged:
				continue

			for vinfo in currentTranscript.IntegratedVariantObjects_CDS_Hits:
				vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
				try :
					if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
						# Variant <-> Fasta Check
						if vinfo.Ref != ghandler.singleSeq(chrName, vinfo.ChrPosition):
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "Sub_REF != Fasta-Seq: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
							#print("Sub_REF != Fasta-Seq: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

						# Variant <-> Transcript with Direction Check
						if currentTranscript.ForwardDirection == TranscriptEnum.FORWARD:
							if vinfo.Ref != currentTranscript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "Sub_REF != CDS-Seq 1: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
								#print(
								#	"Sub_REF != CDS-Seq 1: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

							if vinfo.Ref != currentTranscript.SeqInCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "Sub_REF != CDS-Seq 2: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
								#print(
								#	"Sub_REF != CDS-Seq 2: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

							if vinfo.Alt != currentTranscript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_ALT != CDS-Seq 3: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n" )
								#print(
								#	"Sub_ALT != CDS-Seq 3: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

						elif currentTranscript.ForwardDirection == TranscriptEnum.REVERSE:
							if vinfo.Ref != currentTranscript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition) and not currentTranscript.CDS_Changed:
								# When the used CDS gets changed, then it is of cause not possible to check against it.
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 1.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 1.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 2.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 2.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")
							#print(str(vinfo.ChrPosition) + "\t" + str(currentTranscript.TID) )
							if vinfo.ReverseAlt != currentTranscript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_Alt != CDS-Seq 3.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_Alt != CDS-Seq 3.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 4.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 4.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 5.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 5.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

					elif TranscriptEnum.INSERTION in vinfo.Classification:
						# Variant <-> Fasta Check
						if vinfo.Ref != ghandler.singleSeq(chrName, vinfo.ChrPosition):
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != Fasta-Seq: " + str(currentTranscript.TID) + "\t" + str(
								vinfo.ChrPosition) + "\n" )
							#print("INSERTION_REF != Fasta-Seq: " + str(currentTranscript.TID) + "\t" + str(
							#	vinfo.ChrPosition) + "\n")

						# Variant <-> Transcript with Direction Check
						if currentTranscript.ForwardDirection == TranscriptEnum.FORWARD:
							if vinfo.Ref != currentTranscript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 1: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("INSERTION_REF != CDS-Seq 1: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.Ref != currentTranscript.SeqInCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_REF != CDS-Seq 2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition))

							if vinfo.Alt != currentTranscript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
													vinfo.Changed_CDS_Position + len(vinfo.Alt) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_ALT != CDS-Seq 3: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n" )
								#print("INSERTION_ALT != CDS-Seq 3: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition)+ "\n")

						elif currentTranscript.ForwardDirection == TranscriptEnum.REVERSE:
							if vinfo.Ref != currentTranscript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition)and not currentTranscript.CDS_Changed:
								# When the used CDS gets changed, then it is of cause not possible to check against it.
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 1.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_REF != CDS-Seq 1.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 2.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_REF != CDS-Seq 2.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseAlt != currentTranscript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
													vinfo.Changed_CDS_Position + len(vinfo.Alt) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_Alt != CDS-Seq 3.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_Alt != CDS-Seq 3.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "INSERTION_REF != CDS-Seq 4.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("INSERTION_REF != CDS-Seq 4.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition))

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "INSERTION_REF != CDS-Seq 5.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("INSERTION_REF != CDS-Seq 5.2: " + str(currentTranscript.TID) + "\t" + str(
								#	vinfo.ChrPosition))

					elif TranscriptEnum.DELETION in vinfo.Classification:

						# Variant <-> Fasta Check
						try:
							if vinfo.Ref != ghandler.seq(chrName, vinfo.ChrPosition, vinfo.ChrPosition + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
									"DELETION_REF != Fasta-Seq: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
						except SequenceHandlingError as she:
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
												  str(she.description) + str(
													  currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
							if vinfo.Ref != she.sequence_part:
								#ghandler.seq(chrName, vinfo.ChrPosition,#vinfo.ChrPosition + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
													  "DELETION_REF != Fasta-Seq: " + str(
														  currentTranscript.TID) + "\t" + str(
														  vinfo.ChrPosition) + "\n")

						# Variant <-> Transcript with Direction Check
						if currentTranscript.ForwardDirection == TranscriptEnum.FORWARD:
							if vinfo.Ref != currentTranscript.SeqInCDS(
									vinfo.ChrPosition,
													vinfo.ChrPosition + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 1: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.Ref != currentTranscript.SeqInCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
													vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.Alt != currentTranscript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_ALT != CDS-Seq 3: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

						elif currentTranscript.ForwardDirection == TranscriptEnum.REVERSE:
							if vinfo.Ref != currentTranscript.SeqInCDS(
									vinfo.ChrPosition,
													vinfo.ChrPosition + len(vinfo.Ref) - 1)and not currentTranscript.CDS_Changed:
								# When the used CDS gets changed, then it is of cause not possible to check against it.
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 1.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
													vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 2.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")


							if vinfo.ReverseAlt != currentTranscript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_Alt != CDS-Seq 3.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS(
											vinfo.ChrPosition + (len(vinfo.Ref) - 1),
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 4.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.ReverseRef != currentTranscript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
													vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 5.2: " + str(currentTranscript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

					else:
						LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"No Classification in Complete_Check:\t" + str(currentTranscript.TID) +"\t" + str(vinfo.ChrPosition)+ "\n")
				except :
					LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_CRITICAL_LOG, str(vinfo.ChrPosition) + "\t" + str(currentTranscript.TID))
	def Write_All_Fasta(data_path: str,
						data_name: str,
						date_time: bool,
						Orig_DNA: bool,
						Orig_AA: bool,
						New_DNA: bool,
						New_AA: bool,
						gff3: GFF3_Handler_V3,
						extend_file: bool,
						chrom_name:str):
		"""
		Writes a new fasta file with all CDS sequences from every transcript,
		including possible spliced variant of the transcripts (GFF3 data).
		It is possible to choose, which sequences will be written.
		Tags for the type of every sequence exist:
		ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == New Amino Acid Sequence
		Furthermore: Possible damaged transcripts will be seperated into "<data_name>_damaged.txt".
		Including its information about the gene, utr, exons and all available variant information.
		:param data_path: Path to the folder for the new fasta file.
		:param data_name: Name of the new fasta file.
		:param date_time: The date and time will be included inside the filename, if true.
		:param Orig_DNA:  Writes the original CDS DNA sequence, if true.
		:param Orig_AA: Writes the original CDS A sequence, if true.
		:param New_DNA: Writes the new CDS DNA sequence, if true.
		:param New_AA: Writes the new CDS AA sequence, if true.
		:param gff3: The gff3 data structure with all transcripts and their variants/sequences.
		:return: Nothing.
		"""

		if Orig_DNA or Orig_AA or New_DNA or New_AA:
			if data_name == "":
				data_name = "all_transcripts_data"
			if date_time:
				timestop = str(datetime.now()) + "_"
			else:
				timestop = ""
			if not extend_file:

				normal_data_list_to_write = [
					"#ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == New Amino Acid Sequence\n"]
				Transcript_CDS_damaged_data_list = [
					"#Data in here may be incomplete and incorrect because of deletions, which destroyed parts the splice sides.\n"]
			else:
				normal_data_list_to_write = []
				Transcript_CDS_damaged_data_list = []
			#for chrName in gff3.GetChromosomeNames():
			chrName = chrom_name
			for currentTranscript in gff3.GetChrTranscriptsDict(chrName).values():
				currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)
				try:
					if not currentTranscript.CDS_Exist:
						continue

					if Orig_DNA:
						if currentTranscript.ForwardDirection == TranscriptEnum.FORWARD:
							normal_data_list_to_write.append(
								">" + str(currentTranscript.TID) + "|" + str(Fasta_Enum.ODNA.value) + "\n")
							normal_data_list_to_write.append(str(currentTranscript.Complete_CDS) + "\n")
						elif currentTranscript.ForwardDirection == TranscriptEnum.REVERSE:
							normal_data_list_to_write.append(
								">" + str(currentTranscript.TID) + "|" + str(Fasta_Enum.ODNA.value) + "\n")
							normal_data_list_to_write.append(str(currentTranscript.Rev_CDS) + "\n")
						else:
							LogOrganizer.addToLog(LogEnums.COORDINATOR_FASTA_FILE_ERROR_LOG,"No Direction:Write_All_Fasta: " + str(currentTranscript.TID) + "\n")
					if Orig_AA:
						normal_data_list_to_write.append(
							">" + str(currentTranscript.TID) + "|" + Fasta_Enum.OAA.value + "\n")
						normal_data_list_to_write.append(str(currentTranscript.IV_OriginalTranslation) + "\n")

					if currentTranscript.Transcript_CDS_damaged:
						Transcript_CDS_damaged_data_list.append(
							str(currentTranscript.TID) + " " + str(currentTranscript.Gene_Info_String))
						for utr in currentTranscript.UTR_Description:
							if "\n" in utr:
								Transcript_CDS_damaged_data_list.append(utr)
							else:
								Transcript_CDS_damaged_data_list.append(str(utr) + "\n")
						for exon in currentTranscript.EXON_Descriptin:
							if "\n" in exon:
								Transcript_CDS_damaged_data_list.append(exon)
							else:
								Transcript_CDS_damaged_data_list.append(str(exon) + "\n")
						for vinfo in currentTranscript.IntegratedVariantObjects_CDS_Hits:
							# vinfo = For_Type_Safety.Variant_Information_Storage_Type_Safety(vinfo)
							classificationstring = ""
							class_list = vinfo.Classification
							class_list_length = len(class_list)
							for i in range(0, class_list_length):
								classificationstring += class_list[i].value
								if i != (class_list_length - 1):  # no tab after last entry
									classificationstring += "\t"
							Transcript_CDS_damaged_data_list.append(
								str(vinfo.ChrPosition) + "\t" +
								str(vinfo.Ref) + "\t" +
								str(vinfo.Alt) + "\t" +
								str(classificationstring) + "\n")
						continue
					# Meh more DNA/AA Sequences, thanks to Stop-lost transcription length changes
					normal_data_list_to_write.append(
						">" + str(currentTranscript.TID) + "|" + Fasta_Enum.ouchDNA.value + "\n")
					normal_data_list_to_write.append(str(currentTranscript.uChDNAsequence) + "\n")

					normal_data_list_to_write.append(
						">" + str(currentTranscript.TID) + "|" + Fasta_Enum.ouchAA.value + "\n")
					normal_data_list_to_write.append(str(currentTranscript.uChAAsequence) + "\n")

					if New_DNA:
						normal_data_list_to_write.append(
							">" + str(currentTranscript.TID) + "|" + str(Fasta_Enum.nDNA.value) + "\n")
						normal_data_list_to_write.append(str(currentTranscript.IV_Changed_DNA_CDS_Seq) + "\n")
					if New_AA:
						normal_data_list_to_write.append(
							">" + str(currentTranscript.TID) + "|" + str(Fasta_Enum.nAA.value) + "\n")
						normal_data_list_to_write.append(str(currentTranscript.IV_ChangedTranslation) + "\n")
				except:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG,"Write_All_Fasta" + "\t" + str(currentTranscript.TID) + "\n")

			normal_output = "".join(normal_data_list_to_write)
			if not extend_file:
				New_Fasta_File = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "w")
				New_Fasta_File.write(normal_output)
				New_Fasta_File.close()
			else:
				New_Fasta_File = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "a")
				New_Fasta_File.write(normal_output)
				New_Fasta_File.close()

			normal_output = ""
			normal_data_list_to_write = []
			if not extend_file:
				damaged_output = "".join(Transcript_CDS_damaged_data_list)
				New_Fasta_File_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt", "w")
				New_Fasta_File_damaged.write(damaged_output)
				New_Fasta_File_damaged.close()
			else:
				damaged_output = "".join(Transcript_CDS_damaged_data_list)
				New_Fasta_File_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt", "a")
				New_Fasta_File_damaged.write(damaged_output)
				New_Fasta_File_damaged.close()
			damaged_output = ""
			Transcript_CDS_damaged_data_list = []


	# * for stop.
	genetic_code = {'agg': 'r', 'aga': 'r', 'agc': 's', 'agu': 's',
					'aag': 'k', 'aaa': 'k', 'aac': 'n', 'aau': 'n',
					'acg': 't', 'aca': 't', 'acc': 't', 'acu': 't',
					'aug': 'm', 'aua': 'i', 'auc': 'i', 'auu': 'i',

					'cgg': 'r', 'cga': 'r', 'cgc': 'r', 'cgu': 'r',
					'cag': 'q', 'caa': 'q', 'cac': 'h', 'cau': 'h',
					'ccg': 'p', 'cca': 'p', 'ccc': 'p', 'ccu': 'p',
					'cug': 'l', 'cua': 'l', 'cuc': 'l', 'cuu': 'l',

					'ugg': 'w', 'uga': '*', 'ugc': 'c', 'ugu': 'c',
					'uag': '*', 'uaa': '*', 'uac': 'y', 'uau': 'y',
					'ucg': 's', 'uca': 's', 'ucc': 's', 'ucu': 's',
					'uug': 'l', 'uua': 'l', 'uuc': 'f', 'uuu': 'f',

					'ggg': 'g', 'gga': 'g', 'ggc': 'g', 'ggu': 'g',
					'gag': 'e', 'gaa': 'e', 'gac': 'd', 'gau': 'd',
					'gcg': 'a', 'gca': 'a', 'gcc': 'a', 'gcu': 'a',
					'gug': 'v', 'gua': 'v', 'guc': 'v', 'guu': 'v', }
	stopcodon = "*"
	###
	# Output fasta file
	###
	Orig_DNA = True
	Orig_AA = True
	New_DNA = True
	New_AA = True
	Fasta_Data_Name = ""

	###
	# Output VCF file
	###
	Old_Info = True
	Ref_Codons = True
	Ref_AA = True
	Alt_Codons = True
	Alt_AA = True
	###
	# How many chromosomes will be analysed.
	# 0 stands for all.
	###
	firstxChr = 0
	timeStart = datetime.now()


	vcf_path_and_name = invcf
	gff3_path_and_name = ingff
	fasta_FILE_PATH = infasta
	Output_Data_Path = outpath

	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
	#######################################
	print("read vcf")
	vcf = VCF_HANDLER(vcf_path_and_name, firstxChr)
	print("Done: " + str(datetime.now() - timeStart))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
	#######################################
	print("read gff3")
	#print("Hopefully sorted after seqID")
	gff3 = GFF3_Handler_V3(gff3_path_and_name)
	print("Done: " + str(datetime.now() - timeStart))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
	#######################################
	print("read fa")
	ghandler = GenomeHandler(fasta_FILE_PATH, firstxChr)
	print("Done: " + str(datetime.now() - timeStart))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
	#######################################
	print("Chromosome:Names")
	if len(gff3.GetChromosomeNames()) < 10: print("GFF3: " + str(gff3.GetChromosomeNames()))
	else: print("GFF3: " + str(len(gff3.GetChromosomeNames())) + " contigs/chromosomes.")
	if len(vcf.GetChromosomeNames()) < 10: print("VCF: " + str(vcf.GetChromosomeNames()))
	else: print("GFF3: " + str(len(vcf.GetChromosomeNames())) + " contigs/chromosomes.")
	if len(ghandler.GetChromosomeNames()) < 10: print("fasta: " + str(ghandler.GetChromosomeNames()))
	else: print("fasta:" + str(len(ghandler.GetChromosomeNames())) + " contigs/chromosomes.")
	Shared_Chromosomes = []
	Shared_Chromosomes_FA_GFF3 = []
	#######################################


	for name in gff3.GetChromosomeNames():
		if name in vcf.GetChromosomeNames() and name in ghandler.GetChromosomeNames():
			Shared_Chromosomes.append(name)
		if name in ghandler.GetChromosomeNames():
			Shared_Chromosomes_FA_GFF3.append(name)
	if len(Shared_Chromosomes) < 10: print(print("Shared_Chromosomes(all):" + str(Shared_Chromosomes)))
	else: print("Shared_Chromosomes(all):" + str(len(Shared_Chromosomes)) + " contigs/chromosomes.")

	#print("Shared_Chromosomes, fasta, gff3 (and used):" + str(Shared_Chromosomes_FA_GFF3))
	### For Writing the CDS
	### in all transcripts

	extend_File = False

	for name in Shared_Chromosomes:
		print(name)
		print("Write the CDS and Rev CDS")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

		#chrID = gff3.GetChromosomeID(name)

		#for i in gff3.dictListOfTranscripts[gff3.dictChrNames[name]]:  # only a number, not the object
		#	trancripts[i].CompleteTheCDS(ghandler.seq(name, trancripts[i].StartOfRNA, trancripts[i].EndOfRNA))
		#	if trancripts[i].ForwardDirection == TranscriptEnum.REVERSE:
		#		trancripts[i].ReverseTheCDS()

		for trancript in gff3.GetChrTranscriptsDict(name).values():
			try:
				trancript.CompleteTheCDS(ghandler.seq(name, trancript.StartOfRNA, trancript.EndOfRNA), genetic_code)
			except SequenceHandlingError as she:
				LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG, str(trancript.TID) + "\t"
									  + str(trancript.StartOfRNA) + "\t" + str(trancript.EndOfRNA) + "\t" + str(she.description))
				if she.sequence_part != "":
					trancript.CompleteTheCDS(she.sequence_part, genetic_code)
			if trancript.ForwardDirection == TranscriptEnum.REVERSE:
				trancript.ReverseTheCDS(genetic_code)
		#print("Done: " + str(datetime.now() - timeStart))
		###
		# Connecting all transcripts with all their mutations/variants.
		# Contains handling of trialllele variants, if an original vcf file was used.
		# But it can't handle all possible appearances. VCF preprocessing exist for a reason.
		# The connection triggers a first evaluation of the variants, too.
		###
		transcriptRange = 300 # base pairs before RNA starts and after RNA ends. Just in Case for long InDels.
		variant_transcript_exceeding_warning_int = transcriptRange / 100
		phasingWarning = True
		print("Connect transcripts with variants")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		#print("Search in " + name)

		phasingVariants = [] # should normally be empty
		currentTranscriptID = 0 # initialization, first Round == 0
		firstTranscriptMatch = 0 #
		currentTranscriptList = gff3.GetChrTranscriptList(name)
		transcriptsWithMultiAllelVariants = []
		transcriptsWithoutStop = []


		for variant in vcf.GetChr_VCF_Variant_List(name):
			found = True
			currentTranscriptID = firstTranscriptMatch

			while (True):
				###
				# The Following case: no more available transcripts, because the position of the mutation is after the last transcript
				###
				if currentTranscriptID == len(currentTranscriptList):
					# print("No more transcripts.")
					break
				currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscriptList[currentTranscriptID])
				if (
							currentTranscript.StartOfRNA - transcriptRange) > variant.Position:  # transcriptstart higher than variant.position
					break
				elif (currentTranscript.StartOfRNA - transcriptRange) <= \
						variant.Position \
						<= (currentTranscript.EndOfRNA + transcriptRange):  # between start and end (+ range)

					currentTranscript.Add_Variant_Information(variant)

					if "," in variant.Alternate and currentTranscript not in transcriptsWithMultiAllelVariants:
						if phasingWarning:
							phasingWarning = False
							print('Warning: Still Phases inside VCF. It may not work correctly and can slow down the whole process.')
						# not necessary with vcf-preprocessing - but maybe the data change in the future -> usefull again
						transcriptsWithMultiAllelVariants.append(currentTranscript)
					if "," in variant.Alternate and variant not in phasingVariants:
						phasingVariants.append(variant)
						LogOrganizer.addToLog(LogEnums.COORDINATOR_PHASING_LOG, str(Variant.Chromosome)+
											  str(Variant.Position) +"\t"+ str(Variant.ID) +"\t"+ str(Variant.Reference)
											  + "\t" + str(Variant.Alternate) +"\t"+ str(Variant.Qual) +"\t"+ str(Variant.Filter)
											  +"\t"+  str(Variant.Info))
					if currentTranscript.Lost_Stop:
						transcriptsWithoutStop.append(currentTranscript)
					currentTranscript.ListofVariants.append(variant.Position)
					variant.ListofTranscripts.append(currentTranscript.IndexKey)
					variant.SListofTranscripts.append(currentTranscript.TID)
					if found:
						firstTranscriptMatch = currentTranscriptID
						found = False
					currentTranscriptID += 1
				elif (variant.Position > currentTranscript.EndOfRNA + transcriptRange):
					currentTranscriptID += 1
				else:
					print("Well this should not happen.")
					break
		####
		# For every transcript with multiple allele variants.
		# Creates a copy with all objects and data inside.
		# Takes only the first allele variant for the copy and
		# takes only the second allele variant for the original.
		# The copy gets a new IndexID.
		# Adds new transcript to the existing list in the gff3 handler.
		# Sorts the transcripts per chromosome and position again (like in gff3 handler).
		####
		multiallelwarning = True
		for currentTranscript in transcriptsWithMultiAllelVariants:
			if multiallelwarning:
				print("[multiallelwarning] The dataset contains raw data. This may cause errors.")
				multiallelwarning = False
			currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)
			newTranscript = copy.deepcopy(currentTranscript)
			newTranscript.Remove_Mult_Allel_Entry_In_All_Variant_Information(0)
			currentTranscript.Remove_Mult_Allel_Entry_In_All_Variant_Information(1)
			newTranscript.IndexKey = gff3.GetNextTranscriptIndex()
			currentTranscriptList.append(newTranscript)
			#currentTranscriptList = sorted(currentTranscriptList, key=lambda sTranscript: sTranscript.StartOfRNA)
			#gff3.dictListOfTranscripts[gff3.dictChrNames[name]][newTranscript.TID] = newTranscript
			gff3.AddNewTranscriptToDict(name, newTranscript)


		gff3.updateTransripts(currentTranscriptList,name)
		#gff3.ListOfTranscripts[gff3.dictChrNames[name]] = currentTranscriptList

		#print("Done: " + str(datetime.now() - timeStart))

		kill_transcripts_with_to_many_variants = gff3.GetChrTranscriptsDict(name)
		for trancript in kill_transcripts_with_to_many_variants.values():
			trancript = For_Type_Safety_and_statics.Transcript_Type_Safety(trancript)
			if len(trancript.IntegratedVariantObjects_CDS_Hits) > 500:
				trancript.Transcript_CDS_damaged = True
				print(str(trancript.TID) +"has " + str(len(trancript.IntegratedVariantObjects_CDS_Hits)) + " variations. The Limit is 500. It will be handled as damaged.")


		########
		### Calculate the effect length of mutation effects (version 2)
		###	Creates parts of the changed DNA- and AA-sequence.
		###	The transcripts will be extended, until a stop appears (if there is none).
		########
		print("Calculate the effect length")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		#print(name)
		nr_transcripts_current = 0
		aTranscriptDict = gff3.GetChrTranscriptsDict(
			name)  # self.dictListOfTranscripts = [{}] # List for chromosomes, dict for normal entrys
		for currentTranscript in aTranscriptDict.values():
			currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)
			nr_transcripts_current += 1
			if not currentTranscript.CDS_Exist:
				continue
			if currentTranscript.Transcript_CDS_damaged:
				continue
			#if nr_transcripts_current % 100 == 0:
			#	print(str(nr_transcripts_current) + " ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
			if 'Ma09_t08910.1' == currentTranscript.TID:
				print("bugsearch")
			currentTranscript.Create_IV_Changed_DNA_CDS_Seq(genetic_code, currentTranscript.IntegratedVariantObjects_CDS_Hits,
															stopcodon)
			if currentTranscript.Lost_Stop and stopcodon in currentTranscript.IV_ChangedTranslation:
				#lost stop, but got new somewhere in the transcript
				currentTranscript.Lost_Stop = False
				currentTranscript.prematureStopCodon(stopcodon)
			elif stopcodon not in currentTranscript.IV_ChangedTranslation:
				currentTranscript.Lost_Stop = True
			i = 0
			while currentTranscript.Lost_Stop:
				i+=1
				#print(currentTranscript.TID)
				currentTranscript.Create_IV_ChangedTranslation(genetic_code)
				if stopcodon in currentTranscript.IV_ChangedTranslation:
					# print("new stopcodon already inside transcript" + str(transcriptHier.TID))

					#LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"New Stopcodon is in " + str(currentTranscript.TID) + "\n")
					lastPosInCDS = currentTranscript.LastCDSPosition
					if currentTranscript.Found_New_Stop or stopcodon in currentTranscript.IV_ChangedTranslation:
						break
					currentTranscript.Lost_Stop = False
					currentTranscript.Found_New_Stop = True
					#currentTranscript.Find_New_Stop() todo
					#break #1736 without - 122 with
				# print ("Lost_Stop not completely inside transcript.")
				lastPosInCDS = currentTranscript.LastCDSPosition
				if currentTranscript.ForwardDirection == TranscriptEnum.FORWARD:
					try:
						nextDNA = ghandler.seq(name, lastPosInCDS + 1, lastPosInCDS + 100)
					except SequenceHandlingError as she:
						LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG, str(currentTranscript.TID) + "\t" + she.description)
						if she.sequence_part == "":
							print("transcript " + str(currentTranscript.TID) + " declared as broken(circular?). It reached the end of the contig/chrom.")
							currentTranscript.Transcript_CDS_damaged = True
							break
						else:
							nextDNA = she.sequence_part

				elif currentTranscript.ForwardDirection == TranscriptEnum.REVERSE:
					try:
						nextDNA = ghandler.seq(name, lastPosInCDS - 100, lastPosInCDS - 1)
					except SequenceHandlingError as she:
						LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG,str(currentTranscript.TID) + "\t" + she.description)
						if she.sequence_part == "":
							print("transcript " + str(currentTranscript.TID) + " declared as broken(circular?). It reached the end of the contig/chrom.")
							currentTranscript.Transcript_CDS_damaged = True
							break
						else:
							nextDNA = she.sequence_part
					#nextDNA = For_Type_Safety_and_statics.ReverseSeq(nextDNA) #reverse, because it will be added to the already reversed dna transcript
				else:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Error: No Direction: " + str(currentTranscript.TID)+"\n")
					break
				currentTranscript.Find_New_Stop(nextDNA, genetic_code, stopcodon)

				if len(currentTranscript.IV_Changed_DNA_CDS_Seq) % 3 != 0 \
						or len(currentTranscript.Complete_CDS) % 3 != 0 \
						or len(currentTranscript.Rev_CDS) % 3 != 0:
					#dostuff 	1 -> +2bp
					#			2 -> +1bp
					currentTranscript.checkLastVariants(ghandler, genetic_code, stopcodon)
					# +2 positions means 2 potential more variant effects.
					# repeatly, and check if deletions now have a new effect, because in rev cds it can be .... .... ....
					# and it is possible, that damaged transcripts (rev cds) is still functional because of this, meeeh fuck it
				if i == variant_transcript_exceeding_warning_int +1:
					print(str(currentTranscript.TID) + " exceeds variation range limit of " + str(int(variant_transcript_exceeding_warning_int*100)) + " nucleotides.")
				if i == 20:
					currentTranscript.Transcript_CDS_damaged = True
					print("transcript "+ str(currentTranscript.TID) + " declared as broken, after extending it by: " +str(i*100) + " nucleotides." )
					break
			j = 0
			while currentTranscript.origDNAtoshort:
				#print(currentTranscript.TID)
				j +=1
				if j == 10:
					currentTranscript.Transcript_CDS_damaged = True
					print("transcript " + str(currentTranscript.TID) + " declared as broken, after extending it " + str(j) + " times unsuccessfully." )
					break
				currentTranscript.checkLastVariants(ghandler, genetic_code, stopcodon)
		#print("Done: " + str(datetime.now() - timeStart))
		####################################################
		print("Complete AA sequences")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		Create_Sequences(gff3, Orig_AA, New_AA, genetic_code,name)
		#print("Done: " + str(datetime.now() - timeStart))
		### Find Stop-Codons, if hey are made by a frameshift -> LABEl the mutation
		print("Flag Frameshift caused stops")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		aTranscriptDict = gff3.GetChrTranscriptsDict(name)
		for currentTranscript in aTranscriptDict.values():
			currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)
			if not currentTranscript.CDS_Exist:
				continue
			if currentTranscript.Transcript_CDS_damaged:
				continue
			firstStopPosition = (currentTranscript.IV_ChangedTranslation.find(stopcodon)+1) * 3
			if firstStopPosition == -2:
				continue
			last_frameshifter = ""
			for variant in currentTranscript.IntegratedVariantObjects_CDS_Hits:
				variant = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(variant)
				if TranscriptEnum.STOP_GAINED in variant.Classification \
						or TranscriptEnum.STOP_CHANGED in variant.Classification:  # stop because of variants
					last_frameshifter = ""
					break
				if variant.Changed_CDS_Position < firstStopPosition:
					if (len(variant.Ref) - 1) % 3 != 0 or (
						len(variant.Alt) - 1) % 3 != 0:  # only frameshift and caused no stop
						if TranscriptEnum.STOP_GAINED not in variant.Classification \
								and TranscriptEnum.STOP_CHANGED not in variant.Classification:
							last_frameshifter = variant
					continue
				else:
					break
			if last_frameshifter == "":
				continue
			else:
				last_frameshifter.Classification.append(TranscriptEnum.STOP_CAUSED_IN)
				last_frameshifter.STOP_CAUSED_IN = firstStopPosition - last_frameshifter.Changed_CDS_Position
		#print("Done: " + str(datetime.now() - timeStart))
		#################################################################
		print("Complete Data check.")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		Complete_Check(gff3, ghandler,name)
		#print("Done: " + str(datetime.now() - timeStart))
		#################################################################
		print("Write data")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		"""
		Old, but maybe still usefull functions for further research/programming.
		The functions are now inside old_stuff()
		# Write_New_Special_VCF_File(Output_Data_Path, gff3)
		# Write_All_VCF_Damaged_Transcripts(Output_Data_Path,gff3)
		# Write_VCF_With_Key(Output_Data_Path,gff3)
		# Write_All_VCF_Stop_Lost(Output_Data_Path,gff3)
		# Write_All_VCF_NO_STOP(Output_Data_Path,gff3)
		# Write_All_VCF_NO_START(Output_Data_Path,gff3)
		# Write_All_VCF_To_Many_Stops(Output_Data_Path,gff3)
		"""
		Write_All_VCF_File(Output_Data_Path, gff3, Old_Info, Ref_Codons, Ref_AA, Alt_Codons, Alt_AA,extend_File,name)
		Write_All_Fasta(Output_Data_Path,Fasta_Data_Name,False,Orig_DNA,Orig_AA,New_DNA,New_AA,gff3,extend_File,name)
		print("Done: " + str(datetime.now() - timeStart))
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		extend_File = True
		gff3.freeRAM(name)
		vcf.freeRam(name)
		gc.collect()
		#################################################################
	print("Create log files.")

	LogOrganizer.writeAllLogs(outpath)

	print("Everything is done: " + str(datetime.now() - timeStart))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

	return os.getpid()

	#2232
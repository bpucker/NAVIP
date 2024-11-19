__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"

import copy #to copy objects
from GenomeHandler import GenomeHandler, Fasta_Enum, SequenceHandlingError
from VCF_Handler import VCF_Handler
from VCF_Variant import Variant
from Transcript import TranscriptEnum, For_Type_Safety_and_statics
from GFF3_Handler import GFF3_Handler
from datetime import datetime
from LogOrganizer import LogOrganizer, LogEnums
from SnpEff_HGVS_Converter import SnpEff_HGVS_Converter
import os
import resource
import gc

def navip_main_coordinator(invcf, ingff, infasta, outpath):
	"""
	Coordinates all functions and classes of the main module.
	Also contains a lot of functions, but naming this module "main" or "navip" isn't
	an option, because it would be confusing with the different modules.
	:param invcf: Path and name to a VCF file. Preferred the preprocessed file.
	:param ingff: Path and name to a GFF3 file.
	:param infasta: Path and name to a FASTA file.
	:param outpath: Path to and including the output folder.
	:return: Nothing.
	"""

	def write_all_VCF_file(data_path: str, gff3: GFF3_Handler, old_info: bool,
						   ref_codons: bool, ref_AA: bool,
						   alt_codons: bool, alt_AA: bool,
						   extend_file: bool, chr_name: str):
		"""
		In this function the VCF file with the additional neighborhood awareness information will be written.
		Furthermore, it is possible (and standard) to write a few more information into the NAVIP info:
		TranscriptID, strand direction, classifications, keys with shared effect, ref and alt codon
		with the variant position inside this codon. This information will be written into every new VCF file, too.
		Transcripts without variants inside the CDS, without CDS and those, which are (possible) damaged within the
		splicing structure will be skipped.

		:param data_path: Path into the existing output folder.
		:param gff3: The GFF3 data structure with all transcripts and their variants/sequences.
		:param old_info: The old VCF info column will be placed after the NAVIP info, if true.
		:param ref_codons: The REF DNA codons will be written, if true.
		:param ref_AA: The REF AA will be written, if true.
		:param alt_codons: The ALT codons will be written, if true.
		:param alt_AA: The ALT AA will be written, if true.
		:return: Nothing.
		"""
		# chrom	Pos	ID	Ref	Alt Info
		if not extend_file:
			vcf_file = open(data_path + "All_VCF.vcf", "w")
			vcf_file.write("##Chrom\tPos\tID\tRef\tAlt\tQual\tFilter\tInfo\n")
			vcf_file.write("##NAVIP: All Data Output\n")
			vcf_file.write("##Please note, that the Variant_Position_in_Codon is read from left to right in forward "
						   "and right to left in reverse strand direction.\n")
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
		infoline_parser = SnpEff_HGVS_Converter
		#for chr_name in gff3.get_chromosome_names():

		for transcript in gff3.get_chr_transcripts_dict(chr_name).values():
			transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(transcript)

			if not transcript.CDS_Exist or transcript.Transcript_CDS_damaged:
				#no CDS or transcript exons are damaged
				continue
			if not transcript.IntegratedVariantObjects_CDS_Hits:
				#no variants inside the cds
				continue

			if transcript.ForwardDirection == TranscriptEnum.FORWARD:
				direction = "FOR"
			elif transcript.ForwardDirection == TranscriptEnum.REVERSE:
				direction = "REV"
			else:
				#will never happen
				LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Error:NOT Forward, NOT Reverse:" + str(transcript.TID) + "\n")
				direction = "Error:NOT Forward, NOT Reverse."

			for vinfo in transcript.IntegratedVariantObjects_CDS_Hits:
				vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
				if vinfo.ChrPosition == 5921700 and transcript.TID == "Ma09_t08910.1":
					print("bugsearch")
				if vinfo.ChrPosition == 3917462 and transcript.TID == "Ma00_t01100.1":
					print("bugsearch")
				try:
					snpeff_string = infoline_parser.convert_main(transcript,vinfo)
				except Exception:
					snpeff_string = ""
					LogOrganizer.addToLog(LogEnums.COORDINATOR_VARIANT_LOG,
										  str(vinfo.ChrPosition) + "\t" + str(transcript.TID))
					print('Critical error with variant: ' + str(vinfo.ChrPosition) + '\t' + str(transcript.TID))
				classification_string = ""
				class_list = vinfo.Classification
				class_list_length = len(class_list)
				for i in range(0, class_list_length):
					classification_string += class_list[i].value
					if class_list[i] == TranscriptEnum.STOP_CAUSED_IN:
						classification_string += str(vinfo.STOP_CAUSED_IN)
					if i != (class_list_length - 1):  # no tab after last entry
						classification_string += ","
					elif i == (class_list_length - 1):
						classification_string += "|"
					else:
						LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Bug Write_All_VCF_File: classification list length:" + str(transcript.TID) + "\n")

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
						#seperate entries
						elif i == shared_effect_list_len - 1:
							shared_effect_list += "|"
						#after last entry
						else:
							LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Bug Write_All_VCF_File: shared_effect_list:" + str(transcript.TID) + "\n")

				if len(vinfo.OrigTriplets) % 3 != 0:
					bug = "OrigTriplets: " + str(vinfo.OrigTriplets) \
						  + "\tFrom: " \
						  + transcript.TID \
						  + ": " \
						  + str(vinfo.ChrPosition) \
						  + " " \
						  + classification_string
					LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,bug + "\n")
				elif len(vinfo.ChangedTriplets) % 3 != 0:
					bug = "ChangedTriplets: " + str(vinfo.ChangedTriplets) \
						  + "\tFrom: " \
						  + transcript.TID \
						  + ": " \
						  + str(vinfo.ChrPosition) \
						  + " " \
						  + classification_string
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
				NAVIP_INFO_LIST = ["NAV1=" + str(transcript.TID) + "|", str(direction) + "|", classification_string]
				if shared_effect_list == "":
					NAVIP_INFO_LIST.append("NONE|")
				else:
					NAVIP_INFO_LIST.append(str(shared_effect_list))
				if ref_codons:
					NAVIP_INFO_LIST.append(str(vinfo.OrigTriplets) + "/" + str(vinfo.OrigRaster) + "|")
				if ref_AA:
					NAVIP_INFO_LIST.append(str(vinfo.OrigAmino) + "|")
				NAVIP_INFO_LIST.append(str(vinfo.Unchanged_CDS_Position) + "|")
				if alt_codons:
					NAVIP_INFO_LIST.append(str(vinfo.ChangedTriplets) + "/" + str(vinfo.Changed_Raster) + "|")
				if alt_AA:
					NAVIP_INFO_LIST.append(str(vinfo.NewAmino) + "|")
				NAVIP_INFO_LIST.append(str(vinfo.Changed_CDS_Position) + ";")
				if old_info:
					if "\n" in vinfo.OLD_Info:
						NAVIP_INFO_LIST.append(str(vinfo.OLD_Info))
					else:
						NAVIP_INFO_LIST.append(str(vinfo.OLD_Info) + "\n")

				NAVIP_INFO = "".join(NAVIP_INFO_LIST)

				data_to_write.append(
					str(chr_name) + "\t"
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
		vcf_file.close()

	def create_sequences(gff3: GFF3_Handler, orig_AA: bool, new_AA: bool, genetic_code: dict, chr_name: str):
		"""
		This function calculates the old and new AA CDS to be sure all transcripts are having the AA sequences.
		There are a few possibilities to this point, which can lead to skip a few or all transcripts AA sequences.
		Transcripts without CDS are excluded. Transcripts with a damaged CDS do not get the new AA sequence.
		:param gff3: The GFF3 data structure with all transcripts and their variants/sequences.
		:param orig_AA: Will create the original AA sequence, if true.
		:param new_AA: Will create the new AA sequence, if true.
		:param genetic_code: The dictionary with the genetic code.
		:return: Nothing.
		"""
		#for chr_name in gff3.get_chromosome_names():
		#print(chr_name)
		for transcript in gff3.get_chr_transcripts_dict(chr_name).values():
			try:
				transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(transcript)

				if not transcript.CDS_Exist:
					continue
				if orig_AA:
					transcript.Create_IV_OriginalTranslation(genetic_code)

				if transcript.Transcript_CDS_damaged:
					continue
				if new_AA:
					transcript.Create_IV_ChangedTranslation(genetic_code)
			except Exception:
				LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG,"\t" + str(transcript.TID))

	def complete_check(gff3: GFF3_Handler, ghandler: GenomeHandler, chr_name: str):
		"""
		It's always a good idea to test the data in the end.
		Will print warnings, if there is incorrect data.
		:param gff3: The GFF3 data structure with all transcripts and their variants/sequences.
		:param ghandler: The reference genome sequences for every chromosome is inside.
		:return: Nothing.
		"""
		#for chr_name in shared_chromosomes_FA_GFF3:
		chr_transcript_list = gff3.get_chr_transcripts_dict(chr_name).values()
		for current_transcript in chr_transcript_list:
			current_transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(current_transcript)
			if not current_transcript.CDS_Exist or current_transcript.Transcript_CDS_damaged:
				continue

			for vinfo in current_transcript.IntegratedVariantObjects_CDS_Hits:
				vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
				try :
					if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
						# Variant <-> Fasta Check
						if vinfo.Ref != ghandler.singleSeq(chr_name, vinfo.ChrPosition):
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "Sub_REF != Fasta-Seq: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
							#print("Sub_REF != Fasta-Seq: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

						# Variant <-> Transcript with Direction Check
						if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
							if vinfo.Ref != current_transcript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "Sub_REF != CDS-Seq 1: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
								#print(
								#	"Sub_REF != CDS-Seq 1: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

							if vinfo.Ref != current_transcript.SeqInCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "Sub_REF != CDS-Seq 2: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
								#print(
								#	"Sub_REF != CDS-Seq 2: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

							if vinfo.Alt != current_transcript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_ALT != CDS-Seq 3: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n" )
								#print(
								#	"Sub_ALT != CDS-Seq 3: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

						elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
							if vinfo.Ref != current_transcript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition) and not current_transcript.CDS_Changed:
								# When the used CDS gets changed, then it is of cause not possible to check against it.
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 1.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 1.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 2.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 2.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")
							#print(str(vinfo.ChrPosition) + "\t" + str(current_transcript.TID) )
							if vinfo.ReverseAlt != current_transcript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_Alt != CDS-Seq 3.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_Alt != CDS-Seq 3.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 4.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 4.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"Sub_REF != CDS-Seq 5.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("Sub_REF != CDS-Seq 5.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

					elif TranscriptEnum.INSERTION in vinfo.Classification:
						# Variant <-> Fasta Check
						if vinfo.Ref != ghandler.singleSeq(chr_name, vinfo.ChrPosition):
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != Fasta-Seq: " + str(current_transcript.TID) + "\t" + str(
								vinfo.ChrPosition) + "\n" )
							#print("INSERTION_REF != Fasta-Seq: " + str(current_transcript.TID) + "\t" + str(
							#	vinfo.ChrPosition) + "\n")

						# Variant <-> Transcript with Direction Check
						if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
							if vinfo.Ref != current_transcript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 1: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("INSERTION_REF != CDS-Seq 1: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.Ref != current_transcript.SeqInCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_REF != CDS-Seq 2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition))

							if vinfo.Alt != current_transcript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
													vinfo.Changed_CDS_Position + len(vinfo.Alt) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_ALT != CDS-Seq 3: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n" )
								#print("INSERTION_ALT != CDS-Seq 3: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition)+ "\n")

						elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
							if vinfo.Ref != current_transcript.SeqInCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition)and not current_transcript.CDS_Changed:
								# When the used CDS gets changed, then it is of cause not possible to check against it.
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 1.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_REF != CDS-Seq 1.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_REF != CDS-Seq 2.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_REF != CDS-Seq 2.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseAlt != current_transcript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
													vinfo.Changed_CDS_Position + len(vinfo.Alt) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"INSERTION_Alt != CDS-Seq 3.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n" )
								#print("INSERTION_Alt != CDS-Seq 3.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition) + "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS(
									vinfo.ChrPosition,
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "INSERTION_REF != CDS-Seq 4.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("INSERTION_REF != CDS-Seq 4.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition))

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
									vinfo.Unchanged_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG, "INSERTION_REF != CDS-Seq 5.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition) + "\n")
								#print("INSERTION_REF != CDS-Seq 5.2: " + str(current_transcript.TID) + "\t" + str(
								#	vinfo.ChrPosition))

					elif TranscriptEnum.DELETION in vinfo.Classification:

						# Variant <-> Fasta Check
						try:
							if vinfo.Ref != ghandler.seq(chr_name, vinfo.ChrPosition, vinfo.ChrPosition + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
									"DELETION_REF != Fasta-Seq: " + str(current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
						except SequenceHandlingError as she:
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
												  str(she.description) + str(
													  current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
							if vinfo.Ref != she.sequence_part:
								#ghandler.seq(chr_name, vinfo.ChrPosition,#vinfo.ChrPosition + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
													  "DELETION_REF != Fasta-Seq: " + str(
														  current_transcript.TID) + "\t" + str(
														  vinfo.ChrPosition) + "\n")

						# Variant <-> Transcript with Direction Check
						if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
							if vinfo.Ref != current_transcript.SeqInCDS(
									vinfo.ChrPosition,
													vinfo.ChrPosition + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 1: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.Ref != current_transcript.SeqInCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
													vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.Alt != current_transcript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_ALT != CDS-Seq 3: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

						elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
							if vinfo.Ref != current_transcript.SeqInCDS(
									vinfo.ChrPosition,
													vinfo.ChrPosition + len(vinfo.Ref) - 1)and not current_transcript.CDS_Changed:
								# When the used CDS gets changed, then it is of cause not possible to check against it.
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 1.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
													vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 2.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")


							if vinfo.ReverseAlt != current_transcript.SeqInIV_Changed_DNA_CDS_Seq(
									vinfo.Changed_CDS_Position,
									vinfo.Changed_CDS_Position):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_Alt != CDS-Seq 3.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS(
											vinfo.ChrPosition + (len(vinfo.Ref) - 1),
									vinfo.ChrPosition):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 4.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

							if vinfo.ReverseRef != current_transcript.SeqInRevCDS_OverCDS_Position(
									vinfo.Unchanged_CDS_Position,
													vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
								LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"DELETION_REF != CDS-Seq 5.2: " + str(current_transcript.TID) + "\t" + str(
									vinfo.ChrPosition)+ "\n")

					else:
						LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"No Classification in Complete_Check:\t" + str(current_transcript.TID) +"\t" + str(vinfo.ChrPosition)+ "\n")
				except Exception:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_CRITICAL_LOG, str(vinfo.ChrPosition) + "\t" + str(current_transcript.TID))

	def write_all_fasta(data_path: str, data_name: str, date_time: bool,
						orig_DNA: bool, orig_AA: bool, new_DNA: bool, new_AA: bool,
						gff3: GFF3_Handler, extend_file: bool, chr_name: str):
		"""
		Writes a new fasta file with all CDS sequences from every transcript,
		including possible spliced variant of the transcripts (GFF3 data).
		It is possible to choose, which sequences will be written.
		Tags for the type of every sequence exist:
		ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == New Amino Acid Sequence
		Furthermore: Possible damaged transcripts will be separated into "<data_name>_damaged.txt".
		Including its information about the genes, UTRs, exons, and all available variant information.
		:param data_path: Path to the folder for the new fasta file.
		:param data_name: Name of the new fasta file.
		:param date_time: The date and time will be included inside the filename, if true.
		:param orig_DNA:  Writes the original CDS DNA sequence, if true.
		:param orig_AA: Writes the original CDS A sequence, if true.
		:param new_DNA: Writes the new CDS DNA sequence, if true.
		:param new_AA: Writes the new CDS AA sequence, if true.
		:param gff3: The GFF3 data structure with all transcripts and their variants/sequences.
		:return: Nothing.
		"""

		if orig_DNA or orig_AA or new_DNA or new_AA:
			if data_name == "":
				data_name = "all_transcripts_data"
			if date_time:
				timestop = str(datetime.now()) + "_"
			else:
				timestop = ""
			if not extend_file:

				normal_data_list_to_write = [
					"#ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == New Amino Acid Sequence\n"]
				transcript_CDS_damaged_data_list = [
					"#Data in here may be incomplete and incorrect because of deletions, which destroyed parts the splice sides.\n"]
			else:
				normal_data_list_to_write = []
				transcript_CDS_damaged_data_list = []
			#for chr_name in gff3.get_chromosome_names():

			for current_transcript in gff3.get_chr_transcripts_dict(chr_name).values():
				current_transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(current_transcript)
				try:
					if not current_transcript.CDS_Exist:
						continue

					if orig_DNA:
						if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
							normal_data_list_to_write.append(
								">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.ODNA.value) + "\n")
							normal_data_list_to_write.append(str(current_transcript.Complete_CDS) + "\n")
						elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
							normal_data_list_to_write.append(
								">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.ODNA.value) + "\n")
							normal_data_list_to_write.append(str(current_transcript.Rev_CDS) + "\n")
						else:
							LogOrganizer.addToLog(LogEnums.COORDINATOR_FASTA_FILE_ERROR_LOG,"No Direction:Write_All_Fasta: " + str(current_transcript.TID) + "\n")
					if orig_AA:
						normal_data_list_to_write.append(
							">" + str(current_transcript.TID) + "|" + Fasta_Enum.OAA.value + "\n")
						normal_data_list_to_write.append(str(current_transcript.IV_OriginalTranslation) + "\n")

					if current_transcript.Transcript_CDS_damaged:
						transcript_CDS_damaged_data_list.append(
							str(current_transcript.TID) + " " + str(current_transcript.Gene_Info_String))
						for utr in current_transcript.UTR_Description:
							if "\n" in utr:
								transcript_CDS_damaged_data_list.append(utr)
							else:
								transcript_CDS_damaged_data_list.append(str(utr) + "\n")
						for exon in current_transcript.EXON_Descriptin:
							if "\n" in exon:
								transcript_CDS_damaged_data_list.append(exon)
							else:
								transcript_CDS_damaged_data_list.append(str(exon) + "\n")
						for vinfo in current_transcript.IntegratedVariantObjects_CDS_Hits:
							# vinfo = For_Type_Safety.Variant_Information_Storage_Type_Safety(vinfo)
							classification_string = ""
							class_list = vinfo.Classification
							class_list_length = len(class_list)
							for i in range(0, class_list_length):
								classification_string += class_list[i].value
								if i != (class_list_length - 1):  # no tab after last entry
									classification_string += "\t"
							transcript_CDS_damaged_data_list.append(
								str(vinfo.ChrPosition) + "\t" +
								str(vinfo.Ref) + "\t" +
								str(vinfo.Alt) + "\t" +
								str(classification_string) + "\n")
						continue
					# Meh, more DNA/AA Sequences, thanks to Stop-lost transcription length changes
					normal_data_list_to_write.append(
						">" + str(current_transcript.TID) + "|" + Fasta_Enum.ouchDNA.value + "\n")
					normal_data_list_to_write.append(str(current_transcript.uChDNAsequence) + "\n")

					normal_data_list_to_write.append(
						">" + str(current_transcript.TID) + "|" + Fasta_Enum.ouchAA.value + "\n")
					normal_data_list_to_write.append(str(current_transcript.uChAAsequence) + "\n")

					if new_DNA:
						normal_data_list_to_write.append(
							">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.nDNA.value) + "\n")
						normal_data_list_to_write.append(str(current_transcript.IV_Changed_DNA_CDS_Seq) + "\n")
					if new_AA:
						normal_data_list_to_write.append(
							">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.nAA.value) + "\n")
						normal_data_list_to_write.append(str(current_transcript.IV_ChangedTranslation) + "\n")
				except Exception:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG,"Write_All_Fasta" + "\t" + str(current_transcript.TID) + "\n")

			normal_output = "".join(normal_data_list_to_write)
			if not extend_file:
				new_fasta_file = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "w")
				new_fasta_file.write(normal_output)
				new_fasta_file.close()
			else:
				new_fasta_file = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "a")
				new_fasta_file.write(normal_output)
				new_fasta_file.close()

			if not extend_file:
				damaged_output = "".join(transcript_CDS_damaged_data_list)
				new_fasta_file_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt", "w")
				new_fasta_file_damaged.write(damaged_output)
				new_fasta_file_damaged.close()
			else:
				damaged_output = "".join(transcript_CDS_damaged_data_list)
				new_fasta_file_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt", "a")
				new_fasta_file_damaged.write(damaged_output)
				new_fasta_file_damaged.close()


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
	orig_DNA = True
	orig_AA = True
	new_DNA = True
	new_AA = True
	fasta_data_name = ""

	###
	# Output VCF file
	###
	old_info = True
	ref_codons = True
	ref_AA = True
	alt_codons = True
	alt_AA = True
	###
	# How many chromosomes will be analysed.
	# 0 stands for all.
	###
	firstx_chr = 0
	time_start = datetime.now()


	vcf_path_and_name = invcf
	gff3_path_and_name = ingff
	fasta_file_path = infasta
	output_data_path = outpath

	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
	#######################################
	print("read vcf")
	vcf = VCF_Handler(vcf_path_and_name, firstx_chr)
	print("Done: " + str(datetime.now() - time_start))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
	#######################################
	print("read gff3")
	#print("Hopefully sorted after seqID")
	gff3 = GFF3_Handler(gff3_path_and_name)
	print("Done: " + str(datetime.now() - time_start))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
	#######################################
	print("read fa")
	ghandler = GenomeHandler(fasta_file_path, firstx_chr)
	print("Done: " + str(datetime.now() - time_start))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
	#######################################
	print("Chromosome:Names")
	if len(gff3.get_chromosome_names()) < 10: print("GFF3: " + str(gff3.get_chromosome_names()))
	else: print("GFF3: " + str(len(gff3.get_chromosome_names())) + " contigs/chromosomes.")
	if len(vcf.get_chromosome_names()) < 10: print("VCF: " + str(vcf.get_chromosome_names()))
	else: print("VCF: " + str(len(vcf.get_chromosome_names())) + " contigs/chromosomes.")
	if len(ghandler.GetChromosomeNames()) < 10: print("fasta: " + str(ghandler.GetChromosomeNames()))
	else: print("fasta:" + str(len(ghandler.GetChromosomeNames())) + " contigs/chromosomes.")
	shared_chromosomes = []
	shared_chromosomes_FA_GFF3 = []
	#######################################


	for name in gff3.get_chromosome_names():
		if name in vcf.get_chromosome_names() and name in ghandler.GetChromosomeNames():
			shared_chromosomes.append(name)
		if name in ghandler.GetChromosomeNames():
			shared_chromosomes_FA_GFF3.append(name)
	if len(shared_chromosomes) < 10: print(print("shared_chromosomes(all):" + str(shared_chromosomes)))
	else: print("shared_chromosomes(all):" + str(len(shared_chromosomes)) + " contigs/chromosomes.")

	#print("shared_chromosomes, fasta, gff3 (and used):" + str(shared_chromosomes_FA_GFF3))
	### For Writing the CDS
	### in all transcripts

	extend_file = False

	for name in shared_chromosomes:
		print(name)
		print("Write the CDS and Rev CDS")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

		#chrID = gff3.get_chromosome_ID(name)

		#for i in gff3.dict_list_of_transcripts[gff3.dict_chr_names[name]]:  # only a number, not the object
		#	transcripts[i].CompleteTheCDS(ghandler.seq(name, transcripts[i].StartOfRNA, transcripts[i].EndOfRNA))
		#	if transcripts[i].ForwardDirection == TranscriptEnum.REVERSE:
		#		transcripts[i].ReverseTheCDS()

		for transcript in gff3.get_chr_transcripts_dict(name).values():
			try:
				transcript.CompleteTheCDS(ghandler.seq(name, transcript.StartOfRNA, transcript.EndOfRNA), genetic_code)
			except SequenceHandlingError as she:
				LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG, str(transcript.TID) + "\t"
									  + str(transcript.StartOfRNA) + "\t" + str(transcript.EndOfRNA) + "\t" + str(she.description))
				if she.sequence_part != "":
					transcript.CompleteTheCDS(she.sequence_part, genetic_code)
			if transcript.ForwardDirection == TranscriptEnum.REVERSE:
				transcript.ReverseTheCDS(genetic_code)
		#print("Done: " + str(datetime.now() - time_start))
		###
		# Connecting all transcripts with all their mutations/variants.
		# Contains handling of triallele variants, if an original vcf file was used.
		# But it can't handle all possible appearances. VCF preprocessing exist for a reason.
		# The connection triggers a first evaluation of the variants, too.
		###
		transcript_range = 300 # base pairs before RNA starts and after RNA ends. Just in Case for long InDels.
		variant_transcript_exceeding_warning_int = transcript_range / 100
		phasing_warning = True
		print("Connect transcripts with variants")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		#print("Search in " + name)

		phasing_variants = [] # should normally be empty
		first_transcript_match = 0 #
		current_transcript_list = gff3.get_chr_transcript_list(name)
		transcripts_with_multiallele_variants = []
		transcripts_without_stop = []


		for variant in vcf.get_chr_VCF_variant_list(name):
			found = True
			current_transcript_ID = first_transcript_match

			while True:
				###
				# The Following case: no more available transcripts, because the position of the mutation is after the last transcript
				###
				if current_transcript_ID == len(current_transcript_list):
					# print("No more transcripts.")
					break
				current_transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(current_transcript_list[current_transcript_ID])
				if (
							current_transcript.StartOfRNA - transcript_range) > variant.Position:  # transcript.start higher than variant.position
					break
				elif (current_transcript.StartOfRNA - transcript_range) <= \
						variant.Position \
						<= (current_transcript.EndOfRNA + transcript_range):  # between start and end (+ range)

					current_transcript.Add_Variant_Information(variant)

					if "," in variant.Alternate and current_transcript not in transcripts_with_multiallele_variants:
						if phasing_warning:
							phasing_warning = False
							print('Warning: Still Phases inside VCF. It may not work correctly and can slow down the whole process.')
						# not necessary with vcf-preprocessing - but maybe the data change in the future -> useful again
						transcripts_with_multiallele_variants.append(current_transcript)
					if "," in variant.Alternate and variant not in phasing_variants:
						phasing_variants.append(variant)
						LogOrganizer.addToLog(LogEnums.COORDINATOR_PHASING_LOG, str(Variant.Chromosome)+
											  str(Variant.Position) +"\t"+ str(Variant.ID) +"\t"+ str(Variant.Reference)
											  + "\t" + str(Variant.Alternate) +"\t"+ str(Variant.Qual) +"\t"+ str(Variant.Filter)
											  +"\t"+  str(Variant.Info))
					if current_transcript.Lost_Stop:
						transcripts_without_stop.append(current_transcript)
					current_transcript.ListofVariants.append(variant.Position)
					variant.ListofTranscripts.append(current_transcript.IndexKey)
					variant.SListofTranscripts.append(current_transcript.TID)
					if found:
						first_transcript_match = current_transcript_ID
						found = False
					current_transcript_ID += 1
				elif variant.Position > current_transcript.EndOfRNA + transcript_range:
					current_transcript_ID += 1
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
		multiallele_warning = True
		for current_transcript in transcripts_with_multiallele_variants:
			if multiallele_warning:
				print("[multiallele_warning] The dataset contains raw data. This may cause errors.")
				multiallele_warning = False
			current_transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(current_transcript)
			new_transcript = copy.deepcopy(current_transcript)
			new_transcript.Remove_Mult_Allel_Entry_In_All_Variant_Information(0)
			current_transcript.Remove_Mult_Allel_Entry_In_All_Variant_Information(1)
			new_transcript.IndexKey = gff3.get_next_transcript_index()
			current_transcript_list.append(new_transcript)
			#current_transcript_list = sorted(current_transcript_list, key=lambda s_transcript: s_transcript.StartOfRNA)
			#gff3.dict_list_of_transcripts[gff3.dict_chr_names[name]][new_transcript.TID] = new_transcript
			gff3.add_new_transcript_to_dict(name, new_transcript)


		gff3.update_transcripts(current_transcript_list, name)
		#gff3.list_of_transcripts[gff3.dict_chr_names[name]] = current_transcript_list

		#print("Done: " + str(datetime.now() - time_start))

		kill_transcripts_with_too_many_variants = gff3.get_chr_transcripts_dict(name)
		for transcript in kill_transcripts_with_too_many_variants.values():
			transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(transcript)
			if len(transcript.IntegratedVariantObjects_CDS_Hits) > 500:
				transcript.Transcript_CDS_damaged = True
				print(str(transcript.TID) +"has " + str(len(transcript.IntegratedVariantObjects_CDS_Hits)) + " variations. The Limit is 500. It will be handled as damaged.")


		########
		### Calculate the effect length of mutation effects (version 2)
		###	Creates parts of the changed DNA- and AA-sequence.
		###	The transcripts will be extended, until a stop appears (if there is none).
		########
		print("Calculate the effect length")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		#print(name)
		nr_transcripts_current = 0
		a_transcript_dict = gff3.get_chr_transcripts_dict(
			name)  # self.dictListOfTranscripts = [{}] # List for chromosomes, dict for normal entries
		for current_transcript in a_transcript_dict.values():
			current_transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(current_transcript)
			nr_transcripts_current += 1
			if not current_transcript.CDS_Exist:
				continue
			if current_transcript.Transcript_CDS_damaged:
				continue
			#if nr_transcripts_current % 100 == 0:
			#	print(str(nr_transcripts_current) + " ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
			if 'Ma09_t08910.1' == current_transcript.TID:
				print("bugsearch")
			current_transcript.Create_IV_Changed_DNA_CDS_Seq(genetic_code, current_transcript.IntegratedVariantObjects_CDS_Hits,
															stopcodon)
			if current_transcript.Lost_Stop and stopcodon in current_transcript.IV_ChangedTranslation:
				#lost stop, but got new somewhere in the transcript
				current_transcript.Lost_Stop = False
				current_transcript.prematureStopCodon(stopcodon)
			elif stopcodon not in current_transcript.IV_ChangedTranslation:
				current_transcript.Lost_Stop = True
			i = 0
			while current_transcript.Lost_Stop:
				i+=1
				#print(current_transcript.TID)
				current_transcript.Create_IV_ChangedTranslation(genetic_code)
				if stopcodon in current_transcript.IV_ChangedTranslation:
					# print("new stopcodon already inside transcript" + str(transcript.TID))

					#LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"New Stopcodon is in " + str(current_transcript.TID) + "\n")
					if current_transcript.Found_New_Stop or stopcodon in current_transcript.IV_ChangedTranslation:
						break
					current_transcript.Lost_Stop = False
					current_transcript.Found_New_Stop = True
					#current_transcript.Find_New_Stop() todo
					#break #1736 without - 122 with
				# print ("Lost_Stop not completely inside transcript.")
				last_pos_in_CDS = current_transcript.LastCDSPosition
				if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
					start_position = last_pos_in_CDS + 1; end_position = last_pos_in_CDS + 100
				elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
					start_position = last_pos_in_CDS - 100; end_position = last_pos_in_CDS - 1
					#next_DNA = For_Type_Safety_and_statics.ReverseSeq(next_DNA) #reverse, because it will be added to the already reversed dna transcript
				else:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Error: No Direction: " + str(current_transcript.TID)+"\n")
					break

				try:
					next_DNA = ghandler.seq(name, start_position, end_position)
				except SequenceHandlingError as she:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_TRANSCRIPT_LOG, str(current_transcript.TID) + "\t" + she.description)
					if she.sequence_part == "":
						print("transcript " + str(current_transcript.TID) + " declared as broken(circular?). It reached the end of the contig/chrom.")
						current_transcript.Transcript_CDS_damaged = True
						break
					else:
						next_DNA = she.sequence_part
				current_transcript.Find_New_Stop(next_DNA, genetic_code, stopcodon)

				if len(current_transcript.IV_Changed_DNA_CDS_Seq) % 3 != 0 \
						or len(current_transcript.Complete_CDS) % 3 != 0 \
						or len(current_transcript.Rev_CDS) % 3 != 0:
					#dostuff 	1 -> +2bp
					#			2 -> +1bp
					current_transcript.checkLastVariants(ghandler, genetic_code, stopcodon)
					# +2 positions means 2 potential more variant effects.
					# repeatly, and check if deletions now have a new effect, because in rev cds it can be .... .... ....
					# and it is possible, that damaged transcripts (rev cds) is still functional because of this, meeeh fuck it
				if i == variant_transcript_exceeding_warning_int +1:
					print(str(current_transcript.TID) + " exceeds variation range limit of " + str(int(variant_transcript_exceeding_warning_int*100)) + " nucleotides.")
				if i == 20:
					current_transcript.Transcript_CDS_damaged = True
					print("transcript "+ str(current_transcript.TID) + " declared as broken, after extending it by: " +str(i*100) + " nucleotides." )
					break
			j = 0
			while current_transcript.origDNAtoshort:
				#print(current_transcript.TID)
				j +=1
				if j == 10:
					current_transcript.Transcript_CDS_damaged = True
					print("transcript " + str(current_transcript.TID) + " declared as broken, after extending it " + str(j) + " times unsuccessfully." )
					break
				current_transcript.checkLastVariants(ghandler, genetic_code, stopcodon)
		#print("Done: " + str(datetime.now() - time_start))
		####################################################
		print("Complete AA sequences")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		create_sequences(gff3, orig_AA, new_AA, genetic_code, name)
		#print("Done: " + str(datetime.now() - time_start))
		### Find Stop-Codons, if they are made by a frameshift -> LABEL the mutation
		print("Flag Frameshift caused stops")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		a_transcript_dict = gff3.get_chr_transcripts_dict(name)
		for current_transcript in a_transcript_dict.values():
			current_transcript = For_Type_Safety_and_statics.Transcript_Type_Safety(current_transcript)
			if not current_transcript.CDS_Exist:
				continue
			if current_transcript.Transcript_CDS_damaged:
				continue
			first_stop_position = (current_transcript.IV_ChangedTranslation.find(stopcodon)+1) * 3
			if first_stop_position == -2:
				continue
			last_frameshifter = ""
			for variant in current_transcript.IntegratedVariantObjects_CDS_Hits:
				variant = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(variant)
				if TranscriptEnum.STOP_GAINED in variant.Classification \
						or TranscriptEnum.STOP_CHANGED in variant.Classification:  # stop because of variants
					last_frameshifter = ""
					break
				if variant.Changed_CDS_Position < first_stop_position:
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
				last_frameshifter.STOP_CAUSED_IN = first_stop_position - last_frameshifter.Changed_CDS_Position
		#print("Done: " + str(datetime.now() - time_start))
		#################################################################
		print("Complete Data check.")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		complete_check(gff3, ghandler, name)
		#print("Done: " + str(datetime.now() - time_start))
		#################################################################
		print("Write data")
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		"""
		Old, but maybe still useful functions for further research/programming.
		The functions are now inside old_stuff()
		# Write_New_Special_VCF_File(output_data_path, gff3)
		# Write_All_VCF_Damaged_Transcripts(output_data_path,gff3)
		# Write_VCF_With_Key(output_data_path,gff3)
		# Write_All_VCF_Stop_Lost(output_data_path,gff3)
		# Write_All_VCF_NO_STOP(output_data_path,gff3)
		# Write_All_VCF_NO_START(output_data_path,gff3)
		# Write_All_VCF_Too_Many_Stops(output_data_path,gff3)
		"""
		write_all_VCF_file(output_data_path, gff3, old_info, ref_codons, ref_AA, alt_codons, alt_AA, extend_file, name)
		write_all_fasta(output_data_path, fasta_data_name, False, orig_DNA, orig_AA, new_DNA, new_AA, gff3, extend_file, name)
		print("Done: " + str(datetime.now() - time_start))
		print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
		extend_file = True
		gff3.free_RAM(name)
		vcf.free_RAM(name)
		gc.collect()
		#################################################################
	print("Create log files.")

	LogOrganizer.writeAllLogs(outpath)

	print("Everything is done: " + str(datetime.now() - time_start))
	print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

	return os.getpid()

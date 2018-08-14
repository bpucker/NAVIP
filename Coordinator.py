__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import copy #to copy objects
from GenomHandler import GenomHandler, Fasta_Enum
from VCF_Handler import VCF_HANDLER
from VCF_Variant import Variant, VariantEnum
from Transcript import Transcript, TranscriptEnum, For_Type_Safety_and_statics
from Gff3_Handler import GFF3_Handler_V3
from datetime import datetime
from LogOrganizer import LogOrganizer, LogEnums


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

	def old_stuff():
		"""
		Contains old functions, which may be useful again in the future.
		:return: Nothing.
		"""

		def Write_New_Special_VCF_File(data_path: str, gff3: GFF3_Handler_V3):
			# chrom	Pos	ID	Ref	Alt Info
			vcf_file = open(data_path + "Special_VCF_" + str(datetime.now()) + ".vcf", "w")
			vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tInfo\n")
			for name in gff3.GetChromosomeNames():
				for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
					transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)
					if not transcriptHier.CDS_Exist or transcriptHier.Transcript_CDS_damaged:
						continue
					if transcriptHier.IntegratedVariantObjects_CDS_Hits == []:
						continue

					# remove me maybe, after poster
					'''
					transcriptHier.Create_IV_OriginalTranslation(DictCodeSonne)
					transcriptHier.Create_IV_ChangedTranslation(DictCodeSonne)
	
					orgi_dna = ""
					orig_amino = transcriptHier.IV_OriginalTranslation
					new_dna = transcriptHier.IV_Changed_DNA_CDS_Seq
					new_amino = transcriptHier.IV_ChangedTranslation
					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
						orgi_dna = transcriptHier.Complete_CDS
	
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
						orgi_dna = transcriptHier.Rev_CDS
	
					vcf_file.write("#Transcription_Information: "+ str(transcriptHier.TID) + "  " + direction + " Orig_DNA, Orig_AA, Changed_DNA, Changed_AA\n")
					vcf_file.write(orgi_dna + "\n")
					vcf_file.write(orig_amino + "\n")
					vcf_file.write(new_dna + "\n")
					vcf_file.write(new_amino + "\n")
					'''
					# remove end
					first_vinfo = transcriptHier.IntegratedVariantObjects_CDS_Hits[0]
					For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(first_vinfo)
					variant_List_Of_Interest = []
					markedFirst = False
					markedSecond = False
					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
					else:
						direction = "Error:NOT Forward, NOT Reverse."

					for second_vinfo in transcriptHier.IntegratedVariantObjects_CDS_Hits:
						second_vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(second_vinfo)

						if second_vinfo == first_vinfo:
							continue
						'''
						elif first_vinfo.EndOfOwnEffect == VariantEnum.NO_EndOfOwnEffect and TranscriptEnum.STOP_GAINED not in first_vinfo.Classification:
							markedFirst = True
							markedSecond = True
						elif TranscriptEnum.STOP_GAINED in first_vinfo.Classification and TranscriptEnum.STOP_GAINED in second_vinfo.Classification:
							markedFirst = True
							markedSecond = True
						elif TranscriptEnum.STOP_GAINED in first_vinfo.Classification or TranscriptEnum.STOP_GAINED in second_vinfo.Classification:
							first_vinfo = second_vinfo
	
						elif first_vinfo.EndOfOwnEffect >= second_vinfo.StartOfOwnEffect:
							markedFirst = True
							markedSecond = True
						else:
							first_vinfo = second_vinfo
						'''
						if first_vinfo.EndOfOwnEffect == VariantEnum.NO_EndOfOwnEffect or second_vinfo.EndOfOwnEffect == VariantEnum.NO_EndOfOwnEffect:
							markedFirst = True
							markedSecond = True
						if markedFirst and first_vinfo not in variant_List_Of_Interest:
							variant_List_Of_Interest.append(first_vinfo)
							# chrom	Pos	ID	Ref	Alt Info
							classificationstring = ""
							class_list = first_vinfo.Classification
							class_list_length = len(class_list)
							for i in range(0, class_list_length):
								classificationstring += class_list[i].value
								if i != (class_list_length - 1):  # no tab after last entry
									classificationstring += "\t"

							vcf_file.write(str(name) + "\t" + str(first_vinfo.ChrPosition) + "\t"
										   + str(first_vinfo.ID) + "\t"
										   + str(first_vinfo.Ref) + "\t"
										   + str(first_vinfo.Alt) + "\t"
										   + str(transcriptHier.TID) + ": " + classificationstring + "\t"
										   + direction + "\t"
										   + "REF:" + str(first_vinfo.OrigTriplets) + "\t"
										   + str(first_vinfo.OrigRaster) + "\t"
										   + "ALT:" + str(first_vinfo.ChangedTriplets) + "\t"
										   + str(first_vinfo.Changed_Raster) + "\n")
						if markedSecond and second_vinfo not in variant_List_Of_Interest:
							variant_List_Of_Interest.append(second_vinfo)
							# chrom	Pos	ID	Ref	Alt Info
							classificationstring = ""
							class_list = second_vinfo.Classification
							class_list_length = len(class_list)
							for i in range(0, class_list_length):
								classificationstring += class_list[i].value
								if i != (class_list_length - 1):  # no tab after last entry
									classificationstring += "\t"

							vcf_file.write(str(name) + "\t" + str(second_vinfo.ChrPosition) + "\t"
										   + str(second_vinfo.ID) + "\t"
										   + str(second_vinfo.Ref) + "\t"
										   + str(second_vinfo.Alt) + "\t"
										   + str(transcriptHier.TID) + ": " + classificationstring + "\t"
										   + direction + "\t"
										   + "REF:" + str(second_vinfo.OrigTriplets) + "\t"
										   + str(second_vinfo.OrigRaster) + "\t"
										   + "ALT:" + str(second_vinfo.ChangedTriplets) + "\t"
										   + str(second_vinfo.Changed_Raster) + "\n")
						first_vinfo = second_vinfo
						markedFirst = False
						markedSecond = False
			vcf_file.close()

		def Write_VCF_With_Key(data_path: str, gff3: GFF3_Handler_V3):
			vcf_file = open(data_path + "Special_VCF_" + str(datetime.now()) + ".vcf", "w")
			vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tInfo\n")
			for name in gff3.GetChromosomeNames():
				for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
					transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)
					if not transcriptHier.CDS_Exist or transcriptHier.Transcript_CDS_damaged:
						continue
					if transcriptHier.IntegratedVariantObjects_CDS_Hits == []:
						continue
					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
					else:
						direction = "Error:NOT Forward, NOT Reverse."

					for vinfo in transcriptHier.IntegratedVariantObjects_CDS_Hits:
						vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
						if len(vinfo.SharedEffectsWith) == 0:
							continue

						classificationstring = ""
						shared_effect_list = "Shared Effects with:"
						shared_effect_list_len = len(vinfo.SharedEffectsWith)
						class_list = vinfo.Classification
						class_list_length = len(class_list)
						for i in range(0, class_list_length):
							classificationstring += class_list[i].value
							if i != (class_list_length - 1):  # no tab after last entry
								classificationstring += ","

						for i in range(0, shared_effect_list_len):
							vinfo_effects = vinfo.SharedEffectsWith[i]
							vinfo_effects = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo_effects)
							shared_effect_list += str(vinfo_effects.ChrPosition)
							if i != (shared_effect_list_len - 1):
								shared_effect_list += ","

						vcf_file.write(str(name) + "\t" + str(vinfo.ChrPosition) + "\t"
									   + str(vinfo.ID) + "\t"
									   + str(vinfo.Ref) + "\t"
									   + str(vinfo.Alt) + "\t"
									   + str(transcriptHier.TID) + ": " + str(classificationstring) + "\t"
									   + str(shared_effect_list) + "\t"
									   + direction + "\t"
									   + "REF:" + str(vinfo.OrigTriplets) + "\t"
									   + str(vinfo.OrigRaster) + "\t"
									   + "ALT:" + str(vinfo.ChangedTriplets) + "\t"
									   + str(vinfo.Changed_Raster) + "\n")
			vcf_file.close()

		def Write_All_VCF_Stop_Lost(data_path: str, gff3: GFF3_Handler_V3):
			vcf_file = open(data_path + "Special_VCF_Stop_Lost" + str(datetime.now()) + ".vcf", "w")
			vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tInfo\n")
			for name in gff3.GetChromosomeNames():
				vcf_file.write(str(name) + "\n")
				for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
					transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)
					if not transcriptHier.CDS_Exist or transcriptHier.Transcript_CDS_damaged:
						continue
					if transcriptHier.IntegratedVariantObjects_CDS_Hits == []:
						continue

					transcriptHier.Create_IV_OriginalTranslation(genetic_code)
					transcriptHier.Create_IV_ChangedTranslation(genetic_code)

					orgi_dna = ""
					orig_amino = transcriptHier.IV_OriginalTranslation
					new_dna = transcriptHier.IV_Changed_DNA_CDS_Seq
					new_amino = transcriptHier.IV_ChangedTranslation

					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
						orgi_dna = transcriptHier.Complete_CDS
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
						orgi_dna = transcriptHier.Rev_CDS
					else:
						direction = "Error:NOT Forward, NOT Reverse."

					if orig_amino.count("z") <= 1:
						continue

					# vcf_file.write(str(transcriptHier.TID + "\n"))
					# continue
					vcf_file.write("#Transcription_Information: " + str(
						transcriptHier.TID) + "  " + direction + " Orig_DNA, Orig_AA, Changed_DNA, Changed_AA\n")
					vcf_file.write(orgi_dna + "\n")
					vcf_file.write(orig_amino + "\n")
					vcf_file.write(new_dna + "\n")
					vcf_file.write(new_amino + "\n")

					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"

					for vinfo in transcriptHier.IntegratedVariantObjects_CDS_Hits:
						vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)

						classificationstring = ""
						class_list = vinfo.Classification
						class_list_length = len(class_list)
						for i in range(0, class_list_length):
							classificationstring += class_list[i].value
							if i != (class_list_length - 1):  # no tab after last entry
								classificationstring += ","

						vcf_file.write(str(name) + "\t" + str(vinfo.ChrPosition) + "\t"
									   + str(vinfo.ID) + "\t"
									   + str(vinfo.Ref) + "\t"
									   + str(vinfo.Alt) + "\t"
									   + str(transcriptHier.TID) + ": " + str(classificationstring) + "\t"
									   + direction + "\t"
									   + "REF:" + str(vinfo.OrigTriplets) + "\t"
									   + str(vinfo.OrigRaster) + "\t"
									   + "ALT:" + str(vinfo.ChangedTriplets) + "\t"
									   + str(vinfo.Changed_Raster) + "\n")
			vcf_file.close()

		def Write_All_VCF_NO_STOP(data_path: str, gff3: GFF3_Handler_V3):
			vcf_file = open(data_path + "Special_VCF_NO_STOP" + str(datetime.now()) + ".vcf", "w")
			# vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tInfo\n")
			vcf_file.write("#Transcripts completely without STOP\n")
			for name in gff3.GetChromosomeNames():
				vcf_file.write(str(name)
							   + "\t" + "Transcript/Source"
							   + "\t" + "RNA_Start"
							   + "\t" + "RNA_END"
							   + "\t" + "Gene_Start"
							   + "\t" + "Gene_End" + "\n")
				for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
					transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)
					if not transcriptHier.CDS_Exist:  # or transcriptHier.Transcript_CDS_damaged:
						continue
					# if transcriptHier.IntegratedVariantObjects == []:
					#	continue

					transcriptHier.Create_IV_OriginalTranslation(genetic_code)
					transcriptHier.Create_IV_ChangedTranslation(genetic_code)

					orgi_dna = ""
					orig_amino = transcriptHier.IV_OriginalTranslation
					new_dna = transcriptHier.IV_Changed_DNA_CDS_Seq
					new_amino = transcriptHier.IV_ChangedTranslation

					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
						orgi_dna = transcriptHier.Complete_CDS
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
						orgi_dna = transcriptHier.Rev_CDS
					else:
						direction = "Error:NOT Forward, NOT Reverse."

					if orig_amino.count("z") > 0:
						continue

					# if "hypothetical protein" in transcriptHier.Gene_Info_String:
					#	continue

					UTR_Write = ""
					Exon_Write = ""
					for UTR_Daten in transcriptHier.UTR_Description:
						UTR_Write += UTR_Daten

					for exon_Daten in transcriptHier.EXON_Descriptin:
						Exon_Write += exon_Daten

					vcf_file.write(str(name)
								   + "\t" + str(transcriptHier.TID)
								   + "\t" + str(transcriptHier.StartOfRNA)
								   + "\t" + str(transcriptHier.EndOfRNA)
								   + "\t" + str(transcriptHier.Gene_Start_Position)
								   + "\t" + str(transcriptHier.Gene_End_Position)
								   + "\t" + transcriptHier.Gene_Info_String
								   + UTR_Write
								   + Exon_Write + "\n")

					continue

			vcf_file.close()

		def Write_All_VCF_To_Many_Stops(data_path: str, gff3: GFF3_Handler_V3):
			vcf_file = open(data_path + "Special_VCF_To_Many_Stops" + str(datetime.now()) + ".vcf", "w")
			# vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tInfo\n")
			vcf_file.write("#Transcripts with 2 or more STOPs\n")
			for name in gff3.GetChromosomeNames():
				vcf_file.write(str(name)
							   + "\t" + "Transcript/Source"
							   + "\t" + "RNA_Start"
							   + "\t" + "RNA_END"
							   + "\t" + "Gene_Start"
							   + "\t" + "Gene_End" + "\n")
				for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
					transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)
					if not transcriptHier.CDS_Exist:  # or transcriptHier.Transcript_CDS_damaged:
						continue
					# if transcriptHier.IntegratedVariantObjects == []:
					#	continue

					transcriptHier.Create_IV_OriginalTranslation(genetic_code)
					transcriptHier.Create_IV_ChangedTranslation(genetic_code)

					orgi_dna = ""
					orig_amino = transcriptHier.IV_OriginalTranslation
					new_dna = transcriptHier.IV_Changed_DNA_CDS_Seq
					new_amino = transcriptHier.IV_ChangedTranslation

					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
						orgi_dna = transcriptHier.Complete_CDS
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
						orgi_dna = transcriptHier.Rev_CDS
					else:
						direction = "Error:NOT Forward, NOT Reverse."

					if orig_amino.count("z") <= 1:
						continue
					# if "hypothetical protein" in transcriptHier.Gene_Info_String:
					#	continue

					UTR_Write = ""
					Exon_Write = ""
					for UTR_Daten in transcriptHier.UTR_Description:
						UTR_Write += UTR_Daten

					for exon_Daten in transcriptHier.EXON_Descriptin:
						Exon_Write += exon_Daten

					vcf_file.write(str(name)
								   + "\t" + str(transcriptHier.TID)
								   + "\t" + str(transcriptHier.StartOfRNA)
								   + "\t" + str(transcriptHier.EndOfRNA)
								   + "\t" + str(transcriptHier.Gene_Start_Position)
								   + "\t" + str(transcriptHier.Gene_End_Position)
								   + "\t" + transcriptHier.Gene_Info_String
								   + UTR_Write
								   + Exon_Write + "\n")
					continue
			vcf_file.close()

		def Write_All_VCF_NO_START(data_path: str, gff3: GFF3_Handler_V3):
			vcf_file = open(data_path + "Special_VCF_NO_START" + str(datetime.now()) + ".vcf", "w")
			# vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tInfo\n")
			vcf_file.write("#Transcripts without START in the beginning codon\n")
			for name in gff3.GetChromosomeNames():
				vcf_file.write(str(name)
							   + "\t" + "Transcript/Source"
							   + "\t" + "RNA_Start"
							   + "\t" + "RNA_END"
							   + "\t" + "Gene_Start"
							   + "\t" + "Gene_End" + "\n")
				for transcriptHier in gff3.GetChrTranscriptsDict(name).values():
					transcriptHier = For_Type_Safety_and_statics.Transcript_Type_Safety(transcriptHier)
					if not transcriptHier.CDS_Exist:  # or transcriptHier.Transcript_CDS_damaged:
						continue
					# if transcriptHier.IntegratedVariantObjects == []:
					#	continue

					transcriptHier.Create_IV_OriginalTranslation(genetic_code)
					transcriptHier.Create_IV_ChangedTranslation(genetic_code)

					orgi_dna = ""
					orig_amino = transcriptHier.IV_OriginalTranslation
					new_dna = transcriptHier.IV_Changed_DNA_CDS_Seq
					new_amino = transcriptHier.IV_ChangedTranslation

					if transcriptHier.ForwardDirection == TranscriptEnum.FORWARD:
						direction = "Forward"
						orgi_dna = transcriptHier.Complete_CDS
					elif transcriptHier.ForwardDirection == TranscriptEnum.REVERSE:
						direction = "Reverse"
						orgi_dna = transcriptHier.Rev_CDS
					else:
						direction = "Error:NOT Forward, NOT Reverse."

					if orig_amino[0] == 'm':
						continue
					# if "hypothetical protein" in transcriptHier.Gene_Info_String:UTR_Description
					#	continue

					UTR_Write = ""
					Exon_Write = ""
					for UTR_Daten in transcriptHier.UTR_Description:
						UTR_Write += UTR_Daten

					for exon_Daten in transcriptHier.EXON_Descriptin:
						Exon_Write += exon_Daten

					vcf_file.write(str(name)
								   + "\t" + str(transcriptHier.TID)
								   + "\t" + str(transcriptHier.StartOfRNA)
								   + "\t" + str(transcriptHier.EndOfRNA)
								   + "\t" + str(transcriptHier.Gene_Start_Position)
								   + "\t" + str(transcriptHier.Gene_End_Position)
								   + "\t" + transcriptHier.Gene_Info_String
								   + UTR_Write
								   + Exon_Write + "\n")
					continue
			vcf_file.close()



	def Write_All_VCF_File(data_path: str,
						   gff3: GFF3_Handler_V3,
						   Old_Info: bool,
						   Ref_Codons: bool,
						   Ref_AA: bool,
						   Alt_Codons: bool,
						   Alt_AA: bool):
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
		vcf_file = open(data_path + "All_VCF" + ".vcf", "w")
		vcf_file.write("#All Data Output\n")
		vcf_file.write("#Please note, that the Variant_Position_in_Codon is read from left to right in forward "
					   "and rigth to left in reverse strand direction.\n")
		vcf_file.write("#Chrom\tPos\tID\tRef\tAlt\tQual\tFilter\tInfo\n")
		vcf_file.write("#When complete output is allowed, the Info field is constructed like this:\n")
		vcf_file.write("#Info:TranscriptID|"
					   "Strand_Direction|"
					   "Variant_Classification1,Variant_Classification2,...|"
					   "Shared_EffKey(s)|"
					   "REF_Codon(s);Variant_Position_in_Codon|"
					   "REF_AA|"
					   "old_CDS_Position|"
					   "ALT_Codon(s);Variant_Position_in_Codon|"
					   "ALT_AA|"
					   "new_CDS_Position|"
					   "NAVIP_END|"
					   "<old info field>\n")
		vcf_file.write("#If there are no shared effect keys, the value is:\"NONE\".\n")
		data_to_write = []
		for name in gff3.GetChromosomeNames():
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
					classificationstring = ""
					class_list = vinfo.Classification
					class_list_length = len(class_list)
					for i in range(0, class_list_length):
						classificationstring += class_list[i].value
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
					NAVIP_INFO_LIST = [str(transcriptHier.TID) + "|", str(direction) + "|", classificationstring]
					if shared_effect_list == "":
						NAVIP_INFO_LIST.append("NONE|")
					else:
						NAVIP_INFO_LIST.append(str(shared_effect_list))
					if Ref_Codons:
						NAVIP_INFO_LIST.append(str(vinfo.OrigTriplets) + ";" + str(vinfo.OrigRaster) + "|")
					if Ref_AA:
						NAVIP_INFO_LIST.append(str(vinfo.OrigAmino) + "|")
					NAVIP_INFO_LIST.append(str(vinfo.Unchanged_CDS_Position) + "|")
					if Alt_Codons:
						NAVIP_INFO_LIST.append(str(vinfo.ChangedTriplets) + ";" + str(vinfo.Changed_Raster) + "|")
					if Alt_AA:
						NAVIP_INFO_LIST.append(str(vinfo.NewAmino) + "|")
					NAVIP_INFO_LIST.append(str(vinfo.Changed_CDS_Position) + "|")
					NAVIP_INFO_LIST.append("NAVIP_END|")
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
						+ str(NAVIP_INFO))
			data_to_write = sorted(data_to_write,
								   key=lambda data_line: (int(data_line.split("\t")[1]), str(data_line.split("\t")[7])))
			vcf_file.write("".join(data_to_write))
			data_to_write = []
		vcf_file.close()

	def Create_Sequences(gff3: GFF3_Handler_V3, Orig_AA: bool, New_AA: bool, genetic_code: dict):
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
		for chrName in gff3.GetChromosomeNames():
			print(chrName)
			for currentTranscript in gff3.GetChrTranscriptsDict(chrName).values():
				currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)

				if not currentTranscript.CDS_Exist:
					continue
				if Orig_AA:
					currentTranscript.Create_IV_OriginalTranslation(genetic_code)

				if currentTranscript.Transcript_CDS_damaged:
					continue
				if New_AA:
					currentTranscript.Create_IV_ChangedTranslation(genetic_code)

	def Complete_Check(gff3: GFF3_Handler_V3, ghandler: GenomHandler):
		"""
		Its always a good idea to test the data in the end.
		Will print warnings, if there is incorrect data.
		:param gff3: The gff3 data structure with all transcripts and their variants/sequences.
		:param ghandler: The reference genome sequences for every chromosome is inside.
		:return: Nothing.
		"""
		for chrName in Shared_Chromosomes_FA_GFF3:
			Chr_Transcript_List = gff3.GetChrTranscriptsDict(chrName).values()
			for currentTranscript in Chr_Transcript_List:
				currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)
				if not currentTranscript.CDS_Exist or currentTranscript.Transcript_CDS_damaged:
					continue

				for vinfo in currentTranscript.IntegratedVariantObjects_CDS_Hits:
					vinfo = For_Type_Safety_and_statics.Variant_Information_Storage_Type_Safety(vinfo)
					#try :
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
									vinfo.ChrPosition):
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
									vinfo.ChrPosition):
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

						if vinfo.Ref != ghandler.seq(chrName, vinfo.ChrPosition, vinfo.ChrPosition + len(vinfo.Ref) - 1):
							LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
								"DELETION_REF != Fasta-Seq: " + str(currentTranscript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

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
													vinfo.ChrPosition + len(vinfo.Ref) - 1):
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
						LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,"No Classification in Complete_Check?"+ "\n")
					#except IndexError:
					#	LogOrganizer.addToLog(LogEnums.COORDINATOR_COMPLETE_CHECK_CRITICAL_LOG, str(vinfo.ChrPosition) + "\t" + str(currentTranscript.TID))
	def Write_All_Fasta(data_path: str,
						data_name: str,
						date_time: bool,
						Orig_DNA: bool,
						Orig_AA: bool,
						New_DNA: bool,
						New_AA: bool,
						gff3: GFF3_Handler_V3):
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
			normal_data_list_to_write = [
				"#ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == New Amino Acid Sequence\n"]
			Transcript_CDS_damaged_data_list = [
				"#Data in here may be incomplete and incorrect because of deletions, which destroyed parts the splice sides.\n"]

			for chrName in gff3.GetChromosomeNames():
				for currentTranscript in gff3.GetChrTranscriptsDict(chrName).values():
					currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)

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

					if New_DNA:
						normal_data_list_to_write.append(
							">" + str(currentTranscript.TID) + "|" + str(Fasta_Enum.nDNA.value) + "\n")
						normal_data_list_to_write.append(currentTranscript.IV_Changed_DNA_CDS_Seq + "\n")
					if New_AA:
						normal_data_list_to_write.append(
							">" + str(currentTranscript.TID) + "|" + str(Fasta_Enum.nAA.value) + "\n")
						normal_data_list_to_write.append(str(currentTranscript.IV_ChangedTranslation) + "\n")

			normal_output = "".join(normal_data_list_to_write)
			New_Fasta_File = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "w")
			New_Fasta_File.write(normal_output)
			New_Fasta_File.close()
			normal_output = ""
			normal_data_list_to_write = []

			damaged_output = "".join(Transcript_CDS_damaged_data_list)
			New_Fasta_File_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt", "w")
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

	#######################################
	print("read vcf")
	vcf = VCF_HANDLER(vcf_path_and_name, firstxChr)
	print("Done: " + str(datetime.now() - timeStart))
	#######################################
	print("read gff3")
	#print("Hopefully sorted after seqID")
	gff3 = GFF3_Handler_V3(gff3_path_and_name)
	print("Done: " + str(datetime.now() - timeStart))
	#######################################
	print("read fa")
	ghandler = GenomHandler(fasta_FILE_PATH, firstxChr)
	print("Done: " + str(datetime.now() - timeStart))
	#######################################
	print("Chromosome:Names")
	print("GFF3:" + str(gff3.GetChromosomeNames()))
	print("VCF:" + str(vcf.GetChromosomeNames()))
	print("fasta:" + str(ghandler.GetChromosomeNames()))
	Shared_Chromosomes = []
	Shared_Chromosomes_FA_GFF3 = []
	#######################################

	for name in gff3.GetChromosomeNames():
		if name in vcf.GetChromosomeNames() and name in ghandler.GetChromosomeNames():
			Shared_Chromosomes.append(name)
		if name in ghandler.GetChromosomeNames():
			Shared_Chromosomes_FA_GFF3.append(name)
	print("Shared_Chromosomes(all):" + str(Shared_Chromosomes))
	#print("Shared_Chromosomes, fasta, gff3 (and used):" + str(Shared_Chromosomes_FA_GFF3))

	### For Writing the CDS
	### in all transcripts

	print("Write the CDS and Rev CDS")

	for name in Shared_Chromosomes:
		print(name)
		#chrID = gff3.GetChromosomeID(name)

		#for i in gff3.dictListOfTranscripts[gff3.dictChrNames[name]]:  # only a number, not the object
		#	trancripts[i].CompleteTheCDS(ghandler.seq(name, trancripts[i].StartOfRNA, trancripts[i].EndOfRNA))
		#	if trancripts[i].ForwardDirection == TranscriptEnum.REVERSE:
		#		trancripts[i].ReverseTheCDS()

		for trancript in gff3.GetChrTranscriptsDict(name).values():
			trancript.CompleteTheCDS(ghandler.seq(name, trancript.StartOfRNA, trancript.EndOfRNA))
			if trancript.ForwardDirection == TranscriptEnum.REVERSE:
				trancript.ReverseTheCDS()


	print("Done: " + str(datetime.now() - timeStart))



	###
	# Connecting all transcripts with all their mutations/variants.
	# Contains handling of trialllele variants, if an original vcf file was used.
	# But it can't handle all possible appearances. VCF preprocessing exist for a reason.
	# The connection triggers a first evaluation of the variants, too.
	###
	transcriptRange = 300 # base pairs before RNA starts and after RNA ends. Just in Case for long InDels.
	print("Connect transcripts with variants")
	for name in Shared_Chromosomes:
		print("Search in " + name)

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
						# not necessary with vcf-preprocessing - but maybe the data change in the future -> usefull again
						transcriptsWithMultiAllelVariants.append(currentTranscript)
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
		for currentTranscript in transcriptsWithMultiAllelVariants:
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

	print("Done: " + str(datetime.now() - timeStart))

	########
	### Calculate the effect length of mutation effects (version 2)
	###	Creates parts of the changed DNA- and AA-sequence.
	###	The transcripts will be extended, until a stop appears (if there is none).
	########
	print("Calculate the effect length")

	for name in Shared_Chromosomes:
		print(name)
		aTranscriptDict = gff3.GetChrTranscriptsDict(
			name)  # self.dictListOfTranscripts = [{}] # List for chromosomes, dict for normal entrys
		for currentTranscript in aTranscriptDict.values():
			currentTranscript = For_Type_Safety_and_statics.Transcript_Type_Safety(currentTranscript)

			if not currentTranscript.CDS_Exist:
				continue
			if currentTranscript.Transcript_CDS_damaged:
				continue

			currentTranscript.Create_IV_Changed_DNA_CDS_Seq(genetic_code, currentTranscript.IntegratedVariantObjects_CDS_Hits,
															stopcodon)
			while currentTranscript.Lost_Stop:
				currentTranscript.Create_IV_ChangedTranslation(genetic_code)
				if stopcodon in currentTranscript.IV_ChangedTranslation:
					# print("new stopcodon already inside transcript" + str(transcriptHier.TID))

					#LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"New Stopcodon is in " + str(currentTranscript.TID) + "\n")
					lastPosInCDS = currentTranscript.LastCDSPosition
					currentTranscript.Lost_Stop = False
					currentTranscript.Found_New_Stop = True
					break
				# print ("Lost_Stop not completely inside transcript.")
				lastPosInCDS = currentTranscript.LastCDSPosition
				if currentTranscript.ForwardDirection == TranscriptEnum.FORWARD:
					nextDNA = ghandler.seq(name, lastPosInCDS + 1, lastPosInCDS + 101)
				elif currentTranscript.ForwardDirection == TranscriptEnum.REVERSE:
					nextDNA = ghandler.seq(name, lastPosInCDS - 101, lastPosInCDS - 1)
					#nextDNA = For_Type_Safety_and_statics.ReverseSeq(nextDNA) #reverse, because it will be added to the already reversed dna transcript
				else:
					LogOrganizer.addToLog(LogEnums.COORDINATOR_BUGHUNTING_LOG,"Error: No Direction: " + str(currentTranscript.TID)+"\n")
					break
				currentTranscript.Find_New_Stop(nextDNA, genetic_code, stopcodon)

	print("Done: " + str(datetime.now() - timeStart))
	####################################################
	print("Complete AA sequences (when enabled)")
	Create_Sequences(gff3, Orig_AA, New_AA, genetic_code)
	print("Done: " + str(datetime.now() - timeStart))
	#################################################################
	print("Complete Data check.")
	Complete_Check(gff3, ghandler)
	print("Done: " + str(datetime.now() - timeStart))
	#################################################################

	print("Write all data")
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
	Write_All_VCF_File(Output_Data_Path, gff3, Old_Info, Ref_Codons, Ref_AA, Alt_Codons, Alt_AA)
	Write_All_Fasta(Output_Data_Path,Fasta_Data_Name,False,Orig_DNA,Orig_AA,New_DNA,New_AA,gff3)
	print("Done: " + str(datetime.now() - timeStart))
	#################################################################
	print("Create log files.")

	LogOrganizer.writeAllLogs(outpath)

	print("Everything is done: " + str(datetime.now() - timeStart))



	#2232
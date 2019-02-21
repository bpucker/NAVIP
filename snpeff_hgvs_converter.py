from Transcript import *
from VCF_Variant import *
import math

class snpeff_hgvs_converter():

	def __init__(self):
		pass
	# Info:TranscriptID|Strand_Direction|Variant_Classification1,Variant_Classification2,...|Shared_EffKey(s)|
	# REF_Codon(s);Variant_Position_in_Codon|REF_AA|old_CDS_Position|ALT_Codon(s);Variant_Position_in_Codon|
	# ALT_AA|new_CDS_Position|NAVIP_END|<old info field>
	# Chr1	12584	.	A	C	2961.77	PASS	AT1G01030.1|REV|SUB|NONE|AAC;2|v|357|CAC;2|v|357|NAVIP_END|

	#																			0		1			2
	##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact |
	#	3			4			5			6				7					8		9		10		11.1	11.2
	# Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length |
	#	12.1	12.2		13.1		13.2		14			15
	# CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">


	# 0			1				   2			3			4			5				6			7				8
	# ['C', 'synonymous_variant', 'LOW', 'AT1G01030', 'AT1G01030', 'transcript', 'AT1G01030.1', 'protein_coding', '2/2',
	#	9			10				11			12			13		14
	# 'c.357T>G', 'p.Val119Val', '970/1905', '357/1077', '119/358', '',
	# 0			1					2			3	4			5			6			7			8		9
	#ANN=A|upstream_gene_variant|MODIFIER|AT1G01010|AT1G01010|transcript|AT1G01010.1|protein_coding|	|c.-3408G>A|
	# 10	11	12	13	14	  	15
	# 	|		|	|	|3279|	,A|upstream_gene_variant|MODIFIER|AT1G01010.1|AT1G01010.1-Protein|transcript|TRANSCRIPT_AT1G01010.1-Protein|
	# protein_coding||c.-3408G>A|||||3408|,A|intergenic_region|MODIFIER|CHR_START-AT1G01010|CHR_START-AT1G01010|intergenic_region|
	# CHR_START-AT1G01010|||n.352G>A||||||

	##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
	##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">


	def convert_main(self, transcript:Transcript, vinfo:Variant_Information_Storage) -> str :

		#note to myself: don't forget: info-files are ;-seperated
		aminodict = {}
		aminodict["Ala"] = "A"
		aminodict["A"] = "Ala"
		aminodict["Cys"] = "C"
		aminodict["C"] = "Cys"
		aminodict["Asp"] = "D"
		aminodict["D"] = "Asp"
		aminodict["Glu"] = "E"
		aminodict["E"] = "Glu"
		aminodict["Phe"] = "F"
		aminodict["F"] = "Phe"
		aminodict["Gly"] = "G"
		aminodict["G"] = "Gly"
		aminodict["His"] = "H"
		aminodict["H"] = "His"
		aminodict["Ile"] = "I"
		aminodict["I"] = "Ile"
		aminodict["Lys"] = "K"
		aminodict["K"] = "Lys"
		aminodict["Leu"] = "L"
		aminodict["L"] = "Leu"
		aminodict["Met"] = "M"
		aminodict["M"] = "Met"
		aminodict["Asn"] = "N"
		aminodict["N"] = "Asn"
		aminodict["Pro"] = "P"
		aminodict["P"] = "Pro"
		aminodict["Gln"] = "Q"
		aminodict["Q"] = "Gln"
		aminodict["Arg"] = "R"
		aminodict["R"] = "Arg"
		aminodict["Ser"] = "S"
		aminodict["S"] = "Ser"
		aminodict["Thr"] = "T"
		aminodict["T"] = "Thr"
		aminodict["Val"] = "V"
		aminodict["V"] = "Val"
		aminodict["Trp"] = "W"
		aminodict["W"] = "Trp"
		aminodict["Tyr"] = "Y"
		aminodict["Y"] = "Tyr"
		aminodict["*"] = "*"
		#aminodict["*"] = "stop_gained"
		#aminodict["stop_gained"] = "*"


		Annotation = "" 						#1
		Annotation_Impact = "" 					#2
		Gene_Name = transcript.TID.split(".")[0]#3
		Gene_ID = transcript.TID.split(".")[0] 	#4
		Feature_Type = "transcript" 			#5
		Feature_ID = transcript.TID 			#6
		Transcript_BioType = "protein_coding" 	#7
		Rank = "" 								#8		# ignored
		HGVS_C = "" 							#9		# http://varnomen.hgvs.org/recommendations/DNA/
		HGVS_P = "" 							#10		# http://varnomen.hgvs.org/recommendations/protein/
		cDNA_pos = "" 							#11.1 	# cDNA_pos != CDS_pos, i'm using the matured positions only
		cDNA_length = "" 						#11.2 	# cDNA_length != CDS_length, i'm using the matured positions only
		CDS_pos = vinfo.Unchanged_CDS_Position 	#12.1	# original in snpeff (?); me original
		CDS_length = transcript.uChDNA_length 	#12.2	# original in snpeff (?); me original
		AA_pos = (vinfo.Unchanged_CDS_Position-1)/3	#13.1	# original in snpeff (?); me original
		AA_length = transcript.uChAA_length 	#13.2	# original in snpeff (?); me original
		distance = "" 							#14		# ignored
		errors_warnings = "" 					#15		# maybe ignored

		aa_pos_temp = (vinfo.Unchanged_CDS_Position-1)%3
		if aa_pos_temp == 0:
			AA_pos = (vinfo.Unchanged_CDS_Position-1) / 3 + 1 # aa = 0 is the first, not zero
		elif aa_pos_temp == 1:
			AA_pos = (vinfo.Unchanged_CDS_Position +1) / 3 + 1
		elif aa_pos_temp == 2:
			AA_pos = (vinfo.Unchanged_CDS_Position) / 3 + 1
		AA_pos = int(AA_pos)

		AA_to_next_stop = -1
		new_aa_pos = 0
		aa_pos_temp = (vinfo.Changed_CDS_Position - 1) % 3
		if aa_pos_temp == 0:
			new_aa_pos = (vinfo.Changed_CDS_Position - 1) / 3 + 1  # aa = 0 is the first, not zero
		elif aa_pos_temp == 1:
			new_aa_pos = (vinfo.Changed_CDS_Position + 1) / 3 + 1
		elif aa_pos_temp == 2:
			new_aa_pos = (vinfo.Changed_CDS_Position) / 3 + 1
		new_aa_pos = int(new_aa_pos)
		AA_to_next_stop = transcript.IV_ChangedTranslation[new_aa_pos:].find('*')

		for i,x in enumerate(vinfo.Classification):
			if i < len(vinfo.Classification) -1 :
				Annotation += x.value + ',' # not efficient, but list isn't large
			else:
				Annotation += x.value

		if TranscriptEnum.FRAMESHIFT in vinfo.Classification \
				or TranscriptEnum.FRAMESHIFT_1 in vinfo.Classification \
				or TranscriptEnum.FRAMESHIFT_2 in vinfo.Classification \
				or TranscriptEnum.FRAMESHIFT_1_DEL in vinfo.Classification \
				or TranscriptEnum.FRAMESHIFT_2_DEL in vinfo.Classification \
				or TranscriptEnum.STOP_GAINED in vinfo.Classification \
				or TranscriptEnum.STOP_LOST in vinfo.Classification \
				or TranscriptEnum.START_LOST in vinfo.Classification:
			Annotation_Impact = "HIGH"
		elif TranscriptEnum.DELETION in vinfo.Classification \
				or TranscriptEnum.INSERTION in vinfo.Classification \
				or TranscriptEnum.AA_CHANGE in vinfo.Classification:
			Annotation_Impact = "MODERATE"

		elif TranscriptEnum.STOP_CHANGED  in vinfo.Classification\
				or TranscriptEnum.AA_CHANGE not in vinfo.Classification:
			Annotation_Impact = "LOW"
		else:
			print("Impossible case. Error convert_main.")

		#Annotation_Impact = "MODIFIER" # does not exist? upstream/downstream is not noted


		#HGVS_C
		if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
			# NG_012232.1(NM_004006.1):c.93+1G>T
    		# a substitution of the G nucleotide at c.93+1 (coding DNA reference sequence) by a T
			if TranscriptEnum.FORWARD == transcript.ForwardDirection:
				HGVS_C = "c." + str(vinfo.Unchanged_CDS_Position) + vinfo.Ref + ">" + vinfo.Alt
			else:
				HGVS_C = "c." + str(vinfo.Unchanged_CDS_Position) + vinfo.ReverseRef + ">" + vinfo.ReverseAlt
		elif TranscriptEnum.INSERTION in vinfo.Classification:
			#insertion, duplication, (repeated seq)
			#check for duplication
			if vinfo.Unchanged_CDS_Position >= len(vinfo.Alt) -1: # check if insertion is long enough to be a duplicate
				if transcript.ForwardDirection == TranscriptEnum.FORWARD:
					insert = vinfo.Alt[1:] # first base == refstart, always > 0
				else:
					insert = vinfo.ReverseAlt[1:]

				if len(transcript.uChDNAsequence) >= len(transcript.Complete_CDS):
					refpart = transcript.uChDNAsequence[vinfo.Unchanged_CDS_Position - len(insert):vinfo.Unchanged_CDS_Position +1]
					#print(str(transcript.TID) + "\t" + insert + "\t" + refpart)
					#print("testa1:\t" + vinfo.Alt + "\t" + transcript.uChDNAsequence[vinfo.Unchanged_CDS_Position-1])
					#print("testa2:\t" + vinfo.Alt + "\t" + refpart)
					#print("testa1:\t" + vinfo.ReverseAlt + "\t" + transcript.uChDNAsequence[vinfo.Unchanged_CDS_Position - 1])
					#print("testa2:\t" + vinfo.ReverseAlt + "\t" + refpart)
					"""
						AT5G67350.1	AGC	AGCT
						testa1:	GCTT	C
						testa2:	GCTT	AGCT
						testa1:	AAGC	C
						testa2:	AAGC	AGCT
						c._101insAGC
					"""
				else:
					if TranscriptEnum.FORWARD == transcript.ForwardDirection:
						refpart = transcript.Complete_CDS[vinfo.Unchanged_CDS_Position + 1 - len(insert):vinfo.Unchanged_CDS_Position + 1]
						#print(str(transcript.TID) + "\t" + insert + "\t" + refpart)
						#print("testb1:\t" + vinfo.Alt + "\t" + transcript.Complete_CDS[vinfo.Unchanged_CDS_Position - 1])
						#print("testb2:\t" + vinfo.Alt + "\t" + refpart)
					else:
						refpart = transcript.Rev_CDS[vinfo.Unchanged_CDS_Position + 1 - len(insert):vinfo.Unchanged_CDS_Position + 1]
						#print(str(transcript.TID) + "\t" + insert + "\t" + refpart)
						#print("testc1:\t" + vinfo.Alt + "\t" + transcript.Rev_CDS[vinfo.Unchanged_CDS_Position - 1])
						#print("testc2:\t" + vinfo.Alt + "\t" + refpart)

				HGVS_C = "c."
				if insert == refpart:
					#Format: “prefix”“position(s)_duplicated”“dup”, e.g. g.123_345dup
					vinfo.Classification.append(hgvs_dna.DUP)
					HGVS_C += str(vinfo.Unchanged_CDS_Position - len(insert))
					if len(insert) == 1:
						HGVS_C += 'dup' + insert
					else:
						HGVS_C += '_' + str(vinfo.Unchanged_CDS_Position) + 'dup' + insert
				else:
					#Format: “prefix”“positions_flanking”“ins”“inserted_sequence”, e.g. g.123_124insAGC
					HGVS_C += '_' + str(vinfo.Unchanged_CDS_Position) + 'ins' + insert
				#print(HGVS_C)

		elif TranscriptEnum.DELETION in vinfo.Classification:
			#Format: “prefix”“position(s)_deleted”“del”, e.g. g.123_127del
			HGVS_C = "c."
			if len(vinfo.Ref) == 2: #del of 1 base
				#one nucleotide - NG_012232.1:g.19del
				# a deletion of the T at position g.19 in the sequence AGAATCACA to AGAA_CACA
    			# NOTE: it is allowed to describe the variant as NG_012232.1:g.19delT
				HGVS_C += str(vinfo.Unchanged_CDS_Position +1) + 'del' + vinfo.Ref[1]
			else:
				#NG_012232.1:g.19_21del
    			#a deletion of nucleotides g.19 to g.21 in the sequence AGAATCACA to AGAA___CA
    			#NOTE: it is allowed to describe the variant as NG_012232.1:g.19_21delTCA
				HGVS_C += str(vinfo.Unchanged_CDS_Position +1) + '_' + str(vinfo.Unchanged_CDS_Position + len(vinfo.Ref)) + 'del' + vinfo.Ref[1:]
			#print(HGVS_C)
		else:
			print("Something unexpected happened in the snpeff_hgvs_converter converter.")

		#HGVS_P
		HGVS_P = "p."
		newamino = ''
		oldamino= ''
		for aa in vinfo.NewAmino:
			newamino += aminodict[aa.upper()]
		for aa in vinfo.OrigAmino:
			oldamino += aminodict[aa.upper()]

		if TranscriptEnum.STOP_LOST in vinfo.Classification:
			# Format (C-terminal): “prefix”“Ter_position”“new_amino_acid”“ext”“position_new_termination_site”, e.g. p.Ter110Glnext*17
			# extension, sub, del, ins
			# sub first, because its easy
			ext = ""
			if AA_to_next_stop != -1:
				ext = 'ext*' + str(AA_to_next_stop)
			else:
				ext = 'ext*?'
			HGVS_P += '*' + str(AA_pos) + vinfo.NewAmino + ext

			#if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
			#	HGVS_P += '*' + str(AA_pos) + vinfo.NewAmino + ext
			#elif TranscriptEnum.INSERTION in vinfo.Classification:
			#	HGVS_P += '*' + str(AA_pos) + vinfo.NewAmino + ext
			#elif TranscriptEnum.DELETION in vinfo.Classification:
			#	HGVS_P += '*' + str(AA_pos) + vinfo.NewAmino + ext
		elif TranscriptEnum.STOP_CHANGED in vinfo.Classification and TranscriptEnum.INSERTION in vinfo.Classification and vinfo.NewAmino[0] != '*':
			# insertion in the stop codon, but still a stop codon exist -> extension, but only if the first AA
			# is not a stop codon
			# extensioncheck
			if AA_to_next_stop != -1:
				ext = 'ext*' + str(AA_to_next_stop)
			else:
				ext = 'ext*?'
				print('mhhh')
			HGVS_P += '*' + str(AA_pos) + newamino + ext

		elif TranscriptEnum.SUBSTITUTION in vinfo.Classification or newamino == '*':
			# missense -> AA-Change #LRG_199p1:p.Trp24Cys
			# nonsense -> stop gained #LRG_199p1:p.Trp24Ter (p.Trp24*)
			# silent -> AA is not changed #NP_003997.1:p.Cys188=
			if vinfo.OrigAmino == vinfo.NewAmino: #silent
				HGVS_P += oldamino + str(AA_pos) + '='
			else: #missense and nonsense
				HGVS_P += oldamino + str(AA_pos) + newamino

		elif TranscriptEnum.FRAMESHIFT_1 in vinfo.Classification \
			or TranscriptEnum.FRAMESHIFT_2 in vinfo.Classification \
			   or TranscriptEnum.FRAMESHIFT_1_DEL in vinfo.Classification \
			   or TranscriptEnum.FRAMESHIFT_2_DEL in vinfo.Classification:
			"""
			p.Arg97ProfsTer23 (short p.Arg97fs)
    		a variant with Arg97 as the first amino acid changed, 
    		shifting the reading frame, replacing it for a Pro and terminating at position Ter23.
			"""
			if newamino == '*': # counts as SUB
				pass
			elif '*' in newamino: #no idea what this is, but not always a frameshift
				pass
			else:
				#first changed aa has to be the first aa here
				if oldamino[0] != newamino[0]:
					HGVS_P += oldamino + str(AA_pos) + newamino + 'Ter' + str(vinfo.STOP_CAUSED_IN)
				else:
					#find first changed aa....
					# test every aa, if they are equal, from the variant_position
					i = 0
					for oldaa,newaa in  transcript.uChAAsequence[AA_pos:], transcript.IV_ChangedTranslation[new_aa_pos:]:
						if oldaa == "":
							print("should not happen")
						elif oldaa == newaa:
							i +=1
							continue
						else:
							if vinfo.STOP_CAUSED_IN == -1:
								HGVS_P += oldaa + str(AA_pos + i) + newaa + 'Ter?'
							else:
								HGVS_P += oldaa + str(AA_pos + i) + newaa + 'Ter' + str(vinfo.STOP_CAUSED_IN - i)

		elif TranscriptEnum.INSERTION in vinfo.Classification:
			# extension is already checked ->
			# frameshift, too
			# first duplicationcheck
			# second and last is the normal insertion
			pass
		elif TranscriptEnum.DELETION in vinfo.Classification:
			# deletionstuff
			pass
		else:
			print('what did i forget?')



		snpeff_like_info_string = 'ANN=' + vinfo.Alt + '|' \
								  + Annotation + "|" \
								  + Annotation_Impact + "|" \
								  + Gene_Name + "|" \
								  + Gene_ID + "|" \
								  + Feature_Type + "|" \
								  + Feature_ID + "|" \
								  + Transcript_BioType + "|" \
								  + Rank + "|" \
								  + HGVS_C + "|" \
								  + HGVS_P + "|" \
								  + str(cDNA_pos) + "/" \
								  + str(cDNA_length) + "|" \
								  + str(CDS_pos) + "/" \
								  + str(CDS_length) + "|" \
								  + str(AA_pos) + "/" \
								  + str(AA_length) + "|" \
								  + distance + "|" \
								  + errors_warnings
		return snpeff_like_info_string


@unique
class hgvs_dna(Enum):
	SUB = "Substitution"
	DEL = "Deletion"
	DUP = "Duplication"
	INS = "Insertion"
	INV = "Inversion"
	CON = "Conversion"
	DELIN = "Deletion-Insertion"
	ALL = "Alleles"
	REP = "Repeated sequences"
	COM = "Complex"

@unique
class hgvs_prot(Enum):
	SUB = "Substitution"
	DEL = "Deletion"
	DUP = "Duplication"
	INS = "Insertion"
	DELIN = "Deletion-Insertion"
	ALL = "Alleles"
	REP = "Repeated sequences"
	FS  = "Frame shift"
	EX  = "Extension"
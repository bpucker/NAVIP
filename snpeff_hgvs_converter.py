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
		aminodict['X'] = 'Xaa'
		aminodict['Xaa'] = 'X'
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
			if newamino == '*': # counts as SUB, len == 1
				"""
				p.(Tyr4*)
				the predicted consequence at the protein level of the variant ATGGATGCATACGTCACG.. to 
				ATGGATGCATA\_GTCACG (c.12delC) is a Tyr to translation termination codon. NOTE: the 
				variant is described as a substitution, not as a frame shift (p.Tyr4TerfsTer1)
				"""
				HGVS_P += oldamino + str(AA_pos) + '*'
			#elif '*' in newamino: #no idea what this is, but not always a frameshift
			#	pass
			else:
				#first changed aa has to be the first aa here
				if oldamino[0] != newamino[0]:
					HGVS_P += oldamino + str(AA_pos) + newamino + 'Ter' + str(vinfo.STOP_CAUSED_IN)
				else:
					#find first changed aa....
					# test every aa, if they are equal, from the variant_position
					i = 0
					# python can't do that:
					#for oldaa,newaa in  transcript.uChAAsequence[AA_pos:], transcript.IV_ChangedTranslation[new_aa_pos:]:
					#	if oldaa == "":
					#		print("should not happen")
					#	elif newaa == '*':
					#		HGVS_P += oldamino + str(AA_pos + i) + '*'
					#		break
					#	elif oldaa == newaa:
					#		i +=1
					#		continue
					#	else:
					#		if vinfo.STOP_CAUSED_IN == -1:
					#			HGVS_P += oldaa + str(AA_pos + i) + newaa + '*?'
					#		else:
					#			HGVS_P += oldaa + str(AA_pos + i) + newaa + '*' + str(vinfo.STOP_CAUSED_IN - i)
					#		break
					for j, oldaa in enumerate (transcript.uChAAsequence[AA_pos-1:]):
						newaa = transcript.IV_ChangedTranslation[new_aa_pos-1 + j]
						if oldaa == "":
							print("should not happen")
						elif newaa == '*':
							HGVS_P += oldamino + str(AA_pos + i) + '*'
							break
						elif oldaa == newaa:
							i +=1
							continue
						else:
							if vinfo.STOP_CAUSED_IN == -1:
								HGVS_P += oldaa + str(AA_pos + i) + newaa + '*?'
							else:
								HGVS_P += oldaa + str(AA_pos + i) + newaa + '*' + str(vinfo.STOP_CAUSED_IN - i)
							break

		elif TranscriptEnum.INSERTION in vinfo.Classification:
			# extension is already checked ->
			# frameshift, too
			# first duplicationcheck
			# second and last is the normal insertion

			if ((len (vinfo.Alt) -1) % 3) != 0:
				print('how`?') # well, it should not be a frameshift -> multiple of 3
			else:
				#dup check
				insertAA = vinfo.NewAmino
				if (len(vinfo.Alt) -1) % 3 != 0:
					print("how2?")
				refpartAA = ""

				if len(transcript.uChAAsequence) >= len(transcript.IV_OriginalTranslation):
					if (len(vinfo.Alt) -1) / 3 == len(vinfo.NewAmino):
						#Insertion of complete AA without interfering with other AA
						refpartAA = transcript.uChAAsequence[AA_pos-1 - len(vinfo.NewAmino):AA_pos]
					else:
						#Insertion with interfering with other AA
						# -> the original AA is in vinfo.NewAmino
						# -> do not compare it with itself
						if len(vinfo.NewAmino) == 1:
							print('how?3')
						refpartAA = transcript.uChAAsequence[AA_pos-1 - (len(vinfo.NewAmino) -1) :AA_pos]
				else:
					if (len(vinfo.Alt) -1) / 3 == len(vinfo.NewAmino):
						#Insertion of complete AA without interfering with other AA
						refpartAA = transcript.IV_OriginalTranslation[AA_pos-1 - len(vinfo.NewAmino):AA_pos]
					else:
						#Insertion with interfering with other AA
						# -> the original AA is in vinfo.NewAmino
						# -> do not compare it with itself
						if len(vinfo.NewAmino) == 1:
							print('how?3')
						refpartAA = transcript.IV_OriginalTranslation[AA_pos-1 - (len(vinfo.NewAmino) -1 ):AA_pos ]
				#print(vinfo.NewAmino + "\t" + refpartAA)
				#print(transcript.IV_OriginalTranslation[AA_pos - (len(vinfo.NewAmino) -3 ):AA_pos + 3])

				if len(vinfo.NewAmino) == len(refpartAA):
					dup = True
					#print(refpartAA, vinfo.NewAmino)
					#print(type(refpartAA))
					for i, oa in enumerate(refpartAA):
						if vinfo.NewAmino[i] == oa:
							continue
						else:
							dup = False
					if dup:
						if (len(vinfo.Alt) -1) / 3 == len(vinfo.NewAmino):
							# Insertion of complete AA without interfering with other AA
							if len(vinfo.NewAmino) == 1: # one AA
								"""
								p.Ala3dup (one amino acid)
								a duplication of amino acid Ala3 in the sequence 
								MetGlyAlaArgSerSerHis to MetGlyAlaAlaArgSerSerHis
								"""
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'dup'
							else: # multiple AA
								"""
								p.Ala3_Ser5dup (several amino acids)
								a duplication of amino acids Ala3 to Ser5 in the sequence 
								MetGlyAlaArgSerSerHis to MetGlyAlaArgSerAlaArgSerSerHis
								"""
								HGVS_P += refpartAA[0] + str(AA_pos - (len(vinfo.NewAmino) -1 )) + '_' \
										  + aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'dup'


						else:
							# Insertion with interfering with other AA
							if len(vinfo.NewAmino) == 2:  # one new AA, first AA not changed, because its a dup
								# so it looks like, this can only happen, if vinfo.NewAmino[0] == vinfo.NewAmino[1]
								"""
								p.Ala3dup (one amino acid)
								a duplication of amino acid Ala3 in the sequence 
								MetGlyAlaArgSerSerHis to MetGlyAlaAlaArgSerSerHis
								"""
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'dup'
							else:  # multiple AA
								"""
								p.Ala3_Ser5dup (several amino acids)
								a duplication of amino acids Ala3 to Ser5 in the sequence 
								MetGlyAlaArgSerSerHis to MetGlyAlaArgSerAlaArgSerSerHis
								"""
								#print(refpartAA[0] + str(AA_pos - (len(vinfo.NewAmino) - 1)) + '_' + vinfo.OrigAmino + str(AA_pos) + 'dup' )
								#print(refpartAA + str(AA_pos - (len(vinfo.NewAmino) - 1)) + '_' + vinfo.OrigAmino + str(AA_pos) + 'dup')
								#print(vinfo.ChrPosition)

								HGVS_P += refpartAA[0] + str(AA_pos - (len(vinfo.NewAmino) - 1)) + '_' + \
										  aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'dup'
					else:
						# no dup, but insertion without frameshift
						# delins are possible
						if (len(vinfo.Alt) -1) / 3 == len(vinfo.NewAmino):
							# Insertion of complete AA without interfering with other AA
							"""
								p.His4_Gln5insAla
								the insertion of amino acid Ala between amino acids His4 and Gln5 
								changing MetLysGlyHisGlnGlnCys to MetLysGlyHisAlaGlnGlnCys

								p.Lys2_Gly3insGlnSerLys
								the insertion of amino acids GlnSerLys between amino acids Lys2 and 
								Gly3 changing MetLysGlyHisGlnGlnCys to MetLysGlnSerLysGlyHisGlnGlnCys
							"""
							#if len(vinfo.NewAmino) == 1: # one AA
							HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + '_' + refpartAA[1] \
									  + str(AA_pos+1) + 'ins' + newamino
							#else: # multiple AA
							#	HGVS_P += vinfo.OrigAmino[0] + str(AA_pos) + '_' + refpartAA[AA_pos + 1] \
							#			  + str(AA_pos + 1) + 'ins' + vinfo.NewAmino
						else:
							# looks like, if the first AA is changed -> delin, the ref AA, when the variation is interfering
							# with another AA
							"""
								p.Cys28delinsTrpVal
								a deletion of amino acid Cys28, replaced with TrpVal
							"""
							if vinfo.OrigAmino[0] == vinfo.NewAmino[0]: # no delin
								#print(vinfo.ChrPosition)
								#print(vinfo.OrigAmino[0] + str(AA_pos) + '_')
								#print(refpartAA)
								#print(transcript.IV_OriginalTranslation[AA_pos:])
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + '_' + refpartAA[ 1] \
										  + str(AA_pos + 1) + 'ins' + newamino
							else:
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'delins' + newamino
				else:
					# reasons for a variation to be in this else:
					# insertion in the beginning of the transctript -> insert to long, no complete refAA existent
					#	-> example: insert 5 AA, but its in the start codon
					if len(vinfo.NewAmino) >= AA_pos:
						# no duplication/ refAA possible, because insertion is to big
						"""
							p.Cys28delinsTrpVal
							a deletion of amino acid Cys28, replaced with TrpVal
						"""
						if vinfo.OrigRaster == 2:
							#insertion without interfering with first origAA
							"""
							p.His4_Gln5insAla
								the insertion of amino acid Ala between amino acids His4 
								and Gln5 changing MetLysGlyHisGlnGlnCys to MetLysGlyHisAlaGlnGlnCys
							"""
							if len(vinfo.NewAmino) > 1:
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'ins' + newamino[1:]
							else:
								print("curios effect")
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'ins' + newamino
						else:
							#insertion with interfering with first origAA
							if vinfo.OrigAmino[0] == vinfo.NewAmino[0]:
								#first AA not changed -> insertion after
								HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'ins' + newamino
							else:
								#first AA changed -> SUB or DELINS
								if len(vinfo.NewAmino) == 1 and vinfo.NewAmino == '*':
									#SUB
									#LRG_199p1:p.Trp24Ter (p.Trp24*)
									HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + newamino
								else:
									#DELINS
									"""
										p.Cys28delinsTrpVal
										a deletion of amino acid Cys28, replaced with TrpVal
									"""
									HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'delins' + newamino

					else:
						# stopcodon + stuff at the end of the transcript
						if vinfo.NewAmino == vinfo.OrigAmino:
						#if vinfo.NewAmino[0] == vinfo.OrigAmino[0] and len (vinfo.NewAmino) == 1 and len(vinfo.OrigAmino) == 1:
							# silent sub
							"""
								silent (no change)
								NP_003997.1:p.Cys188=
							"""
							b = ""
							for aa in vinfo.OrigAmino:
								b += aminodict[aa]
							HGVS_P += b + str(AA_pos)  + '='
						else:
							"""
								p.Cys28delinsTrpVal
								a deletion of amino acid Cys28, replaced with TrpVal
								
								p.(Ter315TyrextAsnLysGlyThrTer) (alternatively p.*315TyrextAsnLysGlyThr*)
									a variant in the stop codon (Ter/*) at position 315, 
									changing it to a Tyr-codon (a no-stop variant) and adding a tail of 
									new amino acids to the protein’s C-terminus, ending at a new stop codon (Ter5/*5)
							"""
							if vinfo.NewAmino[0] == vinfo.OrigAmino[0]:
								# no sub, possible here

								if len(vinfo.NewAmino) > 1 and len(vinfo.OrigAmino) > 1:
									#delins
									b = ""
									for aa in vinfo.OrigAmino:
										b += aminodict[aa]
									HGVS_P += b + str(AA_pos +1) + 'delins' + newamino[1:]
								elif len(vinfo.OrigAmino) == 1 and len(vinfo.NewAmino) > 1:
									# here, if the AA 2+ are not identical with the original AA
									HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + 'delins' + newamino
								else: # orig > 1 and newamino == 1 () impossible?
									print('Impossible case triggered.')
							else:
								b = ""
								for aa in vinfo.OrigAmino:
									b += aminodict[aa]
								HGVS_P += b + str(AA_pos) + 'delins' + newamino
		elif TranscriptEnum.DELETION in vinfo.Classification:
			# deletionstuff, no frameshifts possible anymore -> codon-position should always be 2 ???? nope
			# len(vinfo.NewAmino) > 1: true
			if vinfo.OrigRaster == 2:
				#deletion of complete AA, without interfering of another AA (startposition +1 (dna))
				if len(vinfo.OrigAmino) == 2:
					#complete deletion of the second AA
					"""
						one AA
						LRG_199p1:p.Val7del
						a deletion of amino acid Val7 in the reference sequence LRG_199p1
					"""
					HGVS_P += aminodict[vinfo.OrigAmino[1].upper()] + str(AA_pos +1 ) + 'del'
				elif len(vinfo.OrigAmino) == 1:
					#impossible case (or somethi8ng went wrong with the transcript)
					if vinfo.OrigAmino == vinfo.NewAmino:
						#sub
						#NP_003997.1:p.Cys188=
						HGVS_P += oldamino + str(AA_pos) + '='
					else:
						try:
							HGVS_P += oldamino + str(AA_pos) + vinfo.NewAmino[0]
						except IndexError:
							print("Variant error: " + str(vinfo.ChrPosition))
							HGVS_P += 'oldamino' + str(AA_pos) + 'Xaa'
				else:
					"""
						a few AA
						NP_003997.1:p.Lys23_Val25del
						a deletion of amino acids Lys23 to Val25 in reference sequence NP_003997.1
					"""
					HGVS_P += oldamino[1] + str(AA_pos + 1) + '_' + oldamino[len(oldamino)-1] + str(AA_pos + len(oldamino)) +  'del'

			else :
				if vinfo.OrigAmino == "":
					print("Variant error: " + str(vinfo.ChrPosition))
				elif len(vinfo.OrigAmino) == 1:
					#this should not happen here, but who knows
					if vinfo.OrigAmino == vinfo.NewAmino:
						#sub
						#NP_003997.1:p.Cys188=
						HGVS_P += oldamino + str(AA_pos) + '='
					else:
						HGVS_P += oldamino + str(AA_pos) + aminodict[vinfo.NewAmino.upper()]
				elif TranscriptEnum.STOP_GAINED in vinfo.Classification:
					"""
						p.Trp26Ter (p.Trp26*)
						amino acid Trp26 is changed to a stop codon (Ter, *)
						NOTE: this change is not described as a deletion of 
						the C-terminal end of the protein (i.e. p.Trp26_Arg1623del)
					"""
					HGVS_P += aminodict[vinfo.OrigAmino[0]] + str(AA_pos) + aminodict[vinfo.NewAmino.upper()]
				else:
					# always deletions, starting inside the codon, always multiple of 3 bp deleted
					# stop lost, stop changed and stop gained is handled to this point
					# --> only delins (?)
					"""
						p.Cys28_Lys29delinsTrp
						a deletion of amino acids Cys28 and Ly29, replaced with Trp
					"""
					HGVS_P += aminodict[vinfo.OrigAmino[0].upper()] + str(AA_pos) + "_" + aminodict[vinfo.OrigAmino[len(vinfo.OrigAmino)-1].upper()] + str(AA_pos + len(vinfo.OrigAmino)) + 'delins'
		else:
			print('what did i forget?')
		snpeff_like_info_string = 'NAV2=' + vinfo.Alt + '|' \
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
								  + errors_warnings + ";"
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
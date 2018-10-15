__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import Transcript
import sys

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
			try:
				self.gff3_ID_Dict[seqid,gff3_id].append(spline)
			except KeyError:
				self.gff3_ID_Dict[seqid, gff3_id] = [spline]
			if len(parent) == 0:
				#when its a gene entry
				continue
			for pp in parent:
				try:
					self.dict_gff_for_parents[seqid, pp].append(spline)
				except KeyError:
					self.dict_gff_for_parents[seqid, pp] = [spline]

	def createTranscripts(self):

		for key, splinelist in self.gff3_ID_Dict.items():

			for spline in splinelist:
				if spline[2] != "gene":
					continue

				gene_seqid = spline[0]
				gene_gff3_ID = key[1] # because of (seqid,gff3_id)
				gene_info_string = spline[8]

				if gene_seqid not in self.dict_Chr_dict_Transcript:
					self.dict_Chr_dict_Transcript[gene_seqid] = {}


				for spline_child in self.dict_gff_for_parents[gene_seqid,gene_gff3_ID]:
					chr = spline_child[0]
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
						sys.exit()

					# IndexKey: int, TID: str, StartOfRNA: int, EndOfRNA: int, ForwardDirection: TranscriptEnum.REVERSE):
					transcript = Transcript.Transcript(self.GetNextTranscriptIndex(),
											TID,
											start,
											end,
											strand,
											chr)
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




















__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


import Transcript
import sys

class GFF3_Handler:


	def __init__(self, GFF3_data_path: str):
		self.gff3_ID_dict = {}
		self.dict_gff_for_parents = {}
		self.dict_chr_names = {}
		self.dict_chr_dict_transcript = {}
		self.list_of_transcripts = {}
		self.nr_chroms = -1
		self.trans_index = -1

		self.read_GFF3(GFF3_data_path)
		self.create_transcripts()
		self.sort_transcripts()

	def read_GFF3(self, GFF3_data_path: str):
		gff3 = open(GFF3_data_path, 'r')
		#lines = gff3.readlines()
		count_generic_id = 0
		for line in gff3:
			if line.startswith('#') or len(line) <=2:
				continue
			spline = line.split('\t')
			seqid = spline[0]
			if seqid not in self.dict_chr_names.keys():
				self.nr_chroms += 1
				self.dict_chr_names[seqid] = self.nr_chroms
				self.dict_chr_names[self.nr_chroms] = seqid
				self.list_of_transcripts[seqid] = []

			atts = spline[8].split(';')
			gff3_id = ""
			parent = ""
			if "ID=" not in spline[8]:
				#create generic_id: chr_type_start_end
				generic_id = str(spline[0]) + '_' + str(spline[1]) + '_' + str(spline[3]) + '_' + str(spline[4])


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
				count_generic_id +=1
				gff3_id = generic_id
				#print(line)
				#print("Item(s) without ID - is it really gff version 3?(or empty lines....)")
				#sys.exit()
			try:
				self.gff3_ID_dict[seqid,gff3_id].append(spline)
			except KeyError:
				self.gff3_ID_dict[seqid, gff3_id] = [spline]
			if len(parent) == 0:
				#when it's a gene entry
				continue
			for pp in parent:
				try:
					self.dict_gff_for_parents[seqid, pp].append(spline)
				except KeyError:
					self.dict_gff_for_parents[seqid, pp] = [spline]
		gff3.close()

	def create_transcripts(self):

		for key, splinelist in self.gff3_ID_dict.items():

			for spline in splinelist:
				if spline[2] != "gene":
					continue

				gene_seqid = spline[0]
				gene_gff3_ID = key[1] # because of (seqid,gff3_id)
				gene_info_string = spline[8]

				if gene_seqid not in self.dict_chr_dict_transcript:
					self.dict_chr_dict_transcript[gene_seqid] = {}


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
					transcript = Transcript.Transcript(self.get_next_transcript_index(),
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
					self.dict_chr_dict_transcript[gene_seqid][transcript.IndexKey] = transcript
					self.list_of_transcripts[gene_seqid].append(transcript)

	def get_next_transcript_index (self):
		"""
		For creating new transcripts (need for a new ID).
		:return: Next integer ID.
		"""
		self.trans_index +=1
		return self.trans_index

	def get_chromosome_names (self)->list:
		"""
		Returns all names of the chromosomes inside a list.
		:return: List of all chromosome names.
		"""
		count_names = len(self.dict_chr_names)/2
		i = 0
		name_list = []
		while count_names != i:
			name_list.append(self.dict_chr_names[i])
			i +=1
		if not name_list[0]:
			name_list.pop(0)

		return name_list

	def get_chromosome_ID (self, chr_name:str)->int:
		return self.dict_chr_names[chr_name]

	def get_chr_transcripts_dict(self, chr_name: str):
		"""
		Returns a dictionary with all transcripts from the chosen chromosome.
		:param chr_name: Name of the chromosome.
		:return: Dictionary with all transcripts inside this chromosome.
		"""
		return self.dict_chr_dict_transcript[chr_name]

	def free_RAM (self, chr_name: str):
		"""
		Frees the RAM.
		:param chr_name:
		:return:
		"""
		self.dict_chr_dict_transcript[chr_name] = []
		self.list_of_transcripts[chr_name] = []

	def get_chr_transcript_list(self, chr_name:str):
		return self.list_of_transcripts[chr_name]

	def sort_transcripts(self):
		for chr_name in self.list_of_transcripts.keys():
			self.list_of_transcripts[chr_name] = sorted(self.list_of_transcripts[chr_name], key=lambda s_transcript: s_transcript.StartOfRNA)

	def update_transcripts(self, transcript_list: list, chr_name: str):
		#i = self.dict_chr_names[chr_name]
		#self.list_of_transcripts[i] = transcript_list
		self.list_of_transcripts[chr_name] = sorted(transcript_list, key=lambda s_transcript: s_transcript.StartOfRNA)

	def add_new_transcript_to_dict(self, chr_name: str, transcript: Transcript):
		self.dict_chr_dict_transcript[chr_name][transcript.TID] = transcript

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

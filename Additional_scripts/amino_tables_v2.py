from Additional_scripts.variant_parser import nav1_container
from Transcript import TranscriptEnum

class amino_table:
	"""
	This class will create an amino comparison table.
	(20+2) * (20+2)
	20 Amino acids + stop + X for everything else.
	"""
	def __init__(self):
		self.aminolist = ["A", "C", "D" , "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", '*']
		self.aminotable = {}
		#aminotable[original_aminoacid][new_aminoacid]
		for aa in self.aminolist:
			self.aminotable[aa] = {}
			for aa2 in self.aminolist:
				self.aminotable[aa][aa2] = 0

	def add_nav1_container_to_dict(self, nav1_cont: nav1_container):

		if TranscriptEnum.SUBSTITUTION.value not in nav1_cont.C2_Variant_Classifications\
				or TranscriptEnum.INSERTION in nav1_cont.C2_Variant_Classifications\
				or TranscriptEnum.DELETION in nav1_cont.C2_Variant_Classifications:
			raise ValueError("Only Substitutions allowed")

		self.aminotable[nav1_cont.C5_REF_AA.upper()][nav1_cont.C8_ALT_AA.upper()] += 1

	def create_aminotable(self, aminotable_data_file):
		outputtext = ".\t"
		outputtext += "\t".join(self.aminolist) # header line
		for orig_aa in self.aminolist:
			outputtext += "\n" # everything starts in a new line
			outputtext +=  orig_aa
			for new_aa in self.aminolist:
				outputtext += "\t" + str(self.aminotable[orig_aa][new_aa])

		write_data_file = open(aminotable_data_file, 'w')
		write_data_file.write(outputtext)
		write_data_file.close()
from enum import Enum, unique

class Variant:
	"""
	A simple storage class for a vcf variant.
	Contains all its information.
	"""
	def __init__ (self,
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
		Initialization of all information.
		Contains all information from the vcf file and one new ID for variant handling.
		:param chromosome: Chromosome.
		:param position: Position.
		:param ID: VCF ID
		:param usefullID: A real ID in numbers, which is unique.
		:param reference: VCF REF
		:param alternate: VCF ALT
		:param qual: VCF QUAL
		:param filter_: VCF Filter
		:param info: VCF Info
		"""


		#Initialization of variables
		self.Chromosome = chromosome
		self.Position = position
		self.ID = ID
		self.UsefullID = usefullID
		self.Reference = reference
		self.Alternate = alternate
		self.Qual = qual
		self.Filter = filter_
		self.Info = info
		self.ListofTranscripts =[]
		self.SListofTranscripts = []


@unique
class VariantEnum (Enum):
	"""
	Contains a few enums for variant usage only.
	"""
	DEFAULT = -1 # nothing decided (for finding errors - default should not exist in the output)
	NO_RASTER_NOW = -2
	NO_StartOfOwnEffect = -3
	NO_EndOfOwnEffect = -4



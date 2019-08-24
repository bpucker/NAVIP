from enum import Enum


class nav_infoline_enum(Enum):
	A0_NAVID = "NAVID="
	A1_NAVID_START = 6
	A2_NAV1 = "NAV1="
	A3_NAV1_START = 5
	A4_Nav2 = "NAV2="
	A5_NAV2_START = 5


class easy_parser:
	"""
	This class will contain ONE variant. It will create it from one incomming line-string.
	"""

	def __init__(self, vcf_file_line:str):
		if vcf_file_line.startswith('#'):
			raise ValueError("No comment lines allowed: " + vcf_file_line)
		elif  len(vcf_file_line) <= 4:
			raise ValueError("Line to short: " + vcf_file_line)


		##CHROM POS     ID        REF ALT    QUAL FILTER INFO
		##0		1		2			3 4		5		6		7
		spvcf_file_line = vcf_file_line.split("\t")
		# A is here because so the variables are ordered
		self.A0_chrom = spvcf_file_line[0]
		self.A1_position = int(spvcf_file_line[1])
		self.A2_vcf_id = spvcf_file_line[2]
		self.A3_ref = spvcf_file_line[3]
		self.A4_alt = spvcf_file_line[4]
		self.A5_qual = spvcf_file_line[5]
		self.A6_filter = spvcf_file_line[6]
		self.A7_info = spvcf_file_line[7]
		if len(spvcf_file_line) > 8:
			if len(spvcf_file_line) == 9:
				self.A8_after_info = spvcf_file_line[9]
			else:
				self.A8_after_info = "\t".join(spvcf_file_line[9:])
		else:
			self.A8_after_info = ""
		self.B0_navid = ""
		self.B1_nav1 = ""
		self.B2_nav2 = ""
		self.B3_additional_info = ""
		for info in self.A7_info.split(";"):
			if info.startswith(nav_infoline_enum.A0_NAVID.value):
				self.B0_navid = int(info[nav_infoline_enum.A1_NAVID_START.value:])
			elif info.startswith(nav_infoline_enum.A2_NAV1.value):
				self.B1_nav1 = info[nav_infoline_enum.A3_NAV1_START.value:]
			elif info.startswith(nav_infoline_enum.A4_Nav2.value):
				self.B2_nav2 = info[nav_infoline_enum.A5_NAV2_START.value:]
			else:
				self.B3_additional_info += info
		self.C_nav1_container_list = []
		for nav1string in self.B1_nav1.split("@"):
			self.C_nav1_container_list.append(nav1_container(nav1string))
		self.D_nav2_container_list = []
		for nav2string in self.B2_nav2.split("@"):
			self.D_nav2_container_list.append(nav2_container(nav2string))



class nav1_container:

	def __init__(self, one_nav1_string:str):

		##                                              0             1                  2                                               3                    4.0           4.1                  5       6                     7.0     7.1                        8           9
		##Info=<ID=NAV1, Type=String,Number=.,Values=[TranscriptID|Strand_Direction|Variant_Classification1,Variant_Classification2,...|Shared_EffKey(s)|REF_Codon(s)/Variant_Position_in_Codon|REF_AA|old_CDS_Position|ALT_Codon(s)/Variant_Position_in_Codon|ALT_AA|new_CDS_Position]>

		sp_one_nav1_string = one_nav1_string.split("|")
		self.C0_TranscriptID = sp_one_nav1_string[0]
		self.C1_Strand_Direction = sp_one_nav1_string[1]
		self.C2_Variant_Classifications = sp_one_nav1_string[2].split(",")
		self.C3_Shared_EffKey = sp_one_nav1_string[3].split(",")
		self.C4_0_REF_Codons = sp_one_nav1_string[4].split("/")[0]
		self.C4_1_Variant_Position_in_Codon = sp_one_nav1_string[4].split("/")[1]
		self.C5_REF_AA = sp_one_nav1_string[5]
		self.C6_old_CDS_Position = sp_one_nav1_string[6]
		self.C7_0_ALT_Codons = sp_one_nav1_string[7].split("/")[0]
		self.C7_1_Variant_Position_in_Codon = sp_one_nav1_string[7].split("/")[1]
		self.C8_ALT_AA = sp_one_nav1_string[8]
		self.C9_new_CDS_Position = sp_one_nav1_string[9]

class nav2_container:
	def __init__(self,one_nav2_string:str):
		##                                             0    1                2               3        4    5            6             7               8       9    10     11.0     11.1        12.0    12.1            13.0  13.1     14        15
		##Info=<ID=NAV2, Type=String,Number=.,Values=[ALT|ANNOTATIONS|ANNOTATION_IMPACT|GENE_NAME|GENE_ID|FEATURE_TYPE|FEATURE_ID|TRANSCRIPT_BIOTYPE|RANK|HGVS_C|HGVS_P|cDNA_pos/cDNA_length|CDS_pos/CDS_DNA_length|AA_pos/AA_length|distance|errors_warnings]>

		spone_nav2_string = one_nav2_string.split("|")

		self.D00_ALT = spone_nav2_string[0]
		self.D01_ANNOTATIONS = spone_nav2_string[1]
		self.D02_ANNOTATION_IMPACT = spone_nav2_string[2]
		self.D03_GENE_NAME = spone_nav2_string[3]
		self.D04_GENE_ID = spone_nav2_string[4]
		self.D05_FEATURE_TYPE = spone_nav2_string[5]
		self.D06_FEATURE_ID = spone_nav2_string[6]
		self.D07_TRANSCRIPT_BIOTYPE = spone_nav2_string[7]
		self.D08_RANK = spone_nav2_string[8]
		self.D09_HGVS_C = spone_nav2_string[9]
		self.D10_HGVS_P = spone_nav2_string[10]
		self.D11_0_cDNA_pos = spone_nav2_string[11].split("/")[0]
		self.D11_1_cDNA_length = spone_nav2_string[11].split("/")[1]
		self.D12_0_CDS_pos  = spone_nav2_string[12].split("/")[0]
		self.D12_1_CDS_DNA_length  = spone_nav2_string[12].split("/")[1]
		self.D13_0_AA_pos = spone_nav2_string[13].split("/")[0]
		self.D13_1_AA_length = spone_nav2_string[13].split("/")[1]
		self.D14_distance =  spone_nav2_string[14]
		self.D15_errors_warnings =  spone_nav2_string[15]
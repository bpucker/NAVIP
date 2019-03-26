import compensating_indels as cindels
import amino_tables as atables
from datetime import datetime


def sfa2_main(navip_vcf_file_link:str, mod_or_not:bool, outputfolder:str ,formats:str):
	#def find_all_cindels(vcf_file_link:str, mod_or_not: bool, outputfolder: str, max_bp_range: int):
	#cindels.find_all_cindels(navip_vcf_file_link, mod_or_not, outputfolder, max_bp_range )

	#def find_all_cindels_v2(navip_vcf_file_link: str, mod_or_not: bool, outputfolder: str, formats:str):
	time = datetime.now()
	cindels.find_all_cindels_v2(navip_vcf_file_link, mod_or_not, outputfolder, formats)
	print("Finished in: " +str(datetime.now() - time))
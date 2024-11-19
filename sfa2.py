import compensating_indels as cindels
from datetime import datetime

def sfa2_main(navip_vcf_file:str, mod_or_not:bool, output_folder:str, formats:str, max_x_axis_bpr:int):
	#def find_all_cindels_v1(vcf_file_link:str, mod_or_not: bool, output_folder: str, max_bp_range: int):
	#cindels.find_all_cindels_v1(navip_vcf_file, mod_or_not, output_folder, max_bp_range)
	#def find_all_cindels_v2(navip_vcf_file_link: str, mod_or_not: bool, output_folder: str, formats:str):

	print("Starting compensating_indels script")
	starttime = datetime.now()
	try:
		cindels.find_all_cindels_v2(navip_vcf_file, mod_or_not, output_folder, formats, max_x_axis_bpr)
	except Exception:
		import sys; print("Warning: List of cInDels is empty, skipping the analysis...", file=sys.stderr)
	print("Finished in: " + str(datetime.now() - starttime))

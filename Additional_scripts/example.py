import Additional_scripts.only_no_none as onn
import Additional_scripts.ID_Handler as idh
import Additional_scripts.join_variants_together as jvt
import Additional_scripts.variant_parser as vp
import Additional_scripts.amino_tables_v2 as at
from Transcript import TranscriptEnum

def do_stuff(all_vcf_file, vcf_with_id, vcf_with_id_merged, table_out ):
	## existing file
	#all_vcf_file = "/homes/janbaas/NAVIP_prj/data/arabidopsis_navip/arabidopsis_ram/first_run_ara11/All_VCF.vcf"

	## not existing files
	#vcf_with_id = "/homes/janbaas/NAVIP_prj/data/arabidopsis_navip/arabidopsis_ram/first_run_ara11/Test/All_VCF_ID.vcf"
	#vcf_with_id_merged = "/homes/janbaas/NAVIP_prj/data/arabidopsis_navip/arabidopsis_ram/first_run_ara11/Test/All_VCF_id_merged.vcf"
	#table_out = "/homes/janbaas/NAVIP_prj/data/arabidopsis_navip/arabidopsis_ram/first_run_ara11/Test/testtable.txt"


	#create vcf with id tag
	#                      inputfile,    outputfile
	idh.add_stable_id_tag(all_vcf_file,vcf_with_id )

	#create vcf with a single line per position
	#                          inputfile,        outputfile
	jvt.join_variants_together(vcf_with_id, vcf_with_id_merged)


	# code for the outputtable
	# you still have to decide, which variant will go in the table

	vcf = open(vcf_with_id_merged, 'r')
	lines = vcf.readlines()
	vcf.close()

	myaminotable = at.amino_table()

	for line in lines:
		if line.startswith('#'):
			continue
		elif len(line) <= 4:
			continue
		# in case of different aminoacids because of different transcripts the first will be chosen
		nav1_cont = vp.easy_parser(line).C_nav1_container_list[0]


		# only single amino acids changes are possible for the table
		# everything, thats not a substitution will be dismissed
		# feel free to use more definitions
		if TranscriptEnum.SUBSTITUTION.value not in nav1_cont.C2_Variant_Classifications \
				or TranscriptEnum.INSERTION in nav1_cont.C2_Variant_Classifications \
				or TranscriptEnum.DELETION in nav1_cont.C2_Variant_Classifications:
			continue

		myaminotable.add_nav1_container_to_dict(nav1_cont)

	myaminotable.create_aminotable(table_out)

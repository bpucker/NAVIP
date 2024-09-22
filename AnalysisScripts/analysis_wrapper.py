### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.3 ###

__usage__ = """
					python3 analysis_wrapper.py
					--in <VCF_FOLDER>
					--out <OUTPUT_DIRECTORY>
					--fasta <FASTA_FILE>
					--gff <GFF_FILE>
					--navip_clean_vcf <SCRIPT_navip_clean_vcf.py>
					--snpeff <SnpEff_PATH>
					--navip <NAVIP_SCRIPT_PATH>
					--compare_stop_gain_events <SCRIPT_compare_stop_gain_events.py>
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""

import os, sys, glob, subprocess

# --- end of imports --- #

def main( arguments ):
	"""! @brief runs all parts of this script """
	
	input_folder = arguments[ arguments.index( '--in' )+1 ]
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	gff_file = arguments[ arguments.index( '--gff' )+1 ]
	fasta_file = arguments[ arguments.index( '--fasta' )+1 ]
	navip_clean_vcf = arguments[ arguments.index( '--navip_clean_vcf' )+1 ]
	snpeff = arguments[ arguments.index( '--snpeff' )+1 ]
	navip = arguments[ arguments.index( '--navip' )+1 ]
	compare_stop_gain_events = arguments[ arguments.index( '--compare_stop_gain_events' )+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	p = subprocess.Popen( args= "cd " + output_folder, shell=True )
	p.communicate()
	
	vcf_files = glob.glob( input_folder + "*.vcf" )
	for vcf in vcf_files:
		# --- clean VCF --- #
		clean_vcf = output_folder + vcf.split('/')[-1].split('.vcf')[0] + ".clean.vcf"
		if not os.path.isfile( clean_vcf ):
			cmd = "python3 " + navip_clean_vcf + " --in " + vcf + " --out " + clean_vcf
			p = subprocess.Popen( args= cmd, shell=True )
			p.communicate()
		
		# --- run SnpEff --- #
		snpeff_vcf = output_folder + vcf.split('/')[-1].split('.vcf')[0] + ".snpeff.vcf"
		doc_file = output_folder + vcf.split('/')[-1].split('.vcf')[0] + ".snpeff.doc.txt"
		if not os.path.isfile( snpeff_vcf ):
			cmd = "java -Xmx8g -jar " + snpeff + " araport11 " + clean_vcf + " > " + snpeff_vcf + " 2> " + doc_file
			p = subprocess.Popen( args= cmd, shell=True )
			p.communicate()
		
		# --- run NAVIP --- #
		navip_output =  output_folder + vcf.split('/')[-1].split('.vcf')[0] + ".navip/"
		navip_vcf = navip_output + "NAVIP_Main_Output/All_VCF.vcf"
		if not os.path.isfile( navip_vcf ):
			cmd = navip + " -i " + clean_vcf + " -g " + gff_file + " -f " + fasta_file + " -o " + navip_output
			p = subprocess.Popen( args= cmd, shell=True )
			p.communicate()
		
		# --- compare SnpEff and NAVIP --- #
		comparison_output_folder = output_folder + vcf.split('/')[-1].split('.vcf')[0] + "/"
		if not os.path.exists( comparison_output_folder ):
			os.makedirs( comparison_output_folder )
			cmd = "python3 " + compare_stop_gain_events + " --snpeffvcf " + snpeff_vcf + " --navipvcf " + navip_vcf + " --out " + comparison_output_folder
			p = subprocess.Popen( args= cmd, shell=True )
			p.communicate()


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

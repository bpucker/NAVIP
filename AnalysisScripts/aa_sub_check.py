### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.2 ###

__usage__ = """
					python3 aa_sub_check.py
					--snpeffvcf <SnpEff_VCF_OUTPUT_FILE(SELECTION)>
					--navipvcf <NAVIP_VCF_OUTPUT_FILE>
					--out <OUTPUT_FOLDER>
					"""

import os, sys, re

# --- end of imports --- #

def load_aa_substitutions_snpeff( input_vcf_file ):
	"""! @brief load SnpEff stop_gained cases with two SNVs per codon """
	
	aa_substitution_cases = {}
	with open( input_vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				pos =  parts[0] + "_%_" + parts[1]
				if "ANN=" in line:
					blocks = parts[7].split('ANN=')[1]
					if "," in blocks:
						blocks = blocks.split(',')
					else:
						blocks = [ blocks ]
					for block in blocks:
						if 'missense_variant' in block:
							subblocks = block.split('|')	#get block with amino acid substitution details
							substitution = subblocks[10][2:]	#remove "p." prefix
							#print( substitution )
							number = re.findall("\d+", substitution)[0]
							original_aa, new_aa = substitution.split( number )
							#print( original_aa, new_aa )
							aa_substitution_cases.update( { pos: [ original_aa, new_aa ] } )
			line = f.readline()
	return aa_substitution_cases


def load_aa_substitutions_navip( input_vcf_file ):
	"""! @brief load SnpEff stop_gained cases with two SNVs per codon """

	aa_substitution_cases = {}
	with open( input_vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				pos =  parts[0] + "_%_" + parts[1]
				if "NAV2=" in line:
					block = parts[7].split('NAV2=')[1]
					if 'SUB' in block:
						subblocks = block.split('|')	#get block with amino acid substitution details
						substitution = subblocks[10][2:]	#remove "p." prefix
						if "ext" in substitution:
							substitution = substitution.split('ext')[0]
						#print( substitution )
						number = re.findall("\d+", substitution)[0]
						original_aa, new_aa = substitution.split( number )
						#print( original_aa, new_aa )
						aa_substitution_cases.update( { pos: [ original_aa, new_aa ] } )
						aa_substitution_cases.update( { pos: [ original_aa, new_aa ] } )
			line = f.readline()
	return aa_substitution_cases


def identify_matches( snpeff_aas, navip_aas, summary_file ):
	"""! @brief identify matches between SnpEff and NAVIP predictions for potential stop_gained cases with two SNVs per codon """
	
	three2one_mapping = {
											"Ala": "A",  # Alanine
											"Arg": "R",  # Arginine
											"Asn": "N",  # Asparagine
											"Asp": "D",  # Aspartic acid
											"Cys": "C",  # Cysteine
											"Glu": "E",  # Glutamic acid
											"Gln": "Q",  # Glutamine
											"Gly": "G",  # Glycine
											"His": "H",  # Histidine
											"Ile": "I",  # Isoleucine
											"Leu": "L",  # Leucine
											"Lys": "K",  # Lysine
											"Met": "M",  # Methionine
											"Phe": "F",  # Phenylalanine
											"Pro": "P",  # Proline
											"Ser": "S",  # Serine
											"Thr": "T",  # Threonine
											"Trp": "W",  # Tryptophan
											"Tyr": "Y",  # Tyrosine
											"Val": "V"  # Valine
											#"*": "*"	#stop
											}
	matches, differences = {}, {}
	for x in "ARNDCEQGHILKMFPSTWYV":
		for y in "ARNDCEQGHILKMFPSTWYV":
			matches.update( { x+y: 0 } )
			differences.update( { x+y: 0 } )
	for key in list( snpeff_aas.keys() ):
		try:
			navip_effect = navip_aas[ key ]	#get NAVIP effect prediction
			snpeff_effect = snpeff_aas[ key ]	#get SnpEff effect prediction
			if navip_effect[1] == snpeff_effect[1]:	#same effect prediction
				matches[ three2one_mapping[ snpeff_effect[0] ]+three2one_mapping[ snpeff_effect[1] ]  ] += 1
			else:	#different effect predictions
				differences[ three2one_mapping[ navip_effect[0] ]+three2one_mapping[ navip_effect[1] ]  ] += 1
		except KeyError:
			print( key )	#no amino acid change and premature stop codons are not considered in this analysis
	
	with open( summary_file, "w" ) as out:
		out.write( '#AminoAcidSubstitution\tMatches\tDifferences\n' )
		for x in "ARNDCEQGHILKMFPSTWYV":
			for y in "ARNDCEQGHILKMFPSTWYV":
				if x == y:
					out.write( "\t".join( [ x + y, ".", str( differences[ x+y ] ) ] ) + "\n" )
				else:
					out.write( "\t".join( [ x + y, str( matches[ x+y ] ), str( differences[ x+y ] ) ] ) + "\n" )
	
	return matches, differences



def main( arguments ):
	"""! @brief run everything """
	
	snpeff_vcf_file = arguments[ arguments.index('--snpeffvcf')+1 ]
	navip_vcf_file = arguments[ arguments.index('--navipvcf')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- load amino acid substitutions --- #
	snpeff_aas = load_aa_substitutions_snpeff( snpeff_vcf_file )
	navip_aas = load_aa_substitutions_navip( navip_vcf_file )
	
	# --- identify matches and differences --- #
	summary_file = output_folder + "summary.txt"
	matches, differences = identify_matches( snpeff_aas, navip_aas, summary_file )


if '--snpeffvcf' in sys.argv and '--navipvcf' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

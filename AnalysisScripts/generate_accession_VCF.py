### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python3 generate_accession_VCF.py
					--in <INPUT_FILE>
					--out <OUTPUT_FOLDER>
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""

import os, sys, gzip

# --- end of imports --- #

def main( arguments ):
	"""! @brief running everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	with gzip.open( input_file, "r" ) as f:
		line1 = f.readline().decode(encoding='utf-8')	#header1
		line2 = f.readline().decode(encoding='utf-8')	#header2
		line3 = f.readline().decode(encoding='utf-8')	#header3
		line4 = f.readline().decode(encoding='utf-8')	#header4
		line5 = f.readline().decode(encoding='utf-8')	#header5
		line6 = f.readline().decode(encoding='utf-8')	#header6
		line7 = f.readline().decode(encoding='utf-8')	#header7
		line8 = f.readline().decode(encoding='utf-8')	#header8
		line9 = f.readline().decode(encoding='utf-8')	#header9
		accessions = line9.strip().split('\t')[9:]
	
	for accession in accessions:
		print( accession )
		accession_index = line9.index( accession )
		output_file = output_folder + accession + ".vcf"
		with open( output_file, "w" ) as out:
			with gzip.open( input_file, "r" ) as f:
				f.readline()	#header1
				f.readline()	#header2
				f.readline()	#header3
				f.readline()	#header4
				f.readline()	#header5
				f.readline()	#header6
				f.readline()	#header7
				f.readline()	#header8
				f.readline()	#header9
				out.write( line1 + "\t".join( line9.split('\t')[:9] + [ accession ] ) )
				line = f.readline().decode(encoding='utf-8')
				while line:
					parts = line.strip().split('\t')
					if parts[ accession_index ] != "./.:.:.":
						if parts[4] != ".":
							if parts[ accession_index ][:3] in [ "0|1", "1|1" ]:
								out.write( "\t".join( parts[:9] + [ parts[ accession_index ] ] ) + "\n" )
					line = f.readline().decode(encoding='utf-8')


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.21 ###

__usage__ = """
					python3 navip_clean_vcf.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""

import os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief running everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			out.write( f.readline() )
			line = f.readline()
			while line:
				parts = line.split('\t')
				if len( parts[3] ) == 1 and len( parts[4] ) == 1:
					out.write( line )
				line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

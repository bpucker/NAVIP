### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python3 navip_clean_vcf.py
					--fastain <FASTA_INPUT_FILE>
					--fastaout <FASTA_OUTPUT_FILE>
					--gffin <GFF_INPUT_FILE>
					--gffout <GFF_OUTPUT_FILE>
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""

import os, sys

# --- end of imports --- #




def main( arguments ):
	"""! @brief running everything """
	
	fasta_in_file = arguments[ arguments.index('--fastain')+1 ]
	fasta_out_file = arguments[ arguments.index('--fastaout')+1 ]
	gff_in_file = arguments[ arguments.index('--gffin')+1 ]
	gff_out_file = arguments[ arguments.index('--gffout')+1 ]

	with open( fasta_out_file, "w" ) as out:
		with open( fasta_in_file, "r" ) as f:
			counter = 1
			line = f.readline()
			while line:
				if line[0] == ">":
					out.write( '>' + str( counter ) + "\n" )
					counter += 1
				else:
					out.write( line )
				line = f.readline()


	with open( gff_out_file, "w" ) as out:
		with open( gff_in_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					parts[0] = parts[0][-1] + ""
					out.write( "\t".join( parts ) + "\n" )
				line = f.readline()


if '--fastain' in sys.argv and '--fastaout' in sys.argv and '--gffin' in sys.argv and '--gffout' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python3 merge_navip_vs_snpeff_results.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FILE>
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""


import glob, sys, os
import numpy as np

def main( arguments ):
	"""! @brief run everything """
	
	input_folder = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]

	result_files = glob.glob( input_folder + "*/all_genes_with_premature_stop_codons_and_details.txt" )[:200]	#select just 200

	data = {}
	ID_order = []
	number_of_events = []
	for filename in list( sorted( result_files ) ):
		ID = filename.split('/')[-2]
		ID_order.append( ID )
		with open( filename, "r" ) as f:
			counter = 0
			f.readline()	#first header line
			f.readline()	#second header line
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				var_key = parts[0] + "_%_" + parts[1].zfill(7)
				try:
					data[ var_key ].update( { ID: { 'navip': parts[5], 'snpeff': parts[6] } } )
				except KeyError:
					data.update( { var_key: { 	'data': { ID: { 'navip': parts[5], 'snpeff': parts[6] } },
																'chr': parts[0],
																'pos': parts[1],
																'gene': parts[2],
																'ref': parts[3],
																'alt': parts[4]
																			} 
															} )
				counter += 1
				line = f.readline()
			number_of_events.append( counter )

	print( "mean number of events: " + str( np.mean( number_of_events ) ) )
	print( "median number of events: " + str( np.median( number_of_events ) ) )

	all_keys = list( sorted( data.keys() ) )
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "Chr", "Pos", "Gene", "Ref", "Alt" ] + ID_order ) + "\n" )
		for key in all_keys:
			new_line = [ 	data[ key ]['chr'],
									data[ key ]['pos'],
									data[ key ]['gene'],
									data[ key ]['ref'],
									data[ key ]['alt']
								]
			for ID in ID_order:
				try:
					new_line.append( data[ key ]['data'][ ID ]['navip'] + " // " + data[ key ]['data'][ ID ]['snpeff'] )
				except KeyError:
					new_line.append( "." )
			out.write( "\t".join( new_line ) + "\n" )
				
if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

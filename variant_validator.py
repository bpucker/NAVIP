### Boas Pucker ###
### Dakota Howard ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

from operator import itemgetter
import sys, time
# --- end of imports --- #

__usage__ = """
	python variant_validator.py
	--assembly <FULL_PATH_TO_ASSEMBLY_FILE>
	--ref <FULL_PATH_TO_REFERENCE_FILE>
	--invcf <FULL_PATH_TO_INPUT_VCF_FILE>
	--flank <INT, size of query sequence>
	--outvcf <FULL_PATH_TO_OUTPUT_VCF>
	--chr <CHROMOSOME_TO_PROCESS>
	--outerr <FULL_PATH_TO_ERROR_OUTPUT_FILE>
					"""

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split(" ")[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq.upper() } )
					header = line.strip()[1:].split(" ")[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq.upper() } )
	return sequences


def load_variants( vcf_file ):
	"""! @brief load all variants from given VCF file and return sorted list """
	
	variants = []
	variants_split = {}
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if ',' in parts[4]:
					variants.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'ref': parts[3], 'alt': sorted(parts[4].split(',') , key = len, reverse = True ), 'qual' : parts[5], 'filter' : parts[6], 'info' : parts[7] } )
				else:
					variants.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'ref': parts[3], 'alt': parts[4], 'qual' : parts[5], 'filter' : parts[6], 'info' : parts[7] } )
			line = f.readline()
	variants = sorted( variants, key=itemgetter( 'chr', 'pos' ) )
	start = 0
	for idx in range( len( variants ) ):
		if len( variants ) == idx + 1:
			variants_split.update( { variants[idx]['chr'] : ( variants[start:idx + 1] ) } )
		elif variants[idx]['chr'] != variants[idx + 1]['chr']:
			variants_split.update( { variants[idx]['chr'] : ( variants[start:idx] ) } )
			start = idx + 1
	return variants_split


def search_seq_in_assembly( seq, assembly ):
	"""! @brief search seq in assembly """
	
	for contig in assembly.values():
		result = contig.find( seq )
		if result > -1:
			return True
	return False


def check_variants( variants_to_consider, ref_chr, assembly, flank_dist, variants_in_range ):
	"""! @brief construct seq block and check for present; return validation state for all variants """
	
	bad_multi_variants = []
	chr_name = variants_to_consider[ 0 ][ 'chr' ]
	center_pos = variants_to_consider[ 0 ][ 'pos' ]
	seq = ""
	
	#If multiple variants are present
	if isinstance( variants_to_consider[ 0 ][ 'alt' ], list ):
		
		#For each variant present
		for alt_num in range( len( variants_to_consider[ 0 ][ 'alt' ] ) ):
			seq = ""
			num_variants = len( variants_in_range )
			
			#If there are no downstream variants
			if num_variants == 0:
				seq += ref_chr[ variants_to_consider[ 0 ][ 'pos' ] - flank_dist : variants_to_consider[ 0 ][ 'pos' ] - 1]
			 #If there is downstream variants to be added
			else:
				start = variants_to_consider[ 0 ][ 'pos' ] - flank_dist
				for var_num in range( num_variants ):
					seq+= ref_chr[ start : variants_in_range[ var_num ][ 'pos' ] - 1] 
					start = variants_in_range[ var_num ][ 'pos' ]
					seq+= variants_in_range[ var_num ][ 'alt' ]
					
				seq += ref_chr[ start : variants_to_consider[ 0 ][ 'pos' ] - 1]
				
			#Adding focal point seq	
			seq += variants_to_consider[ 0 ]['alt'][ alt_num ]
			match = search_seq_in_assembly( seq, assembly )
			
			#If seq found then break
			if match == True:
				bad_multi_variants = dict(variants_to_consider[ 0 ])
				good_variant = variants_to_consider[ 0 ][ 'alt' ][ alt_num ]
				bad_multi_variants[ 'alt' ].remove(variants_to_consider[ 0 ][ 'alt' ][ alt_num ])
				variants_to_consider[ 0 ][ 'alt' ] = good_variant
				break
		     #If seq not found try reverse
			else:
				seq = revcomp(seq)
				match = search_seq_in_assembly( seq, assembly )
				if match == True:
					bad_multi_variants = dict( variants_to_consider[ 0 ] )
					good_variant = variants_to_consider[ 0 ][ 'alt' ][ alt_num ]					
					bad_multi_variants[ 'alt' ].remove( variants_to_consider[ 0 ][ 'alt' ][ alt_num ])
					variants_to_consider[ 0 ][ 'alt' ] = good_variant
					break
		variants_to_consider[0].update( { 'status': match } )		
		
	#For all cases without multiple variants
	else:
		
		#Add all validated upstream variants to seq in order
		for idx, variant in enumerate( variants_to_consider ): 
			
			#If focal point being added to seq is the first upstream variant
			if idx == 0:
				
				#If there is only one focal point being tested
				if len( variants_to_consider ) == 1: 
					num_variants = len( variants_in_range )
					
					#For the first downstream variant in seq
					if num_variants == 0:
						seq += ref_chr[ variants_to_consider[ 0 ][ 'pos' ] - flank_dist : variants_to_consider[ 0 ][ 'pos' ] - 1]
					 #For the rest of the downstream variants in seq
					else:
						start = variants_to_consider[ 0 ][ 'pos' ] - flank_dist
						for var_num in range( num_variants ):
							seq+= ref_chr[ start : variants_in_range[ var_num ][ 'pos' ] - 1] 
							start = variants_in_range[ var_num ][ 'pos' ]
							seq+= variants_in_range[ var_num ][ 'alt' ]
						seq += ref_chr[ start :variants_to_consider[ 0 ][ 'pos' ] - 1]
						
					seq += variants_to_consider[ 0 ][ 'alt' ]
				  #If there are multiple focal points being tested
				else:
					seq += ref_chr[ center_pos - flank_dist : variant[ 'pos' ] - 1 ] + variant[ 'alt' ]
		    	#If focal point is last upstream variant being added
			elif idx == len( variants_to_consider )-1:
				seq += ref_chr[ variants_to_consider[ idx - 1 ][ 'pos' ] + len( variants_to_consider[ idx - 1 ][ 'alt' ] ) - 1 : variant[ 'pos' ] - 1 ] + variant[ 'alt' ] + ref_chr[ variant[ 'pos' ] + len( variant[ 'alt' ] ) - 1 : center_pos ]	
		     	#For all other upstream variants
			else:
				seq += ref_chr[ variants_to_consider[ idx - 1 ][ 'pos' ] : variant[ 'pos' ] - 1 ] + variant[ 'alt' ]

		match = search_seq_in_assembly( seq, assembly )
		#If false reverse seq
		if match == False:
			seq = revcomp( seq )
			match = search_seq_in_assembly( seq, assembly )
		
		#If true validate all upstream variants
		if match == True:
			for num in range( len( variants_to_consider ) ):
				variants_to_consider[num].update( { 'status' : match } )
				# if multiple SNPs fail together then split and check them individually
		else:   #For each variant to be checked
		 	for idx, variant in enumerate( variants_to_consider ):
				seq = ""
				num_variants = len( variants_in_range )
				
				#If there are no downstream variants to be added to seq
				if num_variants == 0:
					seq += ref_chr[ variants_to_consider[ idx ][ 'pos' ] - flank_dist : variants_to_consider[ idx ][ 'pos' ] - 1]
					#If there are downstream variants to be added to seq
				else:
					start = variants_to_consider[ idx ][ 'pos' ] - flank_dist
					
					for var_num in range( num_variants ):
						seq+= ref_chr[ start: variants_in_range[ var_num ][ 'pos' ] - 1] 
						start = variants_in_range[ var_num ][ 'pos' ]
						seq+= variants_in_range[ var_num ][ 'alt' ]
				
					seq += ref_chr[ start : variants_to_consider[ idx ][ 'pos' ] - 1]
				
				seq += variants_to_consider[ idx ][ 'alt' ]
				match = search_seq_in_assembly( seq, assembly )
				
				#If false reverse seq
				if match == False:
					seq = revcomp(seq)
					match = search_seq_in_assembly( seq, assembly )
				
				#If true add variant to downstream variants list
				if match == True:
					variants_in_range.append(variants_to_consider[idx])
					
				variants_to_consider[idx].update( { 'status': match } )		

	return variants_to_consider, bad_multi_variants


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def variant_validation( variants, output_vcf, flank_dist, ref_seq, assembly, false_multi_variants, chromosome):
	"""! @Collects all downstream variants  """
	
	validated_variants = []
	bad_multi_variants = []
	valid_counter = 0
	invalid_counter = 0
	
	with open( output_vcf, "w", 0 ) as out:

		#If first being added add the header
		if chromosome[-1] == "1":
			out.write( "#Chr\tPos\tID\tRef\tAlt\tQual\tStatus\tInfo\n" )
		i = 0
		
		while i < len( variants ):
			
			#Completion bar
			bar_size = 80
			fill = int(round(bar_size * float(i) / float(len(variants))))
			perc = "%.2f" % (round(100.0 * float(i) / float(len(variants)),2))
			bar = '#' * fill + '-' * (bar_size - fill)
			valid = "%.2f" % (round((((valid_counter + 1) / (invalid_counter + 1 + valid_counter)) * 100),2))
			sys.stdout.write( "validity: " + str( valid ) + "% " + str( bar ) + " " + chromosome + " " + str( perc ) + "% complete\n" )
			sys.stdout.flush()
			
			variants_to_consider = []	#list with all variants needed to construct query sequence
			start_pos = variants[ i ]['pos']-flank_dist	#upstream end of range of variants to consider
			end_pos = variants[ i ]['pos']+flank_dist	#downstream end of range of variants to consider
			variants_in_range = []     #downstream variants needed to construct query sequence
			
			#For validated downstream variants to be added to seq
			if len(validated_variants) != 0:
				n = len(validated_variants)
				pos = validated_variants[ n - 1 ]['pos']

				#Collect all variants within downstream range of flank_dist
				while pos >= start_pos:
					 n = n - 1
					 if n == -1:
							break
					 if validated_variants[ n ]['pos'] >= start_pos:
						variants_in_range.insert ( 0 , validated_variants[ n ] )
					 pos = validated_variants[ n ]['pos']
	
			#break if the end of the variant list is reached
			if i == len( variants ):
				break
		
			variants_to_consider.append( variants[ i ] )	#add current focal point variant

			#If first SNP added has multiple variants check it byitself
			if isinstance(variants[i]['alt'],list):
				pass	
			 #use this case for first variant, because no check of prior variants is required """
			elif len( validated_variants ) == 0:	
				prev_val_variants = 0

				#Used to check multiple variants at a time
				k = 1
				while variants[ i+k ]['pos'] <= end_pos and variants[ i+k ]['chr'] == variants[ i ]['chr']:
					
					#Break if next SNP to be checked has multiple variants
					if isinstance(variants[ i+k ]['alt'],list):
						break
					variants_to_consider.append( variants[ i+k ] )
					k+= 1
					if i+k >= len( variants ):
						break
			 #this case is used for almost all variants (except first one)
			else:
				
				#Used to check multiple variants at a time
				k = 1
				if i+k == len( variants ):
					pass
				else:
					while variants[ i+k ]['pos'] <= end_pos and variants[ i+k ]['chr'] == variants[ i ]['chr']:
						
						#Break if next SNP to be checked has multiple variants					
						if isinstance(variants[ i+k ]['alt'],list):
							break
						variants_to_consider.append( variants[ i+k ] )
						k+= 1
						if i+k >= len( variants) :
							break

			checked_variants = check_variants( variants_to_consider, ref_seq, assembly, flank_dist, variants_in_range )	#construct seq and check
			
			#If checked_variants returned an SNP that had multiple variants place bad variants in a separate file
			if len(checked_variants[1]) != 0:
				bad_multi_variants.append(checked_variants[1])
				
			checked_variants = checked_variants[0] 	
			#Write status of variants to file
			for idx, variant in enumerate(checked_variants):
				if variant['status'] == True:
					validated_variants.append( variant )
					out.write( "\t".join( map( str, [ variant['chr'], variant['pos'], ".", variant['ref'], variant['alt'], variant['qual'], "PASS", variant['info'] ] ) ) + '\n' )
					valid_counter += 1
				else:
					out.write( "\t".join( map( str, [ variant['chr'], variant['pos'], ".", variant['ref'], variant['alt'], variant['qual'], "FAIL", variant['info'] ] ) ) + '\n' )
					invalid_counter += 1
			i+= len( checked_variants ) #Update current position in the VCF file

	#Write to separate file SNP's that were not valid where multiple were present
	with open( false_multi_variants, "w", 0 ) as out:
		if chromosome[-1] == "1":
			out.write( "#Chr\tPos\tID\tRef\tAlt\tQual\tStatus\tInfo\n" )
		for size in range(len(bad_multi_variants)):
			out.write( "\t".join( map( str, [ bad_multi_variants[size]['chr'], bad_multi_variants[size]['pos'], ".", bad_multi_variants[size]['ref'], ",".join(bad_multi_variants[size]['alt']), variant['qual'], "FAIL", variant['info'] ] ) ) + '\n' )

	return [validated_variants]


def main( arguments ):
	"""! @brief run all parts of script """
	
	"""@WARNING: not suitable for investigating species with more than 9 chromosomes! """
	
	assembly_file = arguments[ arguments.index( '--assembly' )+1 ]
	ref_seq_file = arguments[ arguments.index( '--ref' )+1 ]
	input_vcf = arguments[ arguments.index( '--invcf' )+1 ]
	chromosome = arguments[ arguments.index( '--chr' )+1 ]
	output_vcf = arguments[ arguments.index( '--outvcf' )+1 ]
	output_multi_SNP_errors = arguments[ arguments.index( '--outerr' )+1 ]
	
	ref_seq = load_sequences( ref_seq_file )
	assembly = load_sequences( assembly_file )
	variants = load_variants( input_vcf )
	if '--flank' in sys.argv:
		flank_dist = int( arguments[ arguments.index( '--flank' )+1 ] )
	else:
		flank_dist = 30
	
	variant_validation( variants[ chromosome ], output_vcf, flank_dist, ref_seq[ chromosome ], assembly, output_multi_SNP_errors, chromosome )


if __name__ == '__main__':
	if '--assembly' in sys.argv and '--ref' in sys.argv and '--invcf' in sys.argv and '--chr' in sys.argv and '--outvcf' in sys.argv and '--outerr' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

### Boas Pucker ###
### Dakota Howard ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

from operator import itemgetter
import sys, glob, re, os, time, datetime, shutil

# --- end of imports --- #

__usage__ = """
	python variant_validaton_wrapper.py
	--assembly <FULL_PATH_TO_ASSEMBLY_FILE>
	--ref <FULL_PATH_TO_REFERENCE_FILE>
	--vcf <FULL_PATH_TO_INPUT_VCF_FILE>
	--flank <INT, size of query sequence>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	--script <FULL_PATH_TO variant_validator.py>
					"""


def submit_jobs_to_cluster( prefix, assembly_file, ref_seq_file, input_vcf, script, flank_dist, chromosomes ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	para_jobs = 30
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, chromosome in enumerate( chromosomes ):
		ID = "X_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + chromosome + '.sh'
		out_file = prefix + ID  + chromosome+ '.out'
		err_file = prefix + ID + chromosome + '.err'
		
		out_vcf = prefix + chromosome + ".vcf"
		out_err = prefix + chromosome + ".err.txt"
		cmd = " ".join( map( str, [ "python", script, "--assembly", assembly_file, "--ref", ref_seq_file, "--invcf",  input_vcf, "--flank", flank_dist, "--outvcf", out_vcf, "--outerr", out_err, "--chr", chromosome ] ) )
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=1G",
																"-l arch=lx-amd64",
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep(1)
		os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "X_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "X_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )


def final_processing( prefix ):
	"""! @brief blt processing of BLAT results for identification of best hit """
	
	vcfs = sorted( glob.glob( prefix + "*.vcf" ) )
	error_files = sorted( glob.glob( prefix + "*.err.txt" ) )
	
	os.popen( "cat " + " ".join( vcfs ) + " > " + prefix + "final.vcf" )
	os.popen( "cat " + " ".join( error_files ) + " > " + prefix + "final.err.txt" )


def load_sequences(fasta_file):
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


def main( arguments ):
	"""! @brief runs everything """
	
	assembly_file = arguments[ arguments.index( '--assembly' )+1 ]
	ref_seq_file = arguments[ arguments.index( '--ref' )+1 ]
	input_vcf = arguments[ arguments.index( '--vcf' )+1 ]
	prefix = arguments[ arguments.index( '--out' )+1 ]
	script = arguments[ arguments.index( '--script' )+1 ]
	
	if '--flank' in sys.argv:
		flank_dist = int( arguments[ arguments.index( '--flank' )+1 ] )
	else:
		flank_dist = 30
	
	ref = load_sequences( ref_seq_file )
	
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	submit_jobs_to_cluster( prefix, assembly_file, ref_seq_file,  input_vcf, script, flank_dist, sorted( ref.keys() ) )
	final_processing( prefix )


if __name__ == '__main__':
	if '--assembly' in sys.argv and '--ref' in sys.argv and '--vcf' in sys.argv and '--out' in sys.argv and '--script' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	

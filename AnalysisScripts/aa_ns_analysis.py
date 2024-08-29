### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.21 ###

__usage__ = """
					python aa_ns_analysis.py
					--in <NAVIP_OUTPUT_FILE>
					--genes <GENE_IDs_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--gff <GFF_FILE_WITH_ALL_GENE_IDs>
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""

import sys, os, re
from operator import itemgetter
import plotly.io as pio
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import mannwhitneyu
import statistics

# ---- end of imports --- #

def load_gene_IDs( genes_info_file ):
	"""! @brief load gene IDs from input file """
	
	genes = {}
	with open( genes_info_file, "r" ) as f:
		lines = f.read().strip().split('\n')
		for line in lines:
			genes.update( { line: None } )
	return genes


def load_variants_per_gene( navip_vcf_file ):
	"""! @brief load aa_n / aa_s status per variant per gene from given VCF file """
	
	variants_per_gene = {}
	with open( navip_vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				anno_string = parts[7].split('NAV2=')[1].split(';')[0]
				gene_id = anno_string.split('|')[3]
				aa_change = anno_string.split('|')[10][2:]	#extract amino acid change string and remove the 'p.' prefix
				if '*' in aa_change:	#'*' might indicate a following premature stop codon in the protein (needs to be removed)
					if 'ext' in aa_change:	#stop codon lost (extension of protein)
						aa_change = aa_change.split('ext')[0]	#cut down to substitution of stop codon by amino acid (counted later)
					elif aa_change.count('*') > 1:	#if more than one stop codon (stop codon changed)
						aa_change = "*".join( aa_change.split('*')[:2] )	#keep leading asterisk indicating a stop codon
					else:
						aa_change = aa_change.split('*')[0]
				aaS_status = True	#mutation is a synonymous substitution (no change in encoded amino acid)
				if 'ins' in aa_change:	#'ins' = insertion
					aaS_status = False
				elif 'del' in aa_change:	#'del' = deletion
					aaS_status = False
				elif 'dup' in aa_change:	#'dup' = duplication
					aaS_status = False
				elif aa_change == "":	#frameshift causing premature stop
					aaS_status = False
				else:
					try:
						pos = re.findall( "[0-9]+", aa_change )
						if len( pos ) == 1:
							pos = pos[0]
							aa1, aa2 = aa_change.split( pos )
							if aa1 != aa2:
								aaS_status = False
						else:
							print( "ERROR:" + aa_change + "\t" +anno_string )
					except IndexError:
						print( "ERROR: " + aa_change + "\t" + anno_string )
				# --- add information to existing dictionary --- #
				try:
					variants_per_gene[ gene_id ].append( aaS_status )
				except KeyError:
					variants_per_gene.update( { gene_id: [ aaS_status ] } )
			line = f.readline()
	return variants_per_gene


def calculate_aaNaaS_ratio_per_gene( variants_per_gene ):
	"""! @brief calculate aaN/aaS ratio per gene """
	
	aaNaaS_ratio_per_gene = {}
	for gene in list( variants_per_gene.keys() ):
		variants = variants_per_gene[ gene ]
		aaN = float( len( variants ) - sum( variants ) )	#all variants minus synonymous ones
		aaS = sum( variants )
		aaS += 1
		aaN += 1
		aaNaaS_ratio_per_gene.update( { gene: aaN / aaS } )
	return aaNaaS_ratio_per_gene


def plot_aaNaaS_ratios_of_groups( aaNaaS_ratio_per_gene, genes, figure_file, supplied_genes, display_cutoff=10 ):
	"""! @brief plot aaN/aaS ratios of defined genes against all other genes """
	
	group1_values = []
	group2_values = []
	if len( supplied_genes ) > 0:
		for gene in supplied_genes:
			try:
				genes[ gene ]
				try:
					group1_values.append( min( [ aaNaaS_ratio_per_gene[ gene ], display_cutoff ] ) )
				except KeyError:
					group1_values.append( 1 )	#if gene ID not in variant data, gene does not contain variants (1 is used in this case)
			except KeyError:
				try:
					group2_values.append( min( [ aaNaaS_ratio_per_gene[ gene ], display_cutoff ] ) )
				except KeyError:
					group2_values.append( 1 )	#if gene ID not in variant data, gene does not contain variants (1 is used in this case)
	else:
		for gene in list( aaNaaS_ratio_per_gene.keys() ):
			try:
				genes[ gene ]
				group1_values.append( min( [ aaNaaS_ratio_per_gene[ gene ], display_cutoff ] ) )
			except KeyError:
				group2_values.append( min( [ aaNaaS_ratio_per_gene[ gene ], display_cutoff ] ) )
	
	group1_values_df = pd.DataFrame( group1_values, columns=["group1"] )
	group2_values_df = pd.DataFrame( group2_values, columns=["group2"] )
	
	ng1 = len( group1_values )
	ng2 = len( group2_values )
	
	fig = go.Figure()
	
	legend1 = 'genes with stop_gained\n(n=' + str( ng1 ) + ")"
	fig.add_trace(go.Violin( y=group1_values_df['group1'], name=legend1, showlegend=False, meanline=dict(visible=True, color='black') ) )	#, line_color='blue', scalemode='width'
	
	legend2 = 'other genes\n(n=' + str( ng2 ) + ")"
	fig.add_trace(go.Violin( y=group2_values_df['group2'], name=legend2, showlegend=False, meanline=dict(visible=True, color='black') ) )	#, line_color='red', scalemode='width'
	
	#fig.update_traces(meanline_visible=True)	#add line indicating the mean
	fig.update_layout(violingap=0, violinmode='overlay')
	fig.update_yaxes(range=[0, display_cutoff*1.1])
	fig.update_layout(yaxis_title=r"$aa_{N}/aa_{S} \text{ per gene}$")
	
	#fig.update_layout(xaxis=dict(tickfont=dict(size=20)),  #change x-tick fontsize
    #              yaxis=dict(titlefont=dict(size=20)))  #change y-axis labele fontsize


	pio.write_image( fig, figure_file, scale=10 )
	return group1_values, group2_values


def load_gene_IDs_from_gff_file( gff_file, IDtag="ID" ):
	"""! @brief load gene IDs from given GFF3 file """
	
	gene_IDs = []
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					gene = parts[-1].split( IDtag + '=')[1]
					if ";" in gene:
						gene = gene.split(';')[0]
					gene_IDs.append( gene )
			line = f.readline()
	return gene_IDs


def main( arguments ):
	"""! @brief run everything """
	
	navip_vcf_file = arguments[ arguments.index( '--in' )+1 ]	#NAVIP output VCF file for aaN/aaS determination
	genes_info_file = arguments[ arguments.index( '--genes' )+1 ]	#list of genes (one ID per line) to compare against all other genes
	output_folder =  arguments[ arguments.index( '--out' )+1 ]
	
	if not output_folder[-1] == "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--gff' in arguments:	#this option is helpful to have a more complete background with genes not harbouring any sequence variants
		gff_file = arguments[ arguments.index( '--gff' )+1 ]
		supplied_genes = load_gene_IDs_from_gff_file( gff_file )
	else:
		supplied_genes = []
	
	genes = load_gene_IDs( genes_info_file )
	variants_per_gene = load_variants_per_gene( navip_vcf_file )
	
	# --- determine aaN/aaS ratios per gene --- #
	aaNaaS_ratio_per_gene = calculate_aaNaaS_ratio_per_gene( variants_per_gene )
	
	# --- plot ratios of group specified by input file against all other genes --- #
	figure_file = output_folder + "violin_plot.png"
	group1, group2 = plot_aaNaaS_ratios_of_groups( aaNaaS_ratio_per_gene, genes, figure_file, supplied_genes )
	
	# -- run statistical test --- #
	u_value, p_value = mannwhitneyu( group1, group2, alternative='two-sided' )
	print( "U: " + str( u_value ) + "\tp-value: " + str( p_value ) )
	print( "group1: mean=" + str( statistics.mean( group1 ) ) + "\tmedian=" + str( statistics.median( group1 ) ) )
	print( "group2: mean=" + str( statistics.mean( group2 ) ) + "\tmedian=" + str( statistics.median( group2 ) ) )


if '--out' in sys.argv and '--in' in sys.argv and '--genes' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

# Neighborhood-Aware Variant Impact Predictor (NAVIP)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2620396.svg)](https://doi.org/10.5281/zenodo.2620396)

Variant calling is a process to identify differences between strains, accessions, genotypes, or individuals.
Sequence variants can be analyzed to predict their functional effects on encoded proteins based on available structural annotations.
Many existing tools operate on a variant-by-variant basis, and per-variant functional impact predictions are generally considered accurate.
However, challenging cases arise where multiple neighboring variants must be considered simultaneously, especially when these variants influence each other's functional impact.

The Neighborhood-Aware Variant Impact Predictor (NAVIP) addresses this problem by considering all variants that may affect a coding sequence (CDS) during the prediction process.
This comprehensive approach increases the accuracy of predicting functional consequences by considering the broader genomic context surrounding the target variants.
To use NAVIP, users must provide a Variant Call Format (VCF) file, a genomic FASTA file, and an associated Gene Feature Format (GFF3) file as input.
The tool is also freely available on our web server at: https://pbb-tools.de/NAVIP.

Important advice: Please only provide single nucleotide variants (SNVs) and insertions/deletions (InDels). Multiple nucleotide variants (MNVs) are currently not supported.

# Usage

For the main program, there are no strict dependencies other than Linux, [Python 3](https://www.python.org), and [matplotlib](https://matplotlib.org).
For most use cases, it is sufficient to download the source code / clone the git repository and run the script ***runnavip.sh*** with the required arguments, which will guide the user through the entire NAVIP processing pipeline:

```
Usage: runnavip.sh -i <invcf> -g <ingff> -f <infasta> -o <outpath>

where:
   -i | --invcf   Specify the input VCF file
   -g | --ingff   Specify the input GFF file
   -f | --infasta Specify the input FASTA file
   -o | --outpath Specify the output path
   -h | --help    Show this help message
```

The pipeline will generate three folders with the analysis results in the specified output directory: VCF_Preprocessing, NAVIP_Main_Output, and SFA_Output.
For additional options and instructions on how to run individual modules of NAVIP, please have a look at our wiki page: https://github.com/bpucker/NAVIP/wiki.

# Reference (how to cite):

Baasner, J.-S., Rempel, A., Howard, D., Pucker, B. (2024). NAVIP: Unraveling the Influence of Neighboring Small Sequence Variants on Functional Impact Prediction. PLoS Comput Biol 21(2): e1012732. doi: [10.1371/journal.pcbi.1012732](https://doi.org/10.1371/journal.pcbi.1012732).

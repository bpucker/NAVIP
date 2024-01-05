NAVIP has three existing modules:
 1) VCF preprocessing,
 2) the NAVIP main program,
 3) and one simple first analysis of the created data.

You can choose the module with "--mode <module>".
The module shortcuts are "pre", "main" and "sfa" and there is a vcf-format-check "vcfc".


The VCF-Check needs one argument:
 "--invcf <path_with_file>"

VCF preprocessing needs two arguments:
 "--invcf <path_with_file>" and "--outpath <path_to_folder>"

Please be aware, that no new folder will be created.

The NAVIP main program needs four arguments:
 "--invcf <path_with_file>", "--ingff <path_with_file>", "--infasta <path_with_file>" and "--outpath <path_to_folder>"

The best possible output will be available, when the VCF file is 'corrected' by the preprozessing.
However, NAVIP will still be able to deal with most of the 'normal' VCF data and will do its best.

The SFA module needs three arguments:
 "--innavipvcf <path_with_file>", "--innavipfasta <path_with_file>" and "--outpath <path_to_folder>"


Example: VCF preprocessing:
python3 navip.py \
 --mode pre \
 --invcf /.../small_variants.vcf \
 --outpath /.../VCF_Preprocessing/

Example: NAVIP main:
python3 navip.py \
 --mode main \
 --invcf VCF_Preprocessing/first.vcf \
 --ingff /.../Araport11_GFF3_genes_transposons.201606.gff \
 --infasta /.../TAIR10.fa \
 --outpath /.../NAVIP_Main_Output/

Example: SFA:
python3 navip.py \
 --mode sfa \
 --innavipvcf /.../All_VCF.vcf \
 --innavipfasta /.../all_transcripts_data.fa \
 --outpath /.../SFA_Output/

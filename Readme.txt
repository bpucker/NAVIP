NAVIP has three existing modules: VCF preprocessing, the NAVIP main program and one simple first analysis of the created data.
You can choose the module with "--mode <module>".
The module shortcuts are "pre","main" and "sfa" and there is a vcf-format-check "vcfc".

For the VCF-Check: "--mode vcfc --invcf <path_with_file>"

VCF preprocessing needs two more arguments:
 "--invcf <path_with_file>" and "--outpath <path_to_folder>" 
Please be aware, that no new folder will be created. 

The NAVIP main program needs four arguments:
"--invcf <path_with_file>", "--ingff <path_with_file>", "--infasta <path_with_file>" and "--outpath <path_to_folder>"
The best possible output will be available, when the VCF file is 'corrected' by the preprozessing. However, NAVIP will still be able to deal with most of the 'normal' VCF data and will do its best.

The SFA module needs three arguments:
"--innavipvcf <path_with_file>","--innavipfasta <path_with_file>" and "--outpath <path_to_folder>"

Example: VCF preprocessing:
"python3 navip.py --mode pre --invcf /prj/gf-arabseq/project_VariantAnnotation/data/20160806_small_variants.vcf --outpath /prj/gf-arabseq/project_VariantAnnotation/data/VCF_Preprocessing/"

Example: NAVIP main:
"python3 navip.py --mode main --invcf /prj/gf-arabseq/project_VariantAnnotation/data/VCF_Preprocessing/first.vcf --ingff /prj/gf-arabseq/data/Araport11_official_release/Araport11_GFF3_genes_transposons.201606.gff --infasta /prj/gf-arabseq/data/TAIR10/TAIR10.fa --outpath /prj/gf-arabseq/project_VariantAnnotation/data/NAVIP_Main_Output/"

Example: SFA:
"python3 navip.py --mode sfa --innavipvcf /prj/gf-arabseq/project_VariantAnnotation/data/NAVIP_Main_Output/All_VCF.vcf --innavipfasta /prj/gf-arabseq/project_VariantAnnotation/data/NAVIP_Main_Output/all_transcripts_data.fa --outpath /prj/gf-arabseq/project_VariantAnnotation/data/SFA_Output/"

__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"

import sys
import VCF_preprocess,Coordinator, sfa, VCF_Format_Check
###
example_pre = "\"python3 navip.py --mode pre --invcf /prj/gf-arabseq/project_VariantAnnotation/data/20160806_small_variants.vcf" \
              " --outpath /prj/gf-arabseq/project_VariantAnnotation/data/VCF_Preprocessing/\""

example_main = "\"python3 navip.py --mode main --invcf /prj/gf-arabseq/project_VariantAnnotation/data/VCF_Preprocessing/first.vcf" \
               " --ingff /prj/gf-arabseq/data/Araport11_official_release/Araport11_GFF3_genes_transposons.201606.gff" \
               " --infasta /prj/gf-arabseq/data/TAIR10/TAIR10.fa" \
               " --outpath /prj/gf-arabseq/project_VariantAnnotation/data/NAVIP_Main_Output/\""

example_sfa = "\"python3 navip.py --mode sfa --innavipvcf /prj/gf-arabseq/project_VariantAnnotation/data/NAVIP_Main_Output/All_VCF.vcf" \
              " --innavipfasta /prj/gf-arabseq/project_VariantAnnotation/data/NAVIP_Main_Output/all_transcripts_data.fa" \
              " --outpath /prj/gf-arabseq/project_VariantAnnotation/data/SFA_Output/\""

readmetext = "NAVIP has three existing modules: VCF preprocessing, the NAVIP main program and one simple first analysis of the created data.\n" \
             "You can choose the module with \"--mode <module>\".\n" \
             "The module shortcuts are \"pre\",\"main\" and \"sfa\" and there is a vcf-format-check \"vcfc\".\n\n" \
             "For the VCF-Check: \"--mode vcfc --invcf <path_with_file>\"\n\n"\
             "VCF preprocessing needs two more arguments:" \
             "\n \"--invcf <path_with_file>\" and \"--outpath <path_to_folder>\" \n" \
             "Please be aware, that no new folder will be created." \
             " \n\n" \
             "The NAVIP main program needs four arguments:\n" \
             "\"--invcf <path_with_file>\", \"--ingff <path_with_file>\", \"--infasta <path_with_file>\" " \
             "and \"--outpath <path_to_folder>\"\n" \
             "The best possible output will be available, when the VCF file is 'corrected' by the preprozessing. However, NAVIP will still be able to deal with most of the 'normal' VCF data and will do its best.\n" \
             "\n" \
             "The SFA module needs three arguments:\n" \
             "\"--innavipvcf <path_with_file>\",\"--innavipfasta <path_with_file>\" and \"--outpath <path_to_folder>\"\n\n" \
             "Example: VCF preprocessing:\n" + example_pre + "\n\n" \
             "Example: NAVIP main:\n" + example_main + "\n\n"\
             "Example: SFA:\n" + example_sfa + "\n"


#pre == vcf preprocessing
#   input:  one vcf file, path to output folder
#   output: two vcf files
#   desc_i: original vcf file
#   desc_o: triallele variants splitted into two vcf files, including all 'normal' variants each

#main == navip main program
#   input:  one vcf, one gff3, one fasta file, path to output folder
#   output: one vcf, one fasta file, one log file
#   desc_i: preprocessed vcf
#   desc_i: vcf, gff3 and fasta file from the same organism
#   desc_o: vcf file with all original vcf data rows including description from neighbour variant interaction
#   desc_o: fasta file which contains the new and old dna/aa seq of every transcript
#   desc_o: log file contains information about errors or problematic variants/transcripts

#sfa == simple first analysis
#   input:  one vcf and one fasta file, path to output folder
#   output: five vcf amd two txt files + four plots
#   desc_i: vcf file from navip main program, including description from neihgbour variant interaction
#   desc_i: fasta output file from navip main program
#   desc_o: vcf file with interactions only
#   desc_o: vcf file with DVC (double variant codons) only
#   desc_o: vcf file with neutralized indels only
#   desc_o: vcf file with unsync DVC only,
#   desc_o: table with detailed DVC aminoacid change
#   desc_o: list of differences between AA in transcripts
#   desc_o: plots: different view on the transcript length differences


if __name__ == '__main__':

    #teststatement = "--mode main --invcf /prj/gf-arabseq/project_VariantAnnotation/navip_bugsearch/first.vcf --ingff /grp/gf/Alle_temp_Ordner/datenaustausch_jan_sarah/CRIBI_V2.1_extended_20161128.gff3 --infasta /grp/gf/Alle_temp_Ordner/datenaustausch_jan_sarah/Vv12x_CRIBI.fa --outpath /prj/gf-arabseq/project_VariantAnnotation/navip_bugsearch/"
    #teststatement = "--mode main --invcf /vol/tmp/Jens_Jan_Debugging_Share/debugging/Analyse_Influnence_Of_SNP/first.vcf --ingff /vol/tmp/Jens_Jan_Debugging_Share/debugging/CRIBI_V2.1_extended_20161128.gff3 --infasta /vol/tmp/Jens_Jan_Debugging_Share/debugging/Vv12x_CRIBI.fa --outpath /vol/tmp/Jens_Jan_Debugging_Share/de_2"

    #for jens with possible bugs inside
    #jenstest = "--mode pre --invcf /prj/gf-arabseq/project_VariantAnnotation/Bugsearch_for_Jens/SEY/Seyval_variants.vcf " \
    #           "--outpath /prj/gf-arabseq/project_VariantAnnotation/Bugsearch_for_Jens/SEY/"
    jenstest = "--mode main --invcf /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/Bugsearch_for_Jens/PN/testfirst.vcf " \
               "--ingff /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/Bugsearch_for_Jens/test.gff3 " \
               "--infasta /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/Bugsearch_for_Jens/Vv12x_CRIBI.fa " \
               "--outpath /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/Bugsearch_for_Jens/PN/bugsearch/"
    #test = "--mode main --invcf /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/fakeVCFdata.vcf " \
    #           "--ingff /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/fakegff.gff " \
    #           "--infasta /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/fakegenomdata.fa " \
    #           "--outpath /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/Output/"
    #jenstest = "--mode sfa --innavipvcf /prj/gf-arabseq/project_VariantAnnotation/Bugsearch_for_Jens/SEY/second/All_VCF.vcf " \
    #           "--innavipfasta /prj/gf-arabseq/project_VariantAnnotation/Bugsearch_for_Jens/SEY/second/all_transcripts_data.fa " \
    #           "--outpath /prj/gf-arabseq/project_VariantAnnotation/Bugsearch_for_Jens/SEY/second/sfa/"

    sys.argv = jenstest.split(" ")

    if "--mode" in sys.argv:
        args = sys.argv
        if args[args.index("--mode")+1] == "pre" \
                and "--invcf" in args \
                and "--outpath" in args:
            print("Start: VCF preprocessing")
            VCF_preprocess.vcf_preprocessing(args[args.index("--invcf") + 1],
                                               args[args.index("--outpath") + 1])


        elif args[args.index("--mode")+1] == "main" \
                and "--invcf" in args \
                and "--ingff" in args \
                and "--infasta" in args \
                and "--outpath" in args:
            print("Start: NAVIP")
            Coordinator.navip_main_coordinator(args[args.index("--invcf") + 1],
                                               args[args.index("--ingff") + 1],
                                               args[args.index("--infasta") + 1],
                                               args[args.index("--outpath") + 1])


        elif args[args.index("--mode")+1] == "sfa" \
                and "--innavipvcf" in args\
                and "--innavipfasta" in args\
                and "--outpath" in args:
            print("Start: simple first analysis")
            sfa.sfa_main(args[args.index("--innavipvcf") + 1],
                        args[args.index("--innavipfasta") + 1],
                        args[args.index("--outpath") + 1])
        elif args[args.index("--mode")+1] == "vcfc" \
                and "--invcf" in args:
            VCF_Format_Check.VCF_Check(str(args[args.index("--invcf") + 1]))

        else:
            print("Arguments are invalid."
                  "\nPlease look into the Readme.txt or create a new one with \"--read\"")
    elif "--read" in sys.argv:
        readme = open("Readme.txt","w")
        readme.write(readmetext)
        readme.close()
    else:
        print("Arguments are invalid."
              "\nPlease look into the Readme.txt or create a new one with \"--read\"")

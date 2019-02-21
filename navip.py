__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"

import sys
import os
from time import sleep
import VCF_preprocess,Coordinator, VCF_Format_Check
try:
    import sfa
except ModuleNotFoundError:
    print("ModuleNotFoundError")
    print(str(sys.exc_info()[0]))


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

    # janstest = "--mode main --invcf /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/fakeVCFdata.vcf " \
    #           "--ingff /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/fakegff.gff " \
    #           "--infasta /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/fakegenomdata.fa " \
    #           "--outpath /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/Output/"
    # janstest = "--mode sfa --innavipvcf /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/Output/All_VCF.vcf " \
    #           "--innavipfasta /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/Output/all_transcripts_data.fa " \
    #           "--outpath /prj/gf-arabseq/project_VariantAnnotation/members/janbaas/test_Data/Output/sfa/"

    #snpeffTest = "--mode pre --invcf /homes/janbaas/NAVIP_prj/members/bpucker/20180923_SnpEff_on_validated_variants/valid_SNPs_tair10.vcf " \
    #             "--outpath /homes/janbaas/NAVIP_prj/members/Snpeff_VS_NAVIP/2019-02-13/Pre-Module/"

    snpeffTest = "--mode main --invcf /homes/janbaas/NAVIP_prj/members/Snpeff_VS_NAVIP/2019-02-13/Pre-Module/first.vcf " \
                 "--ingff /prj/gf-arabseq/project_VariantAnnotation/members/Snpeff_VS_NAVIP/sources/old/TAIR10_GFF3_genes_IDs.gff " \
                 "--infasta /prj/gf-arabseq/project_VariantAnnotation/members/Snpeff_VS_NAVIP/sources/TAIR10.fa " \
                 "--outpath /homes/janbaas/NAVIP_prj/members/Snpeff_VS_NAVIP/2019-02-13/better_format/ " \
                 "--ow"

    # "--ingff /prj/gf-arabseq/project_VariantAnnotation/members/Snpeff_VS_NAVIP/sources/old/TAIR10_GFF3_genes_IDs.gff " \
    # "--ingff /prj/gf-arabseq/project_VariantAnnotation/members/Snpeff_VS_NAVIP/sources/Araport11_GFF3_genes_transposons.201606_WITHOUT_ORGANELS.gff " \

    sys.argv = snpeffTest.split()
    sleep(2) # time for creating directories. sometimes useful if python is to fast
    if "--mode" in sys.argv:
        args = sys.argv
        overwriting = False
        if "--ow" in sys.argv:
            overwriting = True
        if args[args.index("--mode")+1] == "pre" \
                and "--invcf" in args \
                and "--outpath" in args:
            print("Start: VCF preprocessing")
            if not overwriting:
                try:
                    file = open(args[args.index("--outpath") + 1] + "first.vcf", "r")
                    file.close()
                    sys.exit(args[args.index("--outpath") + 1] + "first.vcf is already existing.")
                except:
                    pass
                try:
                    file = open(args[args.index("--outpath") + 1] + "second.vcf", "r")
                    file.close()
                    sys.exit(args[args.index("--outpath") + 1] + "second.vcf is already existing.")
                except:
                    pass
            try:
                file = open(args[args.index("--invcf")+1],'r')
                file.close()
                file = open(args[args.index("--outpath")+1] + "try-file-for-exceptions", "w")
                file.write("test")
                file.close()
                os.remove(args[args.index("--outpath")+1] + "try-file-for-exceptions")
            except FileNotFoundError:
                e = sys.exc_info()[0]
                sys.exit(e)
            except PermissionError:
                e = sys.exc_info()[0]
                sys.exit(e)
            except :
                e = sys.exc_info()[0]
                sys.exit(e)

            VCF_preprocess.vcf_preprocessing(args[args.index("--invcf") + 1],
                                               args[args.index("--outpath") + 1])


        elif args[args.index("--mode")+1] == "main" \
                and "--invcf" in args \
                and "--ingff" in args \
                and "--infasta" in args \
                and "--outpath" in args:
            print("Start: NAVIP")
            if not overwriting:
                try:
                    file = open(args[args.index("--outpath") + 1] + "all_transcripts_data.fa", "r")
                    file.close()
                    sys.exit(args[args.index("--outpath") + 1] + "all_transcripts_data.fa is already existing.")

                except:
                    pass
                try:
                    file = open(args[args.index("--outpath") + 1] + "all_transcripts_data_damaged.txt", "r")
                    file.close()
                    sys.exit(args[args.index("--outpath") + 1] + "all_transcripts_data_damaged.txt is already existing.")
                except:
                    pass
                try:
                    file = open(args[args.index("--outpath") + 1] + "All_VCF.vcf", "r")
                    file.close()
                    sys.exit(args[args.index("--outpath") + 1] + "All_VCF.vcf is already existing.")
                except:
                    pass
            try:
                file = open(args[args.index("--invcf")+1],'r')
                file.close()
                file = open(args[args.index("--ingff") + 1], 'r')
                file.close()
                file = open(args[args.index("--infasta") + 1], 'r')
                file.close()
                file = open(args[args.index("--outpath")+1] + "try-file-for-exceptions", "w")
                file.write("test")
                file.close()
                os.remove(args[args.index("--outpath")+1] + "try-file-for-exceptions")
            except FileNotFoundError:
                e = sys.exc_info()[0]
                sys.exit(e)
            except PermissionError:
                e = sys.exc_info()[0]
                sys.exit(e)
            except :
                e = sys.exc_info()[0]
                sys.exit(e)

            Coordinator.navip_main_coordinator(args[args.index("--invcf") + 1],
                                               args[args.index("--ingff") + 1],
                                               args[args.index("--infasta") + 1],
                                               args[args.index("--outpath") + 1])


        elif args[args.index("--mode")+1] == "sfa" \
                and "--innavipvcf" in args\
                and "--innavipfasta" in args\
                and "--outpath" in args:
            print("Start: simple first analysis")
            if not overwriting:
                try:
                    file = open(args[args.index("--outpath") + 1] + "no_NONE.vcf", "r")
                    file.close()
                    sys.exit("There already exist some output-data in: " + args[args.index("--outpath") + 1])

                except:
                    pass
            try:
                file = open(args[args.index("--innavipvcf")+1],'r')
                file.close()
                file = open(args[args.index("--innavipfasta") + 1], 'r')
                file.close()
                file = open(args[args.index("--outpath")+1] + "try-file-for-exceptions", "w")
                file.write("test")
                file.close()
                os.remove(args[args.index("--outpath")+1] + "try-file-for-exceptions")
            except FileNotFoundError:
                e = sys.exc_info()[0]
                sys.exit(e)
            except PermissionError:
                e = sys.exc_info()[0]
                sys.exit(e)
            except :
                e = sys.exc_info()[0]
                sys.exit(e)
            #doku for all files, new names, too (better files, too)
            sfa.sfa_main(args[args.index("--innavipvcf") + 1],
                        args[args.index("--innavipfasta") + 1],
                        args[args.index("--outpath") + 1])
        elif args[args.index("--mode")+1] == "vcfc" \
                and "--invcf" in args\
                and "--outpath" in args:
            VCF_Format_Check.VCF_Check(str(args[args.index("--invcf") + 1]),str(args[args.index("--outpath") + 1]))

        else:
            print("Arguments are invalid."
                  "\nPlease look into the Readme.txt.")
    else:
        print("Arguments are invalid."
              "\nPlease look into the Readme.txt.")
    #print(sys.argv)

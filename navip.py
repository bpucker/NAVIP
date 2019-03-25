__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"

import sys
import os
from time import sleep
import VCF_preprocess,Coordinator, VCF_Format_Check
try:
    import sfa2
except ModuleNotFoundError:
    print("ModuleNotFoundError")
    print(str(sys.exc_info()))
except ImportError:
    print("ImportError")
    print(str(sys.exc_info()))




if __name__ == '__main__':
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
                try:
                    file = open(args[args.index("--outpath") + 1] + "try-file-for-exceptions", "r")
                except FileNotFoundError:
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
                and "--outpath" in args:
            print("Start: simple first analysis")
            try:
                file = open(args[args.index("--innavipvcf")+1],'r')
                file.close()
                try:
                    file = open(args[args.index("--outpath") + 1] + "try-file-for-exceptions", "r")
                except FileNotFoundError:
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
            every_cindel_compensation = True
            if "--ecc" in args:
                every_cindel_compensation: False
            picture_formats = "pdf"
            if "--format" in args:
                picture_formats = args[args.index("--format")+1]


            if not overwriting:
                mode_name = ""
                if every_cindel_compensation:
                    mode_name = "cInDels"
                else:
                    mode_name = "only_additive_cInDels"

                all_the_output_file_names = []

                table_outputname1 = "transcripts_isoform_" + mode_name + '.txt'
                table_outputname2 = "transcripts_isoform_" + mode_name + '_TIDs' + '.txt'
                table_outputname3 = "transcripts_unique_" + mode_name + '.txt'
                table_outputname4 = "transcripts_unique_" + mode_name + '_TIDs' + '.txt'
                table_outputname5 = "transcripts_unique_" + mode_name + '_TIDs_detailed' + '.txt'
                table_outputname6 = "transcripts_unique_" + mode_name + '_TIDs_detailed_sorted_by_TID' + '.txt'

                all_the_output_file_names.append(table_outputname1)
                all_the_output_file_names.append(table_outputname2)
                all_the_output_file_names.append(table_outputname3)
                all_the_output_file_names.append(table_outputname4)
                all_the_output_file_names.append(table_outputname5)
                all_the_output_file_names.append(table_outputname6)

                outputfolder = args[args.index("--outpath") + 1]

                probably_supported_formats = ["png", "pdf", "ps", "eps", "svg"]
                for orig_outputname in [table_outputname1,table_outputname3]:
                    for picture_format in picture_formats.split(','):
                        if picture_format.lower() in probably_supported_formats:
                            try:
                                open(outputfolder + orig_outputname.split(".")[0] + "." + picture_format.lower(), 'r')
                                sys.exit("File is existing and overwriting is not active:" + outputfolder + orig_outputname.split(".")[0] + "." + picture_format.lower())
                            except FileNotFoundError:
                                continue
                for file_name in all_the_output_file_names:
                    try:
                        open(outputfolder + file_name, 'r')
                        sys.exit("File is existing and overwriting is not active:" + outputfolder + file_name)
                    except FileNotFoundError:
                        continue

            #doku for all files, new names, too (better files, too)
            # vcf_file_link:str, mod_or_not: bool, outputfolder: str, max_bp_range: int
            sfa2.sfa2_main(args[args.index("--innavipvcf") + 1],
                        every_cindel_compensation,
                        args[args.index("--outpath") + 1],
                           picture_formats)
        elif args[args.index("--mode")+1] == "vcfc" \
                and "--invcf" in args\
                and "--outpath" in args:
            VCF_Format_Check.VCF_Check(str(args[args.index("--invcf") + 1]),str(args[args.index("--outpath") + 1]))

        else:
            print("Arguments are invalid."
                  "\nPlease look into the Readme or the wiki for more information.")
    else:
        print("Arguments are invalid."
              "\nPlease look into the Readme or the wiki for more information.")

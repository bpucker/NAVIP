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
                and "--innavipfasta" in args\
                and "--outpath" in args\
                and "--bpr" in args:
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

            #doku for all files, new names, too (better files, too)
            # vcf_file_link:str, mod_or_not: bool, outputfolder: str, max_bp_range: int
            sfa2.sfa2_main(args[args.index("--innavipvcf") + 1],
                        every_cindel_compensation,
                        args[args.index("--outpath") + 1],
                        int(args[args.index("--bpr") + 1]) )
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

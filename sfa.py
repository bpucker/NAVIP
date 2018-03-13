import matplotlib.pyplot as plt
from datetime import datetime

def use_shared_numbers(old_vcf:str,new_vcf:str,sdvc_vcf:str,InDel_vcf:str, AAC_out_path:str, unsync_dvc_out:str):
    """
    This collection of functions in here creating a few files with analysed data:
        The "no_NONE.vcf" VCF file, only with variants with shared effects.
        The "dvc.vcf" VCF file, only with double variant codon information.
        The "dvc_unsync.vcf", only with double variant codons, which had not the same codon in the original dna, but now.
        The "indel_neutralize.vcf", only with variants, whose indel effects neutralize a frameshift.
        The "AAC_Tab.txt", which is an overview for all DVC
    :param old_vcf: NAVIP VCF file with shared keys.
    :param new_vcf: Path + name for a new output file for shared effects variants only
    :param sdvc_vcf: Path + name for a new output file, for DVC within the same codon in both, original dna and new dna
    :param InDel_vcf: Path + name for a new output file, only containing neutralizing frameshift variants.
    :param AAC_out_path: Path + name for a new output table with an overview table for DVC.
    :param unsync_dvc_out: Path + name for a new output file, containing only dvc with different codons in the original dna
    :return: nothing.
    """

    def no_NONE(path_with_vcf_file:str, output_path_and_name: str ):
        """
        Creates the no_None.vcf file.
        :param path_with_vcf_file: Old navip vcf file.
        :param output_path_and_name:  Output + name of the new vcf file.
        :return: Nothing.
        """
        navip_vcf_file = open(path_with_vcf_file, "r")
        line = navip_vcf_file.readline()
        lines = []
        while line:
            if line.startswith("#"):
                line = navip_vcf_file.readline()
                continue
            elif "NONE" in line:
                line = navip_vcf_file.readline()
                continue
            else:
                lines.append(line)
                line = navip_vcf_file.readline()
        navip_vcf_file.close()

        outie = open(output_path_and_name,"w")
        outie.write("".join(lines))
        outie.close()

    def gene_and_pos_dict(path_with_vcf_file:str)->dict:
        """
        Creates a dictionary with a dictionary for every transcript, containing all its variants.
        :param path_with_vcf_file: The new vcf file, only with variants with shared effects.
        :return: A dictionary, containing a dict for every transcript with all its variants.
        """
        #0      1       2   3   4   5       6       7.0         7.1 7.2                                 7.3     7.4   7.5  7.6 7.7 7.8 7.9
        #Chr1	791050	.	CAG	C	2159.86	PASS	AT1G03230.1|FOR|DEL,Frameshift-2,Amino acid change|791053|GCAGCG;1|aa|GCC;1|a|941|NAVIP_END
        dictdict = {}
        vcf_file = open(path_with_vcf_file, "r")
        line = vcf_file.readline()
        while line:
            tid = str(line.split("\t")[7].split("|")[0])
            pos = int(line.split("\t")[1])
            try:
                dictdict[tid][pos] = line
            except KeyError:
                dictdict[tid] = {}
                dictdict[tid][pos] = line
            line = vcf_file.readline()
        vcf_file.close()
        return dictdict

    def single_double_codon_variant(genedict:dict)-> dict:
        """
        Search and find all DVC variants for chromosome_position_TID
        :param genedict: Dictionary containing all variants per transcript.
        :return: Dictionary, containing all DVC.
        """
        # 0      1       2   3   4   5       6       7.0         7.1 7.2                                 7.3     7.4   7.5  7.6 7.7 7.8 7.9
        # Chr1	791050	.	CAG	C	2159.86	PASS	AT1G03230.1|FOR|DEL,Frameshift-2,Amino acid change|791053|GCAGCG;1|aa|GCC;1|a|941|NAVIP_END
        sortme = []
        for gene in genedict:
            sortme.append(gene)
        sortme = sorted(sortme)
        all_sdcv = {}
        for sortetgene in sortme:
            for pos in genedict[sortetgene]:
                line = genedict[sortetgene][pos]
                if len(line.split("\t")[3]) > 1 or len(line.split("\t")[4]) > 1:
                    continue
                cds = int(line.split("\t")[7].split("|")[8])
                chr = str(line.split("\t")[0])
                pos = str(line.split("\t")[1])
                for shared_id in line.split("\t")[7].split("|")[3].split(","):
                    shared_id = int(shared_id)
                    shared_line = genedict[sortetgene][shared_id]
                    if len(shared_line.split("\t")[3]) > 1 or len(shared_line.split("\t")[4]) > 1:
                        continue
                    elif abs(int(shared_line.split("\t")[7].split("|")[8]) - cds) == 2:
                        if chr +"_"+ pos +"_"+ str(shared_line.split("\t")[0]) in all_sdcv:
                            continue
                        elif chr +"_"+ str(shared_line.split("\t")[0]) +"_"+ pos in all_sdcv:
                            continue
                        else:
                            all_sdcv[chr +"_"+ pos +"_"+ str(shared_line.split("\t")[7].split("|")[0])] = (line,shared_line)

        return all_sdcv

    def write_single_double_variant_codon(sdcv_vcf, unique_sdcv_dict):
        """
        Writing the vcf file with only unique DVC (CHR_POSITION unique).
        :param sdcv_vcf: Path + name for a new output file, for DVC within the same codon in both, original dna and new dna
        :param unique_sdcv_dict: Dictionary, only with CHR_POSITION unique DVC.
        :return:  Nothing.
        """
        sortme = []

        for tuple in unique_sdcv_dict:
            sortme.append(unique_sdcv_dict[tuple][0])

        sortme = sorted(sortme, key = lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
        #key=lambda data_line: (int(data_line.split("\t")[1]), str(data_line.split("\t")[7])))
        outie = open(sdcv_vcf,"w")
        outie.write("".join(sortme))
        outie.close()
        print (len(sortme)/2) # 503

    def indel_neutralize(genedict:dict)->dict:
        """
        Uses the dictionary, with every variant per transcript to filter after indels
        and filter again for indels, whose effects neutralizes together.
        :param genedict: Dict containing all transcripts and its variants.
        :return: Dictionary with frameshift neutral indels.
        """
        # 0      1       2   3   4   5       6       7.0         7.1 7.2                                 7.3     7.4   7.5  7.6 7.7 7.8 7.9
        # Chr1	791050	.	CAG	C	2159.86	PASS	AT1G03230.1|FOR|DEL,Frameshift-2,Amino acid change|791053|GCAGCG;1|aa|GCC;1|a|941|NAVIP_END
        sortme = []
        for gene in genedict:
            sortme.append(gene)
        sortme = sorted(sortme)
        f_neutral_indels = {}
        for sortetgene in sortme:
            for pos in genedict[sortetgene]:
                line = genedict[sortetgene][pos]
                if not "Frameshift" in line:
                    continue
                list_of_shifts = [line]
                int_of_shifts = 0
                classifications = line.split("\t")[7].split("|")[2]
                if "+1" in classifications:
                    int_of_shifts += 1
                elif "+2" in classifications:
                    int_of_shifts +=2
                elif "-1" in classifications:
                    int_of_shifts -=1
                elif "-2" in classifications:
                    int_of_shifts -=2
                else:
                    print("error classifications1")
                for shared_id in line.split("\t")[7].split("|")[3].split(","):
                    shared_id = int(shared_id)
                    shared_line = genedict[sortetgene][shared_id]
                    if not "Frameshift" in shared_line:
                        continue
                    list_of_shifts.append(shared_line)
                    classifications = shared_line.split("\t")[7].split("|")[2]
                    if "+1" in classifications:
                        int_of_shifts += 1
                    elif "+2" in classifications:
                        int_of_shifts += 2
                    elif "-1" in classifications:
                        int_of_shifts -= 1
                    elif "-2" in classifications:
                        int_of_shifts -= 2
                    else:
                        print("error classifications1")
                    if int_of_shifts == 0:
                        f_neutral_indels[sortetgene] = list_of_shifts
        return f_neutral_indels

    def write_indel_neutralize_dict(path:str,indel_neutralize_dict:dict):
        """
        Writes the file, containing all indel, whose neutralized frameshifts.
        :param path: Path and name for the outgoing file.
        :param indel_neutralize_dict: Dictionary containing all speciall neutral frameshift variants.
        :return: NOthing.
        """
        outie = []
        for listi in indel_neutralize_dict:
            for line in indel_neutralize_dict[listi]:
                outie.append(line)
            # sortme.append(uniquq_sdcv_dict[tuple][1])
        outie = sorted(outie, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))
        indel_vcf = open(path,"w")
        indel_vcf.write("".join(outie))
        indel_vcf.close()

    def create_AAC_table(sdvc: str) -> dict:
        """
        Creates a dictionary, containing DVC information.
        From which Original AA ->to which New AA
        :param sdvc: VCF File, which contains only DVC with same codon in both, orig DNA and new DNA
        :return: Dictionary contains the numbers, how the DVC mutated
        """
        AAC_Table = {}
        no_indels = open(sdvc, "r")

        line = no_indels.readline()

        AA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
        AA2 = []
        for a in AA:
            AA2.append(a.lower())
        AA = sorted(AA2)
        for a1 in AA:
            for a2 in AA:
                AAC_Table[a1, a2] = 0

        while line:
            orig_aa = line.split("\t")[7].split("|")[5]
            new__aa = line.split("\t")[7].split("|")[7]
            if (orig_aa, new__aa) in AAC_Table:
                AAC_Table[orig_aa, new__aa] += 1
            else:
                AAC_Table[orig_aa, new__aa] = 1

            line = no_indels.readline()
        no_indels.close()

        return AAC_Table

    def dict_to_table(aac_dict: dict, AAC_out_path:str) -> (str, str):
        """
        Writes the AAC table into a text file.
        :param aac_dict: Dictionary, containing from which original AA to which new AA the variants mutate.
        :param AAC_out_path: Path and file name for the output.
        :return: NOthing.
        """
        sortme = []
        AA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
        AA2 = []
        for a in AA:
            AA2.append(a.lower())
        AA = sorted(AA2)

        for entry in aac_dict:
            sortme.append((entry, aac_dict[entry]))
        sortme = sorted(sortme)
        output = []
        out_round_down = []
        for a in AA:
            output.append("\t" + str(a).upper())
            out_round_down.append("\t" + str(a).upper())

        for left in AA:
            output.append("\n" + str(left).upper())
            for above in AA:
                # count = int(aac_dict[left,above]) + int(aac_dict[above,left])
                count = int(aac_dict[left, above])
                # count = int(aac_dict[above,left])
                output.append("\t" + str(count))

        entrys = []
        sumColumn = 0
        sumRow = 0
        for left in AA:
            out_round_down.append("\n" + str(left).upper())
            count2 = 0
            count3 = 0
            for above in AA:
                # count = int(aac_dict[left,above]) + int(aac_dict[above,left])
                count = int(int(aac_dict[left, above])/2)
                out_round_down.append("\t" + str(count))
                count2 += int(int(aac_dict[above, left])/2)
                count3 += int(int(aac_dict[left, above])/2)
            entrys.append("\t" + str(count2))
            out_round_down.append("\t" + str(count3))
            sumColumn += count2
            sumRow += count3
        out_round_down.append("\n" + "".join(entrys))
        out_round_down.append("\nSum column: " + str(sumColumn))
        out_round_down.append("\nSum row: " + str(sumRow))

        #print(aac_dict)
        write_tab = open(AAC_out_path, "w")
        write_tab.write("Left Original AA, Above New AA\n")
        write_tab.write("".join(out_round_down))
        write_tab.close()

    def all_unsync_dvc(sdvc_vcf:str, unsinc_dvc_out: str):
        """
        Searches all DVC, which does not have the same codon in the origin DNA/AA
        :param sdvc_vcf: Path + file name. VCF contains all DVC data.
        :param unsinc_dvc_out: Output VCF, only contains DVC without same codon usage in origin DNA/AA
        :return: Nothing.
        """
        dvc_unsync = open(sdvc_vcf, "r")
        line = dvc_unsync.readline()
        lines = []
        while (line):
            lines.append(line)
            line = dvc_unsync.readline()
        dvc_unsync.close()
        data_to_write = []

        for i,line in enumerate(lines):
            # 0      1       2   3   4   5       6       7.0         7.1 7.2                                 7.3     7.4   7.5  7.6 7.7 7.8 7.9
            # Chr1	791050	.	CAG	C	2159.86	PASS	AT1G03230.1|FOR|DEL,Frameshift-2,Amino acid change|791053|GCAGCG;1|aa|GCC;1|a|941|NAVIP_END
            if i == len(lines)-1:
                continue
            oaa1 = line.split("\t")[7].split("|")[5]
            naa1 = line.split("\t")[7].split("|")[7]
            tid1 = line.split("\t")[7].split("|")[0]

            oaa2 = lines[i+1].split("\t")[7].split("|")[5]
            naa2 = lines[i+1].split("\t")[7].split("|")[7]
            tid2 = lines[i+1].split("\t")[7].split("|")[0]

            if tid1 != tid2:
                continue
            elif oaa1 != oaa2:
                data_to_write.append("#oaa1 != oaa2\n")
                data_to_write.append(line + lines[i+1])
            elif naa1 != naa2:
                data_to_write.append("#naa1 != naa2\n")
                data_to_write.append(line + lines[i + 1])

        outie = open (unsinc_dvc_out,"w")
        outie.write("".join(data_to_write))
        outie.close()

    def unique_dvc(all_tid_sdcv_dict:dict):
        """
        Creates a dictionary, with all DVC <CHR_Position> unique.
        No more multiple entrys TID.1, TID.2, ....
        :param all_tid_sdcv_dict: Dictionary with all DVC
        :return: Dictionary with unique DVC
        """
        first_tid_dict = {}
        sortme = []
        for key in all_tid_sdcv_dict:
            sortme.append((str(key).split("_")[0],str(key).split("_")[1],key,all_tid_sdcv_dict[key]))
        sortme = sorted(sortme, key = lambda triple: (str(triple[0]), int(triple[1]), int(triple[2].split(".")[1])))

        for entry in sortme:
            key = entry[2].split(".")[0]
            first_tid_dict[key] = entry[3]
        return first_tid_dict

    #old_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/2017-10-17 17:17:56.641179_All_VCF_.vcf"
    #new_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/no_NONE.vcf"
    #sdvc_vcf= "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/dcv.vcf"
    #InDel_vcf="/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/indel_neutralize.vcf"

    #create a vcf file only with vcf with shared ids -> only variants with shared effects
    no_NONE(old_vcf,new_vcf)

    # a dict with all (gene,pos) = line
    genedict = gene_and_pos_dict(new_vcf)


    #dict with {chr_pos_tid} = (line1,line2)
    # line1 and lin2 are together one double_variant_codon
    all_tid_sdcv_dict = single_double_codon_variant(genedict)

    #dict with {chr_pos_tid}, but always the first tid.1 when possible
    first_tid_dvc_dict = unique_dvc(all_tid_sdcv_dict)

    #write the lines from the dict to the file
    write_single_double_variant_codon(sdvc_vcf,first_tid_dvc_dict)

    ####for writing the AAC table
    aac_dict = create_AAC_table(sdvc_vcf)
    # because of 7 AAC without same startcodon -> /2 rounded down
    # (important for use with larger dataset, because 0.5 + 0.5 == 1 but its wrong
    dict_to_table(aac_dict,AAC_out_path)

    #show the mentioned data without same origin AA
    all_unsync_dvc(sdvc_vcf,unsync_dvc_out)

    #create a dict with all indels, which neutralize their effects (+2 -2 ==0)
    indel_neutralize_dict = indel_neutralize(genedict)
    #write the indels into a vcf file
    write_indel_neutralize_dict(InDel_vcf, indel_neutralize_dict)

def main_for_compare_AA_length(navip_fasta_file:str, aa_len_dif_txt:str):
    """
    Managing all data/functions for a aminoacis comparing from
    original CDS (start to stop) to the new AA CDS (start to stop)
    :param navip_fasta_file: Navtip fasta file, including all CDS sequences, with the types aDNA, aAA, nDNA, nAA
    :param aa_len_dif_txt: Path + file name for the new difference file
    :return: Nothing.
    """

    def compare_cds_length_to_stop(navip_fasta_file:str) -> (dict, list):
        """
        Reads the navip fasta file, which contains all CDS sequences in four categories:
            oDNA = originial DNA,
            oAA = original AA,
            nDNA = new DNA,
            nAA = new AA
        and returns a dictionary with TID_type
        and a list of all TIDs (=transcript ID)
        :param navip_fasta_file: NAVIP Fasta file, containing all CDS data.
        :return: dictionary with all sequences, available with TID_TYPE and a list, with all TIDs
        """
        my_fasta_file = open(navip_fasta_file,"r")
        line = my_fasta_file.readline()
        tid_type_dict = {}
        tid_list = []
        while line:
            if line.startswith("#"):
                line = my_fasta_file.readline()
                continue

            if line.startswith(">"):
                tid = line.split("|")[0][1:]
                if tid not in tid_list:
                    tid_list.append(tid)
                type = ""
                if "oDNA" in line:
                    type = "oDNA"
                elif "oAA" in line:
                    type = "oAA"
                elif "nDNA" in line:
                    type = "nDNA"
                elif "nAA" in line:
                    type = "nAA"
                else:
                    print("error" + str(tid))
                    break
                line = my_fasta_file.readline()
                if type == "oDNA" or type == "nDNA":
                    tid_type_dict[tid, type] = len(line)
                else:
                    tid_type_dict[tid, type] = first_stop = line.find("*")
                line = my_fasta_file.readline()
                # print (line)
        tid_list = sorted(tid_list)
        #sortme = sorted(sortme, key=lambda data_line: (str(data_line.split("\t")[0]), int(data_line.split("\t")[1])))

        my_fasta_file.close()
        return (tid_type_dict, tid_list)

    def for_dna_length(tid_type_dict:dict, tid_list:list):
        """
        For comparing the length of the DNA CDS sequence.
        May be used again in the future.
        :param tid_type_dict: Dictionary with tid_type, contains all sequences
        :param tid_list: List with all tids.
        :return: NOthing.
        """

        dif_dict = {}

        for tid in tid_list:
            try:
                old_minus_new = int(tid_type_dict[tid, "oDNA"]) - int(tid_type_dict[tid, "nDNA"])
                try:
                    dif_dict[old_minus_new] = (int(dif_dict[old_minus_new][0]) + 1, dif_dict[old_minus_new][1])
                    dif_dict[old_minus_new][1].append(tid)
                except KeyError:
                    dif_dict[old_minus_new] = (1, [tid])

                '''
                if abs(int(dicti[tid,"oDNA"]) - int(dicti[tid,"nDNA"])) >= 10:
                    print (str(tid) + "\t" + str(dicti[tid,"oDNA"]) + "\t"+ str(dicti[tid,"oAA"])
                           + "\t"+ str(dicti[tid,"nDNA"]) + "\t"+ str(dicti[tid,"nAA"]) + "\n")
                '''
            except KeyError:
                pass
        output = ["Differences between original DNA length and new DNA length"]
        sortme = dif_dict.keys()
        sortme = sorted(sortme)
        for a in sortme:
            output.append((str(a) + "\t" + str(dif_dict[a][0])) + "\n")
            if a > 0:
                print(str(a) + "\t" + str(dif_dict[a][0]))
                if len(dif_dict[a][1]) > 100:
                    continue
                print(str(dif_dict[a][1]))

    def for_aa_stop_length(tid_type_dict:dict, tid_list:list) -> list:
        """
        Calculates the output string for the aa_len_dif_txt file.
        :param tid_type_dict: Dictionary with tid_type, contains all sequences
        :param tid_list: List with all tids.
        :return: Output string for the aa_len_dif_txt file.
        """
        dif_dict = {}
        for tid in tid_list:
            try:
                new_minus_old = - int(tid_type_dict[tid, "oAA"]) + int(tid_type_dict[tid, "nAA"])
                try:
                    dif_dict[new_minus_old] = (int(dif_dict[new_minus_old][0]) + 1, dif_dict[new_minus_old][1])
                    dif_dict[new_minus_old][1].append(tid)
                except KeyError:
                    dif_dict[new_minus_old] = (1, [tid])
            except KeyError:
                pass
        output = ["#Differences between original AA length and new AA length\n"]
        sortme = dif_dict.keys()
        sortme = sorted(sortme)
        for a in sortme:
            output.append((">" + str(a) + "\t" + str(dif_dict[a][0])) + "\n")
            i = 1
            for tid in dif_dict[a][1]:
                output.append(tid + " ")
                if i % 10 == 0:
                    output.append("\n")
                i += 1
            if i % 10 != 1:
                output.append("\n")
            if a > 0:
                if len(dif_dict[a][1]) > 100:
                    continue
        return output

    tid_type_dict, tid_list = compare_cds_length_to_stop(navip_fasta_file)
    length_dif = open(aa_len_dif_txt, "w")
    length_dif.write("".join(for_aa_stop_length(tid_type_dict, tid_list)))
    length_dif.close()

def plots(aa_len_dif:str, outpath: str):
    """
    Reading the difference file from the aminoacids to create four plots.
    But: Most transcripts will have no differences, so there is a probling showing all of them in one plot.
    Every plot shows only a part of this data.
    :param aa_len_dif: Path (including file name) of the length difference data file.
    :param outpath: Path to the existing output folder
    :return: Nothing.
    """
    aa_len_dif_file = open(aa_len_dif, "r")
    line = aa_len_dif_file.readline()
    dif_count_list =[]
    while line :
        if line.startswith(">"):
            dif_count_list.append((int(line.split("\t")[0][1:]),int(line.split("\t")[1])))
        line = aa_len_dif_file.readline()

    aa_len_dif_file.close()
    y = []
    x = []
    #tuple example: (-5,21)
    for tuple in dif_count_list:
        #if tuple[1] < 10 or tuple[1] > 100:
        #    continue
        if -20 <= tuple[0] <= 20 and tuple[0] != 0:
            y.append(tuple[1])
            x.append(tuple[0])




    plt.bar(x,y,1,color="blue")
    plt.title("old length - new length")
    plt.gcf()
    #
    plt.xlabel("Number of Amino Acids")
    plt.ylabel("Number of transcripts")
    plt.savefig(outpath + "20_around_zero.png")
    plt.clf()
    y = []
    x = []
    # tuple example: (-5,21)
    for tuple in dif_count_list:
        # if tuple[1] < 10 or tuple[1] > 100:
        #    continue
        if -20 <= tuple[0] <= 20:
            continue
        y.append(tuple[1])
        x.append(tuple[0])

    plt.bar(x, y, 1, color="blue")
    plt.title("old length - new length")
    plt.gcf()
    #
    plt.xlabel("Number of Amino Acids")
    plt.ylabel("Number of transcripts")
    plt.savefig(outpath + "not_20_around_zero.png")
    plt.clf()


    y = []
    x = []
    # tuple example: (-5,21)
    for tuple in dif_count_list:
        # if tuple[1] < 10 or tuple[1] > 100:
        #    continue
        if  tuple[0] ==0:
            continue
        y.append(tuple[1])
        x.append(tuple[0])

    plt.bar(x, y, 1, color="blue")
    plt.title("old length - new length")
    plt.gcf()
    #
    plt.xlabel("Number of Amino Acids")
    plt.ylabel("Number of transcripts")
    plt.savefig(outpath + "not_zero.png")
    plt.clf()


    y = []
    x = []
    # tuple example: (-5,21)
    for tuple in dif_count_list:
        y.append(tuple[1])
        x.append(tuple[0])

    plt.bar(x, y, 1, color="blue")
    plt.title("old length - new length")
    plt.ylim(0,20)
    #
    plt.xlabel("Number of Amino Acids")
    plt.ylabel("Number of transcripts")
    plt.savefig(outpath + "all.png")
    plt.clf()

def count_stuff_unique(path_to_neutralized_indels:str, no_none_vcf:str, outpath:str):
    """
    (new) Testing area.
    Contains an output for an overview file about:
        Unique,
        AAC = amino acid changes
        StopL = stop lost
        StopG = stop gained
        FraShi = frameshifts
    :param path_to_neutralized_indels: path to the file, which contains only neutralized indels.
    :param no_none_vcf: path to all vcf data, which have effects with other variants
    :return: Nothing.
    """

    def create_dict_and_gap(path_to_neutralized_indels:str)->(dict,list):
        """
        Testing area. Not used for anything now.
        :param path_to_neutralized_indels:
        :return:
        """
        indels = open(path_to_neutralized_indels)
        line = indels.readline()
        tid_indels = {} # {tid} = (count,[])
        uniques = []
        # 0      1       2   3   4   5       6       7.0         7.1 7.2                                 7.3     7.4   7.5  7.6 7.7 7.8 7.9
        # Chr1	791050	.	CAG	C	2159.86	PASS	AT1G03230.1|FOR|DEL,Frameshift-2,Amino acid change|791053|GCAGCG;1|aa|GCC;1|a|941|NAVIP_END
        while line:
            splitedline = line.split("\t")
            if str(splitedline[0]) +"_"+ str(splitedline[1]) +"_"+str(splitedline[7].split("|")[0].split(".")[0]) not in uniques:
                uniques.append(str(splitedline[0]) +"_"+ str(splitedline[1]) +"_"+str(splitedline[7].split("|")[0].split(".")[0]))
            try:
                meh = tid_indels[line.split("\t")[7].split("|")[0]][1]
                meh.append(line)
                tid_indels[line.split("\t")[7].split("|")[0]] = ((tid_indels[line.split("\t")[7].split("|")[0]][0] +1,meh))
                #print(tid_indels[line.split("\t")[7].split("|")[0]][1])
            except KeyError:
                #print(line)
                meh = []
                meh.append(line)
                tid_indels[line.split("\t")[7].split("|")[0]] = (1,meh)
            line = indels.readline()
        indels.close()
        #print("test")
        #count = 0
        #for entry in tid_indels:
        #    print(tid_indels[entry])
        #    count += int(tid_indels[entry][0])
        #print (count)
        #print("uniques:" + str(len(uniques)))

        return (tid_indels,uniques)

    def count_Unique_entry_stop_stoplost__framesthift(new_vcf:str)->(list,list,list,list,list,list):
        """
        Testing area, but still :
            Reads the entire navip vcf, containing all variants with shared keys
            Returns lists with unique entrys.
        :param new_vcf:Navip vcf, containing all variants with shared keys
        :return: unique(only one chr_pos entry),unique_aa_change,unique_stop_lost,unique_stop_gained,unique_frameshift
        """
        unique = []
        unique_aa_change = []
        unique_stop_lost = []
        unique_stop_gained = []
        unique_frameshift = []
        dvc = []

        no_none_vcf = open(new_vcf)
        line = no_none_vcf.readline()
        lines = []
        while line:

            if line.startswith("#"):
                line = no_none_vcf.readline()
                continue
            splitted = line.split("\t")
            id = str(splitted[0]) + "_" + str(splitted[1]) + "_" + str(splitted[7].split("|")[0].split(".")[0])

            if id not in unique:
                unique.append(id)
                lines.append(line)

            if "Frameshift" in line and id not in unique_frameshift:
                unique_frameshift.append(id)
            if "Amino acid change" in line and id not in unique_aa_change:
                unique_aa_change.append(id)
            if "Stop gained" in line and id not in unique_stop_gained:
                unique_stop_gained.append(id)
            if "Stop lost" in line and id not in unique_stop_lost:
                unique_stop_lost.append(id)

            line = no_none_vcf.readline()
        """
        for i,line in enumerate(lines):
            if i == len(line)-1:
                continue
            elif "SUB" not in line:
                continue

            print(line.split("\t")[7].split("|")[8])
            print(lines[i+1].split("\t")[7].split("|"))

            if abs(int(line.split("\t")[7].split("|")[8]) - int(lines[i+1].split("\t")[7].split("|")[8])) == 2:
                if 1 == int(line.split("\t")[7].split("|")[6].split(";")[1]) \
                        or 1 == int(lines[i+1].split("\t")[7].split("|")[6].split(";")[1]):
                    continue
                dvc.append(str(splitted[0]) + "_" + str(splitted[1]) + "_" + str(splitted[7].split("|")[0].split(".")[0]))
        """
        # 0      1       2   3   4   5       6       7.0         7.1 7.2                                 7.3     7.4   7.5  7.6 7.7 7.8 7.9
        # Chr1	791050	.	CAG	C	2159.86	PASS	AT1G03230.1|FOR|DEL,Frameshift-2,Amino acid change|791053|GCAGCG;1|aa|GCC;1|a|941|NAVIP_END
        no_none_vcf.close()
        #print ("Amino acid change: " + str(len(unique_aa_change)))
        #print ("Frameshift: " + str(len(unique_framesthift)))
        #print ("Stop gained: " + str(len(unique_stop_gained)))
        #print ("Stop lost: " + str(len(unique_stop_lost)))
        return(unique,unique_aa_change,unique_stop_lost,unique_stop_gained,unique_frameshift)#,dvc)

    def handle_lists_with_sets(list1:list, remove_from_list_1:list):
        """
        Working with sets and set operations.
        :param list1:
        :param remove_from_list_1:
        :return:
        """
        erg = set(list1) & set(remove_from_list_1)
        return erg



    #tid_indels,uniques = create_dict_and_gap(path_to_neutralized_indels)
    uniquelists = count_Unique_entry_stop_stoplost__framesthift(no_none_vcf)
    description = ["Unique","_AAC","_StopL","_StopG","_FraShi"]

    #print("\t\t"+"\t".join(description))
    outi = "\t" + "\t".join(description)
    for i,x in enumerate(uniquelists):
        outi += "\n"+str(description[i]) + "\t"
        for y in uniquelists:
            outi +=str(len(handle_lists_with_sets(x,y)))+ "\t"
            #if len(handle_lists_with_sets(x,y)) <1000:
            #   outi += "\t"
        #print(outi)

    overview = open(outpath + "overview.txt", "w")
    overview.write(outi)
    overview.close()



def sfa_main(innavipvcf:str, innavipfasta:str, outpath:str):
    """
    Main function for controlling this module.
    Handles the in-going and out-going files and starts all other functions.
    :param innavipvcf:  NAVIP created vcf file. Path + file name.
    :param innavipfasta: NAVIP created fasta file. Path + file name.
    :param outpath: Path to existing output folder.
    :return: Nothing.
    """

    def vcf_stuff():
        """
        Organize the use_shared_numbers() input:
            old_navip_vcf, new_vcf, sdvc_vcf, InDel_vcf, AAC_out_path, unsync_dvc_out
        :return: Nothing.
        """

        # my output vcf (NAVIP)
        #old_navip_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/2017-10-17 17:17:56.641179_All_VCF_.vcf"
        old_navip_vcf = innavipvcf

        # for output of only variants with shared effects
        # new_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/no_NONE.vcf"
        new_vcf = outpath + "no_NONE.vcf"

        # DVC will be sorted after the first TID with a Double-Variant-inside-one-Codon-effect
        # including unsync DVC
        # sdvc_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/dvc.vcf"
        sdvc_vcf = outpath + "dvc.vcf"

        # all indels, which neutralize their effects (+2 -1 -1 == 0)
        # InDel_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/indel_neutralize.vcf"
        InDel_vcf = outpath + "indel_neutralize.vcf"

        # table for AAC overview
        # AAC_out_path = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/AAC_Tab.txt"
        AAC_out_path = outpath + "AAC_Tab.txt"

        # DVC wihtout same origin tripplet (due to frameshifts)
        # unsync_dvc_out = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/dvc_unsync.vcf"
        unsync_dvc_out = outpath + "dvc_unsync.vcf"

        #does all the funny stuff for the 5 output files
        use_shared_numbers(old_navip_vcf, new_vcf, sdvc_vcf, InDel_vcf, AAC_out_path, unsync_dvc_out)

    def aa_len_stuff ():
        """
        Organizing the inout for main_for_compare_AA_length():
            navip_fasta_file, aa_len_dif
        :return: Nothing.
        """

        #name and path to the new datafile
        # > len_dif counted_length
        # all tids, max 10 per row, many rows
        #aa_len_dif = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/AA_length_dif.txt"
        aa_len_dif = outpath + "AA_length_dif.txt"

        #from NAVIP created fasta file with 4x cds from every ttransctipt (original and new DNA/AA)
        #navip_fasta_file = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/2017-10-09 16:32:31.872025_all_transcripts_data.fa"
        navip_fasta_file = innavipfasta

        # organize the stuff for comparing length of AA)
        main_for_compare_AA_length(navip_fasta_file, aa_len_dif)

    def for_plots():
        """
        Starting the plots.
        :return: Nothing.
        """

        # aa_len_dif = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/AA_length_dif.txt"
        aa_len_dif = outpath + "AA_length_dif.txt"

        plots(aa_len_dif,outpath)

    def count_stuff():
        """
        Managing the counting function.
        :return: Nothing.
        """

        # all indels, which neutralize their effects (+2 -1 -1 == 0)
        # InDel_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/indel_neutralize.vcf"
        InDel_vcf = outpath + "indel_neutralize.vcf"

        # for output of only variants with shared effects
        # new_vcf = "/prj/gf-arabseq/project_VariantAnnotation/data/CollisionData/no_NONE.vcf"
        no_none_vcf = outpath + "no_NONE.vcf"

        count_stuff_unique(InDel_vcf, no_none_vcf, outpath)


    time = datetime.now()
    print("##################################################################################")
    print("VCF work start.")
    vcf_stuff()
    print("VCF work is done:" + str(datetime.now() - time))
    print("AA work start.")
    aa_len_stuff ()
    print("AA work is done:" + str(datetime.now() - time))
    print("plots starting")
    for_plots()
    print("plots are done:" + str(datetime.now() - time))
    print("Count things.")
    count_stuff()
    print("Everything done:" + str(datetime.now() - time))
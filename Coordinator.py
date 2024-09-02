__author__ = "Jan-Simon Baasner"
__email__ = "janbaas@cebitec.uni-bielefeld.de"

import copy
from GenomeHandler import Genomehandler, Fasta_Enum, SequenceHandlingError
from VcfHandler import VcfHandler
from VCF_Variant import Variant, VariantEnum
from Transcript import Transcript, TranscriptEnum, ForTypeSafetyAndStatics
from Gff3_Handler import Gff3HandlerV3
from datetime import datetime
from LogOrganizer import LogOrganizer, LogEnums
from snpeffhgvsconverter import SnpeffHgvsConverter
import os
import resource
import gc


def navip_main_coordinator(invcf, ingff, infasta, outpath):
    """
    Coordinates all functions and classes of the main module.
    Also contains a lot of functions, but naming this module "main" or "navip" isn't
    an option, because it would be confusing with the different modules.
    :param invcf: Path and name to a vcf file. Preferred the preprocessed file.
    :param ingff: Path and name to a gff3 file.
    :param infasta: Path and name to a fasta file.
    :param outpath: Path to and including the output folder.
    :return: Nothing.
    """

    def write_all_vcf_file(data_path: str,
                           gff3: Gff3HandlerV3,
                           old_info: bool,
                           ref_codons: bool,
                           ref_aa: bool,
                           alt_codons: bool,
                           alt_aa: bool,
                           extend_file: bool,
                           chr_name: str):
        """
        In this function the VCF-File with the additional neighbourhood awareness information will be written.
        Furthermore, it is possible (and standard) to write a few more Information into the NAVIP info:
        TranscriptID, strand direction, classifications, keys with shared effect, ref and alt codon
        with the variant position inside this codon. This information will be written into every new VCF file, too.
        Transcripts without variants inside the CDS, without CDS and those, which are (possible) damaged within the
        splicing structure will be skipped.

        :param data_path: Path into the existing output folder.
        :param gff3: The gff3 data structure with all transcripts and their variants/sequences.
        :param old_info: The old VCF info column will be placed after the NAVIP info, if this is true.
        :param ref_codons: The REF DNA codons will be written, if true.
        :param ref_aa: The REF AA will be written, if true.
        :param alt_codons: The ALT codons will be written, if true.
        :param alt_aa: The ALT AA will be written, if true.
        :return: Nothing.
        """
        # chrom	Pos	ID	Ref	Alt Info
        if not extend_file:
            vcf_file = open(data_path + "All_VCF" + ".vcf", "w")
            vcf_file.write("##Chrom\tPos\tID\tRef\tAlt\tQual\tFilter\tInfo\n")
            vcf_file.write("##NAVIP: All Data Output\n")
            vcf_file.write("##Please note, that the Variant_Position_in_Codon is read from left to right in forward "
                           "and rigth to left in reverse strand direction.\n")
            vcf_file.write("##If there are no shared effect keys, the value is:\"NONE\".\n")
            vcf_file.write("##Info=<ID=NAV1, Type=String,Number=.,Values=[TranscriptID|"
                           "Strand_Direction|"
                           "Variant_Classification1,Variant_Classification2,...|"
                           "Shared_EffKey(s)|"
                           "REF_Codon(s)/Variant_Position_in_Codon|"
                           "REF_AA|"
                           "old_CDS_Position|"
                           "ALT_Codon(s)/Variant_Position_in_Codon|"
                           "ALT_AA|"
                           "new_CDS_Position]>\n")

            vcf_file.write("##Info=<ID=NAV2, Type=String,Number=.,Values=[ALT|ANNOTATIONS|ANNOTATION_IMPACT|GENE_NAME|"
                           "GENE_ID|FEATURE_TYPE|FEATURE_ID|TRANSCRIPT_BIOTYPE|RANK|HGVS_C|HGVS_P|cDNA_pos/cDNA_length|"
                           "CDS_pos/CDS_DNA_length|AA_pos/AA_length|distance|errors_warnings]>\n")
        else:
            vcf_file = open(data_path + "All_VCF" + ".vcf", "a")
        data_to_write = []
        info_line_parser = SnpeffHgvsConverter
        # for name in gff3.GetChromosomeNames():
        name = chr_name
        for transcript_here in gff3.get_chr_transcripts_dict(name).values():
            transcript_here = ForTypeSafetyAndStatics.transcript_type_safety(transcript_here)

            if not transcript_here.cds_exist or transcript_here.Transcript_CDS_damaged:
                # no CDS or transcript exons are damaged
                continue
            if transcript_here.integrated_variant_objects_cds_hits == []:
                # no variants inside the cds
                continue

            if transcript_here.ForwardDirection == TranscriptEnum.FORWARD:
                direction = "FOR"
            elif transcript_here.ForwardDirection == TranscriptEnum.REVERSE:
                direction = "REV"
            else:
                # will never happen
                LogOrganizer.add_to_log(LogEnums.COORDINATOR_BUGHUNTING_LOG,
                                      "Error:NOT Forward, NOT Reverse:" + str(transcript_here.TID) + "\n")
                direction = "Error:NOT Forward, NOT Reverse."

            for vinfo in transcript_here.integrated_variant_objects_cds_hits:
                vinfo = ForTypeSafetyAndStatics.variant_information_storage_type_safety(vinfo)
                if vinfo.ChrPosition == 5921700 and transcript_here.TID == "Ma09_t08910.1":
                    print("bugsearch")
                if vinfo.ChrPosition == 3917462 and transcript_here.TID == "Ma00_t01100.1":
                    print("bugsearch")
                try:
                    snpeff_string = info_line_parser.convert_main(info_line_parser, transcript_here, vinfo)
                except Exception:
                    snpeff_string = ""
                    LogOrganizer.add_to_log(LogEnums.COORDINATOR_VARIANT_LOG,
                                            str(vinfo.ChrPosition) + "\t" + str(current_transcript.TID))
                    print('Critical error with variant: ' + str(vinfo.ChrPosition) + '\t' + str(transcript_here.TID))
                classificationstring = ""
                class_list = vinfo.Classification
                class_list_length = len(class_list)
                for i in range(0, class_list_length):
                    classificationstring += class_list[i].value
                    if class_list[i] == TranscriptEnum.STOP_CAUSED_IN:
                        classificationstring += str(vinfo.STOP_CAUSED_IN)
                    if i != (class_list_length - 1):  # no tab after last entry
                        classificationstring += ","
                    elif i == (class_list_length - 1):
                        classificationstring += "|"
                    else:
                        LogOrganizer.add_to_log(LogEnums.COORDINATOR_BUGHUNTING_LOG,
                                              "Bug Write_All_VCF_File: classification list length:" + str(
                                                  transcript_here.TID) + "\n")

                shared_effect_list = ""
                if len(vinfo.SharedEffectsWith) > 0:
                    shared_effect_list = ""  # "\tShared Effects with:"
                    shared_effect_list_len = len(vinfo.SharedEffectsWith)
                    # creates the string with all keys
                    for i in range(0, shared_effect_list_len):
                        vinfo_effects = vinfo.SharedEffectsWith[i]
                        vinfo_effects = ForTypeSafetyAndStatics.variant_information_storage_type_safety(
                            vinfo_effects)
                        shared_effect_list += str(vinfo_effects.ChrPosition)
                        if i != (shared_effect_list_len - 1):
                            shared_effect_list += ","
                        # seperate entrys
                        elif i == shared_effect_list_len - 1:
                            shared_effect_list += "|"
                        # after last entry
                        else:
                            LogOrganizer.add_to_log(LogEnums.COORDINATOR_BUGHUNTING_LOG,
                                                  "Bug Write_All_VCF_File: shared_effect_list:" + str(
                                                      transcript_here.TID) + "\n")

                if len(vinfo.OrigTriplets) % 3 != 0:
                    bug = "OrigTriplets: " + str(vinfo.OrigTriplets) \
                          + "\tFrom: " \
                          + transcript_here.TID \
                          + ": " \
                          + str(vinfo.ChrPosition) \
                          + " " \
                          + classificationstring
                    LogOrganizer.add_to_log(LogEnums.COORDINATOR_BUGHUNTING_LOG, bug + "\n")
                elif len(vinfo.ChangedTriplets) % 3 != 0:
                    bug = "ChangedTriplets: " + str(vinfo.ChangedTriplets) \
                          + "\tFrom: " \
                          + transcript_here.TID \
                          + ": " \
                          + str(vinfo.ChrPosition) \
                          + " " \
                          + classificationstring
                    LogOrganizer.add_to_log(LogEnums.COORDINATOR_BUGHUNTING_LOG, bug + "\n")
                # Info:TranscriptID|"
                #   "Strand_Direction|"
                #   "Variant_Classification1,Variant_Classification2,...|"
                #   "Shared_EffKey(s)|"
                #   "REF_Codon(s);Variant_Position_in_Codon|"
                #   "REF_AA|"
                #   "ALT_Codon(s);Variant_Position_in_Codon|"
                #   "ALT_AA|"
                #   "NAVIP_END|"
                #   "<old info field>"
                navip_info = ""
                navip_info_list = ["NAV1=" + str(transcript_here.TID) + "|", str(direction) + "|", classificationstring]
                if shared_effect_list == "":
                    navip_info_list.append("NONE|")
                else:
                    navip_info_list.append(str(shared_effect_list))
                if ref_codons:
                    navip_info_list.append(str(vinfo.OrigTriplets) + "/" + str(vinfo.OrigRaster) + "|")
                if ref_aa:
                    navip_info_list.append(str(vinfo.OrigAmino) + "|")
                navip_info_list.append(str(vinfo.Unchanged_CDS_Position) + "|")
                if alt_codons:
                    navip_info_list.append(str(vinfo.ChangedTriplets) + "/" + str(vinfo.Changed_Raster) + "|")
                if alt_aa:
                    navip_info_list.append(str(vinfo.NewAmino) + "|")
                navip_info_list.append(str(vinfo.Changed_CDS_Position) + ";")
                if old_info:
                    if "\n" in vinfo.OLD_Info:
                        navip_info_list.append(str(vinfo.OLD_Info))
                    else:
                        navip_info_list.append(str(vinfo.OLD_Info) + "\n")

                navip_info = "".join(navip_info_list)
                navip_info_list = []  # clear

                data_to_write.append(
                    str(name) + "\t"
                    + str(vinfo.ChrPosition) + "\t"
                    + str(vinfo.ID) + "\t"
                    + str(vinfo.Ref) + "\t"
                    + str(vinfo.Alt) + "\t"
                    + str(vinfo.Qual) + "\t"
                    + str(vinfo.Filter) + "\t"
                    + snpeff_string
                    + str(navip_info))
        data_to_write = sorted(data_to_write,
                               key=lambda data_line: (int(data_line.split("\t")[1]), str(data_line.split("\t")[7])))
        vcf_file.write("".join(data_to_write))
        data_to_write = []
        vcf_file.close()

    def create_sequences(gff3: Gff3HandlerV3, orig_aa: bool, new_aa: bool, genetic_code: dict, chr_name: str):
        """
        This function calculates the old and new AA CDS to be sure all transcripts are having the AA sequences.
        There are a few possibilities to this point, which can lead to skip a few or all transcripts AA sequences.
        Transcripts without CDS are excluded. Transcripts with a damaged CDS do not get the new AA sequence.
        :param gff3: The gff3 data structure with all transcripts and their variants/sequences.
        :param orig_aa: Will create the original AA sequence, if true.
        :param new_aa: Will create the new AA sequence, if true.
        :param genetic_code: The dictionary with the genetic code.
        :return: Nothing.
        """
        for current_transcript in gff3.get_chr_transcripts_dict(chr_name).values():
            try:
                current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(current_transcript)

                if not current_transcript.cds_exist:
                    continue
                if orig_aa:
                    current_transcript.create_iv_original_translation(genetic_code)

                if current_transcript.Transcript_CDS_damaged:
                    continue
                if new_aa:
                    current_transcript.create_iv_changed_translation(genetic_code)
            except Exception:
                LogOrganizer.add_to_log(LogEnums.COORDINATOR_TRANSCRIPT_LOG, "\t" + str(current_transcript.TID))

    def complete_check(gff3: Gff3HandlerV3, ghandler: Genomehandler, chr_name: str):
        """
        Its always a good idea to test the data in the end.
        Will print warnings, if there is incorrect data.
        :param gff3: The gff3 data structure with all transcripts and their variants/sequences.
        :param ghandler: The reference genome sequences for every chromosome is inside.
        :return: Nothing.
        """
        # for chrName in Shared_Chromosomes_FA_GFF3:
        chr_transcript_list = gff3.get_chr_transcripts_dict(chr_name).values()
        for current_transcript in chr_transcript_list:
            current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(current_transcript)
            if not current_transcript.cds_exist or current_transcript.Transcript_CDS_damaged:
                continue

            for vinfo in current_transcript.integrated_variant_objects_cds_hits:
                vinfo = ForTypeSafetyAndStatics.variant_information_storage_type_safety(vinfo)
                try:
                    if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
                        # Variant <-> Fasta Check
                        if vinfo.Ref != ghandler.singleSeq(chr_name, vinfo.ChrPosition):
                            LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                  "Sub_REF != Fasta-Seq: " + str(current_transcript.TID) + "\t" + str(
                                                      vinfo.ChrPosition) + "\n")
                        # Variant <-> Transcript with Direction Check
                        if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
                            if vinfo.Ref != current_transcript.seq_in_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_REF != CDS-Seq 1: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.Ref != current_transcript.seq_in_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_REF != CDS-Seq 2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.Alt != current_transcript.seq_in_iv_changed_dna_cds_seq(
                                    vinfo.Changed_CDS_Position,
                                    vinfo.Changed_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_ALT != CDS-Seq 3: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                        elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
                            if vinfo.Ref != current_transcript.seq_in_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition) and not current_transcript.CDS_Changed:
                                # When the used CDS gets changed, then it is of cause not possible to check against it.
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_REF != CDS-Seq 1.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_REF != CDS-Seq 2.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseAlt != current_transcript.seq_in_iv_changed_dna_cds_seq(
                                    vinfo.Changed_CDS_Position,
                                    vinfo.Changed_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_Alt != CDS-Seq 3.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_REF != CDS-Seq 4.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "Sub_REF != CDS-Seq 5.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                    elif TranscriptEnum.INSERTION in vinfo.Classification:
                        # Variant <-> Fasta Check
                        if vinfo.Ref != ghandler.singleSeq(chr_name, vinfo.ChrPosition):
                            LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                  "INSERTION_REF != Fasta-Seq: " + str(
                                                      current_transcript.TID) + "\t" + str(
                                                      vinfo.ChrPosition) + "\n")
                        # Variant <-> Transcript with Direction Check
                        if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
                            if vinfo.Ref != current_transcript.seq_in_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_REF != CDS-Seq 1: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.Ref != current_transcript.seq_in_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_REF != CDS-Seq 2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.Alt != current_transcript.seq_in_iv_changed_dna_cds_seq(
                                    vinfo.Changed_CDS_Position,
                                    vinfo.Changed_CDS_Position + len(vinfo.Alt) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_ALT != CDS-Seq 3: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                        elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
                            if vinfo.Ref != current_transcript.seq_in_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition) and not current_transcript.CDS_Changed:
                                # When the used CDS gets changed, then it is of cause not possible to check against it.
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_REF != CDS-Seq 1.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_REF != CDS-Seq 2.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseAlt != current_transcript.seq_in_iv_changed_dna_cds_seq(
                                    vinfo.Changed_CDS_Position,
                                    vinfo.Changed_CDS_Position + len(vinfo.Alt) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_Alt != CDS-Seq 3.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_REF != CDS-Seq 4.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "INSERTION_REF != CDS-Seq 5.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                    elif TranscriptEnum.DELETION in vinfo.Classification:
                        # Variant <-> Fasta Check
                        try:
                            if vinfo.Ref != ghandler.seq(chr_name, vinfo.ChrPosition,
                                                         vinfo.ChrPosition + len(vinfo.Ref) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != Fasta-Seq: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")
                        except SequenceHandlingError as she:
                            LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                    str(she.description) + str(
                                                      current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
                            if vinfo.Ref != she.sequence_part:
                                # ghandler.seq(chrName, vinfo.ChrPosition,#vinfo.ChrPosition + len(vinfo.Ref) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != Fasta-Seq: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                        # Variant <-> Transcript with Direction Check
                        if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
                            if vinfo.Ref != current_transcript.seq_in_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition + len(vinfo.Ref) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != CDS-Seq 1: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                            if vinfo.Ref != current_transcript.seq_in_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != CDS-Seq 2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                            if vinfo.Alt != current_transcript.seq_in_iv_changed_dna_cds_seq(
                                    vinfo.Changed_CDS_Position,
                                    vinfo.Changed_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_ALT != CDS-Seq 3: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                        elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
                            if vinfo.Ref != current_transcript.seq_in_cds(
                                    vinfo.ChrPosition,
                                    vinfo.ChrPosition + len(vinfo.Ref) - 1) and not current_transcript.CDS_Changed:
                                # When the used CDS gets changed, then it is of cause not possible to check against it.
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != CDS-Seq 1.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != CDS-Seq 2.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                            if vinfo.ReverseAlt != current_transcript.seq_in_iv_changed_dna_cds_seq(
                                    vinfo.Changed_CDS_Position,
                                    vinfo.Changed_CDS_Position):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_Alt != CDS-Seq 3.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds(
                                    vinfo.ChrPosition + (len(vinfo.Ref) - 1),
                                    vinfo.ChrPosition):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != CDS-Seq 4.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                            if vinfo.ReverseRef != current_transcript.seq_in_rev_cds_over_cds_position(
                                    vinfo.Unchanged_CDS_Position,
                                    vinfo.Unchanged_CDS_Position + len(vinfo.Ref) - 1):
                                LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                                      "DELETION_REF != CDS-Seq 5.2: " + str(
                                                          current_transcript.TID) + "\t" + str(
                                                          vinfo.ChrPosition) + "\n")

                    else:
                        LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_LOG,
                                              "No Classification in Complete_Check:\t" + str(
                                                  current_transcript.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
                except Exception:
                    LogOrganizer.add_to_log(LogEnums.COORDINATOR_COMPLETE_CHECK_CRITICAL_LOG,
                                            str(vinfo.ChrPosition) + "\t" + str(current_transcript.TID))

    def write_all_fasta(data_path: str,
                        data_name: str,
                        date_time: bool,
                        orig_dna: bool,
                        orig_aa: bool,
                        new_dna: bool,
                        new_aa: bool,
                        gff3: Gff3HandlerV3,
                        extend_file: bool,
                        chrom_name: str):
        """
        Writes a new fasta file with all CDS sequences from every transcript,
        including possible spliced variant of the transcripts (GFF3 data).
        It is possible to choose, which sequences will be written.
        Tags for the type of every sequence exist:
        ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == New Amino Acid Sequence
        Furthermore: Possible damaged transcripts will be seperated into "<data_name>_damaged.txt".
        Including its information about the gene, utr, exons and all available variant information.
        :param data_path: Path to the folder for the new fasta file.
        :param data_name: Name of the new fasta file.
        :param date_time: The date and time will be included inside the filename, if true.
        :param orig_dna:  Writes the original CDS DNA sequence, if true.
        :param orig_aa: Writes the original CDS A sequence, if true.
        :param new_dna: Writes the new CDS DNA sequence, if true.
        :param new_aa: Writes the new CDS AA sequence, if true.
        :param gff3: The gff3 data structure with all transcripts and their variants/sequences.
        :return: Nothing.
        """

        if orig_dna or orig_aa or new_dna or new_aa:
            if data_name == "":
                data_name = "all_transcripts_data"
            if date_time:
                timestop = str(datetime.now()) + "_"
            else:
                timestop = ""
            if not extend_file:
                normal_data_list_to_write = [
                    "#ODNA == Original DNA; OAA == Original Amino Acid Sequence; nDNA == new DNA; nAA == "
                    "New Amino Acid Sequence\n"]
                transcript_cds_damaged_data_list = [
                    "#Data in here may be incomplete and incorrect because of deletions, which destroyed parts "
                    "the splice sides.\n"]
            else:
                normal_data_list_to_write = []
                transcript_cds_damaged_data_list = []
            chr_name = chrom_name
            for current_transcript in gff3.get_chr_transcripts_dict(chr_name).values():
                current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(current_transcript)
                try:
                    if not current_transcript.cds_exist:
                        continue
                    if orig_dna:
                        if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
                            normal_data_list_to_write.append(
                                ">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.ODNA.value) + "\n")
                            normal_data_list_to_write.append(str(current_transcript.Complete_CDS) + "\n")
                        elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
                            normal_data_list_to_write.append(
                                ">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.ODNA.value) + "\n")
                            normal_data_list_to_write.append(str(current_transcript.Rev_CDS) + "\n")
                        else:
                            LogOrganizer.add_to_log(LogEnums.COORDINATOR_FASTA_FILE_ERROR_LOG,
                                                  "No Direction:write_all_fasta: " + str(current_transcript.TID) + "\n")
                    if orig_aa:
                        normal_data_list_to_write.append(
                            ">" + str(current_transcript.TID) + "|" + Fasta_Enum.OAA.value + "\n")
                        normal_data_list_to_write.append(str(current_transcript.IV_OriginalTranslation) + "\n")

                    if current_transcript.Transcript_CDS_damaged:
                        transcript_cds_damaged_data_list.append(
                            str(current_transcript.TID) + " " + str(current_transcript.Gene_Info_String))
                        for utr in current_transcript.UTR_Description:
                            if "\n" in utr:
                                transcript_cds_damaged_data_list.append(utr)
                            else:
                                transcript_cds_damaged_data_list.append(str(utr) + "\n")
                        for exon in current_transcript.EXON_Descriptin:
                            if "\n" in exon:
                                transcript_cds_damaged_data_list.append(exon)
                            else:
                                transcript_cds_damaged_data_list.append(str(exon) + "\n")
                        for vinfo in current_transcript.integrated_variant_objects_cds_hits:
                            # vinfo = For_Type_Safety.Variant_Information_Storage_Type_Safety(vinfo)
                            classificationstring = ""
                            class_list = vinfo.Classification
                            class_list_length = len(class_list)
                            for i in range(0, class_list_length):
                                classificationstring += class_list[i].value
                                if i != (class_list_length - 1):  # no tab after last entry
                                    classificationstring += "\t"
                            transcript_cds_damaged_data_list.append(
                                str(vinfo.ChrPosition) + "\t" +
                                str(vinfo.Ref) + "\t" +
                                str(vinfo.Alt) + "\t" +
                                str(classificationstring) + "\n")
                        continue
                    normal_data_list_to_write.append(
                        ">" + str(current_transcript.TID) + "|" + Fasta_Enum.ouchDNA.value + "\n")
                    normal_data_list_to_write.append(str(current_transcript.uChDNAsequence) + "\n")

                    normal_data_list_to_write.append(
                        ">" + str(current_transcript.TID) + "|" + Fasta_Enum.ouchAA.value + "\n")
                    normal_data_list_to_write.append(str(current_transcript.uChAAsequence) + "\n")

                    if new_dna:
                        normal_data_list_to_write.append(
                            ">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.nDNA.value) + "\n")
                        normal_data_list_to_write.append(str(current_transcript.IV_Changed_DNA_CDS_Seq) + "\n")
                    if new_aa:
                        normal_data_list_to_write.append(
                            ">" + str(current_transcript.TID) + "|" + str(Fasta_Enum.nAA.value) + "\n")
                        normal_data_list_to_write.append(str(current_transcript.iv_changed_translation) + "\n")
                except Exception:
                    LogOrganizer.add_to_log(LogEnums.COORDINATOR_TRANSCRIPT_LOG,
                                          "write_all_fasta" + "\t" + str(current_transcript.TID) + "\n")

            normal_output = "".join(normal_data_list_to_write)
            if not extend_file:
                new_fasta_file = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "w")
                new_fasta_file.write(normal_output)
                new_fasta_file.close()
            else:
                new_fasta_file = open(str(data_path) + str(timestop) + str(data_name) + ".fa", "a")
                new_fasta_file.write(normal_output)
                new_fasta_file.close()

            normal_output = ""
            normal_data_list_to_write = []
            if not extend_file:
                damaged_output = "".join(transcript_cds_damaged_data_list)
                new_fasta_file_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt",
                                              "w")
                new_fasta_file_damaged.write(damaged_output)
                new_fasta_file_damaged.close()
            else:
                damaged_output = "".join(transcript_cds_damaged_data_list)
                new_fasta_file_damaged = open(str(data_path) + str(timestop) + str(data_name) + "_damaged" + ".txt",
                                              "a")
                new_fasta_file_damaged.write(damaged_output)
                new_fasta_file_damaged.close()
            damaged_output = ""
            transcript_cds_damaged_data_list = []

    # * for stop.
    genetic_code = {'agg': 'r', 'aga': 'r', 'agc': 's', 'agu': 's',
                    'aag': 'k', 'aaa': 'k', 'aac': 'n', 'aau': 'n',
                    'acg': 't', 'aca': 't', 'acc': 't', 'acu': 't',
                    'aug': 'm', 'aua': 'i', 'auc': 'i', 'auu': 'i',

                    'cgg': 'r', 'cga': 'r', 'cgc': 'r', 'cgu': 'r',
                    'cag': 'q', 'caa': 'q', 'cac': 'h', 'cau': 'h',
                    'ccg': 'p', 'cca': 'p', 'ccc': 'p', 'ccu': 'p',
                    'cug': 'l', 'cua': 'l', 'cuc': 'l', 'cuu': 'l',

                    'ugg': 'w', 'uga': '*', 'ugc': 'c', 'ugu': 'c',
                    'uag': '*', 'uaa': '*', 'uac': 'y', 'uau': 'y',
                    'ucg': 's', 'uca': 's', 'ucc': 's', 'ucu': 's',
                    'uug': 'l', 'uua': 'l', 'uuc': 'f', 'uuu': 'f',

                    'ggg': 'g', 'gga': 'g', 'ggc': 'g', 'ggu': 'g',
                    'gag': 'e', 'gaa': 'e', 'gac': 'd', 'gau': 'd',
                    'gcg': 'a', 'gca': 'a', 'gcc': 'a', 'gcu': 'a',
                    'gug': 'v', 'gua': 'v', 'guc': 'v', 'guu': 'v', }
    stopcodon = "*"
    ###
    # Output fasta file
    ###
    orig_dna = True
    orig_aa = True
    new_dna = True
    new_aa = True
    fasta_data_name = ""

    ###
    # Output VCF file
    ###
    old_info = True
    ref_codons = True
    ref_aa = True
    alt_codons = True
    alt_aa = True
    ###
    # How many chromosomes will be analysed.
    # 0 stands for all.
    ###
    firstx_chr = 0
    time_start = datetime.now()

    vcf_path_and_name = invcf
    gff3_path_and_name = ingff
    fasta_file_path = infasta
    output_data_path = outpath

    print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    #######################################
    print("read vcf")
    vcf = VcfHandler(vcf_path_and_name, firstx_chr)
    print("Done: " + str(datetime.now() - time_start))
    print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    #######################################
    print("read gff3")
    # print("Hopefully sorted after seqID")
    gff3 = Gff3HandlerV3(gff3_path_and_name)
    print("Done: " + str(datetime.now() - time_start))
    print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    #######################################
    print("read fa")
    ghandler = Genomehandler(fasta_file_path, firstx_chr)
    print("Done: " + str(datetime.now() - time_start))
    print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    #######################################
    print("Chromosome:Names")
    if len(gff3.get_chromosome_names()) < 10:
        print("GFF3: " + str(gff3.get_chromosome_names()))
    else:
        print("GFF3: " + str(len(gff3.get_chromosome_names())) + " contigs/chromosomes.")
    if len(vcf.get_chromosome_names()) < 10:
        print("VCF: " + str(vcf.get_chromosome_names()))
    else:
        print("VCF: " + str(len(vcf.get_chromosome_names())) + " contigs/chromosomes.")
    if len(ghandler.GetChromosomeNames()) < 10:
        print("fasta: " + str(ghandler.GetChromosomeNames()))
    else:
        print("fasta:" + str(len(ghandler.GetChromosomeNames())) + " contigs/chromosomes.")
    shared_chromosomes = []
    shared_chromosomes_fa_gff3 = []
    #######################################

    for name in gff3.get_chromosome_names():
        if name in vcf.get_chromosome_names() and name in ghandler.GetChromosomeNames():
            shared_chromosomes.append(name)
        if name in ghandler.GetChromosomeNames():
            shared_chromosomes_fa_gff3.append(name)
    if len(shared_chromosomes) < 10:
        print(print("Shared_Chromosomes(all):" + str(shared_chromosomes)))
    else:
        print("Shared_Chromosomes(all):" + str(len(shared_chromosomes)) + " contigs/chromosomes.")

    # print("Shared_Chromosomes, fasta, gff3 (and used):" + str(Shared_Chromosomes_FA_GFF3))
    # For Writing the CDS
    # in all transcripts

    extend_file = False

    for name in shared_chromosomes:
        print(name)
        print("Write the CDS and Rev CDS")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

        for transcript in gff3.get_chr_transcripts_dict(name).values():
            try:
                transcript.complete_the_cds(ghandler.seq(name, transcript.StartOfRNA, transcript.EndOfRNA), genetic_code)
            except SequenceHandlingError as she:
                LogOrganizer.add_to_log(LogEnums.COORDINATOR_TRANSCRIPT_LOG, str(transcript.TID) + "\t"
                                        + str(transcript.StartOfRNA) + "\t" + str(transcript.EndOfRNA) + "\t" + str(
                    she.description))
                if she.sequence_part != "":
                    transcript.complete_the_cds(she.sequence_part, genetic_code)
            if transcript.ForwardDirection == TranscriptEnum.REVERSE:
                transcript.reverse_the_cds(genetic_code)
        # print("Done: " + str(datetime.now() - timeStart))
        ###
        # Connecting all transcripts with all their mutations/variants.
        # Contains handling of trialllele variants, if an original vcf file was used.
        # But it can't handle all possible appearances. VCF preprocessing exist for a reason.
        # The connection triggers a first evaluation of the variants, too.
        ###
        transcript_range = 300  # base pairs before RNA starts and after RNA ends. Just in Case for long InDels.
        variant_transcript_exceeding_warning_int = transcript_range / 100
        phasing_warning = True
        print("Connect transcripts with variants")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        # print("Search in " + name)

        phasing_variants = []  # should normally be empty
        current_transcript_id = 0  # initialization, first Round == 0
        first_transcript_match = 0  #
        current_transcript_list = gff3.get_chr_transcript_list(name)
        transcripts_with_multi_allel_variants = []
        transcripts_without_stop = []

        for variant in vcf.get_chr_vcf_variant_list(name):
            found = True
            current_transcript_id = first_transcript_match

            while (True):
                ###
                # The Following case: no more available transcripts, because the position of 
                # the mutation is after the last transcript
                ###
                if current_transcript_id == len(current_transcript_list):
                    # print("No more transcripts.")
                    break
                current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(
                    current_transcript_list[current_transcript_id])
                if (current_transcript.StartOfRNA - transcript_range) > variant.Position:
                    # transcriptstart higher than variant.position
                    break
                elif (current_transcript.StartOfRNA - transcript_range) <= variant.Position \
                        <= (current_transcript.EndOfRNA + transcript_range):  # between start and end (+ range)
                    current_transcript.add_variant_information(variant)
                    if "," in variant.Alternate and current_transcript not in transcripts_with_multi_allel_variants:
                        if phasing_warning:
                            phasing_warning = False
                            print('Warning: Still Phases inside VCF. '
                                  'It may not work correctly and can slow down the whole process.')
                            # not necessary with vcf-preprocessing - but maybe the
                            # data change in the future -> usefull again
                        transcripts_with_multi_allel_variants.append(current_transcript)
                    if "," in variant.Alternate and variant not in phasing_variants:
                        phasing_variants.append(variant)
                        LogOrganizer.add_to_log(LogEnums.COORDINATOR_PHASING_LOG, str(Variant.Chromosome) +
                                                str(Variant.Position) + "\t" + str(Variant.ID) + "\t" + str(
                            Variant.Reference)
                                                + "\t" + str(Variant.Alternate) + "\t" + str(Variant.Qual) + "\t" + str(
                            Variant.Filter)
                                                + "\t" + str(Variant.Info))
                    if current_transcript.lost_stop:
                        transcripts_without_stop.append(current_transcript)
                    current_transcript.ListofVariants.append(variant.Position)
                    variant.ListofTranscripts.append(current_transcript.IndexKey)
                    variant.SListofTranscripts.append(current_transcript.TID)
                    if found:
                        first_transcript_match = current_transcript_id
                        found = False
                    current_transcript_id += 1
                elif (variant.Position > current_transcript.EndOfRNA + transcript_range):
                    current_transcript_id += 1
                else:
                    print("Well this should not happen.")
                    break
        ####
        # For every transcript with multiple allele variants.
        # Creates a copy with all objects and data inside.
        # Takes only the first allele variant for the copy and
        # takes only the second allele variant for the original.
        # The copy gets a new IndexID.
        # Adds new transcript to the existing list in the gff3 handler.
        # Sorts the transcripts per chromosome and position again (like in gff3 handler).
        ####
        multiallelwarning = True
        for current_transcript in transcripts_with_multi_allel_variants:
            if multiallelwarning:
                print("[multiallelwarning] The dataset contains raw data. This may cause errors.")
                multiallelwarning = False
            current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(current_transcript)
            new_transcript = copy.deepcopy(current_transcript)
            new_transcript.remove_mult_allel_entry_in_all_variant_information(0)
            current_transcript.remove_mult_allel_entry_in_all_variant_information(1)
            new_transcript.IndexKey = gff3.get_next_transcript_index()
            current_transcript_list.append(new_transcript)
            # currentTranscriptList = sorted(currentTranscriptList, key=lambda sTranscript: sTranscript.StartOfRNA)
            # gff3.dictListOfTranscripts[gff3.dictChrNames[name]][newTranscript.TID] = newTranscript
            gff3.add_new_transcript_to_dict(name, new_transcript)

        gff3.update_transripts(current_transcript_list, name)
        # gff3.ListOfTranscripts[gff3.dictChrNames[name]] = currentTranscriptList

        # print("Done: " + str(datetime.now() - timeStart))

        kill_transcripts_with_to_many_variants = gff3.get_chr_transcripts_dict(name)
        for transcript in kill_transcripts_with_to_many_variants.values():
            transcript = ForTypeSafetyAndStatics.transcript_type_safety(transcript)
            if len(transcript.integrated_variant_objects_cds_hits) > 500:
                transcript.Transcript_CDS_damaged = True
                print(str(transcript.TID) + "has " + str(len(transcript.integrated_variant_objects_cds_hits)) +
                      " variations. The Limit is 500. It will be handled as damaged.")
        ########
        # Calculate the effect length of mutation effects (version 2)
        # Creates parts of the changed DNA- and AA-sequence.
        # The transcripts will be extended, until a stop appears (if there is none).
        ########
        print("Calculate the effect length")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        # print(name)
        nr_transcripts_current = 0
        a_transcript_dict = gff3.get_chr_transcripts_dict(
            name)  # self.dictListOfTranscripts = [{}] # List for chromosomes, dict for normal entrys
        for current_transcript in a_transcript_dict.values():
            current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(current_transcript)
            nr_transcripts_current += 1
            if not current_transcript.cds_exist:
                continue
            if current_transcript.Transcript_CDS_damaged:
                continue
            current_transcript.create_iv_changed_dna_cds_seq(genetic_code,
                                                             current_transcript.integrated_variant_objects_cds_hits,
                                                             stopcodon)
            if current_transcript.lost_stop and stopcodon in current_transcript.iv_changed_translation:
                current_transcript.lost_stop = False
                current_transcript.premature_stop_codon(stopcodon)
            elif stopcodon not in current_transcript.iv_changed_translation:
                current_transcript.lost_stop = True
            i = 0
            while current_transcript.lost_stop:
                i += 1
                current_transcript.create_iv_changed_translation(genetic_code)
                if stopcodon in current_transcript.iv_changed_translation:
                    last_pos_in_cds = current_transcript.last_cds_position
                    if current_transcript.found_new_stop or stopcodon in current_transcript.iv_changed_translation:
                        break
                    current_transcript.lost_stop = False
                    current_transcript.found_new_stop = True
                last_pos_in_cds = current_transcript.last_cds_position
                if current_transcript.ForwardDirection == TranscriptEnum.FORWARD:
                    try:
                        next_dna = ghandler.seq(name, last_pos_in_cds + 1, last_pos_in_cds + 100)
                    except SequenceHandlingError as she:
                        LogOrganizer.add_to_log(LogEnums.COORDINATOR_TRANSCRIPT_LOG,
                                                str(current_transcript.TID) + "\t" + she.description)
                        if she.sequence_part == "":
                            print("transcript " + str(current_transcript.TID) + " declared as broken(circular?). "
                                                                                "It reached the end of the "
                                                                                "contig/chrom.")
                            current_transcript.Transcript_CDS_damaged = True
                            break
                        else:
                            next_dna = she.sequence_part
                elif current_transcript.ForwardDirection == TranscriptEnum.REVERSE:
                    try:
                        next_dna = ghandler.seq(name, last_pos_in_cds - 100, last_pos_in_cds - 1)
                    except SequenceHandlingError as she:
                        LogOrganizer.add_to_log(LogEnums.COORDINATOR_TRANSCRIPT_LOG,
                                                str(current_transcript.TID) + "\t" + she.description)
                        if she.sequence_part == "":
                            print("transcript " + str(
                                current_transcript.TID) + "declared as broken(circular?). It reached the end of the "
                                                          "contig/chrom.")
                            current_transcript.Transcript_CDS_damaged = True
                            break
                        else:
                            next_dna = she.sequence_part
                else:
                    LogOrganizer.add_to_log(LogEnums.COORDINATOR_BUGHUNTING_LOG,
                                          "Error: No Direction: " + str(current_transcript.TID) + "\n")
                    break
                current_transcript.find_new_stop(next_dna, genetic_code, stopcodon)

                if len(current_transcript.IV_Changed_DNA_CDS_Seq) % 3 != 0 \
                        or len(current_transcript.Complete_CDS) % 3 != 0 \
                        or len(current_transcript.Rev_CDS) % 3 != 0:
                    current_transcript.check_last_variants(ghandler, genetic_code, stopcodon)
                if i == variant_transcript_exceeding_warning_int + 1:
                    print(str(current_transcript.TID) + " exceeds variation range limit of " + str(
                        int(variant_transcript_exceeding_warning_int * 100)) + " nucleotides.")
                if i == 20:
                    current_transcript.Transcript_CDS_damaged = True
                    print("transcript " + str(
                        current_transcript.TID) + " declared as broken, after extending it by: " + str(
                        i * 100) + " nucleotides.")
                    break
            j = 0
            while current_transcript.origDNAtoshort:
                # print(currentTranscript.TID)
                j += 1
                if j == 10:
                    current_transcript.Transcript_CDS_damaged = True
                    print(
                        "transcript " + str(current_transcript.TID) + " declared as broken, after extending it " + str(
                            j) + " times unsuccessfully.")
                    break
                current_transcript.check_last_variants(ghandler, genetic_code, stopcodon)
        ####################################################
        print("Complete AA sequences")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        create_sequences(gff3, orig_aa, new_aa, genetic_code, name)
        # Find Stop-Codons, if hey are made by a frameshift -> LABEl the mutation
        print("Flag Frameshift caused stops")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        a_transcript_dict = gff3.get_chr_transcripts_dict(name)
        for current_transcript in a_transcript_dict.values():
            current_transcript = ForTypeSafetyAndStatics.transcript_type_safety(current_transcript)
            if not current_transcript.cds_exist:
                continue
            if current_transcript.Transcript_CDS_damaged:
                continue
            first_stop_position = (current_transcript.iv_changed_translation.find(stopcodon) + 1) * 3
            if first_stop_position == -2:
                continue
            last_frameshifter = ""
            for variant in current_transcript.integrated_variant_objects_cds_hits:
                variant = ForTypeSafetyAndStatics.variant_information_storage_type_safety(variant)
                if TranscriptEnum.STOP_GAINED in variant.Classification \
                        or TranscriptEnum.STOP_CHANGED in variant.Classification:  # stop because of variants
                    last_frameshifter = ""
                    break
                if variant.Changed_CDS_Position < first_stop_position:
                    if (len(variant.Ref) - 1) % 3 != 0 or (
                            len(variant.Alt) - 1) % 3 != 0:  # only frameshift and caused no stop
                        if TranscriptEnum.STOP_GAINED not in variant.Classification \
                                and TranscriptEnum.STOP_CHANGED not in variant.Classification:
                            last_frameshifter = variant
                    continue
                else:
                    break
            if last_frameshifter == "":
                continue
            else:
                last_frameshifter.Classification.append(TranscriptEnum.STOP_CAUSED_IN)
                last_frameshifter.STOP_CAUSED_IN = first_stop_position - last_frameshifter.Changed_CDS_Position
        # print("Done: " + str(datetime.now() - timeStart))
        #################################################################
        print("Complete Data check.")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        complete_check(gff3, ghandler, name)
        # print("Done: " + str(datetime.now() - timeStart))
        #################################################################
        print("Write data")
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        """
        Old, but maybe still usefull functions for further research/programming.
        The functions are now inside old_stuff()
        # Write_New_Special_VCF_File(Output_Data_Path, gff3)
        # Write_All_VCF_Damaged_Transcripts(Output_Data_Path,gff3)
        # Write_VCF_With_Key(Output_Data_Path,gff3)
        # Write_All_VCF_Stop_Lost(Output_Data_Path,gff3)
        # Write_All_VCF_NO_STOP(Output_Data_Path,gff3)
        # Write_All_VCF_NO_START(Output_Data_Path,gff3)
        # Write_All_VCF_To_Many_Stops(Output_Data_Path,gff3)
        """
        write_all_vcf_file(output_data_path, gff3, old_info, ref_codons, ref_aa, alt_codons, alt_aa, extend_file, name)
        write_all_fasta(output_data_path, fasta_data_name, False, orig_dna, orig_aa, new_dna, new_aa, gff3, extend_file,
                        name)
        print("Done: " + str(datetime.now() - time_start))
        print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        extend_file = True
        gff3.free_ram(name)
        vcf.free_ram(name)
        gc.collect()
    #################################################################
    print("Create log files.")

    LogOrganizer.write_all_logs(outpath)

    print("Everything is done: " + str(datetime.now() - time_start))
    print("ram in mbyte\t" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

    return os.getpid()

# 2232

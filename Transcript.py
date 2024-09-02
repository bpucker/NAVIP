__author__ = "Jan-Simon Baasner"
__email__ = "janbaas@cebitec.uni-bielefeld.de"

from enum import Enum, unique
from VCF_Variant import Variant, VariantEnum
from LogOrganizer import LogEnums, LogOrganizer
from GenomeHandler import Genomehandler, SequenceHandlingError, Fasta_Enum


@unique
class TranscriptEnum(Enum):
    """
    This class supports enums for transcripts and the classification of variants inside these transcripts.
    """
    START_NOT_IN_CDS = "Error - Startposition is not in the CDS"
    END_NOT_IN_CDS = "Error - Endposition is not in the CDS"
    POSITION_NOT_IN_CDS = -1
    FORWARD = True
    REVERSE = False
    UNKNOWN_STRAND_DIRECTION = "No Direction for this transcript available"  # will hopefully never happen

    # Stuff for variants inside the transcript
    INSERTION = "INS"
    DELETION = "DEL"
    SUBSTITUTION = "SUB"
    FRAMESHIFT = "Frameshift"
    FRAMESHIFT_1 = "Frameshift+1"
    FRAMESHIFT_2 = "Frameshift+2"
    FRAMESHIFT_1_DEL = "Frameshift-1"
    FRAMESHIFT_2_DEL = "Frameshift-2"
    STOP_GAINED = "Stop gained"
    STOP_LOST = "Stop lost"
    STOP_CHANGED = "Stop changed"
    STOP_CAUSED_IN = "STOP_CAUSED_IN:".lower()
    AA_CHANGE = "Amino acid change"
    START_LOST = "Start lost"
    MISSENSE_VARIANT = "MISSENSE_VARIANT".lower()  # http://sequenceontology.org/browser/current_svn/term/SO:0001583

    # errors/warnings
    UNKNOWN_NUKLEOTIDE = "UNKNOWN_NUKLEOTIDE".lower()
    UNKNOWN_AMINOACID = "UNKNOWN_AMINOACID".lower()

    # when the transcript structure is damaged
    CDS_START_INVOLVED = "Hits before Start into CDS"
    CDS_STOP_INVOLVED = "Hits CDS Stop and after"
    CDS_INTRON_TO_EXON = "Hits from intron to exon"
    CDS_EXON_TO_INTRON = "Hits from exon to intron"
    CDS_INTRON_DELETION = "CDS_INTRON_DELETION".lower()
    CDS_EXON_DELETION = "CDS_EXON_DELETION".lower()


class Variantinformationstorage:
    """
    A storage class for in-transcript use only. It contains all necessary variant information.
    """

    def __init__(self, chr_position: int, ref: str, alt: str, unchanged_cds_position: int, id: str, qual: str,
                 filter: str, info: str):
        """
        Initialization of all necessary information.
        :param chr_position:  Position of the variant.
        :param ref: VCF REF entry.
        :param alt: VCF ALT entry.
        :param unchanged_cds_position: Position inside the original CDS.
        :param id: VCF ID entry.
        :param qual: VCF QUAL entry.
        :param filter: VCF FILTER entry.
        :param info: VCF INFO entry.
        """
        self.ChrPosition = chr_position
        self.Unchanged_CDS_Position = unchanged_cds_position  # without variant effects
        self.Changed_CDS_Position = 0  # with variant effects
        self.Ref = ref
        self.Alt = alt
        self.ID = id
        self.Qual = qual
        self.Filter = filter
        self.OLD_Info = info

        self.ReverseRef = ForTypeSafetyAndStatics.bio_reverse_seq(ref)
        self.ReverseAlt = ForTypeSafetyAndStatics.bio_reverse_seq(alt)
        if self.Unchanged_CDS_Position != TranscriptEnum.POSITION_NOT_IN_CDS:
            self.OrigRaster = ForTypeSafetyAndStatics.calculate_raster(self.Unchanged_CDS_Position)
        else:
            self.OrigRaster = self.Unchanged_CDS_Position
        self.Classification = []
        self.StartOfOwnEffect = 0
        self.EndOfOwnEffect = 0
        self.OrigTriplets = ""
        self.OrigRevTriplets = ""
        self.OrigAmino = ""
        # after variant effects
        self.Changed_Raster = 0
        self.ChangedTriplets = ""
        self.ChangedRevTriplets = ""
        self.NewAmino = ""

        self.STOP_CAUSED_IN = -1
        self.SharedEffectsWith = []

    def set_unchanged_cds_position(self, unchanged_cds_position: int):
        self.Unchanged_CDS_Position = unchanged_cds_position
        if self.Unchanged_CDS_Position != TranscriptEnum.POSITION_NOT_IN_CDS:
            self.OrigRaster = ForTypeSafetyAndStatics.calculate_raster(self.Unchanged_CDS_Position)
        else:
            self.OrigRaster = self.Unchanged_CDS_Position

    def set_changed_cds_position(self, changed_cds_position: int):
        """
        Sets the value for the changed cds position and calculates the changed raster, too.
        But only, if possible, this position is maybe outside the cds and will note it.
        :param changed_cds_position: New CDS position after all variant effects.
        :return: Nothing.
        """
        self.Changed_CDS_Position = changed_cds_position
        if type(changed_cds_position) == int:
            self.Changed_Raster = ForTypeSafetyAndStatics.calculate_raster(changed_cds_position)
        else:
            self.Changed_Raster = changed_cds_position


class Transcript:
    """
    Contains all information about a single transcript, including all variants in it's range.
    Contains all necessary functions to calculate classification of the variants inside a transcript.
    """

    def __init__(self, index_key: int, tid: str, start_of_rna: int, end_of_rna: int,
                 forward_direction: TranscriptEnum.REVERSE, chr: str):
        """
        Initialization of the transcript-class.
        :param index_key: Unique number, for identifying this transcript: 0...n.
        :param tid: Not unique(tri-allele)! Transcripts identification string from the GFF3-File: For example: "AT5G40340.1"
        :param start_of_rna: Start position inside the genome/dna-data. Not String-Position.
        :param end_of_rna: End position inside genome/dna-data. Not String-Position.
        :param forward_direction: TranscriptEnum.FORWARD or TranscriptEnum.REVERSE for transcript/strand orientation.
        """

        self.IndexKey = index_key
        self.Chr = chr
        self.TID = tid
        self.StartOfRNA = int(start_of_rna)
        self.EndOfRNA = int(end_of_rna)
        self.ListofCDS = []  # [(StartPos, EndPos, Raster)]
        self.ForwardDirection = forward_direction
        self.Complete_CDS = ""
        self.Rev_CDS = ""
        self.ListofVariants = []  # [VariantIDs]
        self.last_cds_position = 0
        self.Gene_Info_String = ""
        self.UTR_Description = []
        self.EXON_Descriptin = []
        self.Gene_Start_Position = 0
        self.Gene_End_Position = 0

        # original not changed data an sequence of the transcript
        self.uChDNAsequence = ""
        self.uChAAsequence = ""
        self.uChDNA_length = 0
        self.uChAA_length = 0

        # for integrated variants action:
        self.integrated_variant_objects_cds_hits = []
        self.IntegratedVariantObjects_NotCDS = []
        self.IV_Changed_DNA_CDS_Seq = ""
        self.IV_Ready = False
        self.IV_NewPositionList = [0]
        self.IV_OriginalTranslation = ""
        self.iv_changed_translation = ""
        self.IV_Count_Stops_in_New_AA = -1
        self.IV_ListOfEffects = []
        self.IV_DNA_length = 0
        self.IV_AA_length = 0

        # state's of transcript
        self.TID_locked = False  # for changing TID because of more than one multi_allel_variants
        self.Transcript_CDS_damaged = False  # exon or intron damaged
        self.cds_exist = False
        self.Rev_CDS_Exist = False
        self.MultiAllelVariants = False
        self.lost_stop = False
        self.found_new_stop = False
        self.CDS_Changed = False
        self.origDNAtoshort = False

    def set_gene_start_position(self, gene_start_position: int):
        """
        Set the start position of the gene.
        :param gene_start_position: Position in the chromosome.
        :return: Nothing.
        """
        self.Gene_Start_Position = gene_start_position

    def set_gene_end_position(self, gene_end_position: int):
        """
        Set the end position of the gene.
        :param gene_end_position: Position in the chromosome.
        :return: Nothing.
        """
        self.Gene_End_Position = gene_end_position

    def set_gene_info_string(self, gene_info_string: str):
        """
        Normally there exist an information string about the classification of the certain gene in the gff3 file.
        Especially hypothetical genes can make some problems, so the information can be stored for explanation,
        if something went wrong.
        :param gene_info_string: String description of the gene, if the gff3 file holds this information.
        :return: Nothing.
        """
        self.Gene_Info_String = gene_info_string

    def add_utr_description(self, utr_description: list):
        """
        Sometimes the gff3 file contains information about the UTR in some genes.
        It can be useful, if the transcript is somehow damaged from variants.
        :param utr_description: String from the gff3 file with the UTR description.
        :return: Nothing.
        """
        self.UTR_Description.append(utr_description)

    def add_exon_descriptin(self, exon_descriptin: list):
        """
        Sometimes the gff3 file contains information about the exon in some genes.
        It can be usefull, if the transcript is somehow damaged from variants.
        :param exon_descriptin: String from the gff3 file with the exon description.
        :return: Nothing
        """
        self.EXON_Descriptin.append(exon_descriptin)

    def change_last_cds_position(self, length: int):
        """
        Updates the CDS Position inside the ListofCDS list.
        Positive values: Longer CDS
        Negative values: smaller CDS
        :return: Nothing.
        """
        self.CDS_Changed = True
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            if length < 0:
                current_last_cds_length = abs(
                    self.ListofCDS[len(self.ListofCDS) - 1][0] - self.ListofCDS[len(self.ListofCDS) - 1][1])
                while current_last_cds_length < abs(length):
                    del (self.ListofCDS[len(self.ListofCDS) - 1])
                    length = -abs(current_last_cds_length + length)
                    current_last_cds_length = abs(
                        self.ListofCDS[len(self.ListofCDS) - 1][0] - self.ListofCDS[len(self.ListofCDS) - 1][1])
            self.ListofCDS[len(self.ListofCDS) - 1][1] = self.last_cds_position + length
            self.last_cds_position += length
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            if length < 0:
                current_last_cds_length = (abs(self.ListofCDS[0][0] - self.ListofCDS[0][1]))
                while current_last_cds_length < abs(length):
                    del (self.ListofCDS[0])
                    length = -abs(current_last_cds_length + length)
                    current_last_cds_length = (abs(self.ListofCDS[0][0] - self.ListofCDS[0][1]))
            self.ListofCDS[0][0] = self.last_cds_position - length
            self.last_cds_position -= length
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No Direction: " + str(self.TID) + "\n")

    def create_iv_original_translation(self, genetic_code: dict) -> bool:
        """
        Create the AminoAcid sequence from the original dna source without any variant effects.
        :param genetic_code: A dictionary which translates lowercase triplet RNA to AA. Example: genetic_code = {'agg':'r'}
        :return: Bool value if the transcription and translation of the DNA was successful.
        """

        if self.ForwardDirection == TranscriptEnum.FORWARD and self.Complete_CDS != "":
            self.IV_OriginalTranslation = ForTypeSafetyAndStatics.translation(
                ForTypeSafetyAndStatics.transcription(self.Complete_CDS), genetic_code)
            return True
        elif self.ForwardDirection == TranscriptEnum.REVERSE and self.Rev_CDS != "":
            self.IV_OriginalTranslation = ForTypeSafetyAndStatics.translation(
                ForTypeSafetyAndStatics.transcription(self.Rev_CDS), genetic_code)
            return True
        else:
            if self.Complete_CDS != "":
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_CDS_CREATION_LOG, "No CDS: " + str(self.TID) + "\n")
                # print("No CDS: " + str(self.TID))
                return False
            elif self.Rev_CDS != "":
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_CDS_CREATION_LOG, "No Rev CDS: " + str(self.TID) + "\n")
                # print("No Rev CDS: " + str(self.TID))
                return False
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_CDS_CREATION_LOG, "No direction?: " + str(self.TID) + "\n")
                # print ("No direction?: " + str(self.TID))
                return False

    def create_iv_changed_translation(self, genetic_code: dict) -> bool:
        """
        Create the AminoAcid sequence with the original dna source AND with variant effects inside the CDS.
        :param genetic_code: A dictionary which translates lowercase triplet RNA to AA. Example: genetic_code = {'agg':'r'}
        :return :Bool value if the transcription and translation of the DNA was successful.
        """

        if self.ForwardDirection == TranscriptEnum.FORWARD:
            self.iv_changed_translation = ForTypeSafetyAndStatics.translation(
                ForTypeSafetyAndStatics.transcription(self.IV_Changed_DNA_CDS_Seq), genetic_code)
            self.IV_AA_length = len(self.iv_changed_translation)
            return True
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            self.iv_changed_translation = ForTypeSafetyAndStatics.translation(
                ForTypeSafetyAndStatics.transcription(self.IV_Changed_DNA_CDS_Seq), genetic_code)
            self.IV_AA_length = len(self.iv_changed_translation)
            return True
        else:
            if self.IV_Changed_DNA_CDS_Seq == "" and self.Complete_CDS == "":
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No CDS: " + str(self.TID) + "\n")
                # print("No CDS: " + str(self.TID))
                return False
            elif self.IV_Changed_DNA_CDS_Seq == "":
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "No Changed_DNA_CDS_Seq, but CDS without changes exist: " + str(
                                            self.TID) + "\n")
                return False
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No direction?: " + str(self.TID) + "\n")
                # print("No direction?: " + str(self.TID))
                return False

    def create_iv_changed_dna_cds_seq(self, genetic_code: dict, integrated_variant_objects: list, stopcodon: str):
        """
        First Step:
            Simple classification for every Variant_Information_Storage object in SUB, DEl and INS.
            While doing this classification, the new CDS position for every variant will be calculated.
            In this step the new DNA sequence of the transcript will be created
        Second Step:
            For every variant the IV_Local_Classification will be used to create the exact changes inside the CDS.
            For every variant the IV_Local_Effect_Length will be used to calculate the length of the effects.
        Third Step.
            For all variants SetKeysToCombinedEffects will be used to add the information about side effects.
            -> After this step every variant contains the information with which other variant it shares an effect.
        :param genetic_code: The genetic code as a dictionary.
        :param integrated_variant_objects: A list from all variants in Variant_Information_Storage format, which all hits the CDS.
        :param stopcodon: Which letter stands for the stop codon.
        :return: Nothing.
        """
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            self.IV_Changed_DNA_CDS_Seq = self.Complete_CDS
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            # self.Rev_CDS = For_Type_Safety_and_statics.BioReverseSeq(self.Complete_CDS)
            self.IV_Changed_DNA_CDS_Seq = self.Rev_CDS

        integrated_variant_objects = sorted(integrated_variant_objects,
                                            key=lambda variant_information: variant_information.Unchanged_CDS_Position)
        # list needs to be ordered after position, so the new cds position can be calculated.
        # because every variants cds position can be changed with previous indels.
        for vinfo in integrated_variant_objects:
            if vinfo.ChrPosition == 5921700 and self.TID == "Ma09_t08910.1":
                print('bugsearch')

            if self.ForwardDirection == TranscriptEnum.FORWARD:
                alt = vinfo.Alt
                ref = vinfo.Ref
            elif self.ForwardDirection == TranscriptEnum.REVERSE:
                alt = vinfo.ReverseAlt
                ref = vinfo.ReverseRef
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Transcript without Direction: " + str(self.TID) + "\n")
                # print("Transcript without Direction: " + str(self.TID))
                continue
            cds_position = vinfo.Unchanged_CDS_Position
            current_additional_position = self.IV_NewPositionList[len(self.IV_NewPositionList) - 1]
            vinfo.set_changed_cds_position(cds_position + current_additional_position)
            if len(ref) == 1 and len(alt) == 1:
                firstkoord = cds_position + current_additional_position - 1
                if firstkoord != 0:
                    first = self.IV_Changed_DNA_CDS_Seq[0:firstkoord]
                else:
                    first = ""
                test = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position - 1]
                if test != ref:
                    LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                            "Error 1: SUB not identical with Ref: " + str(
                                                self.TID) + " " + str(vinfo.ChrPosition) + " " + str(
                                                self.ForwardDirection) + "\n")
                substitution = alt
                second = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position:]
                self.IV_Changed_DNA_CDS_Seq = first + substitution + second
                vinfo.Classification.append(TranscriptEnum.SUBSTITUTION)

                self.IV_NewPositionList.append(current_additional_position)
                self.IV_ListOfEffects.append(
                    (vinfo.ChrPosition, vinfo.Unchanged_CDS_Position, current_additional_position))
            elif len(ref) > 1 and len(alt) == 1:  # del
                dellength = (len(ref) - 1)

                if self.ForwardDirection == TranscriptEnum.REVERSE:
                    if cds_position + current_additional_position + dellength <= 3 \
                            and cds_position + current_additional_position + dellength > 0:
                        first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position - 1]
                        LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_ADDITIONAL_INFO_LOG,
                                                "Note 1: DEL removed start: " + str(self.TID) + "\t" + str(
                                                    vinfo.ChrPosition) + "\n")
                        vinfo.Classification.append(TranscriptEnum.START_LOST)
                    elif cds_position + current_additional_position + dellength <= 0:
                        LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_ADDITIONAL_INFO_LOG,
                                                "Note 1.1: DEL destroyed start exon: " + str(self.TID) + "\t" + str(
                                                    vinfo.ChrPosition) + "\n")
                        vinfo.Classification.append(TranscriptEnum.CDS_START_INVOLVED)
                        first = ""
                    else:
                        first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position - 1]

                    second = self.IV_Changed_DNA_CDS_Seq[cds_position + dellength + current_additional_position - 1:]
                    if second[0] != ref[len(ref) - 1]:
                        LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                                "Error 2: DEL not identical with Ref: " + str(self.TID) + "\t" + str(
                                                    vinfo.ChrPosition) + "\n")
                else:
                    first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position]

                    second = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position + dellength:]
                    if first != "":
                        if first[len(first) - 1] != ref[0]:
                            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                                    "Error 3: DEL not identical with Ref: " + str(
                                                        self.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
                self.IV_Changed_DNA_CDS_Seq = first + second
                self.IV_NewPositionList.append(current_additional_position - (-1 + len(ref)))
                self.IV_ListOfEffects.append(
                    (vinfo.ChrPosition, vinfo.Unchanged_CDS_Position, current_additional_position - (-1 + len(ref))))
                vinfo.Classification.append(TranscriptEnum.DELETION)
            elif len(ref) == 1 and len(alt) > 1:  # insert
                test = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position - 1]
                first = self.IV_Changed_DNA_CDS_Seq[0:cds_position + current_additional_position - 1]
                insert = alt
                second = self.IV_Changed_DNA_CDS_Seq[cds_position + current_additional_position:]
                self.IV_Changed_DNA_CDS_Seq = first + insert + second
                vinfo.Classification.append(TranscriptEnum.INSERTION)

                self.IV_NewPositionList.append(current_additional_position + (-1 + len(alt)))
                self.IV_ListOfEffects.append(
                    (vinfo.ChrPosition, vinfo.Unchanged_CDS_Position, current_additional_position + (-1 + len(alt))))
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Create_IV_IV_Changed_DNA_CDS_Seq - error?" + str(self.TID) + str(
                                            vinfo.ChrPosition) + "\n")
                self.IV_NewPositionList.append(current_additional_position)
                self.IV_ListOfEffects.append(
                    (vinfo.ChrPosition, vinfo.Unchanged_CDS_Position, current_additional_position))
                continue
        for vinfo in integrated_variant_objects:
            self.origDNAtoshort = False
            self.iv_local_classification(vinfo, genetic_code, stopcodon)
            self.iv_local_effect_length(vinfo)
        self.IV_DNA_length = len(self.IV_Changed_DNA_CDS_Seq)
        self.set_keys_to_combined_effects(integrated_variant_objects)

    def set_keys_to_combined_effects(self, integrated_variant_objects: list):
        """
        Every variant, which has an effect to another variant will be calculated inside this function.
        Example:A variant which will cause a frameshift get the keys from every variant after this frameshift.
                Every variant inside this frameshift will get the key from the one it caused it.
        :param integrated_variant_objects: List of all classified variants. They all needs to have the effect length already.
        :return: Nothing.
        """

        if len(integrated_variant_objects) > 0:

            vinfo = ForTypeSafetyAndStatics.variant_information_storage_type_safety(integrated_variant_objects[0])
            list_with_variants = [vinfo]
            marked = False
            in_the_end = False
            for i in range(1, len(integrated_variant_objects)):
                current_vinfo = ForTypeSafetyAndStatics.variant_information_storage_type_safety(
                    integrated_variant_objects[i])
                listed_for_deletion = []
                for old_vinfo in list_with_variants:
                    if old_vinfo.StartOfOwnEffect == current_vinfo.StartOfOwnEffect:
                        marked = True
                    elif old_vinfo.EndOfOwnEffect == VariantEnum.NO_EndOfOwnEffect:
                        if TranscriptEnum.STOP_GAINED in old_vinfo.Classification:
                            if TranscriptEnum.DELETION in old_vinfo.Classification:
                                endeffect_without_stop = old_vinfo.Changed_CDS_Position + (2 - old_vinfo.Changed_Raster)
                            elif TranscriptEnum.INSERTION in old_vinfo.Classification:
                                endeffect_without_stop = old_vinfo.Changed_CDS_Position + (len(old_vinfo.Alt) - 1) + (
                                        2 - old_vinfo.Changed_Raster)
                            elif TranscriptEnum.SUBSTITUTION in old_vinfo.Classification \
                                    or TranscriptEnum.AA_CHANGE in old_vinfo.Classification:
                                endeffect_without_stop = old_vinfo.Changed_CDS_Position + (2 - old_vinfo.Changed_Raster)
                            else:
                                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                                        "Error: SetKeysToCombinedEffects:" + str(self.TID) +
                                                        "\t" + str(old_vinfo.ChrPosition) +
                                                        "\t" + str(current_vinfo.ChrPosition) + "\n")
                                break
                            marked = True
                        else:
                            marked = True
                    elif old_vinfo.EndOfOwnEffect >= current_vinfo.StartOfOwnEffect:
                        marked = True
                    else:
                        listed_for_deletion.append(old_vinfo)
                    if marked:
                        old_vinfo.SharedEffectsWith.append(current_vinfo)
                        current_vinfo.SharedEffectsWith.append(old_vinfo)
                        marked = False
                if len(listed_for_deletion) > 0:
                    for deleteit in listed_for_deletion:
                        if deleteit == []:
                            pass
                        list_with_variants.remove(deleteit)
                if in_the_end:
                    if len(current_vinfo.SharedEffectsWith) > 0:
                        for vinfo_in_shared_effects in current_vinfo.SharedEffectsWith:
                            vinfo_in_shared_effects.SharedEffectsWith.remove(current_vinfo)
                            current_vinfo.SharedEffectsWith.remove(vinfo_in_shared_effects)

                    break
                list_with_variants.append(current_vinfo)

    def raster_zero(self, dna: str, vinfo):
        return dna[vinfo.Unchanged_CDS_Position - 1:vinfo.Unchanged_CDS_Position + 2]

    def raster_one(self, dna: str, vinfo):
        return dna[vinfo.Unchanged_CDS_Position - 2:vinfo.Unchanged_CDS_Position + 1]

    def raster_two(self, dna: str, vinfo):
        return dna[vinfo.Unchanged_CDS_Position - 3:vinfo.Unchanged_CDS_Position]

    def raster_zero_changed(self, dna: str, vinfo):
        return dna[vinfo.Changed_CDS_Position - 1:vinfo.Changed_CDS_Position + 2]

    def raster_one_changed(self, dna: str, vinfo):
        return dna[vinfo.Changed_CDS_Position - 2:vinfo.Changed_CDS_Position + 1]

    def raster_two_changed(self, dna: str, vinfo):
        return dna[vinfo.Changed_CDS_Position - 3:vinfo.Changed_CDS_Position]

    def premature_stop_codon(self, stopcodon):
        stop_start_position_in_changed_cds = 1 + (self.iv_changed_translation.find(stopcodon)) * 3
        stop_end_position_in_changed_cds = 1 + (1 + self.iv_changed_translation.find(stopcodon)) * 3
        last_possible_chr_position = 0

        for tripple in self.IV_ListOfEffects:
            chr_pos = tripple[0]
            cds_pos = tripple[1]
            change = tripple[2]
            if cds_pos + change <= stop_start_position_in_changed_cds:
                last_possible_chr_position = chr_pos
                continue
            elif cds_pos + change <= stop_end_position_in_changed_cds:
                last_possible_chr_position = chr_pos
                continue
            else:
                break

        self.remove_outsider_variants(last_possible_chr_position)

    def remove_outsider_variants(self, chr_position: int):
        marked_for_deletion = []
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            for vinfo in self.integrated_variant_objects_cds_hits:
                if vinfo.ChrPosition > chr_position:
                    marked_for_deletion.append(vinfo)
        else:
            for vinfo in self.integrated_variant_objects_cds_hits:
                if vinfo.ChrPosition < chr_position:
                    marked_for_deletion.append(vinfo)

        for vinfo in marked_for_deletion:
            self.integrated_variant_objects_cds_hits.remove(vinfo)
            self.IntegratedVariantObjects_NotCDS.append(vinfo)

    def iv_local_classification(self, vinfo: Variantinformationstorage, genetic_code: dict, stopcodon: str):
        """
        Adds (and calculates) following information to the Variant_Information_Storage object:
        orig_triplets -> original codons for the position of the variant.
        changed_triplets -> codons for the NEW position in the NEW transcript.
        Classifications -> FRAMESHIFT_1; FRAMESHIFT_2; FRAMESHIFT_1_DEL; FRAMESHIFT_2_DEL;
                        -> STOP_LOST; STOP_GAINED; AA_CHANGE
        :param vinfo: A single Variant_Information_Storage object.
        :param genetic_code: Dictionary of the genetic code.
        :param stopcodon: Character of the stopcodon.
        :return: Nothing.
        """
        chr_pos = vinfo.ChrPosition
        new_raster = vinfo.Changed_Raster
        orig_raster = vinfo.OrigRaster

        cds_position = vinfo.Unchanged_CDS_Position
        cds_position_changed = vinfo.Changed_CDS_Position

        if self.ForwardDirection == TranscriptEnum.FORWARD:
            forward = True
            current_cds = self.Complete_CDS
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            forward = False
            current_cds = self.Rev_CDS
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    "Error: Direction: " + str(self.TID) + " " + str(chr_pos) + "\n")
            return False

        if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
            if orig_raster == 0:
                orig_triplets = self.raster_zero(current_cds, vinfo)
            elif orig_raster == 1:
                orig_triplets = self.raster_one(current_cds, vinfo)
            elif orig_raster == 2:
                orig_triplets = self.raster_two(current_cds, vinfo)
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            orig_raster) + "\n")
                return False
            if new_raster == 0:
                changed_triplets = self.raster_zero_changed(self.IV_Changed_DNA_CDS_Seq, vinfo)
            elif new_raster == 1:
                changed_triplets = self.raster_one_changed(self.IV_Changed_DNA_CDS_Seq, vinfo)
            elif new_raster == 2:
                changed_triplets = self.raster_two_changed(self.IV_Changed_DNA_CDS_Seq, vinfo)
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            new_raster) + "\n")
                return False
            if forward:
                vinfo.OrigTriplets = orig_triplets
                vinfo.OrigAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(orig_triplets), genetic_code)
            else:
                vinfo.OrigTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(orig_triplets)
                vinfo.OrigAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(orig_triplets), genetic_code)
                vinfo.OrigAmino = vinfo.OrigAmino[::-1]
            if forward:
                vinfo.ChangedTriplets = changed_triplets
                vinfo.NewAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(changed_triplets), genetic_code)
            else:
                vinfo.ChangedTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(changed_triplets)
                vinfo.NewAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(changed_triplets), genetic_code)
                vinfo.NewAmino = vinfo.NewAmino[::-1]

        elif TranscriptEnum.INSERTION in vinfo.Classification:
            ###
            #
            # Rasterchange = (len(Alt) -1) % 3
            # Raster	Rasterchange		(Normal without change)	Before	Next
            # 0				0				2						0		2
            # 0				1				2						0		1
            # 0				2				2						0		0
            #
            # 1				0				1						1		1
            # 1				1				1						1		0
            # 1				2				1						1		2
            #
            # 2				0				0						2		0
            # 2				1				0						2		2
            # 2				2				0						2		1
            ###
            raster_change = (len(vinfo.Alt) - 1) % 3

            if raster_change == 1:
                vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_1)
            elif raster_change == 2:
                vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_2)

            if orig_raster == 0:
                # before = ""
                before = current_cds[cds_position - 1:cds_position + 2]
                next = ""
            elif orig_raster == 1:
                # before = current_cds[cds_position-2:cds_position-1]
                before = current_cds[cds_position - 2:cds_position + 1]
                next = ""
            elif orig_raster == 2:
                # before = current_cds[cds_position-3:cds_position-1]
                before = current_cds[cds_position - 3:cds_position]
                next = ""
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            orig_raster) + "\n")
                return False
            if new_raster == 0:
                before2 = ""
            elif new_raster == 1:
                before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed - 2:cds_position_changed - 1]
            elif new_raster == 2:
                before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed - 3:cds_position_changed - 1]
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            new_raster) + "\n")
                return False
            if new_raster == 0 and raster_change == 0:
                next2 = self.IV_Changed_DNA_CDS_Seq[
                        cds_position_changed + len(vinfo.Alt) - 1: cds_position_changed + 2 + len(vinfo.Alt) - 1]
            elif new_raster == 0 and raster_change == 1:
                next2 = self.IV_Changed_DNA_CDS_Seq[
                        cds_position_changed + len(vinfo.Alt) - 1: cds_position_changed + 1 + len(vinfo.Alt) - 1]
            elif new_raster == 0 and raster_change == 2:
                next2 = ""
            elif new_raster == 1 and raster_change == 0:
                next2 = self.IV_Changed_DNA_CDS_Seq[
                        cds_position_changed + len(vinfo.Alt) - 1: cds_position_changed + 1 + len(vinfo.Alt) - 1]
            elif new_raster == 1 and raster_change == 1:
                next2 = ""
            elif new_raster == 1 and raster_change == 2:
                next2 = self.IV_Changed_DNA_CDS_Seq[
                        cds_position_changed + len(vinfo.Alt) - 1: cds_position_changed + 2 + len(vinfo.Alt) - 1]
            elif new_raster == 2 and raster_change == 0:
                next2 = ""
            elif new_raster == 2 and raster_change == 1:
                next2 = self.IV_Changed_DNA_CDS_Seq[
                        cds_position_changed + len(vinfo.Alt) - 1: cds_position_changed + 2 + len(vinfo.Alt) - 1]
            elif new_raster == 2 and raster_change == 2:
                next2 = self.IV_Changed_DNA_CDS_Seq[
                        cds_position_changed + len(vinfo.Alt) - 1: cds_position_changed + 1 + len(vinfo.Alt) - 1]
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            new_raster) + " " + str(raster_change) + "\n")
                return False
            if forward:
                vinfo.OrigTriplets = before + next
                vinfo.OrigAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(vinfo.OrigTriplets), genetic_code)
            else:
                vinfo.OrigTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(before + next)
                vinfo.OrigAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(before + next), genetic_code)
                vinfo.OrigAmino = vinfo.OrigAmino[::-1]

            if forward:
                vinfo.ChangedTriplets = before2 + vinfo.Alt + next2
                vinfo.NewAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(vinfo.ChangedTriplets), genetic_code)
            else:
                vinfo.ChangedTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(
                    before2 + ForTypeSafetyAndStatics.bio_reverse_seq(vinfo.Alt) + next2)
                vinfo.NewAmino = ForTypeSafetyAndStatics.translation(ForTypeSafetyAndStatics.transcription(
                    before2 + ForTypeSafetyAndStatics.bio_reverse_seq(vinfo.Alt) + next2), genetic_code)
                vinfo.NewAmino = vinfo.NewAmino[::-1]

        elif TranscriptEnum.DELETION in vinfo.Classification:

            raster_change = (len(vinfo.Ref) - 1) % 3

            if raster_change == 1:
                vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_1_DEL)
            elif raster_change == 2:
                vinfo.Classification.append(TranscriptEnum.FRAMESHIFT_2_DEL)

            if orig_raster == 0:
                # before = current_cds[cds_position - 1:cds_position] # del position
                # next = current_cds[cds_position + len(vinfo.Ref) -1 :cds_position + len(vinfo.Ref) -1  +2]
                before = ""
                next = ""
            elif orig_raster == 1:
                # before = current_cds[cds_position - 2:cds_position] #del position, too
                # next =  current_cds[cds_position +len(vinfo.Ref) -1 :cds_position +len (vinfo.Ref)-1 +1]
                before = current_cds[cds_position - 2:cds_position - 1]
                next = ""
            elif orig_raster == 2:
                before = current_cds[cds_position - 3:cds_position - 1]
                next = ""
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            orig_raster) + "\n")
                return False

            if orig_raster == 0 and raster_change == 0:
                next = current_cds[cds_position + len(vinfo.Ref) - 1: cds_position + len(vinfo.Ref) - 1 + 2]
            elif orig_raster == 0 and raster_change == 1:
                next = current_cds[cds_position + len(vinfo.Ref) - 1: cds_position + len(vinfo.Ref) - 1 + 1]
            elif orig_raster == 0 and raster_change == 2:
                next = ""
            elif orig_raster == 1 and raster_change == 0:
                next = current_cds[cds_position + len(vinfo.Ref) - 1: cds_position + len(vinfo.Ref) - 1 + 1]
            elif orig_raster == 1 and raster_change == 1:
                next = ""
            elif orig_raster == 1 and raster_change == 2:
                next = current_cds[cds_position + len(vinfo.Ref) - 1: cds_position + len(vinfo.Ref) - 1 + 2]
            elif orig_raster == 2 and raster_change == 0:
                next = ""
            elif orig_raster == 2 and raster_change == 1:
                next = current_cds[cds_position + len(vinfo.Ref) - 1: cds_position + len(vinfo.Ref) - 1 + 2]
            elif orig_raster == 2 and raster_change == 2:
                next = current_cds[cds_position + len(vinfo.Ref) - 1: cds_position + len(vinfo.Ref) - 1 + 1]
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            orig_raster) + " " + raster_change + "\n")
                return False

            if new_raster == 0:
                # before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed - 1: cds_position_changed + 2]
                # next2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed + len(vinfo.Ref) -1:
                # cds_position_changed +len(vinfo.Ref) -1 +2]
                before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed - 1: cds_position_changed + 2]
                next2 = ""
            elif new_raster == 1:
                # before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed-2 : cds_position_changed + 1]
                # next2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed + len(vinfo.Ref) -1 :
                # cds_position_changed + len(vinfo.Ref)-1 +1]
                before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed - 2: cds_position_changed + 1]
                next2 = ""
            elif new_raster == 2:
                # before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed-3: cds_position_changed]
                before2 = self.IV_Changed_DNA_CDS_Seq[cds_position_changed - 3: cds_position_changed]
                next2 = ""
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "(Impossible) Raster Error: " + str(self.TID) + " " + str(chr_pos) + " " + str(
                                            new_raster) + "\n")
                return False

            if forward:
                vinfo.OrigTriplets = before + vinfo.Ref + next
                vinfo.OrigAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(vinfo.OrigTriplets), genetic_code)
            else:
                vinfo.OrigTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(
                    before + ForTypeSafetyAndStatics.bio_reverse_seq(vinfo.Ref) + next)
                vinfo.OrigAmino = ForTypeSafetyAndStatics.translation(ForTypeSafetyAndStatics.transcription(
                    before + ForTypeSafetyAndStatics.bio_reverse_seq(vinfo.Ref) + next), genetic_code)
                vinfo.OrigAmino = vinfo.OrigAmino[::-1]

            if forward:
                vinfo.ChangedTriplets = before2 + next2
                vinfo.NewAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(vinfo.ChangedTriplets), genetic_code)
            else:
                vinfo.ChangedTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(before2 + next2)
                vinfo.NewAmino = ForTypeSafetyAndStatics.translation(
                    ForTypeSafetyAndStatics.transcription(before2 + next2), genetic_code)
                vinfo.NewAmino = vinfo.NewAmino[::-1]

        orig_amino = vinfo.OrigAmino
        new_amino = vinfo.NewAmino
        vinfo.OrigRevTriplets = ForTypeSafetyAndStatics.bio_reverse_seq(vinfo.OrigTriplets)

        if stopcodon in orig_amino and not stopcodon in new_amino:
            vinfo.Classification.append(TranscriptEnum.STOP_LOST)
            self.lost_stop = True
            self.calculate_last_cds_position()
        elif stopcodon in new_amino and not stopcodon in orig_amino:
            vinfo.Classification.append(TranscriptEnum.STOP_GAINED)
        elif stopcodon in new_amino and stopcodon in orig_amino:
            if len(new_amino) == len(orig_amino):
                vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
            elif TranscriptEnum.SUBSTITUTION in vinfo.Classification:
                vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
            elif TranscriptEnum.DELETION in vinfo.Classification:
                vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
            elif TranscriptEnum.INSERTION in vinfo.Classification:
                if vinfo.NewAmino[0] == stopcodon:
                    vinfo.Classification.append(TranscriptEnum.STOP_CHANGED)
                else:
                    vinfo.Classification.append(TranscriptEnum.STOP_LOST)
                    vinfo.Classification.append(TranscriptEnum.STOP_GAINED)
        elif orig_amino != new_amino:
            vinfo.Classification.append(TranscriptEnum.AA_CHANGE)
        if stopcodon in self.iv_changed_translation:
            self.lost_stop = False
        if len(vinfo.ChangedTriplets) % 3 != 0 or len(vinfo.OrigTriplets) % 3 != 0 or len(
                vinfo.OrigRevTriplets) % 3 != 0 or len(vinfo.ChangedRevTriplets) % 3 != 0:
            self.origDNAtoshort = True
        if 'n' in vinfo.ChangedTriplets.lower():
            vinfo.Classification.append(TranscriptEnum.UNKNOWN_NUKLEOTIDE)
        if 'x' in vinfo.NewAmino.lower():
            vinfo.Classification.append(TranscriptEnum.UNKNOWN_AMINOACID)

    def iv_local_effect_length(self, vinfo: Variantinformationstorage):
        """
        Calculates the effect length and add it to a single Variant_Information_Storage (object).
        The starteffect is always a position, because every variant in here is inside the transcript.
        The endeffect can be a position and the NO_EndOfOwnEffect enum value.
        Frameshifts or other high impact effects (Stop) will effect every following variant.
        :param vinfo: A single Variant_Information_Storage object, already containing its classifications.
        :return: False, if anything went wrong.
        """
        starteffect = vinfo.Changed_CDS_Position - vinfo.Changed_Raster
        classification = vinfo.Classification

        if TranscriptEnum.STOP_GAINED in classification \
                or TranscriptEnum.STOP_LOST in classification \
                or TranscriptEnum.FRAMESHIFT in classification \
                or TranscriptEnum.FRAMESHIFT_1 in classification \
                or TranscriptEnum.FRAMESHIFT_2 in classification \
                or TranscriptEnum.FRAMESHIFT_1_DEL in classification \
                or TranscriptEnum.FRAMESHIFT_2_DEL in classification:
            if TranscriptEnum.FRAMESHIFT in classification:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "Still FRAMESHIFT in use." + "\n")
            endeffect = VariantEnum.NO_EndOfOwnEffect
        elif TranscriptEnum.DELETION in classification:
            endeffect = vinfo.Changed_CDS_Position + (2 - vinfo.Changed_Raster)
        elif TranscriptEnum.INSERTION in classification:
            endeffect = vinfo.Changed_CDS_Position + (len(vinfo.Alt) - 1) + (2 - vinfo.Changed_Raster)
        elif TranscriptEnum.SUBSTITUTION in classification or TranscriptEnum.AA_CHANGE in classification:
            endeffect = vinfo.Changed_CDS_Position + (2 - vinfo.Changed_Raster)
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    "Error in IV_Local_Effect_Length: " + str(self.TID) + " " + str(
                                        vinfo.ChrPosition) + "\n")
            return False

        vinfo.StartOfOwnEffect = starteffect
        vinfo.EndOfOwnEffect = endeffect

    def calculate_last_cds_position(self):
        """
        The last position of the CDS in the chromosome may be needed to extend the cds.
        :return: The last position of the CDS in the chromosome. In reverse direction it
                is of course the position on the left side.
        """
        if self.last_cds_position != 0 or len(self.ListofCDS) == 0:
            return False
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            self.last_cds_position = self.ListofCDS[len(self.ListofCDS) - 1][1]
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            self.last_cds_position = self.ListofCDS[0][0]
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No Direction: " + str(self.TID) + "\n")

    def update_last_cds_position(self):
        """
        Updates the CDS Position inside the ListofCDS list.
        :return: Nothing.
        """
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            self.ListofCDS[len(self.ListofCDS) - 1][1] = self.last_cds_position
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            self.ListofCDS[0][0] = self.last_cds_position
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG, "No Direction: " + str(self.TID) + "\n")

    def find_new_stop(self, next_dna: str, genetic_code: dict, stopcodon: str):
        """
        Uses the IV_Changed_DNA_CDS_Seq and add the nextDNA, translates it and search for a new stopcodon.
        :param next_dna: DNA, which will be added to the existing transcript.
        :param genetic_code: Dictionary, used for translation.
        :param stopcodon: Char, which will be used as stop.
        :return: False, if there is no new stop.
        """
        current_raster = ForTypeSafetyAndStatics.calculate_raster(len(self.IV_Changed_DNA_CDS_Seq))
        ###
        # len -1 % 3 --> last position is inside raster x
        ###
        # 0 -> last position + new
        # 1 -> last -1 + last + new
        # 2 -> last -2 + last -1 + last + new -----> only new
        ###
        # new stop position
        # dna to this code -< append changed_cds
        ###
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            if current_raster == 0:
                dna_to_test = self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq) - 1] + next_dna
                old_dna = 1
            elif current_raster == 1:
                dna_to_test = self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq) - 2] + \
                            self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq) - 1] + next_dna
                old_dna = 2
            elif current_raster == 2:
                dna_to_test = next_dna
                old_dna = 0
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Error in find_New_Stop: " + str(self.TID) + " " + str(current_raster) + "\n")
                return False
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            if current_raster == 0:
                dna_to_test = self.IV_Changed_DNA_CDS_Seq[
                                len(self.IV_Changed_DNA_CDS_Seq) - 1] + ForTypeSafetyAndStatics.bio_reverse_seq(next_dna)
                old_dna = 1
            elif current_raster == 1:
                dna_to_test = self.IV_Changed_DNA_CDS_Seq[len(self.IV_Changed_DNA_CDS_Seq) - 2] + \
                            self.IV_Changed_DNA_CDS_Seq[
                                len(self.IV_Changed_DNA_CDS_Seq) - 1] + ForTypeSafetyAndStatics.bio_reverse_seq(next_dna)
                old_dna = 2
            elif current_raster == 2:
                dna_to_test = ForTypeSafetyAndStatics.bio_reverse_seq(next_dna)
                old_dna = 0
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Error in find_New_Stop: " + str(self.TID) + " " + str(current_raster) + "\n")
                return False
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    "Error in find_New_Stop, not transcript direction: " + str(self.TID) + " " + str(
                                        current_raster) + "\n")
            return False
        new_amino = ForTypeSafetyAndStatics.translation(ForTypeSafetyAndStatics.transcription(dna_to_test),
                                                        genetic_code)
        position_in_string = new_amino.find(stopcodon)

        if position_in_string != -1:
            ###
            # dna_position_in_string = 3*Amino_Acid_Position_in_String
            # its alsways the start position -> Position 5 in amino == Stop is in dna_position 15,16,17
            ###
            self.lost_stop = False
            self.found_new_stop = True
            if self.ForwardDirection == TranscriptEnum.FORWARD:
                self.IV_Changed_DNA_CDS_Seq += next_dna[0:(1 + position_in_string) * 3 - old_dna]
                self.change_last_cds_position((1 + position_in_string) * 3 - old_dna)
                self.Complete_CDS += next_dna[0:(1 + position_in_string) * 3 - old_dna]
                self.Rev_CDS += ForTypeSafetyAndStatics.bio_reverse_seq(
                    next_dna[0:(1 + position_in_string) * 3 - old_dna])
            else:
                revdna = ForTypeSafetyAndStatics.bio_reverse_seq(next_dna)[0:(1 + position_in_string) * 3 - old_dna]
                self.IV_Changed_DNA_CDS_Seq += revdna
                self.Rev_CDS += revdna
                self.Complete_CDS = ForTypeSafetyAndStatics.bio_reverse_seq(revdna) + self.Complete_CDS
                self.change_last_cds_position((1 + position_in_string) * 3 - old_dna)
            self.iv_check_for_new_variants(genetic_code, stopcodon)
        else:

            self.Rev_CDS += ForTypeSafetyAndStatics.bio_reverse_seq(next_dna)
            if self.ForwardDirection == TranscriptEnum.FORWARD:
                self.change_last_cds_position(len(next_dna))
                self.IV_Changed_DNA_CDS_Seq += next_dna
                self.Complete_CDS += next_dna
            else:
                self.Complete_CDS = next_dna + self.Complete_CDS
                self.change_last_cds_position(len(next_dna))
                self.IV_Changed_DNA_CDS_Seq += ForTypeSafetyAndStatics.bio_reverse_seq(next_dna)
            return False

    def normalize_variant_classification(self):
        for variant in self.integrated_variant_objects_cds_hits:
            variant = ForTypeSafetyAndStatics.variant_information_storage_type_safety(variant)
            normal_classification = []
            for classification in variant.Classification:
                if classification in normal_classification:
                    continue
                else:
                    normal_classification.append(classification)
            variant.Classification = normal_classification

    def reset_transcript(self):
        self.IV_Changed_DNA_CDS_Seq = ""
        self.IV_NewPositionList = [0]
        self.IV_Changed_DNA_CDS_Seq = []
        self.IV_OriginalTranslation = ""
        self.iv_changed_translation = ""
        self.IV_Count_Stops_in_New_AA = -1
        variantlist = []
        for vinfo in self.integrated_variant_objects_cds_hits:
            vinfo = ForTypeSafetyAndStatics.variant_information_storage_type_safety(vinfo)
            variantlist.append(Variant("nnh",
                                       vinfo.ChrPosition,
                                       vinfo.ID,
                                       -1,
                                       vinfo.Ref,
                                       vinfo.Alt,
                                       vinfo.Qual,
                                       vinfo.Filter,
                                       vinfo.OLD_Info))
        self.integrated_variant_objects_cds_hits = []
        for variant in variantlist:
            self.add_variant_information(variant)

    def iv_check_for_new_variants(self, genetic_code: dict, stopcodon: str):
        """
        Prototype function: It is more a reminder, that it is possible for longer transcripts, to have more variants.
        Not completed.
        :param genetic_code: Dictionary for the genetic code.
        :param stopcodon: Character for the stopcodon.
        :return: Nothing.
        """
        ###
        # need a new function for finding cds position in the new, longer cds.
        # current cds function needs old cds from gff3 file - maybe this can be modified
        # or the last entry can be used for this
        #
        ###
        self.update_last_cds_position()
        rmlist = []
        for vinfo in self.IntegratedVariantObjects_NotCDS:
            vinfo = ForTypeSafetyAndStatics.variant_information_storage_type_safety(vinfo)
            if self.ForwardDirection == TranscriptEnum.FORWARD:
                if self.search_position_in_cds(vinfo.ChrPosition) == TranscriptEnum.POSITION_NOT_IN_CDS:
                    continue
            elif self.ForwardDirection == TranscriptEnum.REVERSE:
                if self.search_position_in_cds_reverse(vinfo.ChrPosition) == TranscriptEnum.POSITION_NOT_IN_CDS:
                    continue  # nnh == not needed here
            variant = Variant("nnh",
                              vinfo.ChrPosition,
                              vinfo.ID,
                              -1,
                              vinfo.Ref,
                              vinfo.Alt,
                              vinfo.Qual,
                              vinfo.Filter,
                              vinfo.OLD_Info)
            rmlist.append(vinfo)
            self.add_variant_information(variant)
        for vinfo in rmlist:
            self.IntegratedVariantObjects_NotCDS.remove(vinfo)
        self.reset_transcript()
        self.create_iv_changed_dna_cds_seq(genetic_code, self.integrated_variant_objects_cds_hits, stopcodon)
        self.normalize_variant_classification()

    def iv_cds_damaging_variant_handler(self, vinfo: Variantinformationstorage):
        """
        Variants, which are classified as deletions do have two different positions:
        Position of the variant and endposition of the deletion.
        If one of them is inside the CDS and one outside the transcript will be marked as damaged.
        The classification of the 'damage' will be made inside this function and added to the
        incomming Variant_Information_Storage object.
        :param vinfo: Variant_Information_Storage object with a transcript 'damaging' deletion
        :return: Nothing.
        """
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            normal_cds = vinfo.Unchanged_CDS_Position
            ###
            # 4 additional cases DELETION (Forward):
            # Del-length + chrPos == Inside CDS :;: BUT ChrPos alone == Outside CDS
            # A1-> before transcript starts, into transcript -> Start lost
            # A2-> starts inside intron, into exon -> Spliceside destroyed
            # Del-length + chrPos == Outside CDS :;: BUT ChrPos alone == Inside CDS
            # B1-> Start: Last-Exon; End: After-Last Exon -> Stop destroyed
            # B2-> Start: In one Exon; End: In one Intron -> Spliceside destroyed
            ###
            # While bughunting, new cases appeared:
            #  -> before transcript starts into intron -> Start + a lot more lost
            #  -> starts inside intron, into next intron -> exon or more lost + Spliceside(s) destroyed
            #  -> starts inside intron, ends after stop -> Stop destroyed + more
            #  -> starts inside exon, ends in next exon -> intron lost + Spliceside(s) destroyed

            second_position = vinfo.ChrPosition + (len(vinfo.Ref) - 1)
            cds_2 = self.search_position_in_cds(second_position)

            first_cds = self.ListofCDS[0]
            last_cds = self.ListofCDS[len(self.ListofCDS) - 1]

            # case A
            if normal_cds == TranscriptEnum.POSITION_NOT_IN_CDS and type(cds_2) == int:

                # case A1
                if first_cds[0] <= second_position <= first_cds[1]:
                    vinfo.Classification.append(TranscriptEnum.CDS_START_INVOLVED)
                # case A2
                else:
                    vinfo.Classification.append(TranscriptEnum.CDS_INTRON_TO_EXON)
            # case B
            elif type(normal_cds) == int and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:

                # case B1
                if last_cds[0] <= vinfo.ChrPosition <= last_cds[1]:
                    vinfo.Classification.append(TranscriptEnum.CDS_STOP_INVOLVED)
                # case B2
                else:
                    vinfo.Classification.append(TranscriptEnum.CDS_EXON_TO_INTRON)
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Error - IV_CDS_Damaging_Variant_Handler?,Forward: " + str(
                                            self.TID) + "\t" + str(vinfo.ChrPosition) + "\n")
        else:
            normal_cds = self.search_position_in_cds_reverse(vinfo.ChrPosition)
            ###
            # 4 additional cases DELETION (Reverse):
            # Del-Length + ChrPos == Inside CDS :;: BUT ChrPos alone == Outside CDS
            # 	A1 -> before transcript starts, into transcript -> Stop lost
            #   A2 -> starts inside intron, into exon -> Spliceside destroyed
            # Del-length + chrPos == Outside CDS :;: BUT ChrPos alone == Inside CDS
            #   B1 -> Start: Last-Exon(Forward Direction); End: After-Last Exon -> Start destroyed
            #   B2 -> Start: In one Exon; End: In one Intron -> Spliceside destroyed
            ###
            second_position = vinfo.ChrPosition + (len(vinfo.Ref) - 1)
            cds_2 = self.search_position_in_cds_reverse(second_position)

            first_cds = self.ListofCDS[0]  # stop in here for reverse transcripts
            last_cds = self.ListofCDS[len(self.ListofCDS) - 1]  # start in here for reverse transcripts

            # case A
            if normal_cds == TranscriptEnum.POSITION_NOT_IN_CDS and type(cds_2) == int:
                # case A1
                if first_cds[0] <= second_position <= first_cds[1]:
                    vinfo.Classification.append(TranscriptEnum.CDS_STOP_INVOLVED)
                # case A2
                else:
                    vinfo.Classification.append(TranscriptEnum.CDS_INTRON_TO_EXON)
            # case B
            elif type(normal_cds) == int and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:
                # case B1
                if last_cds[0] <= vinfo.ChrPosition <= last_cds[1]:
                    vinfo.Classification.append(TranscriptEnum.CDS_START_INVOLVED)
                # case B2
                else:
                    vinfo.Classification.append(TranscriptEnum.CDS_EXON_TO_INTRON)
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Error - IV_CDS_Damaging_Variant_Handler?,Reverse: " + str(
                                            self.TID) + "\t" + str(vinfo.ChrPosition) + "\n")

    def add_variant_information(self, variant: Variant):
        """
        The raw data object VCF_Variant contains all information from the original vcf file.
        But one VCF variant object can have multiple classifications and modifications in different transcripts.
        The solution of this was: A new Variant_Information_Storage object for every variant which is inside
        a transcript.
        --> VCF_Variant goes in
            -> deletion will be checked, if they delete splice sides
                -> IV_CDS_Damaging_Variant_Handler -> classification
            -> store and save data in a Variant_Information_Storage object for this transcript.
                -> list for variants, inside cds
                -> list for variants, outside the cds
        :param variant: VCF_Variant class object, which starting position is between the rna start <-> end (+ range)
        :return: Variant_Information_Storage.
        """

        if self.ForwardDirection == TranscriptEnum.FORWARD:
            cds_position = self.search_position_in_cds(variant.Position)
            new_entry = Variantinformationstorage(variant.Position,
                                                  variant.Reference,
                                                  variant.Alternate,
                                                  cds_position,
                                                  variant.ID,
                                                  variant.Qual,
                                                  variant.Filter,
                                                  variant.Info)
            if len(variant.Reference) > 1:  # DEL
                cds_2 = self.search_position_in_cds(variant.Position + len(variant.Reference) - 1)
                larger_impact_deletion = self.search_position_pair_around_transcript(variant.Position, variant.Position +
                                                                                     len(variant.Reference) - 1)

        else:  # self.ForwardDirection == TranscriptEnum.REVERSE:
            cds_position = self.search_position_in_cds_reverse(variant.Position)

            if len(variant.Reference) > 1:  # DEL
                cds_position = self.search_position_in_cds_reverse(variant.Position + (len(variant.Reference) - 1))
                cds_2 = self.search_position_in_cds_reverse(variant.Position)
                larger_impact_deletion = self.search_position_pair_around_transcript(variant.Position, variant.Position +
                                                                                     len(variant.Reference) - 1)

            new_entry = Variantinformationstorage(variant.Position,
                                                  variant.Reference,
                                                  variant.Alternate,
                                                  cds_position,
                                                  variant.ID,
                                                  variant.Qual,
                                                  variant.Filter,
                                                  variant.Info)
        if len(variant.Reference) > 1:  # DEL
            if larger_impact_deletion[0] != larger_impact_deletion[1]:
                self.Transcript_CDS_damaged = True
                is_it_already_in = False
                for variant_in_list in self.integrated_variant_objects_cds_hits:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.integrated_variant_objects_cds_hits.append(new_entry)
                if larger_impact_deletion[1] % 1 == float(0):
                    new_entry.Classification.append(TranscriptEnum.CDS_INTRON_DELETION)
                elif larger_impact_deletion[1] % 1 == float(0.5):
                    new_entry.Classification.append(TranscriptEnum.CDS_EXON_DELETION)
                else:
                    LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                            "Error 1 :Possible bug in Add_Variant_Information: " + str(
                                                self.TID) + "\t" + str(
                                                variant.Position) + "\n")
            if type(cds_position) == int and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:
                self.Transcript_CDS_damaged = True
                is_it_already_in = False
                for variant_in_list in self.integrated_variant_objects_cds_hits:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.integrated_variant_objects_cds_hits.append(new_entry)
                self.iv_cds_damaging_variant_handler(new_entry)
            elif cds_position == TranscriptEnum.POSITION_NOT_IN_CDS and type(cds_2) == int:
                self.Transcript_CDS_damaged = True
                is_it_already_in = False
                for variant_in_list in self.integrated_variant_objects_cds_hits:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.integrated_variant_objects_cds_hits.append(new_entry)
                self.iv_cds_damaging_variant_handler(new_entry)
            elif cds_position == TranscriptEnum.POSITION_NOT_IN_CDS and cds_2 == TranscriptEnum.POSITION_NOT_IN_CDS:
                is_it_already_in = False
                for variant_in_list in self.IntegratedVariantObjects_NotCDS:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.IntegratedVariantObjects_NotCDS.append(new_entry)
            elif type(cds_position) == int and type(cds_2) == int:
                is_it_already_in = False
                for variant_in_list in self.integrated_variant_objects_cds_hits:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.integrated_variant_objects_cds_hits.append(new_entry)
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Error 2 :Possible bug in Add_Variant_Information: " + str(
                                            self.TID) + "\t" + str(
                                            variant.Position) + "\n")
        else:
            if cds_position == TranscriptEnum.POSITION_NOT_IN_CDS:
                is_it_already_in = False
                for variant_in_list in self.IntegratedVariantObjects_NotCDS:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.IntegratedVariantObjects_NotCDS.append(new_entry)
            elif type(cds_position) == int:
                is_it_already_in = False
                for variant_in_list in self.integrated_variant_objects_cds_hits:
                    current_pos = variant_in_list.ChrPosition
                    if current_pos == variant.Position:
                        is_it_already_in = True
                        break
                if not is_it_already_in:
                    self.integrated_variant_objects_cds_hits.append(new_entry)
                return new_entry
            else:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                        "Error 3 :Possible bug in Add_Variant_Information: " + str(
                                            self.TID) + "\t" + str(variant.Position) + "\n")

    def remove_mult_allel_entry_in_all_variant_information(self, zero_or_one: int):
        """
        When no vcf preprocessing is used there may exist triallel variants.
        Here they will be split into the first and second entry.
        One transcript will only get all first entrys, the other all second entrys.
        Will add an "A" or a "B" to the TID.
        :param zero_or_one: Decides if this transcript will get the first or the second entrys.
        :return: Nothing.
        """
        for multi_allel_variant in self.integrated_variant_objects_cds_hits:
            if "," in multi_allel_variant.Alt:

                self.MultiAllelVariants = True
                multi_allel_variant.Alt = multi_allel_variant.Alt.split(",")[zero_or_one]

                if zero_or_one == 0 and not self.TID_locked:
                    self.TID = str(self.TID) + "B"
                    self.TID_locked = True
                elif zero_or_one == 1 and not self.TID_locked:
                    self.TID = str(self.TID) + "A"
                    self.TID_locked = True
                else:
                    LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                            "Maybe bug in :Remove_Mult_Allel_Entry_In_All_Variant_Information with :"
                                            + self.TID + "\n")

    def add_cds(self, start_position: int, end_position: int, raster: str):
        """
        Adds the CDS from the gff3 file to this transcript, one entry at the time.
        Every entry stands for an exon.
        The raster argument will not be used.
        :param start_position: Startposition of the exon in the chromosome.
        :param end_position: Endposition of the exon in the chromosome.
        :param raster: The gff3 files often contains a raster parameter.
        :return: Nothing.
        """
        self.ListofCDS.append(
            [int(start_position), int(end_position), str(raster)]
        )
        if (self.ListofCDS[0] == []):
            self.ListofCDS.pop(0)

    def complete_the_cds(self, seq_of_transcript: str, genetic_code: dict):
        """
        After the CDS list is complete: This function will calculate the DNA CDS from a given sequence.
        The sequence, which will given to this function has to start at the startposition of the transcript.
        (StartOfRNA argument)
        :param seq_of_transcript: String DNA sequence from the chromosome.
        :return: Nothing.
        """
        completethe_cds = []
        for CDS in self.ListofCDS:
            completethe_cds.append(  # short list of all CDS seqs
                seq_of_transcript[((CDS[0] - self.StartOfRNA)):
                                ((abs(CDS[0] - self.StartOfRNA) + abs(CDS[1] - CDS[0]))) + 1])
        self.Complete_CDS = "".join(completethe_cds)
        self.change_cds_existence(True)
        if (self.ListofCDS == []):
            self.change_cds_existence(False)
        if (self.Complete_CDS == ""):
            self.change_cds_existence(False)
        self.calculate_last_cds_position()
        self.uChDNAsequence = self.Complete_CDS
        self.uChDNA_length = len(self.uChDNAsequence)
        self.uChAAsequence = ForTypeSafetyAndStatics.translation(
            ForTypeSafetyAndStatics.transcription(self.uChDNAsequence), genetic_code)
        self.uChAA_length = len(self.uChAAsequence)

    def reverse_the_cds(self, genetic_code: dict) -> bool:
        """
        This will reverse the existing Complete_CDS for reverse direction transcripts.
        :return: False, if there is no Complete_CDS. True  if it succeded.
        """
        if self.Rev_CDS_Exist:
            return True
        if (self.Complete_CDS == ""):
            return False
        else:
            # A - C and G - T
            # A To c, C to a, G to t, T to g, then uppercase
            self.Rev_CDS = self.Complete_CDS.replace("A", "t")
            self.Rev_CDS = self.Rev_CDS.replace("C", "g")
            self.Rev_CDS = self.Rev_CDS.replace("G", "c")
            self.Rev_CDS = self.Rev_CDS.replace("T", "a")
            self.Rev_CDS = self.Rev_CDS.upper()
            self.Rev_CDS = self.Rev_CDS[::-1]  # reverse #
            self.Rev_CDS_Exist = True
            self.uChDNAsequence = self.Rev_CDS
            self.uChAAsequence = ForTypeSafetyAndStatics.translation(
                ForTypeSafetyAndStatics.transcription(self.uChDNAsequence), genetic_code)
            return True

    def search_position_in_cds(self, position_in_chr: int) -> int:
        """
        Calculates the position inside the CDS.
        Only in forward direction.
        :param position_in_chr: Position inside the chromosome.
        :return: The CDS Position or POSITION_NOT_IN_CDS enum.
        """
        position_in_cds = TranscriptEnum.POSITION_NOT_IN_CDS
        cds_length = 0

        if len(self.ListofCDS) == 0 or position_in_chr < self.ListofCDS[0][0] or position_in_chr > \
                self.ListofCDS[len(self.ListofCDS) - 1][1]:
            return position_in_cds

        for CDS in self.ListofCDS:
            if (CDS[0] <= position_in_chr <= CDS[1]):
                position_in_cds = cds_length + abs(position_in_chr - CDS[0]) + 1
                break  # break, when you found it
            else:
                # +1, because... from position 5 to 10 are not 5 positions, it's 6: P:[5,6,7,8,9,10]
                cds_length += abs(CDS[1] - CDS[0]) + 1  # EndPosition - StartPosition = length

        return position_in_cds

    def search_position_pair_around_transcript(self, startposition_in_chr: int, endposition_in_chr: int):
        """
        CDS are numbered from 0 up. The output of this funtcion is a list with two values from -0.5 to the maximal
        CDS number +0.5 toprescind where the Positions are.
        :param startposition_in_chr: The position of the variant.
        :param endposition_in_chr: The position of the variant + length of deletion.
        :return: List of two floats.
        """
        locked_first = False
        locked_second = False
        result = [float(0.0), float(0.0)]
        for i, CDS in enumerate(self.ListofCDS):
            if startposition_in_chr < CDS[0] and not locked_first:
                result[0] = float(i) - float(0.5)
                locked_first = True
            elif CDS[0] <= startposition_in_chr <= CDS[1] and not locked_first:
                result[0] = float(i)
                locked_first = True
            elif startposition_in_chr > CDS[1] and i - 1 == len(self.ListofCDS) and not locked_first:
                result[0] = float(i) + float(0.5)
                locked_first = True

            if endposition_in_chr < CDS[0] and not locked_second:
                result[1] = float(i) - float(0.5)
                locked_second = True
            elif CDS[0] <= endposition_in_chr <= CDS[1] and not locked_second:
                result[1] = float(i)
                locked_second = True
            elif endposition_in_chr > CDS[1] and i - 1 == len(self.ListofCDS) and not locked_second:
                result[1] = float(i) + float(0.5)
                locked_second = True
        return result

    def search_position_in_cds_reverse(self, position_in_chr: int) -> int:
        """
        Calculates the position inside the CDS.
        Only in reverse direction.
        :param position_in_chr: Position inside the chromosome.
        :return: The CDS Position or POSITION_NOT_IN_CDS enum.
        """
        position_in_cds = TranscriptEnum.POSITION_NOT_IN_CDS
        cds_length = 0

        if len(self.ListofCDS) == 0 or position_in_chr < self.ListofCDS[0][0] or position_in_chr > \
                self.ListofCDS[len(self.ListofCDS) - 1][1]:
            return position_in_cds

        for CDS in self.ListofCDS[::-1]:
            # CDS[0] == Start of CDS
            # CDS[1] == End of CDS
            # CDS is always sorted in direction of forward strand DNA ...
            # for reverse strand DNA we have to reverse the positions, too
            # here in reverse is  CDS[1] the starting Position, not the End,
            # but still have the highest position, because its sorted forward
            if CDS[0] <= position_in_chr <= CDS[1]:
                # End |CDS[0] - PositionInChr| + current length inside CDS
                position_in_cds = cds_length + abs(position_in_chr - CDS[1]) + 1
                break
            else:
                cds_length += abs(CDS[1] - CDS[0]) + 1
        return position_in_cds

    def seq_in_cds(self, start_pos_chr: int, end_pos_chr: int) -> str:
        """
        Calculates a sequence between the start and end position inside the cds.
        Works only in forward strand direction.
        :param start_pos_chr: Position inside the chromosome.
        :param end_pos_chr: Position inside the chromosome.
        :return: The CDS position or START_NOT_IN_CDS/END_NOT_IN_CDS enums.
        """
        start_in_cds = self.search_position_in_cds(start_pos_chr)
        end_in_cds = self.search_position_in_cds(end_pos_chr)

        if (start_in_cds == TranscriptEnum.POSITION_NOT_IN_CDS):
            return TranscriptEnum.START_NOT_IN_CDS
        elif (end_in_cds == TranscriptEnum.POSITION_NOT_IN_CDS):
            return TranscriptEnum.END_NOT_IN_CDS
        else:
            if (start_pos_chr == end_pos_chr):
                return self.Complete_CDS[start_in_cds - 1]
            else:
                return self.Complete_CDS[(start_in_cds - 1): start_in_cds + abs(end_in_cds - start_in_cds)]

    def seq_in_rev_cds(self, start_pos_chr: int, end_pos_chr: int) -> str:
        """
        Calculates a sequence between the start and end position inside the cds.
        Works only in reverse strand direction.
        :param start_pos_chr: Position inside the chromosome.
        :param end_pos_chr: Position inside the chromosome.
        :return: The CDS position or START_NOT_IN_CDS/END_NOT_IN_CDS enums.
        """
        start_in_cds = self.search_position_in_cds_reverse(start_pos_chr)
        end_in_cds = self.search_position_in_cds_reverse(end_pos_chr)

        if (start_in_cds == -1):
            return TranscriptEnum.START_NOT_IN_CDS

        elif (end_in_cds == -1):
            return TranscriptEnum.END_NOT_IN_CDS

        else:
            if (start_pos_chr == end_pos_chr):
                return self.Rev_CDS[start_in_cds - 1]
            else:
                return self.Rev_CDS[start_in_cds - 1: start_in_cds + abs(end_in_cds - start_in_cds)]

    def seq_in_strand_direction(self, start_pos_chr: int, end_pos_chr: int) -> str:
        """
        Calculates a sequence between the start and end position inside the cds.
        Works in both strand direction, using the other two direction functions.
        :param start_pos_chr: Position inside the chromosome.
        :param end_pos_chr: Position inside the chromosome.
        :return: The CDS position or START_NOT_IN_CDS/END_NOT_IN_CDS enums.
        """
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            return self.seq_in_cds(start_pos_chr, end_pos_chr)
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            return self.seq_in_rev_cds(start_pos_chr, end_pos_chr)
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    "Error: No Strand Direction: " + str(self.TID) + "\n")
            # print("Error: No Strand Direction: " + str(self.TID) + "\n")
            return False

    def change_cds_existence(self, does_it_exist: bool):
        """
        Setter method for the bool value CDS_Exist.
        :param does_it_exist:
        :return:
        """
        self.cds_exist = does_it_exist

    def seq_in_cds_over_cds_position(self, start_pos_cds: int, end_pos_cds: int) -> str:
        """
        Calculates the sequence between two CDS positions.
        Works only in forward strand direction.
        :param start_pos_cds: Position inside the CDS.
        :param end_pos_cds: Position inside the CDS.
        :return: Sequence between the two Positions. Or one single base, if start and end position are identical.
        """
        if start_pos_cds == end_pos_cds:
            return self.Complete_CDS[start_pos_cds - 1]
        else:
            return self.Complete_CDS[
                   start_pos_cds - 1:
                   start_pos_cds + abs(start_pos_cds - end_pos_cds)]

    def seq_in_rev_cds_over_cds_position(self, start_pos_cds: int, end_pos_cds: int):
        if start_pos_cds == end_pos_cds:
            return self.Rev_CDS[start_pos_cds - 1]
        else:
            return self.Rev_CDS[
                   start_pos_cds - 1:
                   start_pos_cds + abs(start_pos_cds - end_pos_cds)]

    def seq_in_strand_direction_cds_over_cds_position(self, start_pos_cds: int, end_pos_cds: int):
        """
        Calculates the sequence between two CDS positions.
        Works only in reverse strand direction.
        :param start_pos_cds: Position inside the CDS.
        :param end_pos_cds: Position inside the CDS.
        :return: Sequence between the two Positions. Or one single base, if start and end position are identical.
        """
        if self.ForwardDirection == TranscriptEnum.FORWARD:
            return self.seq_in_cds_over_cds_position(start_pos_cds, end_pos_cds)
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            return self.seq_in_rev_cds_over_cds_position(start_pos_cds, end_pos_cds)

    def seq_in_iv_changed_dna_cds_seq(self, start_pos_new_cds: int, end_pos_new_cds: int) -> str:
        """
        Calculates the sequence between two CDS positions.
        Works in both strand direction, because the IV_Changed_DNA_CDS_Seq
        will only calculated for the current direction.
        :param StartPosCDS: Position inside the CDS.
        :param EndPosCDS: Position inside the CDS.
        :return: Sequence between the two Positions. Or one single base, if start and end position are identical.
        """
        if start_pos_new_cds == end_pos_new_cds:
            return self.IV_Changed_DNA_CDS_Seq[start_pos_new_cds - 1]
        else:
            return self.IV_Changed_DNA_CDS_Seq[
                   start_pos_new_cds - 1:
                   start_pos_new_cds + abs(start_pos_new_cds - end_pos_new_cds)]

    def check_last_variants(self, genomehandler: Genomehandler, genetic_code, stopcodon):
        # +2 positions means 2 potential more variant effects.
        # repeatly, and check if deletions now have a new effect, because in rev cds it can be .... .... ....
        # should be done, because the transcript gets reseted after getting larger cds
        if len(self.IV_Changed_DNA_CDS_Seq) % 3 != 0:
            raster = ForTypeSafetyAndStatics.calculate_raster(
                len(self.IV_Changed_DNA_CDS_Seq))  # 2 is impossible, 1 == +1 base, 0 == +2 base
        elif len(self.Complete_CDS) % 3 != 0:
            raster = ForTypeSafetyAndStatics.calculate_raster(len(self.Complete_CDS))
        elif len(self.Rev_CDS) % 3 != 0:
            raster = ForTypeSafetyAndStatics.calculate_raster(len(self.Rev_CDS))
        else:
            erg = 3
            raster = 2
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    "checkLastVariants\t" + str(self.TID) + "\tIV_Changed_DNA_CDS_Seq\n")
            pass

        erg = 0
        if raster == 0:
            erg = 2
        elif raster == 1:
            erg = 1
        elif raster == 2:
            erg = 3
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    "checkLastVariants\t" + str(self.TID) + "\tIV_Changed_DNA_CDS_Seq\n")
            pass
        further = 1
        if len(self.integrated_variant_objects_cds_hits) > 0:
            if len(self.integrated_variant_objects_cds_hits[len(self.integrated_variant_objects_cds_hits) - 1].Ref) > 1:
                # 1 minimum + complete dellength -1 because of string ---> length
                further = len(
                    self.integrated_variant_objects_cds_hits[len(self.integrated_variant_objects_cds_hits) - 1].Ref)

        if self.ForwardDirection == TranscriptEnum.FORWARD:
            # firstCDS = self.ListofCDS[0]  # stop in here for reverse transcripts
            # lastCDS = self.ListofCDS[len(self.ListofCDS) - 1]  # start in here for reverse transcripts
            old_last_cds_position = self.last_cds_position
            self.change_last_cds_position(further)
            try:
                seq_to_add = genomehandler.seq(self.Chr, old_last_cds_position + 1, self.last_cds_position)
            except SequenceHandlingError as she:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_ADDITIONAL_INFO_LOG,
                                        str(self.TID) + "\t" + str(she.description))
                if she.sequence_part == "":
                    print("Warning: Out of contig/chrom in " + str(self.TID))
                seq_to_add = she.sequence_part
            self.Complete_CDS += seq_to_add
            self.Rev_CDS += ForTypeSafetyAndStatics.bio_reverse_seq(seq_to_add)
            self.IV_Changed_DNA_CDS_Seq += seq_to_add
            self.iv_check_for_new_variants(genetic_code, stopcodon)
        elif self.ForwardDirection == TranscriptEnum.REVERSE:
            old_last_cds_position = self.last_cds_position
            self.change_last_cds_position(further)
            try:
                seq_to_add = genomehandler.seq(self.Chr, self.last_cds_position, old_last_cds_position - 1)
            except SequenceHandlingError as she:
                LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_ADDITIONAL_INFO_LOG,
                                        str(self.TID) + "\t" + str(she.description))
                if she.sequence_part == "":
                    print("Warning: Out of contig/chrom in " + str(self.TID))
                seq_to_add = she.sequence_part
            self.Complete_CDS = seq_to_add + self.Complete_CDS
            self.Rev_CDS += ForTypeSafetyAndStatics.bio_reverse_seq(seq_to_add)
            self.IV_Changed_DNA_CDS_Seq += ForTypeSafetyAndStatics.bio_reverse_seq(seq_to_add)
            self.iv_check_for_new_variants(genetic_code, stopcodon)
        else:
            LogOrganizer.add_to_log(LogEnums.TRANSCRIPT_BUGHUNTING_LOG,
                                    str(self.TID) + "\tproblem with trancript direction\n")
            pass


class ForTypeSafetyAndStatics:
    """
    This class exist for supporting while programming.
    It is more comfortable to use python, when the IDE knows, which type the current variables are.
    Furthermore it has a variety of methods which are used often.
    """

    @staticmethod
    def bp_to_add_because_of_raster(variant_position_inside_cds: int) -> int:

        """
        :param variant_position_inside_cds: self-explanatory.
        :return: 0,1 or 2 for bp to add to complete the triplet.
        """
        """
        Example		Pos in Codon || bp to add
        (4-1) % 3 -> 0				2
        (5-1) % 3 -> 1				1
        (6-3) % 3 -> 2				0
        """
        to_add_dict = {0: 2, 1: 1, 2: 0}
        to_add = (variant_position_inside_cds - 1) % 3
        to_add = to_add_dict[to_add]
        return to_add

    @staticmethod
    def variant_information_storage_type_safety(vinfo: Variantinformationstorage) -> Variantinformationstorage:
        """
        Only here for a trick: Get a Variant_Information_Storage type in, get a Variant_Information_Storage
        out and python and IDE will "know" its type.
        (So you can use all programming functions like usage, refactor and auto extension.)
        :param vinfo: will be returned without change
        :return: Variant_Information_Storage will be returned without change.
        """
        return vinfo

    @staticmethod
    def transcript_type_safety(transcript: Transcript) -> Transcript:
        """
        Only here for a trick: Get a Transcript type in, get a Transcript
        out and python and IDE will "know" its type.
        (So you can use all programming functions like usage, refactor and auto extension.)
        :param transcript: will be returned without change
        :return: Transcript will be returned without change
        """
        return transcript

    @staticmethod
    def bio_reverse_seq(dna: str):
        """
        It will build a biologically reversed dna strand.
        :param dna: string
        :return: biologically reversed dna string
        """
        revdna = dna.upper()
        revdna = revdna.replace("A", "t")
        revdna = revdna.replace("T", "a")
        revdna = revdna.replace("C", "g")
        revdna = revdna.replace("G", "c")
        revdna = revdna.upper()
        revdna = revdna[::-1]
        return revdna

    @staticmethod
    def not_bio_reverse_seq(dna: str):
        revdna = dna.upper()
        revdna = revdna.replace("A", "t")
        revdna = revdna.replace("T", "a")
        revdna = revdna.replace("C", "g")
        revdna = revdna.replace("G", "c")
        revdna = revdna.upper()
        return revdna

    @staticmethod
    def calculate_raster(variant_position_inside_cds: int) -> int:
        """
        A simple calculation of the raster with the currently position inside the CDS.
        Zero == Start of the triplet
        One == In the middle of the triplet
        Two == Last Position of the triplet
        :param variant_position_inside_cds: self-explanatory.
        :return: 0,1 or 2 for current position inside triplet.
        """
        return (variant_position_inside_cds - 1) % 3

    @staticmethod
    def translation(rna: str, genetic_code: dict):
        """
        Translating the RNA to AA with the given genetic code.
        Will not translate incomplete triplets, like the last 2 bp left.
        "x" for every match, which is not inside the genetic code.
        :param rna: string, lower or upper case
        :param genetic_code: dictionary, all in lowercase
        :return: amino acid string, lower case
        """
        rna = rna.lower()
        amino = ""
        rnalen = len(rna)
        position = 0

        while ((position + 3) <= rnalen):
            try:
                amino += genetic_code[rna[position:(position + 3)]]
            except:
                amino += 'x'
            position += 3
        return amino

    @staticmethod
    def transcription(dna: str):
        """
        Simple replacement from 't' to 'u'.
        :param dna: string, lower or upper case
        :return: string in lowercase
        """
        return (str(dna).lower()).replace("t", "u").upper()

# 1825

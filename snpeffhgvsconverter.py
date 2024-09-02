from Transcript import *
from VCF_Variant import *


class SnpeffHgvsConverter():

    def __init__(self):
        pass

    def convert_main(self, transcript: Transcript, vinfo: Variantinformationstorage) -> str:
        aminodict = {}
        aminodict["Ala"] = "A"
        aminodict["A"] = "Ala"
        aminodict["Cys"] = "C"
        aminodict["C"] = "Cys"
        aminodict["Asp"] = "D"
        aminodict["D"] = "Asp"
        aminodict["Glu"] = "E"
        aminodict["E"] = "Glu"
        aminodict["Phe"] = "F"
        aminodict["F"] = "Phe"
        aminodict["Gly"] = "G"
        aminodict["G"] = "Gly"
        aminodict["His"] = "H"
        aminodict["H"] = "His"
        aminodict["Ile"] = "I"
        aminodict["I"] = "Ile"
        aminodict["Lys"] = "K"
        aminodict["K"] = "Lys"
        aminodict["Leu"] = "L"
        aminodict["L"] = "Leu"
        aminodict["Met"] = "M"
        aminodict["M"] = "Met"
        aminodict["Asn"] = "N"
        aminodict["N"] = "Asn"
        aminodict["Pro"] = "P"
        aminodict["P"] = "Pro"
        aminodict["Gln"] = "Q"
        aminodict["Q"] = "Gln"
        aminodict["Arg"] = "R"
        aminodict["R"] = "Arg"
        aminodict["Ser"] = "S"
        aminodict["S"] = "Ser"
        aminodict["Thr"] = "T"
        aminodict["T"] = "Thr"
        aminodict["Val"] = "V"
        aminodict["V"] = "Val"
        aminodict["Trp"] = "W"
        aminodict["W"] = "Trp"
        aminodict["Tyr"] = "Y"
        aminodict["Y"] = "Tyr"
        aminodict["*"] = "*"
        aminodict['X'] = 'Xaa'
        aminodict['Xaa'] = 'X'

        annotation = ""
        annotation_impact = ""
        gene_name = transcript.TID.split(".")[0]
        gene_id = transcript.TID.split(".")[0]
        feature_type = "transcript"
        feature_id = transcript.TID
        transcript_bio_type = "protein_coding"
        rank = ""
        hgvs_c = ""
        hgvs_p = ""
        c_dna_pos = ""
        c_dna_length = ""
        cds_pos = vinfo.Unchanged_CDS_Position
        cds_length = transcript.uChDNA_length
        aa_pos = (vinfo.Unchanged_CDS_Position - 1) / 3
        aa_length = transcript.uChAA_length
        distance = ""
        errors_warnings = ""
        aa_pos_temp = (vinfo.Unchanged_CDS_Position - 1) % 3
        if aa_pos_temp == 0:
            aa_pos = (vinfo.Unchanged_CDS_Position - 1 + 3) / 3  # aa = 0 is the first, not zero
        elif aa_pos_temp == 1:
            aa_pos = (vinfo.Unchanged_CDS_Position - 1 + 2) / 3
        elif aa_pos_temp == 2:
            aa_pos = (vinfo.Unchanged_CDS_Position - 1 + 1) / 3
        aa_pos = int(aa_pos)
        aa_to_next_stop = -1
        new_aa_pos = 0
        aa_pos_temp = (vinfo.Changed_CDS_Position - 1) % 3
        if aa_pos_temp == 0:
            new_aa_pos = (vinfo.Changed_CDS_Position - 1) / 3 + 1  # aa = 0 is the first, not zero
        elif aa_pos_temp == 1:
            new_aa_pos = (vinfo.Changed_CDS_Position + 1) / 3 + 1
        elif aa_pos_temp == 2:
            new_aa_pos = (vinfo.Changed_CDS_Position) / 3 + 1
        new_aa_pos = int(new_aa_pos)
        aa_to_next_stop = transcript.iv_changed_translation[new_aa_pos:].find('*')

        for i, x in enumerate(vinfo.Classification):
            if i < len(vinfo.Classification) - 1:
                annotation += x.value + ','
            # not efficient, but list isn't large
            else:
                annotation += x.value

        if TranscriptEnum.FRAMESHIFT in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_1 in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_2 in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_1_DEL in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_2_DEL in vinfo.Classification \
                or TranscriptEnum.STOP_GAINED in vinfo.Classification \
                or TranscriptEnum.STOP_LOST in vinfo.Classification \
                or TranscriptEnum.START_LOST in vinfo.Classification:
            annotation_impact = "HIGH"
        elif TranscriptEnum.DELETION in vinfo.Classification \
                or TranscriptEnum.INSERTION in vinfo.Classification \
                or TranscriptEnum.AA_CHANGE in vinfo.Classification:
            annotation_impact = "MODERATE"

        elif TranscriptEnum.STOP_CHANGED in vinfo.Classification \
                or TranscriptEnum.AA_CHANGE not in vinfo.Classification:
            annotation_impact = "LOW"
        else:
            print("Impossible case. Error convert_main.")
        # HGVS_C
        if TranscriptEnum.SUBSTITUTION in vinfo.Classification:
            if TranscriptEnum.FORWARD == transcript.ForwardDirection:
                hgvs_c = "c." + str(vinfo.Unchanged_CDS_Position) + vinfo.Ref + ">" + vinfo.Alt
            else:
                hgvs_c = "c." + str(vinfo.Unchanged_CDS_Position) + vinfo.ReverseRef + ">" + vinfo.ReverseAlt
        elif TranscriptEnum.INSERTION in vinfo.Classification:
            if vinfo.Unchanged_CDS_Position >= len(
                    vinfo.Alt) - 1:
                if transcript.ForwardDirection == TranscriptEnum.FORWARD:
                    insert = vinfo.Alt[1:]
                else:
                    insert = vinfo.ReverseAlt[:len(vinfo.ReverseAlt) - 1]
                if len(transcript.uChDNAsequence) >= len(transcript.Complete_CDS):
                    refpart = transcript.uChDNAsequence[
                              vinfo.Unchanged_CDS_Position - len(insert):vinfo.Unchanged_CDS_Position + 1]
                else:
                    if TranscriptEnum.FORWARD == transcript.ForwardDirection:
                        refpart = transcript.Complete_CDS[
                                  vinfo.Unchanged_CDS_Position + 1 - len(insert):vinfo.Unchanged_CDS_Position + 1]
                    else:
                        refpart = transcript.Rev_CDS[
                                  vinfo.Unchanged_CDS_Position + 1 - len(insert):vinfo.Unchanged_CDS_Position + 1]
                hgvs_c = "c."
                if insert == refpart:
                    vinfo.Classification.append(HgvsDna.DUP)
                    hgvs_c += str(vinfo.Unchanged_CDS_Position - len(insert))
                    if len(insert) == 1:
                        hgvs_c += 'dup' + insert
                    else:
                        hgvs_c += '_' + str(vinfo.Unchanged_CDS_Position) + 'dup' + insert
                else:
                    hgvs_c += str(vinfo.Unchanged_CDS_Position - 1) + '_' + str(
                        vinfo.Unchanged_CDS_Position) + 'ins' + insert
        elif TranscriptEnum.DELETION in vinfo.Classification:
            hgvs_c = "c."
            if len(vinfo.Ref) == 2:
                hgvs_c += str(vinfo.Unchanged_CDS_Position + 1) + 'del' + vinfo.Ref[1]
            else:
                hgvs_c += str(vinfo.Unchanged_CDS_Position + 1) + '_' + str(
                    vinfo.Unchanged_CDS_Position + len(vinfo.Ref)) + 'del' + vinfo.Ref[1:]
        else:
            LogOrganizer.add_to_log(LogEnums.CONVERTER_LOG,
                                    "No Classification in: " + str(transcript.TID) + "\t" + str(vinfo.ChrPosition))
        # HGVS_P
        hgvs_p = "p."
        new_amino = ''
        old_amino = ''
        for aa in vinfo.NewAmino:
            new_amino += aminodict[aa.upper()]
        for aa in vinfo.OrigAmino:
            old_amino += aminodict[aa.upper()]

        if TranscriptEnum.STOP_LOST in vinfo.Classification:
            if aa_to_next_stop != -1:
                ext = 'ext*' + str(aa_to_next_stop)
            else:
                ext = 'ext*?'
            hgvs_p += '*' + str(aa_pos) + new_amino + ext
        elif TranscriptEnum.STOP_CHANGED in vinfo.Classification and TranscriptEnum.INSERTION in vinfo.Classification and \
                vinfo.NewAmino[0] != '*':
            if aa_to_next_stop != -1:
                ext = 'ext*' + str(aa_to_next_stop)
            else:
                ext = 'ext*?'
            hgvs_p += '*' + str(aa_pos) + new_amino + ext

        elif TranscriptEnum.SUBSTITUTION in vinfo.Classification or new_amino == '*':
            if vinfo.OrigAmino == vinfo.NewAmino:  # silent
                hgvs_p += old_amino + str(aa_pos) + old_amino
            else:
                hgvs_p += old_amino + str(aa_pos) + new_amino
        elif TranscriptEnum.FRAMESHIFT_1 in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_2 in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_1_DEL in vinfo.Classification \
                or TranscriptEnum.FRAMESHIFT_2_DEL in vinfo.Classification:
            if new_amino == '*':
                hgvs_p += old_amino + str(aa_pos) + '*'
            else:
                # first changed aa has to be the first aa here
                if vinfo.OrigAmino[0] != vinfo.NewAmino[0]:
                    if vinfo.STOP_CAUSED_IN == -1:
                        hgvs_p += old_amino + str(aa_pos) + new_amino + '*' + '?'
                    else:
                        hgvs_p += old_amino + str(aa_pos) + new_amino + '*' + str(vinfo.STOP_CAUSED_IN)
                else:
                    i = 0
                    for j, old_aa in enumerate(transcript.uChAAsequence[aa_pos - 1:]):
                        new_aa = transcript.iv_changed_translation[new_aa_pos - 1 + j]
                        if old_aa == "":
                            print("should not happen")
                        elif new_aa == '*':
                            hgvs_p += old_amino + str(aa_pos + i) + '*'
                            break
                        elif old_aa == new_aa:
                            i += 1
                            continue
                        else:
                            # find better names, here it is to switch single aa-code to 3 letter aa-code
                            abc1 = ""
                            abc2 = ""
                            for abcabc in new_aa:
                                abc1 += aminodict[abcabc.upper()]
                            for abcabc in old_aa:
                                abc2 += aminodict[abcabc.upper()]
                            if vinfo.STOP_CAUSED_IN == -1:
                                hgvs_p += abc2 + str(aa_pos + i) + abc1 + '*?'
                            else:
                                hgvs_p += abc2 + str(aa_pos + i) + abc1 + '*' + str(vinfo.STOP_CAUSED_IN - i)
                            break
        elif TranscriptEnum.INSERTION in vinfo.Classification:
            if ((len(vinfo.Alt) - 1) % 3) != 0:
                print('how`?')
            else:
                # dup check
                insert_aa = vinfo.NewAmino
                if (len(vinfo.Alt) - 1) % 3 != 0:
                    print("how2?")
                refpart_aa = ""
                if len(transcript.uChAAsequence) >= len(transcript.IV_OriginalTranslation):
                    if (len(vinfo.Alt) - 1) / 3 == len(vinfo.NewAmino):
                        # Insertion of complete AA without interfering with other AA
                        refpart_aa = transcript.uChAAsequence[aa_pos - 1 - len(vinfo.NewAmino):aa_pos]
                    else:
                        # Insertion with interfering with other AA
                        # -> the original AA is in vinfo.NewAmino
                        # -> do not compare it with itself
                        if len(vinfo.NewAmino) == 1:
                            print('how?3')
                        refpart_aa = transcript.uChAAsequence[aa_pos - 1 - (len(vinfo.NewAmino) - 1):aa_pos]
                else:
                    if (len(vinfo.Alt) - 1) / 3 == len(vinfo.NewAmino):
                        # Insertion of complete AA without interfering with other AA
                        refpart_aa = transcript.IV_OriginalTranslation[aa_pos - 1 - len(vinfo.NewAmino):aa_pos]
                    else:
                        # Insertion with interfering with other AA
                        # -> the original AA is in vinfo.NewAmino
                        # -> do not compare it with itself
                        if len(vinfo.NewAmino) == 1:
                            print('how?3')
                        refpart_aa = transcript.IV_OriginalTranslation[aa_pos - 1 - (len(vinfo.NewAmino) - 1):aa_pos]
                if len(vinfo.NewAmino) == len(refpart_aa):
                    dup = True
                    for i, oa in enumerate(refpart_aa):
                        if vinfo.NewAmino[i] == oa:
                            continue
                        else:
                            dup = False
                    if dup:
                        if (len(vinfo.Alt) - 1) / 3 == len(vinfo.NewAmino):
                            # Insertion of complete AA without interfering with other AA
                            if len(vinfo.NewAmino) == 1:  # one AA
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'dup'
                            else:  # multiple AA
                                hgvs_p += aminodict[refpart_aa[0].upper()] + str(
                                    aa_pos - (len(vinfo.NewAmino) - 1)) + '_' \
                                          + aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'dup'
                        else:
                            # Insertion with interfering with other AA
                            if len(vinfo.NewAmino) == 2:  # one new AA, first AA not changed, because its a dup
                                # so it looks like, this can only happen, if vinfo.NewAmino[0] == vinfo.NewAmino[1]
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'dup'
                            else:  # multiple AA
                                hgvs_p += aminodict[refpart_aa[0].upper()] + str(
                                    aa_pos - (len(vinfo.NewAmino) - 1)) + '_' + \
                                          aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'dup'
                    else:
                        # no dup, but insertion without frameshift
                        # delins are possible
                        if (len(vinfo.Alt) - 1) / 3 == len(vinfo.NewAmino):
                            # Insertion of complete AA without interfering with other AA
                            hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + '_' + aminodict[
                                refpart_aa[1].upper()] \
                                      + str(aa_pos + 1) + 'ins' + new_amino
                        else:
                            if vinfo.OrigAmino[0] == vinfo.NewAmino[0]:  # no delin
                                # print(vinfo.ChrPosition)
                                # print(vinfo.OrigAmino[0] + str(AA_pos) + '_')
                                # print(refpartAA)
                                # print(transcript.IV_OriginalTranslation[AA_pos:])
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + '_' + aminodict[
                                    refpart_aa[1].upper()] \
                                          + str(aa_pos + 1) + 'ins' + new_amino
                            else:
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'delins' + new_amino
                else:
                    if len(vinfo.NewAmino) >= aa_pos:
                        # no duplication/ refAA possible, because insertion is to big
                        if vinfo.OrigRaster == 2:
                            # insertion without interfering with first origAA
                            if len(vinfo.NewAmino) > 1:
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'ins' + new_amino[1:]
                            else:
                                print("curios effect")
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'ins' + new_amino
                        else:
                            # insertion with interfering with first origAA
                            if vinfo.OrigAmino[0] == vinfo.NewAmino[0]:
                                # first AA not changed -> insertion after
                                hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'ins' + new_amino
                            else:
                                # first AA changed -> SUB or DELINS
                                if len(vinfo.NewAmino) == 1 and vinfo.NewAmino == '*':
                                    hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + new_amino
                                else:
                                    # DELINS
                                    hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'delins' + new_amino
                    else:
                        # stopcodon + stuff at the end of the transcript
                        if vinfo.NewAmino == vinfo.OrigAmino:
                            b = ""
                            for aa in vinfo.OrigAmino:
                                b += aminodict[aa.upper()]
                            hgvs_p += b + str(aa_pos) + b
                        else:
                            if vinfo.NewAmino[0] == vinfo.OrigAmino[0]:
                                # no sub, possible here
                                if len(vinfo.NewAmino) > 1 and len(vinfo.OrigAmino) > 1:
                                    # delins
                                    b = ""
                                    for aa in vinfo.OrigAmino:
                                        b += aminodict[aa.upper()]
                                    hgvs_p += b + str(aa_pos + 1) + 'delins' + new_amino[1:]
                                elif len(vinfo.OrigAmino) == 1 and len(vinfo.NewAmino) > 1:
                                    # here, if the AA 2+ are not identical with the original AA
                                    hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + 'delins' + new_amino
                                else:  # orig > 1 and newamino == 1 () impossible?
                                    print('Impossible case triggered.')
                            else:
                                b = ""
                                for aa in vinfo.OrigAmino:
                                    b += aminodict[aa.upper()]
                                hgvs_p += b + str(aa_pos) + 'delins' + new_amino
        elif TranscriptEnum.DELETION in vinfo.Classification:
            # deletionstuff, no frameshifts possible anymore -> codon-position should always be 2 ???? nope
            # len(vinfo.NewAmino) > 1: true
            if vinfo.OrigRaster == 2:
                # deletion of complete AA, without interfering of another AA (startposition +1 (dna))
                if len(vinfo.OrigAmino) == 2:
                    # complete deletion of the second AA
                    hgvs_p += aminodict[vinfo.OrigAmino[1].upper()] + str(aa_pos + 1) + 'del'
                elif len(vinfo.OrigAmino) == 1:
                    # impossible case (or somethi8ng went wrong with the transcript)
                    if vinfo.OrigAmino == vinfo.NewAmino:
                        # sub
                        # NP_003997.1:p.Cys188=
                        hgvs_p += old_amino + str(aa_pos) + old_amino
                    else:
                        try:
                            hgvs_p += old_amino + str(aa_pos) + aminodict[vinfo.NewAmino[0].upper(9)]
                        except IndexError:
                            print("Variant error: " + str(vinfo.ChrPosition))
                            hgvs_p += 'oldamino' + str(aa_pos) + 'Xaa'
                else:
                    hgvs_p += aminodict[vinfo.OrigAmino[1].upper()] + str(aa_pos + 1) + '_' + \
                              aminodict[vinfo.OrigAmino[len(vinfo.OrigAmino) - 1].upper()] + str(
                        aa_pos + 1 + len(vinfo.OrigAmino) - 1) + 'del'

            else:
                if vinfo.OrigAmino == "":
                    print("Variant error: " + str(vinfo.ChrPosition))
                elif len(vinfo.OrigAmino) == 1:
                    # this should not happen here, but who knows
                    if vinfo.OrigAmino == vinfo.NewAmino:
                        # sub
                        # NP_003997.1:p.Cys188=
                        hgvs_p += old_amino + str(aa_pos) + old_amino
                    else:
                        hgvs_p += old_amino + str(aa_pos) + aminodict[vinfo.NewAmino.upper()]
                elif TranscriptEnum.STOP_GAINED in vinfo.Classification:
                    hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + aminodict[vinfo.NewAmino.upper()]
                else:
                    hgvs_p += aminodict[vinfo.OrigAmino[0].upper()] + str(aa_pos) + "_" + aminodict[
                        vinfo.OrigAmino[len(vinfo.OrigAmino) - 1].upper()] + str(
                        aa_pos + len(vinfo.OrigAmino)) + 'delins'
        elif TranscriptEnum.UNKNOWN_AMINOACID in vinfo.Classification:
            # this here should be a substitution of X into X (because bad chromosome data)
            hgvs_p += old_amino + str(aa_pos) + new_amino
        else:
            LogOrganizer.add_to_log(LogEnums.CONVERTER_LOG,
                                    "No Classification in: " + str(transcript.TID) + "\t" + str(vinfo.ChrPosition))
        snpeff_like_info_string = 'NAV2=' + vinfo.Alt + '|' \
                                  + annotation + "|" \
                                  + annotation_impact + "|" \
                                  + gene_name + "|" \
                                  + gene_id + "|" \
                                  + feature_type + "|" \
                                  + feature_id + "|" \
                                  + transcript_bio_type + "|" \
                                  + rank + "|" \
                                  + hgvs_c + "|" \
                                  + hgvs_p + "|" \
                                  + str(c_dna_pos) + "/" \
                                  + str(c_dna_length) + "|" \
                                  + str(cds_pos) + "/" \
                                  + str(cds_length) + "|" \
                                  + str(aa_pos) + "/" \
                                  + str(aa_length) + "|" \
                                  + distance + "|" \
                                  + errors_warnings + ";"
        return snpeff_like_info_string


@unique
class HgvsDna(Enum):
    SUB = "Substitution"
    DEL = "Deletion"
    DUP = "Duplication"
    INS = "Insertion"
    INV = "Inversion"
    CON = "Conversion"
    DELIN = "Deletion-Insertion"
    ALL = "Alleles"
    REP = "Repeated sequences"
    COM = "Complex"


@unique
class HgvsProt(Enum):
    SUB = "Substitution"
    DEL = "Deletion"
    DUP = "Duplication"
    INS = "Insertion"
    DELIN = "Deletion-Insertion"
    ALL = "Alleles"
    REP = "Repeated sequences"
    FS = "Frame shift"
    EX = "Extension"

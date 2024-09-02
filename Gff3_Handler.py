__author__ = "Jan-Simon Baasner"
__email__ = "janbaas@cebitec.uni-bielefeld.de"

import Transcript
import sys


class Gff3HandlerV3:
    def __init__(self, gff3_data_path: str):
        self.gff3_ID_Dict = {}
        self.dict_gff_for_parents = {}
        self.dict_Chr_Names = {}
        self.dict_Chr_dict_Transcript = {}
        self.List_Of_Transcripts = {}
        self.nr_chroms = -1
        self.transIndex = -1

        self.read_gff3(gff3_data_path)
        self.create_transcripts()
        self.sort_transcripts()

    def read_gff3(self, gff3_data_path: str):
        gff3 = open(gff3_data_path, 'r')
        lines = gff3.readlines()
        gff3.close()
        count_generic_id = 0
        for line in lines:
            if line.startswith('#') or len(line) <= 2:
                continue
            spline = line.split('\t')
            seqid = spline[0]
            if seqid not in self.dict_Chr_Names.keys():
                self.nr_chroms += 1
                self.dict_Chr_Names[seqid] = self.nr_chroms
                self.dict_Chr_Names[self.nr_chroms] = seqid
                self.List_Of_Transcripts[seqid] = []

            atts = spline[8].split(';')
            gff3_id = ""
            parent = ""
            if "ID=" not in spline[8]:
                generic_id = str(spline[0]) + '_' + str(spline[1]) + '_' + str(spline[3]) + '_' + str(spline[4])
            for att in atts:
                if 'ID=' in att:
                    gff3_id = att[3:].replace("\n", "")
                elif 'Parent=' in att:
                    parent = att[7:].replace("\n", "")
                    if ',' in parent:
                        parent = parent.split(',')
                    else:
                        parent = [parent]
            if gff3_id == "":
                count_generic_id += 1
                gff3_id = generic_id
            try:
                self.gff3_ID_Dict[seqid, gff3_id].append(spline)
            except KeyError:
                self.gff3_ID_Dict[seqid, gff3_id] = [spline]
            if len(parent) == 0:
                # when its a gene entry
                continue
            for pp in parent:
                try:
                    self.dict_gff_for_parents[seqid, pp].append(spline)
                except KeyError:
                    self.dict_gff_for_parents[seqid, pp] = [spline]

    def create_transcripts(self):
        for key, splinelist in self.gff3_ID_Dict.items():
            for spline in splinelist:
                if spline[2] != "gene":
                    continue
                gene_seqid = spline[0]
                gene_gff3_id = key[1]
                gene_info_string = spline[8]
                if gene_seqid not in self.dict_Chr_dict_Transcript:
                    self.dict_Chr_dict_Transcript[gene_seqid] = {}
                for spline_child in self.dict_gff_for_parents[gene_seqid, gene_gff3_id]:
                    chr = spline_child[0]
                    gfftype = spline_child[2]
                    if 'mRNA' not in gfftype:
                        continue
                    start = int(spline_child[3])
                    end = int(spline_child[4])
                    if spline_child[6] == '+':
                        strand = Transcript.TranscriptEnum.FORWARD
                    elif spline_child[6] == '-':
                        strand = Transcript.TranscriptEnum.REVERSE
                    else:
                        strand = Transcript.TranscriptEnum.UNKNOWN_STRAND_DIRECTION
                    phase = spline_child[7]
                    atts = spline_child[8].split(';')
                    tid = ""
                    for att in atts:
                        if 'ID=' in att:
                            tid = att[3:]
                            break
                    if tid == "":
                        print("No ID:\n" + spline_child)
                        sys.exit()
                    transcript = Transcript.Transcript(self.get_next_transcript_index(), tid, start, end, strand, chr)
                    cdslist = []
                    for rna_child_spline in self.dict_gff_for_parents[gene_seqid, tid]:
                        if "CDS" in rna_child_spline[2]:
                            cdslist.append((int(rna_child_spline[3]), int(rna_child_spline[4]), rna_child_spline[7]))
                        elif "exon" in rna_child_spline[2]:
                            transcript.add_exon_descriptin("\t".join(rna_child_spline))
                        elif "utr" in rna_child_spline[2]:
                            transcript.add_utr_description("\t".join(rna_child_spline))
                    cdslist = sorted(cdslist)
                    for cds in cdslist:
                        transcript.add_cds(cds[0], cds[1], cds[2])
                    transcript.set_gene_info_string(gene_info_string)
                    self.dict_Chr_dict_Transcript[gene_seqid][transcript.IndexKey] = transcript
                    self.List_Of_Transcripts[gene_seqid].append(transcript)

    def get_next_transcript_index(self):
        """
		For creating new transcripts (need for a new ID).
		:return: Next integer ID.
		"""
        self.transIndex += 1
        return self.transIndex

    def get_chromosome_names(self) -> list:
        """
		Returns all names of the chromosomes inside a list.
		:return: List of all chromosome names.
		"""
        count_names = len(self.dict_Chr_Names) / 2
        i = 0
        name_list = []
        while count_names != i:
            name_list.append(self.dict_Chr_Names[i])
            i += 1
        if name_list[0] == []:
            name_list.pop(0)

        return name_list

    def get_chromosome_id(self, chr_name: str) -> int:
        return self.dict_Chr_Names[chr_name]

    def get_chr_transcripts_dict(self, chr_name: str):
        """
        Returns a dictionary with all transcripts from the choosen chromosome.
        :param chr_name: Name of the chromosome.
        :return: Dictionary with all transcripts inside this chromosome.
        """
        return self.dict_Chr_dict_Transcript[chr_name]

    def free_ram(self, chr_name: str):
        """
        Frees the RAM.
        :param chr_name:
        :return:
        """
        self.dict_Chr_dict_Transcript[chr_name] = []
        self.List_Of_Transcripts[chr_name] = []

    def get_chr_transcript_list(self, chr_name: str):
        return self.List_Of_Transcripts[chr_name]

    def sort_transcripts(self):
        for chrName in self.List_Of_Transcripts.keys():
            self.List_Of_Transcripts[chrName] = sorted(self.List_Of_Transcripts[chrName],
                                                       key=lambda s_transcript: s_transcript.StartOfRNA)

    def update_transripts(self, transcript_list: list, chr_name: str):
        # i = self.dict_Chr_Names[ChrName]
        # self.List_Of_Transcripts[i] = transcriptList
        self.List_Of_Transcripts[chr_name] = sorted(transcript_list, key=lambda s_transcript: s_transcript.StartOfRNA)

    def add_new_transcript_to_dict(self, chr_name: str, transcript: Transcript):
        self.dict_Chr_dict_Transcript[chr_name][transcript.TID] = transcript

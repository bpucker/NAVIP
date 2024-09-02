__author__ = "Jan-Simon Baasner"
__email__ = "janbaas@cebitec.uni-bielefeld.de"

from enum import Enum, unique


class Variant:
    """
    A simple storage class for a vcf variant.
    Contains all its information.
    """

    def __init__(self,
                 chromosome: str,
                 position: int,
                 id: int,  # from vcf-file
                 usefull_id: int,  # not from vcf-file (count variants)
                 reference: str,
                 alternate: str,
                 qual: str,
                 filter_: str,
                 info: str):
        """
        Initialization of all information.
        Contains all information from the vcf file and one new ID for variant handling.
        :param chromosome: Chromosome.
        :param position: Position.
        :param id: VCF ID
        :param usefull_id: A real ID in numbers, which is unique.
        :param reference: VCF REF
        :param alternate: VCF ALT
        :param qual: VCF QUAL
        :param filter_: VCF Filter
        :param info: VCF Info
        """

        # Initialization of variables
        self.Chromosome = chromosome
        self.Position = position
        self.ID = id
        self.UsefullID = usefull_id
        self.Reference = reference
        self.Alternate = alternate
        self.Qual = qual
        self.Filter = filter_
        self.Info = info
        self.ListofTranscripts = []
        self.SListofTranscripts = []


class VariantNavip(Variant):

    def __init__(self, variation: Variant):
        super().__init__(variation.Chromosome, variation.Position, variation.ID, variation.UsefullID,
                         variation.Reference, variation.Alternate, variation.Qual, variation.Filter, variation.Info)
        navip_line = self.Info.split("|NAVIP_END|")[0].split("|")
        self.TranskriptID = navip_line[0]
        self.Strand_Direction = navip_line[1]
        self.Variant_Classification = navip_line[2]
        self.Shared_Eff_Keys = navip_line[3]
        self.Ref_Codon = navip_line[4].split(";")[0]
        self.Ref_in_Codon_Position = navip_line[4].split(";")[1]
        self.Ref_AA = navip_line[5]
        self.Ref_Codon_CDS_Position = navip_line[6]
        self.Alt_Codon = navip_line[7].split(";")[0]
        self.Alt_in_Codon_Position = navip_line[7].split(";")[1]
        self.Alt_AA = navip_line[8]
        self.Alt_Codon_CDS_Position = navip_line[9]


@unique
class VariantEnum(Enum):
    """
    Contains a few enums for variant usage only.
    """
    DEFAULT = -1  # nothing decided (for finding errors - default should not exist in the output)
    NO_RASTER_NOW = -2
    NO_StartOfOwnEffect = -3
    NO_EndOfOwnEffect = -4

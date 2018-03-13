__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"



from enum import Enum, unique

class GenomHandler():
    """
    This class reads the given fasta file, stores the chromosome data
    and contains a few functions to make the sequences available.
    """

    #JustXChromosomes --> 0 == all
    def __init__(self, Fasta_File_Path: str, JustXChromosomes: int):
        """
        The given fasta file will be read and the names and sequences of the
        chromosomes will be stored.
        It is possible to read only x chromosomes in the fasta file.
        :param Fasta_File_Path: Path to the fasta file including the filename.
        :param JustXChromosomes: The number of the first X chromosomes to read. 0 for all chromosomes available.
        """

        self.originChromosomeData = [[]]
        self.originChromosomeNames = []
        self.dictChrName = {}
        countChr = 0
        with open(Fasta_File_Path, "r") as DataFile:
            line = DataFile.readline()
            self.originChromosomeNames.append(((line.split(" ")[0])[1:]).replace("\n",""))
            line = DataFile.readline()
            lines = []
            while (line):
                if (line.startswith(">")):
                    self.originChromosomeData.append("".join(lines).replace("\n", ""))
                    lines = []
                    if (self.originChromosomeData[0] == []):
                        self.originChromosomeData.pop(0)

                    countChr += 1
                    if (countChr >= JustXChromosomes and (JustXChromosomes != 0)):
                        break

                    self.originChromosomeNames.append(((line.split(" ")[0])[1:]).replace("\n",""))

                elif line.startswith(";"):
                    print("Comment is not used: " + line)
                else:
                    lines.append(line)
                # print lines
                line = DataFile.readline()
        if countChr < JustXChromosomes and (JustXChromosomes != 0):
            finallines = "".join(lines).replace("\n", "")
            self.originChromosomeData.append(finallines)
        elif JustXChromosomes == 0:
            finallines = "".join(lines).replace("\n", "")
            self.originChromosomeData.append(finallines)
        # originChromosomeData.append(lines)
        DataFile.close()
        someint = 0
        for name in self.originChromosomeNames:
            self.dictChrName[name] = someint
            self.dictChrName[someint] = name
            someint += 1

    def GetChromosomeNames(self)-> list:
        """
        Returns the names of all chromosomes.
        :return: A list containing all chromosome names.
        """
        countNames = len(self.dictChrName) / 2
        i = 0
        NameList = []
        while countNames != i:
            NameList.append(self.dictChrName[i])
            i += 1
        if NameList[0] == []:
            NameList.pop(0)

        return NameList

    #position in chromosome - not in string (stringposition -1)
    def singleSeq(self,ChrName: str, Position: int ) -> str :
        """
        Returns a single base from the chromosome.
        Note: The position in the chromosome is the position-1 inside the string.
        This will be considered inside this function.
        :param ChrName: Name of the wanted chromosome.
        :param Position: Position inside the chromosome (not inside the string!).
        :return: One single Base.
        """
        return (self.originChromosomeData[self.dictChrName[ChrName]])[Position-1]

    # position in chromosome - not in string (stringposition -1 for start)
    def seq(self,ChrName: str, StartPosition: int , EndPosition: int) -> str :
        """
        Returns a sequence or a single base from the wanted chromosome and position.
        The output will be a single base, if the start und end position are the same.
        :param ChrName: The wanted chromosome.
        :param StartPosition: The start position inside the chromosome (not inside the string!)
        :param EndPosition: The end position inside the chromosome (not inside the string!)
        :return: A sequence of bases or a single base.
        """
        try:
            a = self.dictChrName[ChrName]
            if StartPosition != EndPosition:
                return (self.originChromosomeData[self.dictChrName[ChrName]])[StartPosition-1:EndPosition]
            else:
                return self.singleSeq(ChrName, StartPosition)
        except KeyError as e:
            return []

    def GetChromosome (self, ChrName:str)  -> str:
        """
        This function will deliver the DNA sequence of the wanted chromosome.
        :param ChrName: Name of the chromosome.
        :return: The complete sequence.
        """
        return self.originChromosomeData[self.dictChrName[ChrName]]

@unique
class Fasta_Enum(Enum):
    """
    This class contains énums for transcript description and usage only.
    """
    ODNA = "oDNA"   #original DNA
    OAA = "oAA"     #original AA sequence
    nDNA = "nDNA"   #new DNA
    nAA = "nAA"     #new AA sequence

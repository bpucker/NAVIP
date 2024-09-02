__author__ = "Jan-Simon Baasner"
__email__ = "janbaas@cebitec.uni-bielefeld.de"

from enum import Enum, unique


@unique
class LogEnums(Enum):
    VCF_FORMAT_CHECK_LOG = "VCF_FORMAT_CHECK_LOG".lower()
    TRANSCRIPT_CDS_CREATION_LOG = "TRANSCRIPT_CDS_CREATION_LOG".lower()
    TRANSCRIPT_BUGHUNTING_LOG = "TRANSCRIPT_BUGHUNTING_LOG".lower()
    TRANSCRIPT_ADDITIONAL_INFO_LOG = "TRANSCRIPT_ADDITIONAL_INFO_LOG".lower()
    COORDINATOR_COMPLETE_CHECK_LOG = "COORDINATOR_COMPLETE_CHECK_LOG".lower()
    COORDINATOR_COMPLETE_CHECK_CRITICAL_LOG = "COORDINATOR_COMPLETE_CHECK_CRITICAL_LOG".lower()
    COORDINATOR_FASTA_FILE_ERROR_LOG = "COORDINATOR_FASTA_FILE_ERROR_LOG".lower()
    COORDINATOR_BUGHUNTING_LOG = "COORDINATOR_BUGHUNTING_LOG".lower()
    COORDINATOR_PHASING_LOG = "COORDINATOR_PHASING_LOG".lower()
    COORDINATOR_VARIANT_LOG = "COORDINATOR_VARIANT_LOG".lower()
    COORDINATOR_TRANSCRIPT_LOG = "COORDINATOR_TRANSCRIPT_LOG".lower()
    CONVERTER_LOG = "CONVERTER_LOG".lower()


class LogOrganizer:
    log = {}

    @staticmethod
    def add_to_log(log_name: LogEnums, text: str):
        if "\n" not in text:
            text += "\n"
        if str(log_name.value) in LogOrganizer.log:
            LogOrganizer.log[str(log_name.value)].append(text)
        else:
            LogOrganizer.log[str(log_name.value)] = [text]

    @staticmethod
    def write_all_logs(outpath: str):
        for key in LogOrganizer.log.keys():
            logdata = open(str(outpath) + str(key) + ".log", "w")
            logdata.write("".join(LogOrganizer.log[key]))
            logdata.close()

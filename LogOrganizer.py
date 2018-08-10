__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


from enum import Enum, unique


@unique
class LogEnums (Enum):
	VCF_Format_Check_log = "VCF_Format_Check_log"

class LogOrganizer:

	log = {}

	@staticmethod
	def addToLog(LogName: LogEnums, text: str):
		if str(LogName.value) in LogOrganizer.log:
			LogOrganizer.log[str(LogName.value)].append(text)
		else:
			LogOrganizer.log[str(LogName.value)] = text

	@staticmethod
	def writeLog(LogName: LogEnums, outpath:str):
		logdata = open(str(outpath)+str(LogName.value)+".log","w")
		logdata.write("".join(LogOrganizer.log[str(LogName.value)]))
		logdata.close()



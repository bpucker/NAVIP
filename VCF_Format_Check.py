__author__  = "Jan-Simon Baasner"
__email__   = "janbaas@cebitec.uni-bielefeld.de"


from datetime import datetime
import LogOrganizer
import sys


def VCF_Check(vcf_path_amd_file: str, outpath:str):
	print ("Checking for compatibily.")
	time = datetime.now()
	compatible = True

	vcf = open(vcf_path_amd_file, "r")
	line = vcf.readline()

	#first init
	chr = ""
	pos = 0  # very important
	id = ""
	Ref = ""
	Alternate = ""
	Qual = ""
	Filter = ""
	Info = ""

	lid = -1
	intlist = []
	while line:
		lid += 1
		if line.startswith("#"):
			line = vcf.readline()
			continue

		splitline = line.split("\t")
		try:
			# Thats how it works inside the Handler
			# Now[0],  # Chr
			# int(Now[1]),  # Pos
			# Now[2],  # ID
			# Now[3],  # Ref
			# Now[4],  # Alternate
			# Now[5],  # Qual
			# Now[6],  # Filter
			# Now[7]  # Info

			# test for important cases
			chr = str(splitline[0])
			pos = int(splitline[1]) #very important
			intlist.append((chr,pos))
			id =  str(splitline[2])
			Ref = str(splitline[3])
			Alternate = str(splitline[4])
			Qual = str(splitline[5])
			Filter = str(splitline[6])
			Info = str(splitline[7])
		except:
			print("Incompatible in line " + str(lid))
			print("Line: " + str(line))
			print("Don't forget the '#' for info-lines. All others: STR, INT, STR, STR, STR, STR, STR, STR.")
			print("Use tabs between the entrys - entrys without 'values' simply have a point: '.'")
			e = sys.exc_info()[0]
			print(e)
			compatible = False
			break
		line = vcf.readline()

	for i, tuple in enumerate(intlist):
		chr = tuple[0]
		pos = tuple[1]
		if i == 0:
			continue
		if intlist[i-1][1] > pos:
			if intlist[i - 1][0] == chr:
				logdata = "Position: " + str(intlist[i-1]) + " has to be after position: " + str(intlist[i]) + "\n"
				LogOrganizer.LogOrganizer.addToLog(LogOrganizer.LogEnums.VCF_Format_Check_log, logdata)
				compatible = False
		if intlist[i-1][1] == pos:
			if intlist[i-1][0] == chr:
				logdata = "There is more than one entry in: " + str(pos) + " the program still (probably) works, but will only consider one entry.\n"
				LogOrganizer.LogOrganizer.addToLog(LogOrganizer.LogEnums.VCF_Format_Check_log, logdata)
				compatible = False



	print("Check done in: " + str(datetime.now() -time))
	if compatible:
		print("Looks like it is compatible.")
	else:
		print("Looks like it is incompatible.")
	vcf.close()
	LogOrganizer.LogOrganizer.writeLog(LogOrganizer.LogEnums.VCF_Format_Check_log, outpath)






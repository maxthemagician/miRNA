import sys
import os.path

# substrings mirnas starting at position 2
def extractMiRNAs(writer, inputFile, i):
	inputFile = open(inputFile, "r")
	for line in inputFile:
		if line[0] == ">":
			writer.write(line)
		else:
			writer.write(line[1:i+1] + "\n")
	inputFile.close()


if __name__ == '__main__':
			
	
	# create mature-miRNA-substrings of length 5-8
	for i in range(5, 9):
		fileCut = open( os.path.join(sys.argv[1][0:-3]+"_"+str(i)+".fa"), "w")
		extractMiRNAs(fileCut, sys.argv[1], i)
		fileCut.close()
	print("Done")
	
	

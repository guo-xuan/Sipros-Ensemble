 # -*- coding: utf-8 -*-

 
import getopt, sys 
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator

	
FT2FileName = sys.argv[1] 	
FT2file     = open(FT2FileName)
outputpath  = sys.argv[2]
filenum     = int (sys.argv[3]) 

(pathRoot, pathExt) = os.path.splitext(FT2FileName)

FT2subfiles = []

for i in range (filenum) :
	#os.popen("mkdir "+outputpath+str(i))
	(FT2FileNameRoot, FT2FileNameExt) = os.path.splitext(FT2FileName)
	curFT2subfile = open(outputpath+os.sep+os.path.basename(FT2FileNameRoot)+"_subfile"+str(i)+"." + pathExt, "w")
	FT2subfiles.append(curFT2subfile)


scannum = 0
fileid  = -1
	
for eachline in FT2file :
	eachline = eachline.replace('\r','') #remove end line symbol
	eachline = eachline.replace('\n','') #remove end line symbol
	if (eachline.startswith("H")):
		continue
	if (eachline.startswith("S")):
		scannum = scannum + 1
		fileid  = scannum % filenum
	FT2subfiles[fileid].write(eachline+"\n")


	
# close files	
for i in range (filenum) :
	FT2subfiles[i].close()
		
FT2file.close()

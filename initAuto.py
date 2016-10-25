
#! /usr/bin/python

import os, time
import sys
import fileinput
from glob import glob
import subprocess

fUser = os.getenv("USER")
sampledPars 		= "/afs/ipp-garching.mpg.de/home/l/lscyboz/mcMerge/"
InputFileNameGen 	= "tT_matchbox_NLO.run"
SettingsFolder    	= "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/"
GenericInputFileLO 	= "tT_matchbox_LO.in"
GenericInputFileNLO 	= "tT_matchbox_NLO.in"
WorkFolder		="/afs/ipp-garching.mpg.de/home/l/lscyboz/"
SetupFileNameGen    	= "setupfile.in"

GoSamLO			= SettingsFolder+"gosamtT.rc"
GoSamNLO		= SettingsFolder+"gosamtTNLO.rc"

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

def printSetupLinesInSubmitFileRivet(file):

#    setupLines = ["export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/ \n",
#                  "source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh \n",
#                  "asetup 20.8.2 \n"
#                  ]

    setupLines = ["source /afs/ipp-garching.mpg.de/home/l/lscyboz/Herwig-7.0.3/bin/activate\n",
		  "source /afs/ipp-garching.mpg.de/home/l/lscyboz/Herwig-7.0.3/src/Rivet-2.4.0/rivetenv.sh\n",
		  "export RIVET_ANALYSIS_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/RivetCustomAnalyses/:$RIVET_ANALYSIS_PATH\n"]

    file.writelines(setupLines)

def initRun():

    os.system("chmod a+x initRun.sh")
    os.system("./initRun.sh")
    
def createInputFile(order, Ecm, scale, PDF, shower):


    if order=="LO":

      inputFileName=WorkFolder+"GenericLO/"+"Herwig_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+".in"
      GenericfileLO=open(SettingsFolder+GenericInputFileLO,'r')
      fileLO=open(inputFileName,'w')

      for line in GenericfileLO:
    	if line.find("set EventHandler:LuminosityFunction")!=-1:
	  fileLO.write("set EventHandler:LuminosityFunction:Energy "+Ecm+"*GeV"+"\n")
	elif line.find("set Factory:ScaleChoice")!=-1:
	  if scale!="TopPairMassScale" and scale!="TopPairMTScale":
	    fileLO.write("set Factory:ScaleChoice Scales/FixedScale\n")
	    fileLO.write("set Scales/FixedScale:FixedScale "+scale+"*GeV"+"\n")
	  else:
	    fileLO.write("set Factory:ScaleChoice Scales/"+scale+"\n")
	elif line.find("set myPDFset:PDFName")!=-1:
    	  fileLO.write("set myPDFset:PDFName "+PDF+"lo68cl\n")
	elif line.find("read Matchbox/LO-DefaultShower.in")!=-1:
	  if order=="LO" and shower=="default": fileLO.write("read Matchbox/LO-DefaultShower.in\n")
	  elif order=="LO" and shower=="dipole": fileLO.write("read Matchbox/LO-DipoleShower.in\n")
	  else: print "Wrong shower setting\n"
	elif line.find("read Matchbox/FiveFlavourScheme")!=-1:
	  if shower=="default": fileLO.write("read Matchbox/FiveFlavourScheme.in\n")
	  elif shower=="dipole": fileLO.write("read Matchbox/FiveFlavourNoBMassScheme.in\n")
	  else: print "Wrong shower setting\n"
	elif line.find("saverun")!=-1: fileLO.write("saverun tT_matchbox_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+" EventGenerator")
	else: fileLO.write(line)	

    if order=="NLO":

      inputFileName=WorkFolder+"Generic/"+"Herwig_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+".in"
      GenericfileNLO=open(SettingsFolder+GenericInputFileNLO,'r')
      fileNLO=open(inputFileName,'w')
      
      for line in GenericfileNLO:
        if line.find("set EventHandler:LuminosityFunction")!=-1:
          fileNLO.write("set EventHandler:LuminosityFunction:Energy "+Ecm+"*GeV"+"\n")
        elif line.find("set Factory:ScaleChoice")!=-1:
          if scale!="TopPairMassScale" and scale!="TopPairMTScale":
            fileNLO.write("set Factory:ScaleChoice Scales/FixedScale\n")
            fileNLO.write("set Scales/FixedScale:FixedScale "+scale+"*GeV"+"\n")
          else:
            fileNLO.write("set Factory:ScaleChoice Scales/"+scale+"\n")
        elif line.find("set myPDFset:PDFName")!=-1:
          fileNLO.write("set myPDFset:PDFName "+PDF+"nlo68cl\n")
        elif line.find("read Matchbox/LO-DefaultShower.in")!=-1:
          if order=="NLO" and shower=="default": fileNLO.write("read Matchbox/MCatNLO-DefaultShower.in\n")
          elif order=="NLO" and shower=="dipole": fileNLO.write("read Matchbox/MCatNLO-DipoleShower.in\n")
          else: print "Wrong shower setting\n"
	elif line.find("read Matchbox/FiveFlavourScheme")!=-1:
          if shower=="default": fileNLO.write("read Matchbox/FiveFlavourScheme.in\n")
          elif shower=="dipole": fileNLO.write("read Matchbox/FiveFlavourNoBMassScheme.in\n")
          else: print "Wrong shower setting\n"
        elif line.find("saverun")!=-1: fileNLO.write("saverun tT_matchbox_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+" EventGenerator")
        else: fileNLO.write(line)

def SubmitHerwigJob(inputfile):

    outputfile	     = inputfile.split(".in")[0]+".run"

    submitFileNameSH = os.getcwd()+"/Submit_"+inputfile+".sh"

    if not os.path.exists(outputfile):

        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        if outputfile.find("_LO")!=-1: InputFolder=WorkFolder+"GenericLO/"
        elif outputfile.find("_NLO")!=-1: InputFolder=WorkFolder+"Generic/"

        codeLines2 = []
#	codeLines2.append("mkdir -p "+tmpFolder)
        codeLines2.append("cd "+InputFolder)

        codeLines2.append("Herwig read "+InputFolder+inputfile)

        for codeLine in codeLines2:
          submitfile2.write(codeLine+" \n")

        submitfile2.write("rm "+ submitFileNameSH + " \n")
        submitfile2.close()

        cmd = "chmod a+x " + submitFileNameSH
        os.system(cmd)
        cmd = "qsub -e /dev/null -o /dev/null "+ submitFileNameSH
        os.system(cmd)       
        return True

    else:
      return False

optionsFile = open("options.in", 'r')
options = optionsFile.read().split("\n")

os.system("cp "+GoSamLO+" "+WorkFolder+"GenericLO/")
os.system("cp "+GoSamNLO+" "+WorkFolder+"Generic/")

for orders in options[0].split("\t"):

  order = orders
  for energies in options[1].split("\t"):

    Ecm=energies
    for scales in options[2].split("\t"):

	scale=scales
	for pdfs in options[3].split("\t"):

	  pdf=pdfs
	  for showers in options[4].split("\t"):

		shower=showers
		createInputFile(order, Ecm, scale, pdf, shower)
		SubmitHerwigJob("Herwig_"+order+"_"+Ecm+"_"+scale+"_"+pdf+"_"+shower+".in")

#initRun()
#for i in range(100):
#for i in range(nRuns):
#	spec='%03.0f' % (i,)
#	if not os.path.exists(sampledPars+spec):
#		os.system("mkdir -p "+sampledPars+spec)
#	os.system("cp "+InputFolder+InputFileNameGen.split(".run")[0]+".in "+sampledPars)
#	SubmitHerwigJob(nEvPerFile, i)

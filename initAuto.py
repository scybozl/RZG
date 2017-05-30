
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
GenericInputFileNLO 	= "tT_matchbox_NLO_71.in"
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
    
def createInputFile(order, Ecm, scale, PDF, shower, matching, topmass, dileptons, CFactor):


    if order=="LO":

      gosamFilename=WorkFolder+"GenericLO/gosamtT"+topmass+".rc"
      GenericGosamFile=open(GoSamLO, 'r')
      gosamFile=open(gosamFilename, 'w')

      for line in GenericGosamFile:
	if line.find("mT=172.5")!=-1:
	  gosamFile.write("                        GF=0.0000116637, mT="+topmass.split("*GeV")[0]+", wT=1.4426\n")
	else: gosamFile.write(line)
      

      if dileptons==False: inputFileName=WorkFolder+"GenericLO/"+"Herwig_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+".in"
      else: inputFileName=WorkFolder+"GenericLO/"+"Herwig_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+"_dilepton.in"

      GenericfileLO=open(SettingsFolder+GenericInputFileLO,'r')
      fileLO=open(inputFileName,'w')

      for line in GenericfileLO:
    	if line.find("set EventHandler:LuminosityFunction")!=-1:
	  fileLO.write("set EventHandler:LuminosityFunction:Energy "+Ecm+"*GeV"+"\n")
	elif line.find("set Factory:ScaleChoice")!=-1:
	  if scale!="TopPairMassScale" and scale!="TopPairMTScale":
	    fileLO.write("set Factory:ScaleChoice Scales/FixedScale\n")
	    fileLO.write("set Scales/FixedScale:FixedScale "+scale+"*GeV"+"\n")
	  elif scale=="TopPairMassScale" or scale=="TopPairMTScale":
	    fileLO.write("set Factory:ScaleChoice Scales/"+scale+"\n")
	elif line.find("set myPDFset:PDFName")!=-1:
    	  fileLO.write("set myPDFset:PDFName "+PDF+"lo68cl\n")
	elif line.find("read Matchbox/LO-DefaultShower.in")!=-1:
	  if order=="LO" and shower=="default": 
		fileLO.write("read Matchbox/LO-DefaultShower.in\n")
	  	fileLO.write("set /Herwig/Shower/GtoQQbarSplitFn:AngularOrdered Yes\n")
	  elif order=="LO" and shower=="dipole": fileLO.write("read Matchbox/LO-DipoleShower.in\n")
	  else: print "Wrong shower setting\n"
	elif line.find("read Matchbox/FiveFlavourScheme")!=-1:
	  if shower=="default": fileLO.write("read Matchbox/FiveFlavourScheme.in\n")
	  elif shower=="dipole": fileLO.write("read Matchbox/FiveFlavourNoBMassScheme.in\n")
	  else: print "Wrong shower setting\n"
	elif line.find("set t:NominalMass")!=-1:
	  fileLO.write("set t:NominalMass "+topmass+"*GeV\n")
	  fileLO.write("set t:HardProcessMass "+topmass+"*GeV\n")
	elif line.find("gosamtT")!=-1:
	  fileLO.write("set Amplitudes/GoSam:SetupInFilename gosamtT"+topmass+".rc")
	elif line.find("saverun")!=-1 and dileptons==False: fileLO.write("saverun tT_matchbox_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+" EventGenerator")
	elif line.find("saverun")!=-1 and dileptons==True: fileLO.write("saverun tT_matchbox_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+"_dilepton EventGenerator")
	else: fileLO.write(line)	

    if order=="NLO":

      gosamFilename=WorkFolder+"Generic71/gosamtTNLO"+topmass+".rc"
      GenericGosamFile=open(GoSamNLO, 'r')
      gosamFile=open(gosamFilename, 'w')

      for line in GenericGosamFile:
        if line.find("mT=172.5")!=-1:
          gosamFile.write("                        GF=0.0000116637, mT="+topmass.split("*GeV")[0]+", wT=1.4426\n")
        else: gosamFile.write(line)


      if dileptons==False: inputFileName=WorkFolder+"Generic71/"+"Herwig_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+"_C"+str(CFactor)+".in"
      else: inputFileName=WorkFolder+"Generic71/"+"Herwig_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+"_C"+str(CFactor)+"_dilepton.in"

      GenericfileNLO=open(SettingsFolder+GenericInputFileNLO,'r')
      fileNLO=open(inputFileName,'w')
      
      for line in GenericfileNLO:
	if dileptons==True and line.find("set t:Width")!=-1:
	   fileNLO.write("do t:SelectDecayModes t->nu_mu,mu+,b; t->nu_e,e+,b; t->nu_tau,tau+,b;\n")
	   fileNLO.write("create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n")
	   fileNLO.write("insert /Herwig/Generators/EventGenerator:EventHandler:PostDecayHandlers 0 /Herwig/Generators/BRReweighter\n")
        if line.find("set EventHandler:LuminosityFunction")!=-1:
           fileNLO.write("set EventHandler:LuminosityFunction:Energy "+Ecm+"*GeV"+"\n")
        elif line.find("set Factory:ScaleChoice")!=-1:
           if scale!="TopPairMassScale" and scale!="TopPairMTScale" and scale!="TopPairMTScaleHalf" and scale!="MyScale" and scale!="ETScale":
             fileNLO.write("set Factory:ScaleChoice Scales/FixedScale\n")
             fileNLO.write("set Scales/FixedScale:FixedScale "+scale+"*GeV"+"\n")
           elif scale=="TopPairMassScale" or scale=="TopPairMTScale":
             fileNLO.write("set Factory:ScaleChoice Scales/"+scale+"\n")
	   elif scale=="TopPairMTScaleHalf":
	     fileNLO.write("set Factory:ScaleChoice Scales/TopPairMTScale\n")
	     fileNLO.write("read Matchbox/MuDown.in\n")
	   elif scale=="MyScale":
	     fileNLO.write("library MyScale.so\n");
	     fileNLO.write("create Herwig::MyScale MyScale\n")
	     fileNLO.write("set Factory:ScaleChoice MyScale\n")
	   elif scale=="ETScale":
             fileNLO.write("library ETScale.so\n");
             fileNLO.write("create Herwig::MyScale ETScale\n")
             fileNLO.write("set Factory:ScaleChoice ETScale\n")

	elif line.find("read Matchbox/MuUp.in")!=-1:

		fileNLO.write("cd /Herwig/MatrixElements/Matchbox\n")
		fileNLO.write("set Factory:RenormalizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set Factory:FactorizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set MEMatching:RenormalizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set MEMatching:FactorizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set /Herwig/DipoleShower/DipoleShowerHandler:RenormalizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set /Herwig/DipoleShower/DipoleShowerHandler:FactorizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set /Herwig/Shower/ShowerHandler:RenormalizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set /Herwig/Shower/ShowerHandler:FactorizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set /Herwig/Shower/PowhegShowerHandler:RenormalizationScaleFactor "+str(CFactor)+"\n")
		fileNLO.write("set /Herwig/Shower/PowhegShowerHandler:FactorizationScaleFactor "+str(CFactor)+"\n")

        elif line.find("set myPDFset:PDFName")!=-1:
 	  if PDF=="MSTW2008nnlo": 
 		fileNLO.write("set myPDFset:PDFName "+PDF+"68cl\n")
 		fileNLO.write("cd /Herwig/Couplings\n")

		fileNLO.write("set NLOAlphaS:input_scale 91.1876*GeV\n")
		fileNLO.write("set NLOAlphaS:input_alpha_s 0.11707\n")
		fileNLO.write("set NLOAlphaS:QuarkMasses 0, 0, 0, 1.4, 4.75, 1e+10\n")
		fileNLO.write("set NLOAlphaS:max_active_flavours 5\n")
	  elif PDF=="NNPDF":
		fileNLO.write("set myPDFset:PDFName "+PDF+"30_nlo_as_0118\n")
                fileNLO.write("cd /Herwig/Couplings\n")

                fileNLO.write("set NLOAlphaS:input_scale 91.199997*GeV\n")
                fileNLO.write("set NLOAlphaS:input_alpha_s 0.118\n")
                fileNLO.write("set NLOAlphaS:QuarkMasses 0, 0, 0, 1.275, 4.18, 173.07\n")
                fileNLO.write("set NLOAlphaS:max_active_flavours 5\n")

          else: fileNLO.write("set myPDFset:PDFName "+PDF+"nlo68cl\n")

        elif line.find("read Matchbox/LO-DefaultShower.in")!=-1:
           if order=="NLO" and shower=="default": 
 		if matching=="": fileNLO.write("read Matchbox/MCatNLO-DefaultShower.in\n")
 		elif matching=="POWHEG": fileNLO.write("read Matchbox/Powheg-DefaultShower.in\n")
           elif order=="NLO" and shower=="dipole":
 		if matching=="": fileNLO.write("read Matchbox/MCatNLO-DipoleShower.in\n")
 		elif matching=="POWHEG": fileNLO.write("read Matchbox/Powheg-DipoleShower.in\n")
           else: print "Wrong shower setting\n"
        elif line.find("read Matchbox/FiveFlavourScheme")!=-1:
           if shower=="default": fileNLO.write("read Matchbox/FiveFlavourScheme.in\n")
           elif shower=="dipole": fileNLO.write("read Matchbox/FiveFlavourNoBMassScheme.in\n")
           else: print "Wrong shower setting\n"
 	elif line.find("set t:NominalMass")!=-1:
           fileNLO.write("set t:NominalMass "+topmass+"*GeV\n")
           fileNLO.write("set t:HardProcessMass "+topmass+"*GeV\n")
 	elif line.find("gosamtT")!=-1:
           fileNLO.write("set Amplitudes/GoSam:SetupInFilename gosamtTNLO"+topmass+".rc") 
	elif line.find("saverun")!=-1 and dileptons==False: fileNLO.write("saverun tT_matchbox_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+"_C"+str(CFactor)+" EventGenerator")
        elif line.find("saverun")!=-1 and dileptons==True: fileNLO.write("saverun tT_matchbox_"+order+"_"+Ecm+"_"+scale+"_"+PDF+"_"+shower+matching+"_"+topmass+"_C"+str(CFactor)+"_dilepton EventGenerator")
        else: fileNLO.write(line)

def SubmitHerwigJob(inputfile):

    outputfile	     = inputfile.split(".in")[0]+".run"

    submitFileNameSH = os.getcwd()+"/Submit_"+inputfile+".sh"

    if not os.path.exists(outputfile):

        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        if outputfile.find("_LO")!=-1: InputFolder=WorkFolder+"GenericLO/"
        elif outputfile.find("_NLO")!=-1: InputFolder=WorkFolder+"Generic71/"

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
#        cmd = "qsub "+ submitFileNameSH
#        os.system(cmd)       
        return True

    else:
      return False

optionsFile = open("options.in", 'r')
options = optionsFile.read().split("\n")

os.system("cp "+GoSamLO+" "+WorkFolder+"GenericLO/")
os.system("cp "+GoSamNLO+" "+WorkFolder+"Generic71/")

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
	     for matchings in options[5].split("\t"):

		matching=matchings
		for topmasses in options[6].split("\t"):

		  topmass=topmasses

		  if scale=="ETScale":
		    CFactors = [0.25, 0.5, 1.0, 2.0] 

		  else: CFactors = [1.0]

		  for CFactor in CFactors:

		    createInputFile(order, Ecm, scale, pdf, shower, matching, topmass, False,CFactor)
		    SubmitHerwigJob("Herwig_"+order+"_"+Ecm+"_"+scale+"_"+pdf+"_"+shower+matching+"_"+topmass+"_C"+str(CFactor)+".in")

		    if Ecm=="13000":
                      createInputFile(order, Ecm, scale, pdf, shower, matching, topmass, True,CFactor)
		      SubmitHerwigJob("Herwig_"+order+"_"+Ecm+"_"+scale+"_"+pdf+"_"+shower+matching+"_"+topmass+"_C"+str(CFactor)+"_dilepton.in")

#initRun()
#for i in range(100):
#for i in range(nRuns):
#	spec='%03.0f' % (i,)
#	if not os.path.exists(sampledPars+spec):
#		os.system("mkdir -p "+sampledPars+spec)
#	os.system("cp "+InputFolder+InputFileNameGen.split(".run")[0]+".in "+sampledPars)
#	SubmitHerwigJob(nEvPerFile, i)

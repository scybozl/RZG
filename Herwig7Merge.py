
#! /usr/bin/python

import os, time
import sys
import fileinput
from glob import glob
import subprocess

nEvPerFile = 5000
nRuns = 400

fUser = os.getenv("USER")
sampledPars = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Pozzorini/"
InputFileNameGen = "Pozzorini.run"
SettingsFolder    = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/"
InputFolder="/afs/ipp-garching.mpg.de/home/l/lscyboz/Generic/"
SetupFileNameGen    = "setupfile_noshower.in"

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
    

def SubmitHerwigJob(nEvents, seed):


    specStr          = '%03.0f' % (seed,)
#    tmpFolder        = sampledPars+specStr+"/"
    OutputFile       = "seed_"+specStr+".hepmc"
    OutputFolder     = sampledPars+specStr+"/"
    OutputYoda        = OutputFolder+"seed_"+specStr+".yoda"
    tmp		     = "/tmp/lscyboz/"+specStr+"/"

#    os.system("mkdir -p "+OutputFolder)

    OutputFileFinal  = OutputFolder+OutputFile

    submitFileNameSH = os.getcwd()+"/Submit_"+specStr+".sh"
    redo = False
    if os.path.exists(OutputYoda):
	f=open(OutputYoda,'r')
	for i, line in enumerate(f):
	  if "END YODA_COUNTER" in line:
		number=i
	f=open(OutputYoda,'r')
	string=f.readlines()[number-1]
	nE=float(string.split("\t")[2].split("\n")[0])
	if nE!=nEvents:
	  redo = True

    if not os.path.exists(OutputYoda) or redo==True:

        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        codeLines2 = []
#	codeLines2.append("mkdir -p "+tmpFolder)
	codeLines2.append("cd "+InputFolder)
	codeLines2.append("mkdir -p "+tmp)
	codeLines2.append("cp "+SettingsFolder+SetupFileNameGen+" "+OutputFolder)
	codeLines2.append("echo 'set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed "+str(seed)+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Analysis/HepMCFile:Filename "+tmp+OutputFile+"' >> "+OutputFolder+SetupFileNameGen)

	codeLines2.append("Herwig run "+InputFileNameGen+" -N "+str(nEvents)+" -x "+OutputFolder+SetupFileNameGen)
	codeLines2.append("rivet -a bB4l_comparison "+tmp+OutputFile+" -H "+OutputYoda+" -x 251.659") 
	codeLines2.append("rivet -a bB4l_comparison "+tmp+OutputFile+" -H "+OutputYoda.split(".yoda")[0]+"_unnorm.yoda")

#        codeLines2.append("cp "+tmpFolder+OutputFile+" "+OutputFileFinal)
#	codeLines2.append("cp "+tmpFolder+SetupFileNameGen+" "+OutputFolder)
#        codeLines2.append("rm -rf "+tmpFolder)

        for codeLine in codeLines2:
            submitfile2.write(codeLine+" \n")

        submitfile2.write("rm "+ submitFileNameSH + " \n")
	submitfile2.write("rm -r "+ tmp + " \n")
        submitfile2.close()

        cmd = "chmod a+x " + submitFileNameSH
        os.system(cmd)
        cmd = "qsub "+ submitFileNameSH
        os.system(cmd)
        
        return True

    else:
        return False

#initRun()
#for i in range(100):
for i in range(nRuns):
	spec='%03.0f' % (i,)
	if not os.path.exists(sampledPars+spec):
		os.system("mkdir -p "+sampledPars+spec)
	os.system("cp "+InputFolder+InputFileNameGen.split(".run")[0]+".in "+sampledPars)
	SubmitHerwigJob(nEvPerFile, i)

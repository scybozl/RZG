
#! /usr/bin/python

import os, time
import sys
import fileinput
from glob import glob
import subprocess

fUser = os.getenv("USER")
sampledPars = "/afs/ipp-garching.mpg.de/home/l/lscyboz/mc/"

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
    

def SubmitHerwigJob(nEvents, seed, alphaSMZ, lambdaQCD, Qmin, pT0min, b):

    InputFolder="/afs/ipp-garching.mpg.de/home/l/lscyboz/Generic/"
    InputFileNameGen = "tT_matchbox_LO.run"
    SettingsFolder    = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/"
    SetupFileNameGen    = "setupfile.in"

    specStr          = '%03.0f' % (seed,)
#    tmpFolder        = sampledPars+specStr+"/"
    OutputFile       = "seed_"+specStr+".hepmc"
    OutputFolder     = sampledPars+specStr+"/"
    OutputYoda        = OutputFolder+"seed_"+specStr+".yoda"

#    os.system("mkdir -p "+OutputFolder)

    OutputFileFinal  = OutputFolder+OutputFile

    submitFileNameSH = os.getcwd()+"/Submit_"+specStr+".sh"

    if not os.path.exists(OutputYoda):

        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        codeLines2 = []
#	codeLines2.append("mkdir -p "+tmpFolder)
	codeLines2.append("cd "+InputFolder)
	codeLines2.append("cp "+SettingsFolder+SetupFileNameGen+" "+OutputFolder)
	codeLines2.append("echo 'set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed "+str(seed)+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Analysis/HepMCFile:Filename "+OutputFileFinal+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:AlphaMZ "+alphaSMZ+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:LambdaQCD "+lambdaQCD+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:Qmin "+Qmin+"' >> "+OutputFolder+SetupFileNameGen)
#	codeLines2.append("echo 'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 "+pT0min+"' >> "+OutputFolder+SetupFileNameGen)
#	codeLines2.append("echo 'set /Herwig/UnderlyingEvent/MPIHandler:Power "+b+"' >> "+OutputFolder+SetupFileNameGen)
	
	codeLines2.append("echo 'set /Herwig/Hadronization/ClusterFissioner:PSplitLight "+pT0min+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Hadronization/ClusterFissioner:ClPowLight "+b+"' >> "+OutputFolder+SetupFileNameGen)

	codeLines2.append("Herwig run "+InputFileNameGen+" -N "+str(nEvents)+" -x "+OutputFolder+SetupFileNameGen)
	codeLines2.append("rivet -a ATLAS_2014_I1304688_custom -a ATLAS_2012_I1094568 -a ATLAS_2013_I1243871 "+OutputFileFinal+" -H "+OutputYoda) 

#        codeLines2.append("cp "+tmpFolder+OutputFile+" "+OutputFileFinal)
#	codeLines2.append("cp "+tmpFolder+SetupFileNameGen+" "+OutputFolder)
#        codeLines2.append("rm -rf "+tmpFolder)

        for codeLine in codeLines2:
            submitfile2.write(codeLine+" \n")

        submitfile2.write("rm "+ submitFileNameSH + " \n")
	submitfile2.write("rm "+ OutputFileFinal + " \n")
        submitfile2.close()

        cmd = "chmod a+x " + submitFileNameSH
        os.system(cmd)
        cmd = "qsub -l h_rt=08:00:00 -m as -M scyboz@mpp.mpg.de " + submitFileNameSH
        os.system(cmd)
        
        return True

    else:
        return False

#initRun()
#for i in range(100):
for subdir in SubDirPath(sampledPars):
	params=open(subdir+"/used_params",'r')
	for line in params:
		if 'alphaSMZ' in line:
		  alphaSMZ=line.split()[1]
		if 'lambdaQCD' in line:
		  lambdaQCD=line.split()[1]
		if 'Qmin' in line:
		  Qmin=line.split()[1]
		if 'Psplit' in line:
		  pT0min=line.split()[1]
		if 'Clpow' in line:
		  b=line.split()[1]
	i=int(subdir.split("mc/")[1])
	SubmitHerwigJob(2000, i, alphaSMZ, lambdaQCD, Qmin, pT0min, b)

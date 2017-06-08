#! /usr/bin/python

import os, time
import sys
import fileinput
import glob
import subprocess

fUser = os.getenv("USER")
nEvPerFile = 50000
nRuns = 500
newMerge = True
newControl = True
ControlIndex = ""
EnergyIndex = ""

fUser = os.getenv("USER")
WorkFolder        = "/afs/ipp-garching.mpg.de/home/l/lscyboz/"
SettingsFolder    = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/"
SetupFileNameGen    = "setupfiledipole.in"
sampledPars	  = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Draco/Output_dipole/"
InputFolder	  = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Draco/"

flag=False


def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

def printSetupLinesInSubmitFileRivet(file):

#    setupLines = ["export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/ \n",
#                  "source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh \n",
#                  "asetup 20.8.2 \n"
#                  ]

    setupLines = ["source /afs/ipp-garching.mpg.de/home/l/lscyboz/setup.sh\n",
		  "source /afs/ipp-garching.mpg.de/home/l/lscyboz/Herwig-7.1.0/bin/activate\n",
		  "export RIVET_ANALYSIS_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/RivetCustomAnalyses/:$RIVET_ANALYSIS_PATH\n"]

    file.writelines(setupLines)

def initRun():

    os.system("chmod a+x initRun.sh")
    os.system("./initRun.sh")
    

def SubmitHerwigJob(nEvents, seed, InputFileNameGen):


    specStr          = '%03.0f' % (seed,)
#    tmpFolder        = sampledPars+specStr+"/"
    OutputFile       = "seed_"+specStr+".hepmc"
    OutputFolder     = sampledPars+specStr+"/"
    OutputYoda        = OutputFolder+"seed_"+specStr+"_"
    tmp              = "$TMPDIR/lscyboz/"+specStr+"/"

#    os.system("mkdir -p "+OutputFolder)

    OutputFileFinal  = OutputFolder+OutputFile

    submitFileNameSH = WorkFolder+"Submit_tTShower_"+specStr+".sh"


    if redo==True:

        flag=True

        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        codeLines2 = []
#       codeLines2.append("mkdir -p "+tmpFolder)
        codeLines2.append("cd "+InputFolder)
        codeLines2.append("mkdir -p "+tmp)
        codeLines2.append("cp "+SettingsFolder+SetupFileNameGen+" "+OutputFolder)
        codeLines2.append("echo 'set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed "+str(seed)+"' >> "+OutputFolder+SetupFileNameGen)
        codeLines2.append("echo \"set /Herwig/Analysis/HepMCFile:Filename "+tmp+OutputFile+"\" >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/EventHandlers/EventHandler:HadronizationHandler NULL' >> "+OutputFolder+SetupFileNameGen)

        codeLines2.append("Herwig run "+InputFileNameGen+" -N "+str(nEvents)+" -x "+OutputFolder+SetupFileNameGen)
        codeLines2.append("rivet -a MC_MARKUS13TEV_inclusive "+tmp+OutputFile+" -H "+OutputYoda+"unnorm.yoda")

        for codeLine in codeLines2:
            submitfile2.write(codeLine+" \n")

	submitfile2.write("rm -r "+ tmp + " \n")
        submitfile2.close()

        cmd = "chmod a+x " + submitFileNameSH
        os.system(cmd)
        cmd = "qsub -l h_rt=05:30:00 "+ submitFileNameSH
 #       os.system(cmd)

        return True

    else:
        return False


## Options file for systematic generation: the user should set the settings required for the different runs there

os.system("export RIVET_ANALYSIS_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/RivetCustomAnalyses/:$RIVET_ANALYSIS_PATH")
redo=True
## Name tag for the run

## Submit the job to Herwig
for i in range(500,500+nRuns):
   spec='%03.0f' % (i,)
   if not os.path.exists(sampledPars+spec):
      os.system("mkdir -p "+sampledPars+spec)
   os.system("cp "+InputFolder+"tTShower_dipole.in "+sampledPars)
   if (i+1)%100==0: print "Processing run #"+str(i)
   SubmitHerwigJob(nEvPerFile, i, "tTShower_dipole.run")

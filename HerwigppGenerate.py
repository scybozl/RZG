
#! /usr/bin/python

import os, time
import sys
import fileinput
from glob import glob
import subprocess

fUser = os.getenv("USER")

def printSetupLinesInSubmitFileRivet(file):

    setupLines = ["export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/ \n",
                  "source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh \n",
                  "asetup 20.8.2 \n"
                  ]

    file.writelines(setupLines)



def SubmitHerwigJob(nEvents, seed):

    CommonOpts       = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/herwigpp.common.py"
    TopOpts          = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/herwigpp.top.py"
    StandalonePars   = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/standalonePars.in"

    OutputFile       = "seed_"+str(seed)+".hepmc"

    OutputFolder     = "/afs/ipp-garching.mpg.de/home/l/lscyboz/JobOutput/"+str(seed)

    os.system("mkdir -p "+OutputFolder)

    OutputFileFinal  = OutputFolder+"/"+OutputFile

    submitFileNamePY = os.getcwd()+"/RunHerwig_"+str(seed)+".py"
    submitFileNameSH = os.getcwd()+"/Submit_"+str(seed)+".sh"

    if not os.path.exists(OutputFileFinal):

        submitfile = open(submitFileNamePY, "w")
        codeLines = []
        codeLines.append("from AthenaCommon.AppMgr import ServiceMgr")
        codeLines.append("ServiceMgr.MessageSvc.OutputLevel = DEBUG")
        codeLines.append("theApp.EvtMax = "+str(nEvents))
        codeLines.append("from AthenaServices.AthenaServicesConf import AtRndmGenSvc")
        codeLines.append("ServiceMgr += AtRndmGenSvc()")
        codeLines.append("from AthenaCommon.AlgSequence import AlgSequence")
        codeLines.append("job=AlgSequence()")
        codeLines.append("from Herwigpp_i.Herwigpp_iConf import Herwigpp")
        codeLines.append("job += Herwigpp()")
        codeLines.append("cmds = []")
        codeLines.append("include('herwigpp.common.py')")
        codeLines.append("include('herwigpp.top.py')")
        codeLines.append("f=open('standalonePars.in', 'r')")
        codeLines.append("for line in f:")
        codeLines.append("\tcmds.append(line.split('\\n')[0])")
	codeLines.append("cmds.append('set /Herwig/Generators/LHCGenerator:RandomNumberGenerator:Seed "+str(seed)+"')")
        codeLines.append("cmds.append('set /Herwig/Analysis/HepMCFile:Filename "+OutputFile+"')")
        codeLines.append("job.Herwigpp.Commands += cmds")
        
        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        codeLines2 = []
        codeLines2.append("mkdir -p /tmp/$USER/$PBS_JOBID")
        codeLines2.append("cd /tmp/$USER/$PBS_JOBID")
	codeLines2.append("cp "+submitFileNamePY+" .")
        codeLines2.append("cp "+CommonOpts+" . ")
        codeLines2.append("cp "+TopOpts+" . ")
        codeLines2.append("cp "+StandalonePars+" . ")
        codeLines2.append("athena.py "+"RunHerwig"+submitFileNamePY.split("RunHerwig")[1])
        codeLines2.append("cp -r "+OutputFile+" "+OutputFileFinal)
	codeLines2.append("cp "+TopOpts+" "+OutputFolder)
        codeLines2.append("cd ../")
        codeLines2.append("rm -rf /tmp/$USER/$PBS_JOBID")

        for codeLine in codeLines:
            submitfile.write(codeLine+" \n")

        submitfile.close()

        for codeLine in codeLines2:
            submitfile2.write(codeLine+" \n")

        submitfile2.write("rm "+ submitFileNamePY + " \n")
        submitfile2.write("rm "+ submitFileNameSH + " \n")
        submitfile2.close()

        cmd = "chmod a+x " + submitFileNameSH
        os.system(cmd)
        cmd = "qsub -m as -M scyboz@mpp.mpg.de -o /dev/null -e /dev/null " + submitFileNameSH
        os.system(cmd)
        
        return True

    else:
        return False


for i in range(400):
SubmitHerwigJob(2000, i)

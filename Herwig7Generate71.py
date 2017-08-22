#! /usr/bin/python

import os, time
import sys
import fileinput
import glob
import subprocess

fUser = os.getenv("USER")
nEvPerFile = 10000
nRuns = 1000
newMerge = False
newControl = False
ControlIndex = ""
EnergyIndex = ""

fUser = os.getenv("USER")
SettingsFolder    = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/"
SetupFileNameGen    = "setupfile.in"
WorkFolder        = "/afs/ipp-garching.mpg.de/home/l/lscyboz/"

pars 		= "/afs/ipp-garching.mpg.de/home/l/lscyboz/scan/"

flag=False


def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

def printSetupLinesInSubmitFileRivet(file):

#    setupLines = ["export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/ \n",
#                  "source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh \n",
#                  "asetup 20.8.2 \n"
#                  ]

    setupLines = ["#!/bin/bash -l\n",
		  "# Standard output and error:\n",
		  "#SBATCH -o ./tjob.out.%j\n",
		  "#SBATCH -e ./tjob.err.%j\n",
		  "# Initial working directory:\n",
		  "#SBATCH -D ./\n",
		  "# Job Name:\n",
		  "#SBATCH -J test_slurm\n",
	 	  "# Queue (Partition):\n",
		  "#SBATCH --partition=standard\n",
#		  "#SBATCH --nodes=1\n",
#		  "#SBATCH --ntasks-per-node=32\n",
		  "source /afs/ipp-garching.mpg.de/home/l/lscyboz/setup.sh\n",
		  "source /afs/ipp-garching.mpg.de/home/l/lscyboz/Herwig-7.1.0/bin/activate\n",
		  "export RIVET_ANALYSIS_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/RivetCustomAnalyses/:$RIVET_ANALYSIS_PATH\n"]

    file.writelines(setupLines)

def initRun():

    os.system("chmod a+x initRun.sh")
    os.system("./initRun.sh")
    

def SubmitHerwigJob(nEvents, seed, alphaSMZ, InputFileNameGen, index):

    specStr          = '%03.0f' % (seed,)
#    tmpFolder        = sampledPars+specStr+"/"
    OutputFile       = "seed_"+specStr+".hepmc"
    OutputFolder     = sampledPars+specStr+"/"
    OutputYoda        = OutputFolder+"seed_"+specStr+"_"
    tmp              = "$TMPDIR/lscyboz/"+settings+"_"+specStr+"/"

#    os.system("mkdir -p "+OutputFolder)

    OutputFileFinal  = OutputFolder+OutputFile

    submitFileNameSH = WorkFolder+"Submit_"+settings+"_"+specStr+".sh"

    flag = False
    redo = False
    for norms in options[index+1].split("\t"):
      if os.path.exists(OutputYoda+norms+".yoda"):
        f=open(OutputYoda+norms+".yoda",'r')
        for i, line in enumerate(f):
          if "END YODA_COUNTER" in line:
                number=i
        f=open(OutputYoda+norms+".yoda",'r')
        string=f.readlines()[number-1]
        nE=float(string.split("\t")[2].split("\n")[0])
        if nE!=nEvents:
          redo = True
      else: redo = True
    os.chdir(OutputFolder)
    if len(glob.glob('*.yoda')) <= len(options[index+1].split("\t")):
        redo = True

    if not os.path.exists(OutputYoda+options[index+1].split("\t")[0]+".yoda") or redo==True:

        flag=True

        submitfile2 = open(submitFileNameSH, "w")
        printSetupLinesInSubmitFileRivet(submitfile2)

        codeLines2 = []
#       codeLines2.append("mkdir -p "+tmpFolder)
        codeLines2.append("cd "+InputFolder)
        codeLines2.append("mkdir -p "+tmp)
	codeLines2.append("cp "+SettingsFolder+SetupFileNameGen+" "+OutputFolder)
        codeLines2.append("echo 'set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed "+str(5000+seed)+"' >> "+OutputFolder+SetupFileNameGen)
        codeLines2.append("echo \"set /Herwig/Analysis/HepMCFile:Filename "+tmp+OutputFile+"\" >> "+OutputFolder+SetupFileNameGen)

	if float(alphaSMZ)>=0.145:
	  codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:Qmin 1.200' >> "+OutputFolder+SetupFileNameGen)
	  codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:ThresholdOption Current' >> "+OutputFolder+SetupFileNameGen)
	  codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:NumberOfLoops 3' >> "+OutputFolder+SetupFileNameGen)

	if(InputFileNameGen.find("dipole")!=-1):
	  codeLines2.append("echo 'set /Herwig/DipoleShower/NLOAlphaS:input_alpha_s "+alphaSMZ+"' >> "+OutputFolder+SetupFileNameGen)
	else: codeLines2.append("echo 'set /Herwig/Shower/AlphaQCD:AlphaMZ "+alphaSMZ+"' >> "+OutputFolder+SetupFileNameGen)

        codeLines2.append("Herwig run "+InputFileNameGen+" -N "+str(nEvents)+" -x "+OutputFolder+SetupFileNameGen)

        analyses=""
        for routines in options[index].split("\t"):
                analyses += " -a "+routines
        for norms in options[index+1].split("\t"):
                codeLines2.append("rivet"+analyses+" "+tmp+OutputFile+" -H "+OutputYoda+norms+".yoda -x "+norms)
        codeLines2.append("rivet"+analyses+" "+tmp+OutputFile+" -H "+OutputYoda+"unnorm.yoda")

#        codeLines2.append("cp "+tmpFolder+OutputFile+" "+OutputFileFinal)
#       codeLines2.append("cp "+tmpFolder+SetupFileNameGen+" "+OutputFolder)
#        codeLines2.append("rm -rf "+tmpFolder)

        for codeLine in codeLines2:
            submitfile2.write(codeLine+" \n")

        submitfile2.write("rm "+ submitFileNameSH + " \n")
        submitfile2.write("rm -r "+ tmp + " \n")
        submitfile2.close()

        cmd = "chmod a+x " + submitFileNameSH
        os.system(cmd)
        cmd = "sbatch "+ submitFileNameSH
        os.system(cmd)

        return True

    else:
        return False


## Options file for systematic generation: the user should set the settings required for the different runs there

optionsFile = open("options_dilep.in", 'r')
options = optionsFile.read().split("\n")
os.system("export RIVET_ANALYSIS_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/RivetCustomAnalyses/:$RIVET_ANALYSIS_PATH")


#initRun()
#for i in range(100):
for subdir in SubDirPath(pars):
	params=open(subdir+"/params.dat",'r')
	for line in params:
		if 'alphaSMZ' in line:
		  alphaSMZ=line.split()[1]
		if 'CFactor' in line:
		  CFactor=line.split()[1]
		if 'pTmin' in line:
		  pTmin=line.split()[1]
		if 'ClMaxLight' in line:
		  ClMaxLight=line.split()[1]
		if 'PSplitLight' in line:
		  PSplitLight=line.split()[1]
	
	## Loop through all possible combinations

	## LO, NLO
	for orders in options[0].split("\t"):

	  order = orders
	  ## CM energy
	  for e,energies in enumerate(options[1].split("\t")):

	    Ecm=energies
	    ## Renormalization and factorization scale choices
	    for scales in options[2].split("\t"):

	        scale=scales
	        ## PDF choices (only indicate the PDF name, the order is chosen
	        ## automatically as the order of the process!)
	        for pdfs in options[3].split("\t"):

	          pdf=pdfs
	          ## Shower (default or dipole)
	          for showers in options[4].split("\t"):

	             shower=showers
		     ## Matching (MCatNLO or POWHEG)
		     for matchings in options[5].split("\t"):

		      matching=matchings
		      for topmasses in options[6].split("\t"):
			
			topmass=topmasses
	                ## Name tag for the run
#	                settings=order+"_"+Ecm+"_"+scale+"_"+pdf+"_"+shower+matching+"_"+topmass+"_C"+CFactor+"_"+alphaSMZ
	                settings=order+"_"+Ecm+"_"+scale+"_"+pdf+"_"+shower+matching+"_"+topmass+"_C"+CFactor+"_dilepton_"+alphaSMZ

	                sampledPars = "/afs/ipp-garching.mpg.de/home/l/lscyboz/MC_Herwig_"+settings+"/"

	                ## To choose the Rivet routine according to the cm-energy, look into the
	                ## options file at the right placee
	                index=3*e+8

	                ## If no control (number of runs, right number of events...)
	                ## is needed, just control if one of the final yoda files exists.
	                breakLoop = False or Ecm!="13000" or alphaSMZ!="0.132" or CFactor!="0.5"# or shower!="dipole" 
	                if os.path.exists(sampledPars+"MC_Herwig_"+settings+"_"+options[index+1].split("\t")[0]+".yoda") and newControl == False or breakLoop: break
	                if order=="LO":
	                          InputFolder="/afs/ipp-garching.mpg.de/home/l/lscyboz/GenericLO/"
	                elif order=="NLO":
	                          InputFolder="/afs/ipp-garching.mpg.de/home/l/lscyboz/Generic71/"

	                ## IF some yoda files were generated a second time, re-run the yoda merging
	                flag=False
	                print "Now processing "+settings+"...\n"


	                ## Submit the job to Herwig
	                for i in range(1000,1000+nRuns):
	                        spec='%03.0f' % (i,)
	                        if not os.path.exists(sampledPars+spec):
	                          os.system("mkdir -p "+sampledPars+spec)
	                        os.system("cp "+InputFolder+"Herwig_"+settings.split("_"+alphaSMZ)[0]+".in "+sampledPars)
	                        if (i+1)%100==0: print "Processing run #"+str(i)
	                        SubmitHerwigJob(nEvPerFile, i, alphaSMZ, "tT_matchbox_"+settings.split("_"+alphaSMZ)[0]+".run", index)
				if (i+1)%100==0:
					while True:
		                          os.system('squeue -u lscyboz > file')
                		          strn=open('file', 'r')
#           		                  if len(strn) <= 800: break
                                          if sum(1 for line in strn)<495: break
                                          time.sleep(15)
                                          print "."
                                        print "\n"

			## As long as there are processed jobs in the queue, wait
	                while True:
	                  os.system('squeue -u lscyboz > file')
	                  strn=open('file', 'r')
#	                  if len(strn) <= 800: break
	                  if sum(1 for line in strn)<495: break
	                  time.sleep(15)
	                  print "."
	                print "\n"

	                ## Yoda-merge the files from the different runs 
	                for norms in options[index+1].split("\t"):
	                  if not os.path.exists(sampledPars+"MC_Herwig_"+settings+"_"+norms+".yoda") or flag==True or newMerge==True:
	                        print "Yoda-merging "+settings+" at xs="+norms+" pb"
	                        os.system("yodamerge "+sampledPars+"*/*"+norms+".yoda -o "+sampledPars+"MC_Herwig_"+settings+"_"+norms+".yoda")
	                if not os.path.exists(sampledPars+"MC_Herwig_"+settings+"_unnorm.yoda") or flag==True or newMerge==True:
	                  print "Yoda-merging "+settings+" at generated cross-section"
	                  os.system("yodamerge "+sampledPars+"*/*"+norms+".yoda -o "+sampledPars+"MC_Herwig_"+settings+"_unnorm.yoda")

			os.system("cp "+sampledPars+"MC_Herwig_"+settings+"*.yoda "+subdir)


#! /usr/bin/python

import os, time
import sys
import fileinput
import glob
import subprocess

nEvPerFile = 5000
nRuns = 400
newMerge = True


fUser = os.getenv("USER")
SettingsFolder    = "/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings/"
SetupFileNameGen    = "setupfile.in"
WorkFolder	  = "/afs/ipp-garching.mpg.de/home/l/lscyboz/"

flag=False

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

def SubmitHerwigJob(nEvents, seed, InputFileNameGen, index):


    specStr          = '%03.0f' % (seed,)
#    tmpFolder        = sampledPars+specStr+"/"
    OutputFile       = "seed_"+specStr+".hepmc"
    OutputFolder     = sampledPars+specStr+"/"
    OutputYoda        = OutputFolder+"seed_"+specStr+"_"
    tmp		     = "/tmp/lscyboz/"+settings+"_"+specStr+"/"

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
#	codeLines2.append("mkdir -p "+tmpFolder)
	codeLines2.append("cd "+InputFolder)
	codeLines2.append("mkdir -p "+tmp)
	codeLines2.append("cp "+SettingsFolder+SetupFileNameGen+" "+OutputFolder)
	codeLines2.append("echo 'set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed "+str(seed)+"' >> "+OutputFolder+SetupFileNameGen)
	codeLines2.append("echo 'set /Herwig/Analysis/HepMCFile:Filename "+tmp+OutputFile+"' >> "+OutputFolder+SetupFileNameGen)

	codeLines2.append("Herwig run "+InputFileNameGen+" -N "+str(nEvents)+" -x "+OutputFolder+SetupFileNameGen)

	analyses=""
	for routines in options[index].split("\t"):
		analyses += " -a "+routines
        for norms in options[index+1].split("\t"):
                codeLines2.append("rivet"+analyses+" "+tmp+OutputFile+" -H "+OutputYoda+norms+".yoda -x "+norms) 
	codeLines2.append("rivet"+analyses+" "+tmp+OutputFile+" -H "+OutputYoda+"unnorm.yoda")

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
        cmd = "qsub -e /dev/null -o /dev/null "+ submitFileNameSH
        os.system(cmd)
        
        return True

    else:
        return False


## Options file for systematic generation: the user should set the settings required for the different runs there

optionsFile = open("options.in", 'r')
options = optionsFile.read().split("\n")
os.system("source /afs/ipp-garching.mpg.de/home/l/lscyboz/Herwig-7.0.3/src/Rivet-2.4.0/rivetenv.sh")
os.system("export RIVET_ANALYSIS_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/RivetCustomAnalyses/:$RIVET_ANALYSIS_PATH")


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
		## Name tag for the run
		settings=order+"_"+Ecm+"_"+scale+"_"+pdf+"_"+shower
		if order=="LO":
                          InputFolder="/afs/ipp-garching.mpg.de/home/l/lscyboz/GenericLO/"
                elif order=="NLO":
                          InputFolder="/afs/ipp-garching.mpg.de/home/l/lscyboz/Generic/"

		## To choose the Rivet routine according to the cm-energy, look into the
		## options file at the right placee
		index=3*e+6

		## IF some yoda files were generated a second time, re-run the yoda merging
		flag=False
		print "Now processing "+settings+"...\n"

		## Submit the job to Herwig
		for i in range(nRuns):
			spec='%03.0f' % (i,)
			sampledPars = "/afs/ipp-garching.mpg.de/home/l/lscyboz/MC_Herwig_"+settings+"/"
			if not os.path.exists(sampledPars+spec):
			  os.system("mkdir -p "+sampledPars+spec)
			os.system("cp "+InputFolder+"Herwig_"+settings+".in "+sampledPars)
			if (i+1)%100==0: print "Processing run #"+str(i)
			SubmitHerwigJob(nEvPerFile, i, "tT_matchbox_"+settings+".run", index)

		## As long as there are jobs in the queue, wait
		while True:
                  os.system('qstat -u lscyboz > file')
                  strn=open('file', 'r').read()
		  if strn=='': break
                  #if sum(1 for line in strn)<402: break
                  time.sleep(5)
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

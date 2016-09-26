#==============================================================
#--------------------------------------------------------------
# Private Application Configuration option
from AthenaCommon.AppMgr import ServiceMgr
ServiceMgr.MessageSvc.OutputLevel = DEBUG

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------
# Number of events to be processed (default is 10)
theApp.EvtMax = 10
#--------------------------------------------------------------
# Algorithms Private Options
#--------------------------------------------------------------
from AthenaServices.AthenaServicesConf import AtRndmGenSvc
ServiceMgr += AtRndmGenSvc()
## ServiceMgr.AtRndmGenSvc.Seeds = ["PYTHIA 4789899 989240512",
##                                  "PYTHIA_INIT 820021 2347532"]
# ServiceManager.AtRndmGenSvc.ReadFromFile = true;

from AthenaCommon.AlgSequence import AlgSequence
job=AlgSequence()
from Herwigpp_i.Herwigpp_iConf import Herwigpp
job += Herwigpp()

#herwigpp = job.Herwigpp
cmds = []

include("herwigpp.common.py")
include("herwigpp.top.py")

f=open("standalonePars.in", 'r')
for line in f:
	cmds.append(line.split("\n")[0])

cmds.append("set /Herwig/Analysis/HepMCFile:Filename seed_0.hepmc")

job.Herwigpp.Commands += cmds

#sj#from TruthExamples.TruthExamplesConf import DumpMC
#sj#job += DumpMC()


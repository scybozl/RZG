import os
sampledPars = "/afs/ipp-garching.mpg.de/home/l/lscyboz/mc/"

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])


for subdir in SubDirPath(sampledPars):
	  if os.path.exists(subdir+"/seed_"+subdir.split("mc/")[1]+".yoda"):
		os.system("./yoda2aida "+subdir+"/seed_"+subdir.split("mc/")[1]+".yoda "+subdir+"/out"+".aida")

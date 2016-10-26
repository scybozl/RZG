import os
sampledPars = "/afs/ipp-garching.mpg.de/home/l/lscyboz/tuning_TopPairMassScale_MMHT2014_default/mc/"

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])


for subdir in SubDirPath(sampledPars):
    os.system("./yoda2aida "+subdir+"/out.yoda "+subdir+"/out.aida")

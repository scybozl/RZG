import os
import numpy as np
import sys

alphaSMZ = [0.125, 0.140]
ClMaxLight = [0.5, 3.5]
PSplitLight = [0.5, 3.5]

nsteps = 4 
nstepsClMaxLight = 3
nstepsPSplit = 3

if os.path.isdir("scan"):
  sys.exit("Scan directory already exists... Aborting.")

os.system("mkdir scan")
os.chdir("scan")

index=0

for i in np.linspace(alphaSMZ[0], alphaSMZ[1], nsteps):
  for j in np.linspace(ClMaxLight[0], ClMaxLight[1], nstepsClMaxLight):
    for k in np.linspace(PSplitLight[0], PSplitLight[1], nstepsPSplit):

	os.system("mkdir "+str(index))
	os.chdir(str(index))
	parFile = open("params.dat", 'w')

	parFile.write("alphaSMZ\t"+str(i)+"\n")
	parFile.write("ClMaxLight\t"+str(j)+"\n")
	parFile.write("PSplitLight\t"+str(k))

	parFile.close()

	index += 1
	os.chdir("..")

import os
import numpy as np
import sys

alphaSMZ = [0.125, 0.140]
ClMaxLight = [0.5, 2.0]
PSplitLight = [0.5, 2.0]

nsteps = 2 

if os.path.isdir("scan"):
  sys.exit("Scan directory already exists... Aborting.")

os.system("mkdir scan")
os.chdir("scan")

index=0

for i in np.linspace(alphaSMZ[0], alphaSMZ[1], nsteps):
  for j in np.linspace(ClMaxLight[0], ClMaxLight[1], nsteps):
    for k in np.linspace(PSplitLight[0], PSplitLight[1], nsteps):

	os.system("mkdir "+str(index))
	os.chdir(str(index))
	parFile = open("params.dat", 'w')

	parFile.write("alphaSMZ "+str(i)+"\n")
	parFile.write("ClMaxLight "+str(j)+"\n")
	parFile.write("PSplitLight "+str(k))

	parFile.close()

	index += 1
	os.chdir("..")

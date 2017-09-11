import os
import numpy as np
import sys

alphaSMZ = [0.127, 0.137]
CFactor = [0.25, 0.5, 1.0, 2.0]

nstepsalphaSMZ = 3

if os.path.isdir("scan"):
  sys.exit("Scan directory already exists... Aborting.")

os.system("mkdir scan")
os.chdir("scan")

index=0

for i in np.linspace(alphaSMZ[0], alphaSMZ[1], nstepsalphaSMZ):
    for k in CFactor:

	os.system("mkdir "+str(index))
	os.chdir(str(index))
	parFile = open("params.dat", 'w')

	parFile.write("alphaSMZ\t"+str(i)+"\n")
	parFile.write("CFactor\t"+str(k))

	parFile.close()

	index += 1
	os.chdir("..")

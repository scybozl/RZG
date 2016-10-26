import os

d="/afs/ipp-garching.mpg.de/home/l/lscyboz/mc"
root="/afs/ipp-garching.mpg.de/home/l/lscyboz/"

for files in os.listdir(d+"/000"):

  if files.find("8000")!=-1: continue
  if files.find("used_params")!=-1: continue
  if files.find("unnorm")!=-1: continue
  
  parent="tuning_"+files.split("7000_")[1].split("_1.150000e-01")[0]+"/"


  os.chdir(root)
  os.system("mkdir -p "+parent+"mc/")
  os.chdir(parent+"mc")

  os.system("mkdir "+parent+"ref/")

  for i in range(10):
	
	files1=files.split("1.150000e-01")[0]+"*"
	print files1

	files2="MC_Herwig_NLO_8000_"+files1.split("7000_")[1]

	dest=root+parent+"mc/00"+str(i)

	os.system("mkdir 00"+str(i))
	os.system("cp "+d+"/00"+str(i)+"/"+files1+"177.yoda "+dest)
	os.system("cp "+d+"/00"+str(i)+"/"+files2+"253.yoda "+dest)
	os.system("yodamerge "+dest+"/*177.yoda "+dest+"/*253.yoda -o "+dest+"/out.yoda")

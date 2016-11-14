import os

d="/afs/ipp-garching.mpg.de/home/l/lscyboz/mc_31.10.16"
root="/afs/ipp-garching.mpg.de/home/l/lscyboz/"
rivet=d+"/ref/*.aida"


os.system("source "+root+"professor_setup.sh")
os.system("source "+root+"rivet.sh")


for files in os.listdir(d+"/000"):

  if files.find("8000")!=-1: continue
  if files.find("used_params")!=-1: continue
  if files.find("unnorm")!=-1: continue
  
  parent="tuning_"+files.split("7000_")[1].split("_1.150000e-01")[0]+"/"


  os.chdir(root)
  os.system("mkdir -p "+parent+"mc/")
  os.system("mkdir -p "+parent+"ref/")
  os.system("cp "+rivet+" "+parent+"ref/")

  os.chdir(d)
  os.system("cp hypercube.params "+root+parent)
  os.system("cp "+d+"/myweights "+root+parent)

  for i in range(10):
	
	os.chdir(root+parent+"mc")

	files1=files.split("1.150000e-01")[0]+"*"
	print files1

	files2="MC_Herwig_NLO_8000_"+files1.split("7000_")[1]
	params="used_params"
 
	dest=root+parent+"mc/00"+str(i)

	os.system("mkdir "+dest)
	os.system("cp "+d+"/00"+str(i)+"/"+files1+"177.yoda "+dest)
	os.system("cp "+d+"/00"+str(i)+"/"+files2+"253.yoda "+dest)
	os.system("cp "+d+"/00"+str(i)+"/"+params+" "+dest)
	os.system("yodamerge "+dest+"/*177.yoda "+dest+"/*253.yoda -o "+dest+"/out.yoda")

	os.chdir(root)
	os.system("./yoda2aida "+dest+"/out.yoda "+dest+"/out.aida")

  os.chdir(root+parent)
  os.system("prof-envelopes --datadir .")
#  os.chdir("envelopes")
#  os.system("make-plots *.dat")

  os.system("prof-runcombs --datadir . -c 0:1 -o runcombs.dat")
  os.system("prof-interpolate --datadir . --weights myweights --runs runcombs.dat --ipol quadratic")
  os.system("prof-tune --datadir . --runs runcombs.dat --weights myweights --ipol quadratic")
  

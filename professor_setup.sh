# Setup environment for Professor usage.
# Just source this file in your current shell session via
# source /afs/cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/setup.sh

# Load an up-to-date compiler version
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5-gcc43-opt/setup.sh

# Load a current Rivet installation for make-plots and rivet-config.
#source /afs/cern.ch/sw/lcg/external/MCGenerators/rivet/1.6.0/x86_64-slc5-gcc43-opt/rivetenv.sh
source /afs/ipp-garching.mpg.de/home/l/lscyboz/Herwig-7.0.3/src/Rivet-2.4.0/rivetenv.sh

# setup up-to-date Python version
#PATH=/afs/cern.ch/sw/lcg/external/Python/2.6.6/x86_64-slc5-gcc43-opt/bin:${PATH}
# setup Professor binaries
PATH=/afs/cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/bin:${PATH}
export PATH
# setup 3rd party python modules
PYTHONPATH=/usr/lib64/python2.6/site-packages/:${PYTHONPATH}
# setup Professor imports
PYTHONPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/lib/python2.6/site-packages:${PYTHONPATH}
export PYTHONPATH

export PYTHONPATH=/usr/lib/python2.6/site-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/usr/lib/python2.6/site-packages/:$LD_LIBRARY_PATH

# Include ROOT for Minuit access even if ROOTSYS is set to get the Minuit version PyMinuit was compiled with.
#export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/lib:/afs/cern.ch/sw/lcg/external/Python/2.6.5/x86_64-slc5-gcc43-opt/lib:/afs/cern.ch/sw/lcg/app/releases/ROOT/5.27.02/x86_64-slc5-gcc43-opt/root/lib:${LD_LIBRARY_PATH}

export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.27.02/x86_64-slc5-gcc43-opt/root/lib:${LD_LIBRARY_PATH}

export LD_LIBRARY_PATH=/afs.cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/lib:${LD_LIBRARY_PATH}

export LD_LIBRARY_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/PyMinuit/lib64/python2.6/site-packages/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/Minuit-1.7.9/lib/:$LD_LIBRARY_PATH
export PYTHONPATH=/afs/ipp-garching.mpg.de/home/l/lscyboz/PyMinuit/lib64/python2.6/site-packages/:$PYTHONPATH

# Source the tab completion script if possible
if (complete &> /dev/null); then
    test -e "/afs/cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/prof-completion" && source "/afs/cern.ch/sw/lcg/external/MCGenerators/professor/1.3.0/x86_64-slc5-gcc43-opt/prof-completion"
fi

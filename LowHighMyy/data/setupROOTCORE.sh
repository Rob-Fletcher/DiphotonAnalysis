export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

export MYROOTCOREDIR=$(pwd)

localSetupPyAMI
localSetupFAX
# set here your favourite XrootD access point
#export STORAGEPREFIX=root://gridftp-a1-1.mi.infn.it:1094/
localSetupDQ2Client testing-SL6 --rucioVersion testing-SL6 --skipConfirm
localSetupPandaClient
# If you want to try PoD
#localSetupPoD --skipConfirm PoD-3.16p1-python2.7-x86_64-slc6-gcc47-boost1.55

rcSetup

voms-proxy-init -voms atlas

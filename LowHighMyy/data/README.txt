# First download all the packages
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/HiggsPhys/Run2/HGamma/xAOD/LowHighMyy/tags/LowHighMyy-00-01-20 LowHighMyy
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/HiggsPhys/Run2/HGamma/xAOD/HGamAnalysisFramework/tags/HGamAnalysisFramework-00-02-26 HGamAnalysisFramework
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/SUSYTools-00-06-15 SUSYTools
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/ElectronPhotonID/PhotonVertexSelection/tags/PhotonVertexSelection-00-01-01 PhotonVertexSelection
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/ElectronPhotonID/IsolationCorrections/tags/IsolationCorrections-00-00-35 IsolationCorrections
svn co svn+ssh://svn.cern.ch/reps/atlasoff/Reconstruction/RecoTools/IsolationTool/tags/IsolationTool-00-12-00 IsolationTool
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/ElectronPhotonID/ElectronPhotonSelectorTools/tags/ElectronPhotonSelectorTools-00-02-57 ElectronPhotonSelectorTools
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/MCTruthClassifier/tags/MCTruthClassifier-00-01-41 MCTruthClassifier
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/AnalysisCommon/IsolationSelection/tags/IsolationSelection-00-02-01 IsolationSelection

# to make grid work with release 2.3.31
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/D3PDTools/EventLoopGrid/tags/EventLoopGrid-00-00-44 EventLoopGrid

# Only to run with prooflite
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/D3PDTools/RootCore/tags/RootCore-00-04-26 RootCore

# Do the setup
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

localSetupPyAMI
localSetupFAX
localSetupDQ2Client testing-SL6 --rucioVersion testing-SL6 --skipConfirm
localSetupPandaClient

rcSetup Base,2.3.31

voms-proxy-init -voms atlas

# find packages and compile
rc find_packages
rc compile

# Run the analysis (only xAOD output)
./RootCoreBin/bin/x86_64-slc6-gcc48-opt/runDiphotonAnalysisExotic LowHighMyy/data/DiphotonAnalysis.cfg rootfile

# Run the analysis (xAOD + NTUP output)
# to add a dataset to the dataset directory add it in LowHighMyy/data/datasets.txt or do a new datasets file
# If the dataset is derived it should also find the non-derived number of events and attach it as metadata
./RootCoreBin/python/LowHighMyy/create_dataset.py --input-file file_with_list_of_datasets.txt (default is in LowHighMyy/data/datasets.txt)
# Run the analysis (see the options first)
./RootCoreBin/python/LowHighMyy/Diphotonrun.py -h
./RootCoreBin/python/LowHighMyy/Diphotonrun.py

# Grid with a different macro for now because of an non-resolved issue with the original macro, will re-merge the two when resolved
# To run on grid. No need to run create_dataset.py before, the --datasets options also accepts wildcards (be careful!)
./RootCoreBin/python/LowHighMyy/Grid_Diphotonrun.py -h
./RootCoreBin/python/LowHighMyy/Grid_Diphotonrun.py

# to redo all the setup after the first time
source LowHighMyy/data/setupROOTCORE.sh


+++ NOTE: weights, xs and initial events
For the initial number of events (before derivation) cutbookkeeper is used. Note that a number *per file* will be saved, it will be ok to renormalize the sample BUT if you want the total initial number you would have to sum this value for each file. A TParameter (weight_sum_before) is saved with the sum of all the initial events pre-derivation and will be merged, there is the right number for the complete dataset. If the output is not merged however you would have one per file...

To have cross section and filter efficiency for MC sample you have to use either ./create_datasets.py with the --pyAMI option (it only works locally, not on grid) or add a line in file LowHighMyy/data/cross_sections_13TeV.txt so it will be used by the susytool xs tool.

if all works fine you should use the branch "prel_weight" to sum all the sample to have the right weight (pileup, MC) and re-normalization (xs, filter efficiency).

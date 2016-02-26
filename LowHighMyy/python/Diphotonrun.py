#!/usr/bin/env python

#
#	GGTestrun.py
#	Author: Simone Michele Mazza
#	Mail: simone.mazza@mi.infn.it
#
#	Steering macro for a diphoton decaying graviton
#

import ROOT
# Workaround to fix threadlock issues with GUI
ROOT.PyConfig.StartGuiThread = False
## suppress ROOT command line, use python optparse
ROOT.PyConfig.IgnoreCommandLineOptions = True

import logging
logging.basicConfig(level=logging.DEBUG)
from optparse import OptionParser, OptionGroup, OptionValueError
import socket
import time
import fnmatch, re
import datetime
from itertools import izip
import shutil

def set_additional_metadata(sh_input, sh_output_hist, sh_output_data, user_metadata=None):
    for sample_hist in sh_output_hist:
        weight_before = sample_hist.readHist('weight_sum_before').GetVal()
        weight_selected = sample_hist.readHist('weight_sum_selected').GetVal()

        sample_data = sh_output_data.get(sample_hist.name())
        sample_original = sh_input.get(sample_hist.name())
        cross_section_ami = sample_original.getMetaDouble("nc_xs")
        nevents = sample_original.getMetaDouble("nc_nevt")
        filt_eff = sample_original.getMetaDouble("filt_eff")

        for s in (sample_hist, sample_data):
            s.setMetaDouble('weight_sum_before', weight_before)
            s.setMetaDouble('weight_sum_selected', weight_selected)
            s.setMetaDouble('xsection_ami', cross_section_ami)
            s.setMetaDouble('filt_eff', filt_eff)
            s.setMetaDouble('nc_nevt', nevents)

logging.info('loading packages')
shutil.copyfile(ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C'), 'load_packages.C')
ROOT.gROOT.ProcessLine('.x load_packages.C')
logging.info('loaded packages')

parser = OptionParser(epilog=
"""Examples:

  ./GggTestrun.py -w
""")
og_eventloop = OptionGroup(parser, 'EventLoop options', 'Configure how to run')

og_eventloop.add_option('--submitDir', help='dir to store the output', default='submit_dir')
og_eventloop.add_option('--driver', help='select where to run', choices=('direct', 'prooflite', 'proof', 'local', 'LSF', 'condor', 'grid','prun'), default='direct')
og_eventloop.add_option('--proof-server', default='localhost', help='set the server for proof, default if "localhost", you can use for example "pod://"')
og_eventloop.add_option('--output-stream', default='NTUP_output', help='stream output for ntuple (default=NTUP_output)')
og_eventloop.add_option('-w', '--overwrite', action='store_true', default=False, help='overwrite previous submitDir')
og_eventloop.add_option('-s', '--syst', action='store_true', default=False, help='Run with systematics')
og_eventloop.add_option('--all', action='store_true', default=False, help='Save all events')
og_eventloop.add_option('--all-pres', action='store_true', default=False, help='Save all events after preselection')
og_eventloop.add_option('--no-corr-isol', action='store_true', default=False, help='Do not Pt re-correct isolation')
og_eventloop.add_option('--all-photons', action='store_true', default=False, help='Save all photons')
og_eventloop.add_option('-p', action='store_true', default=False, help='save performance tree')
og_eventloop.add_option('--Higgs', action='store_true', default=False, help='Run higgs analysis (Default: Exotic analysis)')
og_eventloop.add_option('--test-bkg', action='store_true', default=False, help='Run searching for jet faking photons')
og_eventloop.add_option('--workers', default='', help='set the number of workers for proof')
og_eventloop.add_option('--sample_type', default='mc', help='analyze data/mc')
og_eventloop.add_option('--datasets', help='dataset to run, it accept wildcard',
                        default='mc15_13TeV.302463.Pythia8EvtGen_A14NNPDF23LO_Ggammagamma_01_3000.merge.DAOD_HIGG1D1.e3985_s2608_r6765_r6282_p2421')
og_eventloop.add_option('--no-dq2', action='store_true', default=False,
                        help="don't use dq2 tools (for example if you don't have a valid certificate)")
og_eventloop.add_option('--datasetsDir', help='dir to look for datasets', default='datasets')
og_eventloop.add_option('--tree-name', help='tree name (default=CollectionTree)', default='CollectionTree')
og_eventloop.add_option('--remote', help='remote xrootd address', default='local')
og_eventloop.add_option('--gridDirect', action='store_true', default=False, help='Do not use gridDirect')
og_eventloop.add_option('--nevents', type=int, help='number of events to process for all the datasets', default=0)
og_eventloop.add_option('--singleDir', help='folder with the input dataset', default='')
og_eventloop.add_option('--conf_file', help='Configuration file for HgammaFramework',
                        default='LowHighMyy/DiphotonAnalysis.cfg')
og_eventloop.add_option('--no-xAOD', action='store_true', default=False, help='Disable the xAOD output-stream')

parser.add_option_group(og_eventloop)

(options, args) = parser.parse_args()

import atexit
@atexit.register
def quite_exit():
    logging.info("quite exiting")
    ROOT.gSystem.Exit(0)

if options.overwrite:
    import shutil
    import os
    if os.path.isdir(os.path.join(os.getcwd(), options.submitDir)):
        logging.info("removing directory %s" % options.submitDir)
        shutil.rmtree(options.submitDir, True)
        # just retry
        if os.path.isdir(os.path.join(os.getcwd(), options.submitDir)):
            logging.warning('removing directory failed: retry to remove directory')
            shutil.rmtree(options.submitDir, False)
    else:
        logging.warning('no directory %s to delete, ignoring' % options.submitDir)

# Set up the job for xAOD access:
ROOT.xAOD.Init().ignore()

# create a new sample handler to describe the data files we use

logging.info('search dataset %s in local folder', options.datasets)
sh_all = ROOT.SH.SampleHandler()
sh = ROOT.SH.SampleHandler()

if options.singleDir != "":
	ROOT.SH.scanDir (sh_all, options.singleDir);
else:
	try:
	    sh_all.load(options.datasetsDir)
	except:
	    logging.error('cannot find datasets folder, run create_dataset first')
	    print """
		**************************************************
		it seems your dataset is not found
		did you run ./create_dataset.py ?
		if yes have a look to the datasets/ directory
		**************************************************
		"""
	    raise

datasets = []

# read file lists
for dataset in [options.datasets]:
    if '.txt' in dataset:
        logging.debug("find txt file as dataset")
        for l in open(dataset):
            datasets.append(l.strip().rstrip('/'))
    else:
        datasets.append(dataset)

for dataset in datasets:
	regex_name = re.compile(fnmatch.translate(dataset))

	for sample in sh_all:
		name = sample.name()
		m = regex_name.search(name)
		if m:
			logging.debug('sample match: %s (%s)', sample, sample.name())
			sh.add(sample)

if len(sh) == 0:
    print "sample handler is empty, no dataset to process"
    exit()

# setting metadata
logging.info('setting metadata')

# set the name of the tree in our files
sh.setMetaString("nc_tree", options.tree_name)

#sh.printContent()
logging.debug("sample size: %d", len(sh))

hostname = socket.gethostname()
# hostname = socket.gethostbyaddr(socket.gethostname())[0]
#if False:
if 'mi.infn.it' in hostname and not options.no_dq2 and options.driver!='grid' and options.driver!='prun' and options.gridDirect:
    logging.info('using GridDirect')
    if options.remote != 'local':
        if options.remote == 'naples':
            logging.info("taking data from Naples")
            ROOT.SH.makeGridDirect(sh, 'INFN-NAPOLI-ATLAS_LOCALGROUPDISK', 'srm://t2-dpm-01.na.infn.it/', 'root://t2-dpm-01.na.infn.it:1094//', True)
        elif options.remote == 'rome':
            logging.info("taking data from Rome")
            ROOT.SH.makeGridDirect(sh, 'INFN-ROMA1_LOCALGROUPDISK', 'srm://grid-cert-03.roma1.infn.it/', 'root://grid-cert-03.roma1.infn.it:1094//', True)
        elif options.remote == 'frascati-localgroupdisk':
            logging.info("taking data from Frascati")
            ROOT.SH.makeGridDirect(sh, 'INFN-FRASCATI_LOCALGROUPDISK', 'srm://atlasse.lnf.infn.it/', 'root://atlasse.lnf.infn.it:1094//', True)
        elif options.remote == 'frascati-datadisk':
            logging.info("taking data from Frascati")
            ROOT.SH.makeGridDirect(sh, 'INFN-FRASCATI_DATADISK', 'srm://atlasse.lnf.infn.it/', 'root://atlasse.lnf.infn.it:1094//', True)
        elif options.remote == 'taiwan-lcg2':
            logging.info("taking data from Taiwan")
            ROOT.SH.makeGridDirect(sh, 'TAIWAN-LCG2_DATADISK', 'srm://f-dpm001.grid.sinica.edu.tw/', 'root://f-dpm000.grid.sinica.edu.tw:11000//', True)
        elif options.remote == 'in2p3':
            logging.info("taking data from in2p3")
            ROOT.SH.makeGridDirect(sh, 'IN2P3-CC_DATADISK', 'srm://ccsrm.in2p3.fr/', 'root://ccxrdatlas.in2p3.fr:1094//', True)
        else:
            logging.info("ERROR: xrootd site not yet implemented")
        raise
    else:
        try:
            logging.info("using INFN-MILANO-ATLASC_LOCALGROUPDISK on /gpfs/storage_1/")
            ROOT.SH.makeGridDirect(sh, 'INFN-MILANO-ATLASC_LOCALGROUPDISK', 'srm://t2cmcondor.mi.infn.it/', '/gpfs/storage_1/', True)  # for new data copy
            sh.printContent()
        except Exception as ex:
            if True: #"a sample with that name already exists" in str(ex):
                logging.info("using INFN-MILANO-ATLASC_PHYS-SM on /gpfs/storage_2/")
                ROOT.SH.makeGridDirect(sh, 'INFN-MILANO-ATLASC_PHYS-SM', 'srm://t2cmcondor.mi.infn.it/', '/gpfs/storage_2/', True)
                sh.printContent()
                # for mc and old data copy
                #ROOT.SH.makeGridDirect(sh, 'INFN-MILANO-ATLASC_PERF-EGAMMA', 'srm://t2cmcondor.mi.infn.it/', '/gpfs/storage_2/', False)

            else:
                logging.error('no dq2 setup')
                print """
			    **************************************************
			    did you setup dq2?
			    localSetupDQ2Client --skipConfirm
			    voms-proxy-init -voms atlas --valid 48:00
			    **************************************************
			    """
                raise

logging.info("sample handler has %d selected samples:", len(sh))
for isample, sample in enumerate(sh):
    logging.info("    %s", sample.name())
    if isample > 10:
        logging.info("    <other %d samples>", len(sh) - 10)
        break

#logging.info("first file: %s", sh[0].fileName(0))
sh.printContent()

logging.info('creating new job')
job = ROOT.EL.Job()

job.sampleHandler(sh)

# setting job options
if options.nevents != 0:
	logging.info('processing only %d events', options.nevents)
	job.options().setDouble(ROOT.EL.Job.optMaxEvents, options.nevents)

if options.p == True:
	job.options().setDouble(ROOT.EL.Job.optPerfTree, 1)

#job.options().setDouble(ROOT.EL.Job.optCacheSize, 10 * 1024 * 1024)

# add our algorithm to the job
logging.info('creating algorithms')
main_alg = ROOT.DiphotonAnalysis()
main_alg.mem = [] # use this to prevent ownwership problems

# ********** Configuration of the algorithm ***********

# creating pileup weight tool
main_alg.outputStreamName = options.output_stream

# Configuration file
conf = ROOT.HG.Config()
conf.addFile(options.conf_file)
# Set the configuration file
main_alg.setConfig(conf)

# using ntuple service (not working with xAOD)
#logging.info('adding ntuple service')
output = ROOT.EL.OutputStream(options.output_stream)
job.outputAdd(output)
ntuple = ROOT.EL.NTupleSvc(options.output_stream)
job.algsAdd(ntuple)

# setting algorithm variables
# TODO: armonize this and configuration files from framework

# common cuts
main_alg.min_nvertex = 1
main_alg.min_ntracks = 3

if options.Higgs == True:
    # cuts for Higgs
    main_alg.AnalysisBranch=1;
    main_alg.isolation_cut = 6.*1E3;
    main_alg.isolation_track_cut = 2.6*1E3;
    main_alg.leading_rel_cut_pt = 0.4; #Relative cut E_T/m_gg
    main_alg.subleading_rel_cut_pt = 0.3;
else:
    # cuts for exotics
    main_alg.AnalysisBranch=11
    main_alg.isolation_cut = 10.*1E3
    main_alg.leading_min_pt = 55.*1E3
    main_alg.subleading_min_pt = 55.*1E3

# configurations
main_alg.correct_isolation = not options.no_corr_isol
main_alg.do_systematics = options.syst
main_alg.save_all_events = options.all
main_alg.save_all_preselection = options.all_pres
main_alg.save_all_photons = options.all_photons
main_alg.test_bkg = False
main_alg.no_xaod = options.no_xAOD
if options.test_bkg == True:
    main_alg.test_bkg = True
    main_alg.all = True

#******************************************************

# adding the algorithm
logging.info('adding algorithms')
job.algsAdd(main_alg)

# configure the driver
logging.info('creating driver')
driver = None
if (options.driver == 'direct'):
    logging.info('running on direct')
    driver = ROOT.EL.DirectDriver()
elif options.driver == 'local':
    logging.info('running Local driver')
    driver = ROOT.EL.LocalDriver()
elif (options.driver == 'prooflite'):
    logging.info('running on prooflite')
    driver = ROOT.EL.ProofDriver()
    if(options.workers != ''):
        driver.numWorkers = abs(int(options.workers))
        logging.info('using %s workers', options.workers)
elif (options.driver == 'proof'):
    logging.info('running on proof')
    driver = ROOT.EL.ProofDriver()
    logging.info("master: %s", options.proof_server)
    driver.proofMaster = options.proof_server
    if(options.workers != ''):
        driver.numWorkers = abs(int(options.workers))
        logging.info('using %s workers', options.workers)
    driver.makeParOptions = ""
    # driver.returnFiles = False
elif (options.driver == 'LSF'):
    logging.info('running on LSF')
    driver = ROOT.EL.LSFDriver()
    driver.shellInit = 'export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase;source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh;localSetupROOT'
elif (options.driver == 'condor'):
    logging.info('running on Condor')
    driver = ROOT.EL.CondorDriver()
    driver.shellInit =  "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase; source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh || exit $?; source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalROOTSetup.sh --skipConfirm || exit $?;"
    driver.nFilesPerJob = 20
    job.options().setDouble(ROOT.EL.Job.optGridNFilesPerJob, 20)
    job.options().setString(ROOT.EL.Job.optCondorConf, 'Requirements = (OpSysMajorVer == 6 && Machine != "proof-08.mi.infn.it")')
elif (options.driver == 'grid'):
    logging.info('running on grid')
    grid_sample_name = r"user.%nickname%.EXOT10.%in:name[2]%"
    grid_sample_name += datetime.datetime.now().strftime('.%Y%m%d%H%M%S')
    logging.debug("grid sample name template: %s" % grid_sample_name)
    driver = ROOT.EL.GridDriver()
    driver.outputSampleName = grid_sample_name
    driver.athenaTag = '17.2.7.4.1,64,slc5'
    #driver.nc_mergeOutput = False
    driver.options().setString("nc_mergeOutput", "false");
    if options.remote == 'naples':
    	logging.info('submitting on Napoli site')
    	driver.site = "ANALY_INFN-NAPOLI"
    elif options.remote == 'rome':
    	logging.info('submitting on Roma site')
    	driver.site = "ANALY_INFN-ROMA1"
    elif options.remote == 'milan':
	    logging.info('submitting on Milano site')
	    driver.site = "ANALY_INFN-MILANO-ATLASC"
elif (options.driver == 'prun'):
    logging.info('running on grid with prun')
    grid_sample_name = "user.%nickname%.%in:name[2]%.%in:name[6]%"
    grid_sample_name += datetime.datetime.now().strftime('.%Y%m%d%H%M%S')
    logging.debug("grid sample name template: %s" % grid_sample_name)
    driver = ROOT.EL.PrunDriver()
    driver.options().setString("nc_mergeOutput", "false")
    job.options().setString(ROOT.EL.Job.optSubmitFlags, "--mergeScript=")
    driver.options().setString("nc_outputSampleName", grid_sample_name)
    driver.athenaTag = '17.6.0,slc6'


# process the job using the driver
logging.info('submit job')
submitDir = options.submitDir+"_"+datetime.datetime.now().strftime('%Y%m%d')+"_"+datetime.datetime.now().strftime('%H%M%S')
driver.submitOnly(job, submitDir)

logging.info('job terminated')
if options.driver != 'prun' and options.driver != 'grid':
    sh_output_hist_all = ROOT.SH.SampleHandler()
    sh_output_data_all = ROOT.SH.SampleHandler()
    sh_output_hist_all.load(submitDir + '/hist')
    sh_output_data_all.load(submitDir + '/output-NTUP_output')

    sh_output_hist = ROOT.SH.SampleHandler()
    sh_output_data = ROOT.SH.SampleHandler()

    for s in sh:
        sh_output_hist.add(sh_output_hist_all.get(s.name()))
        sh_output_data.add(sh_output_data_all.get(s.name()))

    if len(sh) != len(sh_output_hist):
        logging.error('number of output histogram sample != number of input sample')
    if len(sh) != len(sh_output_data):
        logging.error('number of output data sample != number of input sample')
    logging.info('found %d output sample', len(sh_output_data))

    logging.info('setting additional metadata to output')

    set_additional_metadata(sh, sh_output_hist, sh_output_data)
    sh_output_data.setMetaString('nc_tree', options.tree_name)

    sh_output_hist.save(submitDir + '/hist')
    sh_output_data.save(submitDir + '/output-NTUP_output')

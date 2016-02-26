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

logging.info('loading packages')
shutil.copyfile(ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C'), 'load_packages.C')
ROOT.gROOT.ProcessLine('.x load_packages.C')
logging.info('loaded packages')

parser = OptionParser(epilog=
"""Examples:

  ./GggTestrun.py -w
""")
og_eventloop = OptionGroup(parser, 'EventLoop options', 'Configure how to run')
og_eventloop.add_option('--datasets', help='dataset to run, it accept wildcard',
                        default='mc15_13TeV.303740.Sherpa_CT10_2DP20_myy_4000_4500.merge.AOD.e4331_s2608_s2183_r6869_r6282')
og_eventloop.add_option('--all', action='store_true', default=False, help='Save all events')
og_eventloop.add_option('--no-corr-isol', action='store_true', default=False, help='Do not Pt re-correct isolation')
og_eventloop.add_option('--all-pres', action='store_true', default=False, help='Save all events after preselection')
og_eventloop.add_option('--test-bkg', action='store_true', default=False, help='Run searching for jet faking photons')
og_eventloop.add_option('--Higgs', action='store_true', default=False, help='Run higgs analysis (Default: Exotic analysis)')
og_eventloop.add_option('--no-xAOD', action='store_true', default=False, help='Disable the xAOD output-stream')
og_eventloop.add_option('-s', '--syst', action='store_true', default=False, help='Run with systematics')

parser.add_option_group(og_eventloop)

(options, args) = parser.parse_args()

# Set up the job for xAOD access:
ROOT.xAOD.Init().ignore()

# create a new sample handler to describe the data files we use
logging.info('search dataset')
sh = ROOT.SH.SampleHandler()
ROOT.SH.scanDQ2(sh, options.datasets)
sh.setMetaString("nc_tree", "CollectionTree")

logging.debug("sample size: %d", len(sh))

sh.printContent()

logging.info('creating new job')
job = ROOT.EL.Job()

job.sampleHandler(sh)
#job.options().setDouble(ROOT.EL.Job.optMaxEvents, 100)

# add our algorithm to the job
logging.info('creating algorithms')
#main_alg = ROOT.IsolAnalysis_el()
main_alg = ROOT.DiphotonAnalysis()
main_alg.mem = [] # use this to prevent ownwership problems
main_alg.outputStreamName = "NTUP_output"

output = ROOT.EL.OutputStream("NTUP_output")
job.outputAdd(output)
ntuple = ROOT.EL.NTupleSvc("NTUP_output")
job.algsAdd(ntuple)

#job.options().setDouble(ROOT.EL.Job.optMaxEvents, 100)

# setting algorithm variables
# Configuration file
conf = ROOT.HG.Config()
conf.addFile("LowHighMyy/DiphotonAnalysis.cfg")
# Set the configuration file
main_alg.setConfig(conf)

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
    main_alg.correct_isolation = False;
else:
    # cuts for exotics
    main_alg.AnalysisBranch=11
    main_alg.isolation_cut = 10.*1E3
    main_alg.leading_min_pt = 55.*1E3
    main_alg.subleading_min_pt = 55.*1E3
    main_alg.correct_isolation = False

# configurations
main_alg.correct_isolation = not options.no_corr_isol
main_alg.do_systematics = options.syst
main_alg.save_all_events = options.all
main_alg.save_all_preselection = options.all_pres
main_alg.save_all_photons = False
main_alg.test_bkg = options.test_bkg
main_alg.no_xaod = options.no_xAOD

# adding the algorithm
logging.info('adding algorithms')
job.algsAdd(main_alg)

# configure the driver
logging.info('creating driver')
driver = None

logging.info('running on grid with prun')

grid_sample_name = "user.%nickname%.%in:name[2]%.%in:name[3]%"
grid_sample_name += "."+datetime.datetime.now().strftime('%Y%m%d')+"_"+datetime.datetime.now().strftime('%H%M%S')

logging.debug("grid sample name template: %s" % grid_sample_name)
driver = ROOT.EL.PrunDriver()
#driver = ROOT.EL.GridDriver()
#driver = ROOT.EL.DirectDriver()
driver.options().setString("nc_mergeOutput", "false")
job.options().setString(ROOT.EL.Job.optSubmitFlags, "--mergeScript=")
driver.options().setString("nc_outputSampleName", grid_sample_name)
driver.athenaTag = '17.6.0,slc6'

# process the job using the driver
logging.info('submit job')
driver.submitOnly(job, "submit_dir_"+datetime.datetime.now().strftime('%Y%m%d%H%M%S'))

logging.info('job submitted')

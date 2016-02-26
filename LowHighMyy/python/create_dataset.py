#!/usr/bin/env python

import ROOT
# Workaround to fix threadlock issues with GUI
ROOT.PyConfig.StartGuiThread = False

from optparse import OptionParser
import logging
logging.basicConfig(level=logging.INFO)
import socket

try:
    import colorer
except ImportError:
    print "mmm... you don't have the colorer module, well... no colors for you"

# XS*BR for signal samples
def xsbr_sample(mass, k):
	if mass == 500:
		if k == 0.1:
			return 8.50085e-09
	if mass == 1000:
		if k == 0.01:
			return 1.95286e-12
		if k == 0.03:
			return 1.7547e-11
		if k == 0.05:
			return 4.90385e-11
		if k == 0.1:
			return 1.94796e-10
	if mass == 1500:
		if k == 0.1:
			return 1.5171e-11
	if mass == 2000:
		if k == 0.05:
			return 4.73145e-13
		if k == 0.1:
			return 1.87731e-12
	if mass == 2500:
		if k == 0.1:
			return 2.85205e-13
	if mass == 3000:
		if k == 0.1:
			return 4.68533e-14
	return 1

def add_mymetadata(sh, pyAMI=True):
    import re
    import pyAMI
    import pyAMI.client
    import pyAMI.atlas.api as AtlasAPI

    amiClient = pyAMI.client.Client('atlas')
    AtlasAPI.init()
    regex_signal = re.compile(r"mc[0-9]+_[0-9]+TeV\.[0-9]+\.(?P<generator>.+?)_.*?_Ggammagamma_(?P<k>[0-9]{2,3})_(?P<mass>[0-9]{3,4}).recon")
    regex_p = re.compile(r"_p[0-9]{4}")
    regex_DAOD = re.compile(r"\.DAOD_[A-Z,0-9]{4,8}\.")

    for sample in sh:
        name = sample.name()
        print(name + ':')
        original_totalEvents = 1
        x_sec = 1
        genfilt_eff = 1
        totalEvents = 1
        genfilt_eff = float(pyAMI.atlas.api.get_dataset_info(amiClient, name)[0]['approx_GenFiltEff'])

        if "mc" in name:
            try:
                x_sec = float(pyAMI.atlas.api.get_dataset_info(amiClient, name)[0]['crossSection'])
            except:
                continue
            if "DAOD" in name:
                original_name = regex_DAOD.sub('.AOD.', name)
                original_name = regex_p.sub('', original_name)
                original_totalEvents = int(pyAMI.atlas.api.get_dataset_info(amiClient, original_name)[0]['totalEvents'])
                print(original_name + ' non-derived')
            genfilt_eff = float(pyAMI.atlas.api.get_dataset_info(amiClient, name)[0]['approx_GenFiltEff'])
            totalEvents = int(pyAMI.atlas.api.get_dataset_info(amiClient, name)[0]['totalEvents'])
            print('     xs: ',x_sec, ' filter eff: ',genfilt_eff,' nevents: ',totalEvents, ' non-derived nevents: ', original_totalEvents)

        sample.setMetaDouble("filter_eff", genfilt_eff)
        sample.setMetaDouble("total_events", totalEvents)
        sample.setMetaDouble("non_derived_total_events", original_totalEvents)
        sample.setMetaDouble("xs", x_sec)

        m = regex_signal.search(name)
        if not m:
            continue
        generator = m.group('generator')
        mass = float(m.group('mass'))
        k = float("0."+m.group('k'))*10.
        xsbr = xsbr_sample(mass,k)

        if xsbr == 0:
            logging.warning("BR not found for dataset: '%s'", name)

        print('    gen: ' + generator + ' mass: ',mass,' k: ',k, ' xsbr: ',xsbr)

        sample.setMetaString("generator", generator)
        sample.setMetaDouble("mass", mass)
        sample.setMetaDouble("k", k)
        sample.setMetaDouble("xsbr", xsbr)

def create_datasets_dq2(datasets):
    sh_all = ROOT.SH.SampleHandler()
    logging.info('scanning dq2 for %d queries', len(datasets))
    try:
        for n, d in enumerate(datasets):
            print('%d/%d %s' % (n + 1, len(datasets), d))
            ROOT.SH.scanDQ2(sh_all, d)
    except:
        logging.error('no dq2 setup')
        print """
        **************************************************
        did you setup dq2?
        localSetupDQ2Client --skipConfirm
        voms-proxy-init -voms atlas
        **************************************************
    """
        raise
    logging.info("%d different datasets found scanning with dq2", len(sh_all))
    return sh_all


def create_datasets_EOS(path):
    sh_all = ROOT.SH.SampleHandler()
    for p in path:
        logging.info("scanning EOS dir %s", p)
        sh = ROOT.SH.SampleHandler()
        eos = ROOT.SH.DiskListEOS(p, "root://eosatlas/" + p)
        ROOT.SH.scanDir(sh, eos)
        sh_all.add(sh)
    return sh_all


def create_datasets_dir(datasets, path):
    sh_all = ROOT.SH.SampleHandler()
    for p in path:
        logging.info('scanning dir %s', p)
        ROOT.SH.scanDir(sh_all, p)
    return sh_all


def read_txt_files(filename):
    for d in open(filename):
        if '#' in d:
            continue
        s = d.strip()
        if not s:
            continue
        yield s

def generate_dataset_list(only_MC, pattern = None, infile=""):
    f = open(infile, 'r')
    datalist = f.read().splitlines()
    datasets = []
    for data in datalist:
        if "#" in data or data =="":
            continue
        datasets.append(data)
    for dataset in datasets:
        print dataset
    return datasets

def create_datasets(method, only_MC=False, pattern=None, infile=""):
    """
    method is one of AUTO, DQ2, EOS, DIR
    """
    host = socket.gethostbyaddr(socket.gethostname())[0]
    sh_all = None

    if method == "AUTO":
         if 'cern.ch' in host:
             method = 'EOS'
         elif 'mi.infn.it' in host or 'in2p3.fr' in host:
             method = 'DQ2'
         else:
             raise ValueError("AUTO method works only for cern, Milano and in2p3")

    #logging.info('method to create datasets is %s', method)

    if method == 'EOS':
        eos_path = None
        if 'cern.ch' in host:
            eos_path = ("/eos/atlas/atlasgroupdisk/phys-higgs/HSG1/MC12/p1344/",)
            if not only_MC:
                eos_path.append("/eos/atlas/atlasgroupdisk/phys-higgs/HSG1/data12_8TeV/p1341/")
        if not eos_path:
            logging.error("EOS is only for CERN")
            raise ValueError()
        sh_all = create_datasets_EOS(eos_path)
    elif method == 'DQ2':
        datasets = generate_dataset_list(only_MC, pattern, infile)
        sh_all = create_datasets_dq2(datasets)
    elif method == 'DIR':
        path = None
        if 'mi.infn.it' in host:
            path = ("/gpfs/storage_2/atlas/atlasgroupdisk/phys-sm/mc12_8TeV/NTUP_PHOTON/e1386_s1499_s1504_r3658_r3549_p1344/", )
        if not path:
            logging.error("directory scan is only for Milano")
            raise ValueError()
        sh_all = create_datasets_dir(path)
    else:
        raise ValueError("method %s not supported" % method)

    if not sh_all:
        raise ValueError("cannot create dataset")

    return sh_all


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('--method', choices=("EOS", "DIR", "DQ2", "AUTO"), default="AUTO",
                      help="source to create dataset: EOS, DIR, DQ2, default: AUTO, depending on the site")
    parser.add_option('--pyAMI', action='store_true', default=False,
                      help="Use pyAMI (extract automatic metadata")
    parser.add_option('--only-MC', action='store_true', default=False, help='works only for DQ2')
    parser.add_option('--pattern', default=None, help='pattern to filter datasets')
    parser.add_option('--output-directory', default="datasets", help="default=datasets")
    parser.add_option('--input-file', default="LowHighMyy/data/datasets.txt", help="default=LowHighMyy/data/datasets.txt")
    (options, args) = parser.parse_args()

    logging.info('loading packages')
    ROOT.gROOT.ProcessLine('.x $ROOTCOREDIR/scripts/load_packages.C')
    try:
        # check loading libraries
        ROOT.SH
    except AttributeError:
        logging.error('no RootCore setup')
        print """
    **************************************************
    did you setup the environment?
    source ../RootCore/scripts/setup.sh
    **************************************************
    """

    sh_all = create_datasets(options.method, options.only_MC, options.pattern, options.input_file)
    if options.pyAMI:
    	add_mymetadata(sh_all)
    logging.info('saving SH')
    sh_all.save(options.output_directory)

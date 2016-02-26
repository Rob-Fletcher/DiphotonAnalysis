#include "HGamAnalysisFramework/RunUtils.h"
#include "HGamAnalysisFramework/Config.h"
#include <SampleHandler/SampleLocal.h>
#include <EventLoop/ProofDriver.h>
#include <EventLoop/CondorDriver.h>
#include <EventLoopGrid/GridDriver.h>
#include <EventLoopGrid/PrunDriver.h>

// See header file for documentation

namespace HG {

  TString getHelp(TString progName)
  {
    TString help="\n  HELP\n    "+progName+" [CONFIG-FILES] [root files] [KEY: VALUE]\n\n";
    help+="    [CONFIG-FILES] TEnv text files with settings used both by algortihm and job submission.\n";
    help+="    [root files] if argument is of the form *.root* these will be used as input\n";
    help+="    Some basic config KEYs, which follow the same format as the CONFIG-FILE, are:\n";
    help+="      InputFile:     specifies an input ROOT file. Can be used multiple times\n";
    help+="      InputFileList: specifies a text file containing a list of ROOT files to run over\n";
    help+="      GridDS:        specifies a grid data sample to run over\n";
    help+="      OutputDir:     specifies ouput directory. DATE is replaced by date+time [default: "+progName+"_DATE]\n";
    help+="      SampleName:    specifies sample name [default: sample]\n";
    help+="      BaseConfig:    overrides the default base configuration file (calibration smearing etc)\n";
    help+="\n  EXAMPLE:\n";
    help+="    "+progName+" InputFileList: ~/myFilelists/ZH_files_on_EOS.list SampleName: ZH OutputDir: ~/data/output/HDM_ZH BaseConfig: HGamAnalysisFramework/HgammaConfig.cfg\n";
    help+="  or\n    "+progName+" HDM_ZH.cfg\n";
    help+="  in the latter case, the above options can all be specified in the config file (HDM_ZH.cfg).\n";
    return help;
  }
    
  HG::StrV parseArguments(Config *conf, int argc, char **argv)
  {
    StrV files;
    for (int argi=1; argi<argc; ++argi) {
      TString arg=argv[argi];
      
      // 1. Check if argument is a configuration file. If so add it!
      if ( arg.EndsWith(".cfg") || arg.EndsWith(".config" )) {
        conf->addFile(arg); continue;
      }
      
      // 2. If the argument contains ".root", i.e. *.root*, add it to files to run over
      if (arg.Contains(".root")) {
        files.push_back(arg); continue;
      }
      
      // 3. If the arguemnt ends with colon, add a new argument
      if (arg.EndsWith(":")) {
        TString key(arg); key.ReplaceAll(":","");
        conf->setValue(key,argv[++argi]); continue;
      }
      
      // if we get here, the arguemnt cannot be understood
      HG::fatal("Cannot interpret argument: "+arg+getHelp(argv[0]));
    }

    return files;
  }
  
  void runJob(HgammaAnalysis *alg, int argc, char** argv) {
    /*
     *  Fetch the program name and create a help
     */
    TString help = getHelp(argv[0]);
    
    /*
     *   1. Let's read the configuration
     */
    
    Config conf;
    StrV files = parseArguments(&conf, argc, argv);

    // what's the sample name ?
    TString progName = argv[0];
    TString sampleName = conf.getStr("SampleName","sample");
    TString submitDir  = conf.getStr("OutputDir",progName+"_DATE");
    if (submitDir.Contains("DATE")) {
      TDatime now = TDatime();
      submitDir.ReplaceAll("DATE",Form("%d.%.2d.%.2d_%.2d.%.2d.%.2d",
                                       now.GetYear(),now.GetMonth(),now.GetDay(),
                                       now.GetHour(),now.GetMinute(),now.GetSecond()));
      conf.setValue("OutputDir", submitDir.Data());
    }
    
    // print wether the DC14 is activated (in Makefile.RootCore)
#ifdef __DC14__
    printf("  Code setup assumes Rel19 input (DC14). That is, you should be using AnaBase,2.1.X\n");
#else
    printf("  Code setup assumes Rel20 input. That is, you should be using AnaBase,2.3.X\n");
#endif
    
    printf("\n");
    printf("  %20s  %s\n", "SampleName:", sampleName.Data());
    printf("  %20s  %s\n", "OutputDir:", submitDir.Data());

    if ( HG::fileExist(submitDir) )
      fatal("Output directory "+submitDir+" already exist.\n"+
            "  Rerun after deleting it, or specifying a new one, like below.\n"+
            "    OutputDir: NEWDIR");
    
    // Add the confiuration to the algorithm!
    alg->setConfig(conf);
    
    /******************
     *
     *  Handling of the input
     *
     */

    // create a new sample handler to describe the data files we use
    SH::SampleHandler sh;
    
    // this is the basic description of our job
    EL::Job job;

    // Set to branch access mode
    if (conf.getBool("xAODBranchAccessMode", true))
      job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_branch);
    else
      job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_class);

    // Add our algorithm to the job
    job.algsAdd(alg);
    
    /*
     *  1. Individual input files
     */

    // add to the input file list if specified
    if (conf.isDefined("InputFile"))
      for (auto file : conf.getStrV("InputFile")) files.push_back(file);
    
    if (files.size()) {
      // If one or several root files are specified as arugment, run over these
      std::auto_ptr<SH::SampleLocal> sample (new SH::SampleLocal (sampleName.Data()));
      printf("\n  %s\n","InputFiles:");
      for ( TString file : files ) {
        printf("    %s\n",file.Data());
        sample->add (file.Data());
      }
      printf("\n");

      // add the files to the sample handler
      sh.add (sample.release());

    } else if (conf.isDefined("InputFaxDS")) {
      SH::addGrid(sh, conf.getStr("InputFaxDS").Data());
    } else if (conf.isDefined("GridDS")) {

      /*
       *  2. Grid submission
       */

      // Construct sample to run on
      for (TString sample: conf.getStrV("GridDS"))
        SH::addGrid(sh, sample.Data());

      // add the grid dataset here
      sh.setMetaString("nc_tree", "CollectionTree");

      // Set the file pattern, if specified in the config
      if (conf.isDefined("ng_grid_filter"))
        sh.setMetaString("nc_grid_filter", conf.getStr("nc_grid_filter").Data());

      job.sampleHandler(sh);

      if (!conf.isDefined("OutputDS"))
         fatal("To submit to the grid, you MUST define an OutputDS of the form: user.<UserName>.<UniqueString>");

      
      // Check for an OutputDS
      TString outDS = conf.getStr("OutputDS");
      
      printf("  %20s  %s\n", "OutputDS:", outDS.Data());
      
      EL::PrunDriver driver;
      driver.options()->setString("nc_outputSampleName", outDS.Data());
      
      for (auto opt : {"nc_nFiles", "nc_nFilesPerJob", "nc_nJobs"})
        if (conf.isDefined(opt))
          driver.options()->setDouble(opt, conf.getInt(opt));

      for (auto opt : {"nc_excludedSite", "nc_EventLoop_SubmitFlags", "nc_mergeOutput"})
        if (conf.isDefined(opt))
          driver.options()->setString(opt, conf.getStr(opt).Data());

      if (conf.getBool("UsexAODMerge", false)) {
        if (conf.isDefined("nc_EventLoop_SubmitFlags"))
          fatal("You can NOT specify UsexAODMerge and nc_EventLoop_SubmitFlags at the same time. Exiting.");

        driver.options()->setString(EL::Job::optSubmitFlags,
            "--mergeScript=__panda_rootCoreWorkDir/HGamAnalysisFramework/scripts/xaodmerge \%OUT \%IN");
      }
      
      if ( conf.getBool("SubmitAndDownload",false) ) {
        // this will submit and automatically download the dataset
        // once the job is finished.
        driver.submit(job, submitDir.Data());
      } else {
        // this will only submit the job
        // user will have to download using dq2-get
        driver.submitOnly(job, submitDir.Data());
      }      
      
      // after submitting, we're done.
      return;
      
    } else if (conf.isDefined("InputFileList")) {

      /*
       *  3. Input file sumbission
       */
      TString fileList=conf.getStr("InputFileList");
      printf("  %20s  %s\n", "InputFileList:", fileList.Data());
      if (!HG::fileExist(fileList))
	fileList = gSystem->ExpandPathName(("$ROOTCOREBIN/data/"+fileList).Data());
      if (!HG::fileExist(fileList))
        HG::fatal("The input file-list specified: "+conf.getStr("InputFileList")+" doesn't exist.");
      SH::readFileList(sh, sampleName.Data(), fileList.Data());

      printf("\n  %20s  %s\n\n","InputFileList:",fileList.Data());
      // print out the files
      sh.print();
    } else {
      fatal("No input specified!"+help);
    }

    // set the name of the tree in our files
    // in the xAOD the TTree containing the EDM containers is "CollectionTree"
    // what does "nc_tree" mean??
    sh.setMetaString ("nc_tree", "CollectionTree");
    job.sampleHandler(sh);

    // Determine if more than one CPU should be used
    Int_t numWorkers = 1;

    // First check whether ROOTCORE is set to use multiple CPUs
    // Since it's buggy on lxplus, force user to specify in config file for now
    // TString var = gSystem->ExpandPathName("$ROOTCORE_NCPUS");
    // if (var != "") numWorkers = atoi(var.Data());

    // Check if a different number is specified in configuration
    numWorkers = conf.getNum("NumberOfProofWorkers", numWorkers);

    // Set number of events to be processed
    if (conf.isDefined("NumEvents"))
      job.options()->setDouble(EL::Job::optMaxEvents,conf.getInt("NumEvents"));

    // Run analysis
    if (numWorkers != 1 ) {
      printf("\nUsing PROOF lite on %d worker nodes.\n\n", numWorkers);
      EL::ProofDriver driver;
      driver.numWorkers = numWorkers;
      driver.submit(job, submitDir.Data());
    } else {
      EL::DirectDriver driver;
      driver.submit(job, submitDir.Data());
    }
    //-----------------------------------------------------------------
    
    printf("\n  Output directory: %s\n",submitDir.Data());
    printf("  xAOD ouptut in:   %s/data-MxAOD\n\n",submitDir.Data());
    
  } // runJob
  
  
} // namespace HG

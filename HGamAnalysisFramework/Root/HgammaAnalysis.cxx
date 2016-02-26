#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "PATInterfaces/SystematicVariation.h"

#include <HGamAnalysisFramework/HgammaAnalysis.h>
#include <HGamAnalysisFramework/HgammaUtils.h>
#include <HGamAnalysisFramework/HGamVariables.h>
#include <HGamAnalysisFramework/HGamCategories.h>

#include "PhotonVertexSelection/PhotonVertexHelpers.h"
#include "PhotonVertexSelection/PhotonPointingTool.h"

// #include "AsgTools/SgTEventMeta.h"
// #include "xAODMetaData/FileMetaData.h"

#include "TTree.h"
#include "TBranch.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HgammaAnalysis)

HgammaAnalysis :: HgammaAnalysis (const char *name)
: m_event(nullptr)
, m_store(nullptr)
, m_histoStore(nullptr)
, m_name(name)
, m_vertexTool(nullptr)
, m_photonHandler(nullptr)
, m_electronHandler(nullptr)
, m_jetHandler(nullptr)
, m_muonHandler(nullptr)
, m_eventHandler(nullptr)
, m_truthHandler(nullptr)
, m_overlapHandler(nullptr)
, m_etmissHandler(nullptr)
, m_isInit(false)
, m_isAOD(false)
, m_isMC(false)
, m_isData(false)
{
  // Must have no pointer initialization, for CINT
}

EL::StatusCode HgammaAnalysis :: setupJob (EL::Job& job)
{
  job.useXAOD ();
  
  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init(m_name.Data()).ignore(); // call before opening first file
  
  // tell EventLoop about our output ntuple:
  EL::OutputStream out("MxAOD", "xAODNoMeta");
  job.outputAdd(out);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HgammaAnalysis :: createOutput() {
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HgammaAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  m_histoStore = new HistogramStore();

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (MetaData == nullptr)
    HG::fatal("Couldn't find MetaData TTree in event, is this a proper xAOD file? Exiting.");

  m_isAOD  = MetaData->GetBranch("StreamAOD");
  m_isMAOD = !MetaData->GetBranch("TriggerMenu");
  m_isDAOD = !m_isAOD && !m_isMAOD;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  // file name
  Info("changeInput", "Processing file \"%s\"",wk()->inputFile()->GetName());
  Info("changeInput", "This file has %lli entries", wk()->xaodEvent()->getEntries());

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  m_eventCounter = 0;

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  
  HG::VarHandler::getInstance()->setEventAndStore(event(), store());

  m_isMC = eventInfo()->eventType(xAOD::EventInfo::IS_SIMULATION);
  m_isData = !m_isMC;

  // asg::SgTEventMeta meta(asg::SgTEventMeta::InputStore);
  // const xAOD::FileMetaData *metad = nullptr;
  // if (meta.retrieve(metad, "MetaData").isFailure())
  //   HG::fatal("Couldn't retrieve MetaData");
  // std::string val;
  // metad->value(xAOD::FileMetaData::conditionsTag, val);

  // TEnv uses value from first file it's specified in.
  // If specified, read in additional configuration
  if (m_config.isDefined("Include"))
    for (TString cfg : m_config.getStrV("Include"))
      m_config.addFile(cfg);

  // Fill unspecified values from default config, specified here.
  if (!m_config.isDefined("BaseConfig")) {
     HG::fatal("You must specify a base configuration file, option: BaseConfig. Exiting.");
  } else {
     m_config.addFile(m_config.getStr("BaseConfig"));
  }

  // Currently a hack, for passing whether it's an MxAOD to EventHandler
  if (m_isMAOD)
    m_config.setValue("IsMxAOD", "YES");
  else
    m_config.setValue("IsMxAOD", "NO");

  // Print configuration database, if requested
  if (m_config.getBool("HgammaAnalysis.PrintConfig", true)) {
    Info("initialize()", "Printing full configuration:");
    m_config.printDB();
    Info("initialize()", " ");
  }

   // sample name
  TString sampleName = wk()->metaData()->castString("sample_name");
  Info("initialize()", "Sample name = %s", sampleName.Data());
  
  // Vertex selection tool
  CP::PhotonPointingTool *pointTool = new CP::PhotonPointingTool("PointingTool");
  if (pointTool->initialize().isFailure())
    HG::fatal("Failed vertex init");

  ToolHandle<CP::IPhotonPointingTool> tpoint(pointTool);

  m_vertexTool = new CP::PhotonVertexSelectionTool("PhotonVertexSelectionTool");
  CP_CHECK(m_name, m_vertexTool->setProperty("PhotonPointingTool", tpoint));

  if (m_vertexTool->initialize().isFailure())
    HG::fatal("Failed vertex init");

  m_photonHandler = new HG::PhotonHandler("PhotonHandler", m_event, m_store);
  m_photonHandler->initialize(m_config);

  m_electronHandler = new HG::ElectronHandler("ElectronHandler", m_event, m_store);
  m_electronHandler->initialize(m_config);
  
  m_jetHandler = new HG::JetHandler("JetHandler", m_event, m_store);
  m_jetHandler->setAOD(m_isAOD);
  m_jetHandler->initialize(m_config);

  m_muonHandler = new HG::MuonHandler("MuonHandler", m_event, m_store);
  m_muonHandler->initialize(m_config);
  
  m_eventHandler = new HG::EventHandler(m_event, m_store);
  m_eventHandler->initialize(m_config);

  m_truthHandler = new HG::TruthHandler(m_event, m_store);
  m_truthHandler->initialize(m_config);

  m_overlapHandler = new HG::OverlapRemovalHandler();
  m_overlapHandler->initialize(m_config);

  m_etmissHandler = new HG::ETmissHandler("ETmissHandler", m_event, m_store);
  m_etmissHandler->initialize(m_config);
  
  // Check for HgammaAnalysis specific configs
  m_doTwoGoodPhotonsCut = m_config.getBool("HgammaAnalysis.CheckTwoGoodPhotons", true);

  m_doRelPtCut = m_config.getBool("HgammaAnalysis.CheckRelativePtCuts", true);
  m_relPtCut1  = m_config.getNum ("HgammaAnalysis.RelPtFractionFirst" , 0.35);
  m_relPtCut2  = m_config.getNum ("HgammaAnalysis.RelPtFractionSecond", 0.25);

  m_doVertex   = m_config.getBool("HgammaAnalysis.SelectVertex", true);
  m_doHardPV   = m_config.getBool("HgammaAnalysis.UseHardestVertex", false);

  m_doMyyCut   = m_config.getBool("HgammaAnalysis.CheckMyyWindowCut", true);
  m_myyLow     = m_config.getNum("HgammaAnalysis.LowMyyGeV",105.0)*HG::GeV;
  m_myyHigh    = m_config.getNum("HgammaAnalysis.HighMyyGeV",160.0)*HG::GeV;
  
  m_doJetClean  = m_config.getBool("HgammaAnalysis.CheckJetEventCleaning", false);
  m_jetCleanJvt = m_config.getNum("JetHandler.Selection.JVT", -1.0);
  m_jetCleanPt  = 20.0*HG::GeV;
  if (m_jetCleanJvt < 0.0)
    m_jetCleanPt = 50.0*HG::GeV;
  if (m_doJetClean && m_config.getStr("JetHandler.Selection.CutLevel", "") != "LooseBad")
    HG::fatal("Currently you must clean jets with LooseBad to check CheckJetEventCleaning.");


  // Set up trigger matching map
  m_doTrigMatch      = m_config.getBool("EventHandler.CheckTriggerMatching", false);
  m_requiredTriggers = m_config.getStrV("EventHandler.RequiredTriggers");
  for (auto trig: m_requiredTriggers) {
    m_trigMatch[trig] = TrigType::Undefined;
    TString temp = m_config.getStr("EventHandler.TriggerMatchType."+trig, "");
    if (temp == "DiPhoton"      ) m_trigMatch[trig] = TrigType::DiPhoton;
    if (temp == "DiMuon"        ) m_trigMatch[trig] = TrigType::DiMuon;
    if (temp == "DiElectron"    ) m_trigMatch[trig] = TrigType::DiElectron;
    if (temp == "SinglePhoton"  ) m_trigMatch[trig] = TrigType::SinglePhoton;
    if (temp == "SingleMuon"    ) m_trigMatch[trig] = TrigType::SingleMuon;
    if (temp == "SingleElectron") m_trigMatch[trig] = TrigType::SingleElectron;
  }
  
  // Get list of systematic uncertainties
  // Must be done after all helper tools defined
  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  auto recommendedSystematics = registry.recommendedSystematics();

  m_sysList.push_back(CP::SystematicSet());
  for (auto sys: recommendedSystematics) {
    if (sys.name().find("continuous") != std::string::npos) {
      TString sysname = sys.name();
      sysname.ReplaceAll("__continuous", "");

      m_sysList.push_back(CP::SystematicSet());
      m_sysList.back().insert(CP::SystematicVariation(sysname.Data(), 1));

      m_sysList.push_back(CP::SystematicSet());
      m_sysList.back().insert(CP::SystematicVariation(sysname.Data(), -1));
    } else {
      m_sysList.push_back(CP::SystematicSet());
      m_sysList.back().insert(sys);
    }
  }

  // Need this after tools initialized, so that all systematic histograms are made
  // histInitialize();
  createOutput();

  // move to the user class? This is pretty standard
  TFile *file = wk()->getOutputFile ("MxAOD");
  if(!m_event->writeTo(file).isSuccess()){
    Error("initialize()", "Failed to write event to output file!");
    return EL::StatusCode::FAILURE;
  }

  // register all histograms
  for (auto *histo : m_histoStore->getListOfHistograms()) {
    wk()->addOutput(histo);
  }

  m_isInit = true;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if (!m_isInit) 
    HG::fatal("HgammaAnalysis was not initiialized. Did you forget to call HgammaAnalysis::initialize() ?");

  if(m_eventCounter==0) m_startTime = time(nullptr); //time in seconds

  static int progressInterval = config()->getInt("OutputMessage.ProcessedEventsInterval",1000);
  if ( m_eventCounter && m_eventCounter % progressInterval == 0 ) {
    Info("execute()","%i events processed so far  <<<===",
         static_cast< int >(m_eventCounter));
    Info("execute()","Processing rate = %.3f Hz",
         float(m_eventCounter)/(time(nullptr)-m_startTime));
  }
  m_eventCounter++;

  // This function will print the errors, no checking is required
  CP_CHECK(m_name, applySystematicVariation(CP::SystematicSet()));

  // Clear containers which point to objects from previous event
  HG::VarHandler::getInstance()->clearContainers();

  if (m_doVertex) selectVertex();

  setWeightInitial();

  return EL::StatusCode::SUCCESS;
}




EL::StatusCode HgammaAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HgammaAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  Info("finalize()","Finished processing %i events", m_eventCounter);
  double nSecs = time(nullptr)-m_startTime;
  Info("finalize()","Total time elapsed: %dh %dm %ds",
       int(nSecs)/3600,(int(nSecs)%3600)/60,int(nSecs)%60);
  Info("finalize()","Processing rate = %.3f Hz", float(m_eventCounter)/(time(nullptr)-m_startTime));

  SafeDelete(m_vertexTool);
  SafeDelete(m_photonHandler);
  SafeDelete(m_electronHandler);
  SafeDelete(m_jetHandler);
  SafeDelete(m_muonHandler);
  SafeDelete(m_histoStore);
  SafeDelete(m_eventHandler);
  SafeDelete(m_truthHandler);
  SafeDelete(m_overlapHandler);
  SafeDelete(m_etmissHandler);

  TFile *file = wk()->getOutputFile ("MxAOD");

  if(!m_event->finishWritingTo( file ).isSuccess() ) {
    Error("finalize()","Failed to finish writing event to output file!");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}

/// Configures all handlers for systematic variation according to the specified systematic set
CP::SystematicCode HgammaAnalysis::applySystematicVariation(const CP::SystematicSet &sys)
{
  static const char *METHOD = "HgammaAnalysis::applySystematicVariation";
  CP_CHECK( METHOD, m_photonHandler              ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_muonHandler                ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_electronHandler            ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_jetHandler                 ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_etmissHandler              ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, HG::VarHandler::getInstance()->applySystematicVariation(sys) );

  setWeightInitial();

  return CP::SystematicCode::Ok;
}



bool HgammaAnalysis::pass(const xAOD::PhotonContainer   *photons,
                          const xAOD::ElectronContainer *electrons,
                          const xAOD::MuonContainer     *muons,
                          const xAOD::JetContainer      *jets)
{
  if (var::isPassed.exists())
    return var::isPassed();

  if (m_doTwoGoodPhotonsCut) {
    if (photons == nullptr) return false;
    if (!passTwoGoodPhotonsCut(*photons)) return false;
  }

  if (m_doRelPtCut) {
    if (photons == nullptr) return false;
    if (!passRelativePtCuts(*photons)) return false;
  }

  if (m_doMyyCut) {
    if (photons == nullptr) return false;
    if (!passMyyWindowCut(*photons)) return false;
  }

  if (m_doJetClean &&
      !passJetCleaning())
    return false;

  if (m_doTrigMatch &&
      !passTriggerMatch(photons, electrons, muons, jets))
    return false;

  return true;
}

// Return the generator Higgs mass, in GeV, from config specified by
//   GeneratorHiggsMass.MCCHANNELNUMBER
// if not defined, the code is aborted
double HgammaAnalysis::getGeneratorHiggsMass(int mcID)
{
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_higgsMass.find(mcID) == m_higgsMass.end())
    m_higgsMass[mcID] = config()->getNum(Form("GeneratorHiggsMass.%d",mcID));
  return m_higgsMass[mcID];
}

// Return the generator efficiency from config specified by
//   GeneratorEfficiency.MCCHANNELNUMBER
// if not defined, 1.0 is returned
double HgammaAnalysis::getGeneratorEfficiency(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_genEffs.find(mcID) == m_genEffs.end())
    m_genEffs[mcID] = config()->getNum(Form("GeneratorEfficiency.%d",mcID),1.0);
  return m_genEffs[mcID];
}

// Return the kFactor from config specified by
//   GeneratorEfficiency.MCCHANNELNUMBER
// if not defined, 1.0 is returned
double HgammaAnalysis::getKFactor(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_kFactors.find(mcID) == m_kFactors.end())
    m_kFactors[mcID] = config()->getNum(Form("kFactor.%d",mcID),1.0);
  return m_kFactors[mcID];
}

// Return the cross section, in pb, from config specified by
//   CrossSection.MCCHANNELNUMBER
// if not defined, the code is aborted
double HgammaAnalysis::getCrossSection(int mcID)
{
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_crossSections.find(mcID) == m_crossSections.end())
    m_crossSections[mcID] = config()->getNum(Form("CrossSection.%d", mcID));
  return m_crossSections[mcID];
}

TString HgammaAnalysis::getMCSampleName(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_mcNames.find(mcID) == m_mcNames.end())
    m_mcNames[mcID] = config()->getStr(Form("SampleName.%d",mcID));
  return m_mcNames[mcID];
}

int HgammaAnalysis::getNtotalEvents(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_nTotEvents.find(mcID) == m_nTotEvents.end())
    m_nTotEvents[mcID] = config()->getInt(Form("TotalNEvents.%d",mcID));
  return m_nTotEvents[mcID];
}

TH1F* HgammaAnalysis::getCutFlowHistogram(int mcID, TString suffix) {
  // access the initial number of weighed events                                                                                                                                                    
  TString cutFlowName(Form("CutFlow_MC%d%s",mcID,suffix.Data()));
  bool hasMCname = config()->isDefined(Form("SampleName.%d",mcID));
  if (hasMCname) cutFlowName = Form("CutFlow_%s%s",
				    getMCSampleName(mcID).Data(),suffix.Data());
  else Warning("","SampleName.%d not specfied in config!",mcID);
  TH1F *cflowHist = (TH1F*)wk()->inputFile()->Get(cutFlowName);
  if (cflowHist==nullptr)
    HG::fatal("Cannot access cut-flow histogram "+cutFlowName+" in input file");
  return cflowHist;
}

// intial sum of events, including pileup weights
double HgammaAnalysis::getIntialSumOfWeights(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_NevtsInitial.find(mcID) == m_NevtsInitial.end()) {
    // Hard-coding to bin number 3 = ALLEVTS
    m_NevtsInitial[mcID] = getCutFlowHistogram(mcID,"_noDalitz_weighted")->GetBinContent(3);
  }
  return m_NevtsInitial[mcID];
}

double HgammaAnalysis::lumiXsecWeight(double intLumi, int mcID, bool printFirst) {
  if (intLumi < 0) intLumi = config()->getNum("IntegratedLuminosity_fbInv",1.0);
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_weightXsec.find(mcID) == m_weightXsec.end()) {
    double sigma      = getCrossSection(mcID);
    double gen_eff    = getGeneratorEfficiency(mcID);
    double kFactor    = getKFactor(mcID);
    double sumInitial = getIntialSumOfWeights(mcID);
    
    // Hard-coding to bin number 1,2
    double NxAOD      = getCutFlowHistogram(mcID,"_weighted")->GetBinContent(1);
    double NDxAOD     = getCutFlowHistogram(mcID,"_weighted")->GetBinContent(2);

    int Ntot = config()->isDefined(Form("TotalNEvents.%d",mcID)) ? getNtotalEvents(mcID) : -1;
    double skim_eff = NDxAOD / NxAOD;

    m_weightXsec[mcID] = intLumi * 1e3 * sigma * gen_eff * skim_eff * kFactor / sumInitial;
    if (printFirst) {
      printf("\nMC sample %d: %s\n",mcID,config()->getStr(Form("SampleName.%d",mcID),"").Data());      
      printf("  Cross section:                %10.4e pb\n",sigma);
      if (gen_eff!=1.0) printf("  Generator efficiency:         %10.4e\n",gen_eff);
      if (kFactor!=1.0) printf("  k-factor:                     %10.2f\n",kFactor);
      if (Ntot) printf("  N events in AMI:              %10d\n",Ntot);
      printf("  sum w in xAOD:                %10.2f\n",NxAOD);
      printf("  sum w in DxAOD:               %10.2f\n",NDxAOD);
      if (skim_eff!=1.0) 
	printf("  DxAOD efficiency:             %10.2f%%\n",skim_eff*100);
      printf("  Sum of inital event weights:  %10.2f\n\n",sumInitial);
      // L * sigma * eff * kFactor / Nevts
      printf("  Integrated lumi.:             %10.4f fb-1\n",intLumi);
      printf("  N exp. events for analysis:   %10.2e\n",intLumi*1e3*sigma*gen_eff*skim_eff*kFactor);
      printf("  Cross section event weight:   %10.4e\n\n",m_weightXsec[mcID]);
    }
  }
  return m_weightXsec[mcID];
}

void HgammaAnalysis::selectVertex()
{

  // If the event doesn't contain PVs, can't correct anything
  if (!m_event->contains<xAOD::VertexContainer>("PrimaryVertices")) {
    static bool first=true;
    if (first && not m_isMAOD) Warning("selectVertex","No PrimaryVertices container.%s",
                                       " No PV correction can be applied!!");
    first=false;
    return;
  }

  m_photonHandler->setVertexCorrected(false);
  m_jetHandler->setVertexCorrected(false);
  m_electronHandler->setVertexCorrected(false);
  m_muonHandler->setVertexCorrected(false);

  const xAOD::Vertex *vertex = nullptr;
  if (m_doHardPV) {
    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Couldn't retrieve PrimaryVertices, exiting!");

    vertex = xAOD::PVHelpers::getHardestVertex(vertices);
  } else {
    xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer();
    xAOD::PhotonContainer presel  = photonHandler()->applyPreSelection(photons);

    // If there aren't two photons, just use the hardest vertex
    if (presel.size() < 2) {
      const xAOD::VertexContainer *vertices = nullptr;
      if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
        HG::fatal("Couldn't retrieve PrimaryVertices, exiting!");

      vertex = xAOD::PVHelpers::getHardestVertex(vertices);
    } else {
      // Only use the two leading photons
      presel.resize(2);

      // Get the pointed vertex
      m_vertexTool->getVertex(presel, vertex).ignore();
    }
  }

  if (vertex == nullptr) {
    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Couldn't retrieve PrimaryVertices, exiting!");

    vertex = xAOD::PVHelpers::getHardestVertex(vertices);
  }

  // Only set the vertex to use for kinematic corrections if it's not a nullptr
  if (vertex != nullptr) {
    ConstDataVector<xAOD::VertexContainer> *vcont = new ConstDataVector<xAOD::VertexContainer>(SG::VIEW_ELEMENTS);
    vcont->push_back(vertex);

    if (m_store->record(vcont, "HGamVertices").isFailure())
      HG::fatal("Couldn't add HGamVertices to TStore, exiting.");
  }

}



void HgammaAnalysis::setWeightInitial()
{
  // If already set, don't need to do it again
  // Also allows the code to run on MxAODs
  if (var::weightInitial.exists()) return;

  // Determine the initial event weight
  double weight = 1.0;

  if (isMC()) {
    // First MC generator weights
    weight *= eventHandler()->mcWeight();

    // Pileup weight
    weight *= eventHandler()->pileupWeight();

    // z-vertex weight
  }

  var::weightInitial.setValue(weight);
}


/// Get initial event weight: MC, pileup, z-vertex
double HgammaAnalysis::weightInitial()
{
  if (!var::weightInitial.exists()) {
    HG::fatal("Initial event weight not found, did you call HgammaAnalysis::execute() ?");
    return -1.0; // should never get called
  }

  return var::weightInitial();
}



/// Get event weight: initial weight * leading two photon SFs
double HgammaAnalysis::weight()
{
  if (!var::weight.exists()) {
    HG::fatal("You must call setSelectedObjects before retrieving the weight!");
    return -1.0; // should never get called
  }

  return var::weight();
}



/// Get category weight: weight * objects used for category selection
double HgammaAnalysis::weightCategory()
{
  if (!var::weightCategory.exists()) {
    HG::fatal("You must call setSelectedObjects before retrieving the category weight!");
    return -1.0; // should never get called
  }

  return var::weightCategory();
}



/// Set selected collections
void HgammaAnalysis::setSelectedTruthObjects(const xAOD::TruthParticleContainer *photons  ,
                                             const xAOD::TruthParticleContainer *electrons,
                                             const xAOD::TruthParticleContainer *muons    ,
                                             const xAOD::JetContainer           *jets     )
{
  HG::VarHandler::getInstance()->setTruthContainers(photons, electrons, muons, jets);

  setWeightInitial();

  if (!var::weight.exists())
    var::weight.setValue(weightInitial());

  if (!var::category.exists()) {
    var::category.setValue(0);
    var::weightCategory.setValue(weightInitial());
  }
}



/// Set selected collections
void HgammaAnalysis::setSelectedObjects(const xAOD::PhotonContainer   *photons  ,
                                        const xAOD::ElectronContainer *electrons,
                                        const xAOD::MuonContainer     *muons    ,
                                        const xAOD::JetContainer      *jets     )
{
  HG::VarHandler::getInstance()->setContainers(photons, electrons, muons, jets);

  if (!var::weight.exists()) {
    // Determine total weight (leading photonSFs)
    static SG::AuxElement::Accessor<float> scaleFactor("scaleFactor");
    float myweight = weightInitial();
    if (photons != nullptr) {
      for (size_t i = 0; i < photons->size() && i < 2; ++i)
        myweight *= scaleFactor(*photons->at(i));
    }

    var::weight.setValue(myweight);
  }

  if (!var::category.exists()) {
    // Determine the category and weight
    std::pair<int, float> category = HG::getCategoryAndWeight(photons, electrons, muons, jets);

    var::category.setValue(category.first);
    var::weightCategory.setValue(category.second);
  }
}



/// Checks if event level jet cleaning cut is passed
bool HgammaAnalysis::passJetCleaning()
{
  static SG::AuxElement::ConstAccessor<char>  isClean("isClean");
  static SG::AuxElement::ConstAccessor<float> Jvt("Jvt");

  xAOD::JetContainer jets = m_jetHandler->getCorrectedContainer(); 
  for (auto jet: jets) {
    if (m_jetCleanJvt > 0.0) {
      if (Jvt(*jet) > m_jetCleanJvt &&
          jet->pt() > m_jetCleanPt  &&
          !isClean(*jet))
        return false;
    } else {
      if (jet->pt() > m_jetCleanPt  &&
          !isClean(*jet))
        return false;
    }
  }
  return true;
}



//______________________________________________________________________________
bool HgammaAnalysis::passTwoGoodPhotonsCut(const xAOD::PhotonContainer &photons)
{
  if (photons.size() < 2)
    return false;

  xAOD::PhotonContainer leading = photons;
  leading.resize(2);

  xAOD::PhotonContainer sel = photonHandler()->applySelection(leading);
  if (sel.size() < 2)
    return false;

  return true;
}



/// Checks if relative pT cuts for photons are passed
bool HgammaAnalysis::passRelativePtCuts(const xAOD::PhotonContainer &photons)
{
  // If there aren't two photons, the cut fails
  if (photons.size() < 2) return false;

  // Assume Higgs mass from two leading photons
  double myy = (photons[0]->p4() + photons[1]->p4()).M();

  // Check if relative pT cuts are satisfied
  if (photons[0]->pt()/myy < m_relPtCut1) return false;
  if (photons[1]->pt()/myy < m_relPtCut2) return false;

  return true;
}

/// Checks if myy is in the required window
bool HgammaAnalysis::passMyyWindowCut(const xAOD::PhotonContainer &photons)
{
  // If there aren't two photons, the cut fails
  if (photons.size() < 2) return false;
  double myy = (photons[0]->p4() + photons[1]->p4()).M();
  return m_myyLow <= myy && myy < m_myyHigh;
}


bool HgammaAnalysis::passTriggerMatch(const xAOD::PhotonContainer   *photons,
                                      const xAOD::ElectronContainer *electrons,
                                      const xAOD::MuonContainer     *muons,
                                      const xAOD::JetContainer      *jets)
{
  // Check whether at least one passing trigger is matched to selected objects
  for (auto trig: m_requiredTriggers) {
    if (m_eventHandler->passTrigger(trig) &&
        passTriggerMatch(trig.Data(), photons, electrons, muons, jets))
      return true;
  }
  
  return false;
}



//______________________________________________________________________________
bool HgammaAnalysis::passTriggerMatch(const TString &trig,
                                      const xAOD::PhotonContainer   *photons,
                                      const xAOD::ElectronContainer *electrons,
                                      const xAOD::MuonContainer     *muons,
                                      const xAOD::JetContainer      */*jets*/)
{
  switch (m_trigMatch[trig]) {
    case TrigType::Undefined:
      return true;
    case TrigType::DiPhoton:
      return photons && photons->size() > 1 &&
             m_eventHandler->passTriggerMatch_DiPhoton(trig,
                                                       *photons->at(0),
                                                       *photons->at(1));
    case TrigType::DiMuon:
      return muons && muons->size() > 1 &&
             m_eventHandler->passTriggerMatch_DiMuon(trig,
                                                     *muons->at(0),
                                                     *muons->at(1));
    case TrigType::DiElectron:
      return electrons && electrons->size() > 1 &&
             m_eventHandler->passTriggerMatch_DiElectron(trig,
                                                         *electrons->at(0),
                                                         *electrons->at(1));
    case TrigType::SinglePhoton:
      return photons && 
             ( (photons->size() > 0 && m_eventHandler->passTriggerMatch_SinglePhoton(trig, *photons->at(0))) ||
               (photons->size() > 1 && m_eventHandler->passTriggerMatch_SinglePhoton(trig, *photons->at(1))) );
    case TrigType::SingleMuon:
      return muons && 
             ( (muons->size() > 0 && m_eventHandler->passTriggerMatch_SingleMuon(trig, *muons->at(0))) ||
               (muons->size() > 1 && m_eventHandler->passTriggerMatch_SingleMuon(trig, *muons->at(1))) );
    case TrigType::SingleElectron:
      return electrons && 
             ( (electrons->size() > 0 && m_eventHandler->passTriggerMatch_SingleElectron(trig, *electrons->at(0))) ||
               (electrons->size() > 1 && m_eventHandler->passTriggerMatch_SingleElectron(trig, *electrons->at(1))) );
  }

  // If option isn't recognized above, default to failing match
  return false;
}



enum HgammaAnalysis::TrigType HgammaAnalysis::getTriggerType(TString Trigger)
{
  if (m_trigMatch.count(Trigger)>0)
    return m_trigMatch[Trigger]; 
  else
    return TrigType::Undefined;
}

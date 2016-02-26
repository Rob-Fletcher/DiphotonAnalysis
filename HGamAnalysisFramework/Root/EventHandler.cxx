#include "HGamAnalysisFramework/EventHandler.h"

#include "HGamAnalysisFramework/HGamVariables.h"
#include "HGamAnalysisFramework/JetHandler.h"

#include "PhotonVertexSelection/PhotonVertexHelpers.h"

#include "xAODEventShape/EventShape.h"

#ifndef __DC14__
#include "TrigEgammaMatchingTool/TrigEgammaMatchingTool.h"
#endif

namespace HG {

  SG::AuxElement::Decorator<float> EventHandler::myy("myy");
  SG::AuxElement::Decorator<double> EventHandler::PileupWeight("PileupWeight");

  //______________________________________________________________________________
  EventHandler::EventHandler(xAOD::TEvent *event, xAOD::TStore *store)
  : m_event(event)
  , m_store(store)
  , m_grl(nullptr)
  , m_pileupRW(nullptr)
  , m_configTool(nullptr)
  , m_trigDecTool(nullptr)
  , m_trigMuonMatchTool(nullptr)
  , m_trigElectronMatchTool(nullptr)
  , m_checkGRL(false)
  , m_checkTile(false)
  , m_checkLAr(false)
  , m_checkCore(false)
  , m_checkTrig(false)
  { }

  //______________________________________________________________________________
  EventHandler::~EventHandler()
  {
    SafeDelete(m_grl);
    SafeDelete(m_pileupRW);
    SafeDelete(m_configTool);
    SafeDelete(m_trigDecTool);
    SafeDelete(m_trigMuonMatchTool);
#ifndef __DC14__
    SafeDelete(m_trigElectronMatchTool);
#endif

    for (auto dec: m_trigDec) SafeDelete(dec.second);
    for (auto dec: m_trigAcc) SafeDelete(dec.second);
  }

  //______________________________________________________________________________
  EL::StatusCode EventHandler::initialize(Config &config)
  {
    const char *APP_NAME = "HG::EventHandler";

    m_is50ns = config.getBool("Is50ns", false);

    // GRL selection
    std::vector<std::string> vecGRL;
    if (m_is50ns) {
      for (auto grl: config.getStrV("EventHandler.GRL50ns"))
        vecGRL.push_back(PathResolverFindCalibFile(grl.Data()));
    } else {
      for (auto grl: config.getStrV("EventHandler.GRL"))
        vecGRL.push_back(PathResolverFindCalibFile(grl.Data()));
    }

    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    CHECK(m_grl->setProperty("GoodRunsListVec", vecGRL));

    if (m_grl->initialize().isFailure())
      fatal("Failed to initialize GRL tool");

    // Pileup weighting
    std::vector<std::string> confFiles;
    std::vector<std::string> lcalcFiles;
    if (m_is50ns) {
      for (TString val: config.getStrV("EventHandler.PRW.ConfigFiles50ns"))
        confFiles.push_back(val.Data());

      for (TString val: config.getStrV("EventHandler.PRW.LumiCalcFiles50ns"))
        lcalcFiles.push_back(val.Data());
    } else {
      for (TString val: config.getStrV("EventHandler.PRW.ConfigFiles"))
        confFiles.push_back(val.Data());

      for (TString val: config.getStrV("EventHandler.PRW.LumiCalcFiles"))
        lcalcFiles.push_back(val.Data());
    }

    m_prwSF     = config.getNum("EventHandler.PRW.DataScaleFactor", 0.862069);
    int defChan = config.getNum("EventHandler.PRW.DefaultChannel" , 314000  );

    m_pileupRW = new CP::PileupReweightingTool("Pileup");
    CHECK(m_pileupRW->setProperty("ConfigFiles"    , confFiles ));
    CHECK(m_pileupRW->setProperty("LumiCalcFiles"  , lcalcFiles));
    CHECK(m_pileupRW->setProperty("DataScaleFactor", m_prwSF   ));
    CHECK(m_pileupRW->setProperty("DefaultChannel" , defChan   ));
    if (m_pileupRW->initialize().isFailure())
      fatal("Failed to initialize PRW tool");

    const xAOD::EventInfo *eventInfo = nullptr;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }
    m_isMC    = eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION);
    m_isMxAOD = config.getBool("IsMxAOD", false);
    m_jvt     = config.getNum("JetHandler.Selection.JVT", 0.941);

    //if (config.getBool("EventHandler.CheckDalitz"  , true)) m_checkDalitz = true;
    if (config.getBool("EventHandler.CheckGRL"     , true)) m_checkGRL    = true;
    if (config.getBool("EventHandler.CheckTile"    , true)) m_checkTile   = true;
    if (config.getBool("EventHandler.CheckLAr"     , true)) m_checkLAr    = true;
    if (config.getBool("EventHandler.CheckCore"    , true)) m_checkCore   = true;
    if (config.getBool("EventHandler.CheckVertex"  , true)) m_checkVertex = true;
    if (config.getBool("EventHandler.CheckTriggers", true)) m_checkTrig   = true;

    m_truthPtclName = config.getStr("TruthParticles.ContainerName","TruthParticle");

    //if (not m_isMC) m_checkDalitz = false;
    if (m_isMC) m_checkGRL = false;

    m_requiredTriggers = config.getStrV("EventHandler.RequiredTriggers");

    // Trigger decision tool
    if (not m_isMxAOD) {
      m_configTool = new TrigConf::xAODConfigTool("xAODConfigTool");
      ToolHandle<TrigConf::ITrigConfigTool> configHandle(m_configTool);
      if (configHandle->initialize().isFailure()) {
        fatal("Failed to initialize trigger config handle");
      }

      m_trigDecTool = new Trig::TrigDecisionTool("TrigDecisionTool");
      CHECK(m_trigDecTool->setProperty("ConfigTool"     , configHandle  ));
      CHECK(m_trigDecTool->setProperty("TrigDecisionKey","xTrigDecision"));

      if (m_trigDecTool->initialize().isFailure()) {
        fatal("Failed to initialize trigger decision tool");
      }

      // Set up trigger matching map
      for (auto trig: m_requiredTriggers) {
        m_trigMatch[trig] = config.getStr("EventHandler.TriggerMatch."+trig, "");
      }

      // Trigger matching for muons
      m_trigMuonMatchTool = new Trig::TrigMuonMatching("TrigMuonMatching");
      CHECK(m_trigMuonMatchTool->setProperty("TriggerTool", ToolHandle<Trig::TrigDecisionTool>(m_trigDecTool)));
      if (m_trigMuonMatchTool->initialize().isFailure()) {
        fatal("Failed to initialize TrigMuonMatchingTool");
      }

#ifndef __DC14__
      m_trigElectronMatchTool = new Trig::TrigEgammaMatchingTool("TrigEgammaMatchingTool");
#endif

    }

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  float EventHandler::mcWeight()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
      float mcweight = 1.0;

      const std::vector<float> weights = eventInfo->mcEventWeights();
      if (weights.size() > 0) mcweight = weights[0];

      return mcweight;
    }

    return 1.0;
  }

  //______________________________________________________________________________
  float EventHandler::pileupWeight()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
      PileupWeight(*eventInfo) = m_pileupRW->getCombinedWeight(*eventInfo);
      return PileupWeight(*eventInfo);
    }

    return 1.0;
  }

  //______________________________________________________________________________
  bool EventHandler::pass()
  {
    if (var::isPassedBasic.exists())
      return var::isPassedBasic();
    
    // if (m_checkDalitz && isDalitz()) return false;
    if (m_checkTrig   && !passTriggers()) return false;
    return passDQ();
  }
  
  bool EventHandler::passDQ()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("Cannot access EventInfo");
    
    if (m_checkGRL    && !passGRL   (eventInfo)) return false;
    if (m_checkLAr    && !passLAr   (eventInfo)) return false;
    if (m_checkTile   && !passTile  (eventInfo)) return false;
    if (m_checkCore   && !passCore  (eventInfo)) return false;
    if (m_checkVertex && !passVertex(eventInfo)) return false;
    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::isDalitz()
  {
    if (m_isMC) {
      if (var::isDalitzEvent.exists())
        return var::isDalitzEvent();
    
      const xAOD::TruthParticleContainer *truthParticles = nullptr;
      if (m_event->retrieve(truthParticles, m_truthPtclName.Data()).isFailure())
        HG::fatal("Can't access TruthParticleContainer");

      return HG::isDalitz(truthParticles);
    }
    
    // By default (for data) return false
    return false;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggers()
  {
    // Require at least one trigger to be passed
    for (auto trig: m_requiredTriggers) {
      if (passTrigger(trig.Data())) return true;
    }
    return false;
  }

  //______________________________________________________________________________
  StrV EventHandler::getPassedTriggers()
  {
    StrV passedTrigs;
    for (auto trig: m_requiredTriggers) {
      if (passTrigger(trig))
        passedTrigs.push_back(trig);
    }
    return passedTrigs;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_SinglePhoton(const TString &trig,
                                                     const xAOD::Photon &ph)
  {
    return m_trigElectronMatchTool->matchHLT(&ph, trig.Data());
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_DiPhoton(const TString &trig,
                                               const xAOD::Photon &photon1,
                                               const xAOD::Photon &photon2)
  {
    if (passTriggerMatch_SinglePhoton(trig, photon1) &&
        passTriggerMatch_SinglePhoton(trig, photon2))
      return true;

    // Return false if both photons weren't matched to a trigger object
    return false;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_DiMuon(const TString &trig,
                                             const xAOD::Muon &muon1,
                                             const xAOD::Muon &muon2)
  {
    std::pair<bool, bool> isMatchedMuon1, isMatchedMuon2;
    isMatchedMuon1 = std::make_pair(false, false);
    isMatchedMuon2 = std::make_pair(false, false);

    m_trigMuonMatchTool->matchDimuon(&muon1, &muon2, trig.Data(), isMatchedMuon1, isMatchedMuon2);

    return isMatchedMuon1.first && isMatchedMuon2.first;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_SingleMuon(const TString &trig,
                                                 const xAOD::Muon &muon)
  {
    return m_trigMuonMatchTool->match(&muon, trig.Data());
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_SingleElectron(const TString &trig,
                                                     const xAOD::Electron &el)
  {
#ifndef __DC14__
    return m_trigElectronMatchTool->matchHLT(&el, trig.Data());
#else
    return false;
#endif
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_DiElectron(const TString &trig,
                                                 const xAOD::Electron &el1,
                                                 const xAOD::Electron &el2)
  {
    return (passTriggerMatch_SingleElectron(trig, el1) &&
	    passTriggerMatch_SingleElectron(trig, el2));
  }

  //______________________________________________________________________________
  float EventHandler::selectedVertexZ()
  {
    if (var::selectedVertexZ.exists())
      return var::selectedVertexZ();

    if (!m_store->contains<ConstDataVector<xAOD::VertexContainer> >("HGamVertices"))
      return -999;

    ConstDataVector<xAOD::VertexContainer> *hgamvertices = nullptr;
    if (m_store->retrieve(hgamvertices, "HGamVertices").isFailure())
      return -999;

    if (hgamvertices && hgamvertices->size() > 0) {
      var::selectedVertexZ.setValue((*hgamvertices)[0]->z());
      return var::selectedVertexZ();
    }

    return -999;
  }

  //______________________________________________________________________________
  float EventHandler::hardestVertexZ()
  {
    if (var::hardestVertexZ.exists())
      return var::hardestVertexZ();

    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Cannot access PrimaryVertices");

    const xAOD::Vertex *hardest = xAOD::PVHelpers::getHardestVertex(vertices);

    if (hardest) {
      var::hardestVertexZ.setValue(hardest->z());
      return var::hardestVertexZ();
    }

    return -999;
  }

  //______________________________________________________________________________
  float EventHandler::eventShapeDensity()
  {
    if (var::eventShapeDensity.exists())
      return var::eventShapeDensity();

    const xAOD::EventShape *eventShape = nullptr;
    if (m_event->retrieve(eventShape, "Kt4EMTopoEventShape").isFailure())
      HG::fatal("Couldn't retrieve Kt4EMTopoEventShape from TEvent");

    double rho;
    eventShape->getDensity(xAOD::EventShape::Density, rho);

    var::eventShapeDensity.setValue(rho);

    return rho;
  }

  //______________________________________________________________________________
  float EventHandler::mu()
  {
    if (var::mu.exists())
      return var::mu();

    const xAOD::EventInfo *eventInfo = nullptr;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("Cannot access EventInfo in mu()");

    double mu = eventInfo->averageInteractionsPerCrossing();
    if (not m_isMC)
      mu = m_pileupRW->getLumiBlockMu(*eventInfo)*m_prwSF;

    var::mu.setValue(mu);

    return mu;
  }

  //______________________________________________________________________________
  int EventHandler::numberOfPrimaryVertices()
  {
    if (var::numberOfPrimaryVertices.exists())
      return var::numberOfPrimaryVertices();

    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Cannot access PrimaryVertices");

    int NPV = 0;
    for (auto vertex: *vertices) {
      if (vertex->vertexType() == xAOD::VxType::PriVtx ||
          vertex->vertexType() == xAOD::VxType::PileUp)
        NPV++;
    }

    var::numberOfPrimaryVertices.setValue(NPV);

    return NPV;
  }


  //______________________________________________________________________________
  xAOD::Photon* EventHandler::getClosestHLTObject(const TString &trig, const xAOD::Photon &pho)
  {
#ifndef __DC14__
    if (m_trigElectronMatchTool->matchHLT(&pho, trig.Data()))
    {
      return (xAOD::Photon*) m_trigElectronMatchTool->closestHLTObject(&pho, trig.Data());
    }
#endif
    return NULL;
  }

  //______________________________________________________________________________
  xAOD::Electron* EventHandler::getClosestHLTObject(const TString &trig, const xAOD::Electron &el)
  {
#ifndef __DC14__
    if (m_trigElectronMatchTool->matchHLT(&el, trig.Data()))
    {
      return (xAOD::Electron*) m_trigElectronMatchTool->closestHLTObject(&el, trig.Data());
    }
#endif
    return NULL;
  }

  //______________________________________________________________________________
  xAOD::Muon* EventHandler::getClosestHLTObject(const TString &/*trig*/, const xAOD::Muon &/*mu*/)
  {
    // This function intentionally left blank.  No Muon tool method to retrieve closest HLT object
    return NULL;
  }

  //______________________________________________________________________________
  bool EventHandler::passGRL(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        !m_grl->passRunLB(*eventInfo)) {
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passTile(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error) {
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passLAr(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error) {
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passCore(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)) {
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passVertex(const xAOD::EventInfo */*eventInfo*/)
  {
    // Retrieve PV collection from TEvent
    const xAOD::VertexContainer* vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure()) {
      HG::fatal("Couldn't retrieve PrimaryVertices, exiting.");
      return false;
    }

    for (auto vertex: *vertices)
      if (vertex->vertexType() == xAOD::VxType::VertexType::PriVtx ||
          vertex->vertexType() == xAOD::VxType::VertexType::PileUp)
        return true;

    return false;
  }

  //______________________________________________________________________________
  bool EventHandler::passEventClean(const xAOD::JetContainer *allJets)
  {
    if (var::isPassedEventClean.exists())
      return var::isPassedEventClean();

    // Check whether the event passes cleaning or not
    bool isEventClean = true;
    static SG::AuxElement::ConstAccessor<char> isClean("isClean");
    for (auto jet: *allJets) {
      if (jet->pt() > 20*HG::GeV             &&
          JetHandler::passJVTCut(jet, m_jvt) &&
          not isClean(*jet)                  )
        isEventClean = false;
    }
    
    var::isPassedEventClean.setValue(isEventClean);

    return isEventClean;
  }

  //______________________________________________________________________________
  EL::StatusCode EventHandler::writeVars(TString eventName)
  {
    if (eventName == "")
      eventName = HG::VarHandler::getInstance()->getEventInfoName();

    const xAOD::EventInfo *eventInfo = 0;
    if (m_store->retrieve(eventInfo, eventName.Data()).isFailure()) {
      if (m_event->retrieve(eventInfo, eventName.Data()).isFailure()) {
        HG::fatal("EventHandler::write() cannot access " + eventName + ". Exiting.");
      }
    }

    // create containers
    xAOD::EventInfo *output       = new xAOD::EventInfo();
    xAOD::AuxInfoBase * outputAux = new xAOD::AuxInfoBase();
    output->setStore(outputAux);

    *output = *eventInfo;

    // record event info (yes, can be done before setting the actual values)
    if (!m_event->record(output, eventName.Data())) { return EL::StatusCode::FAILURE; }

    eventName += "Aux.";
    if (!m_event->record(outputAux, eventName.Data())) { return EL::StatusCode::FAILURE; }

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  EL::StatusCode EventHandler::writeEventInfo()
  {
    if (m_event->copy("EventInfo").isFailure())
      Warning("EventHandler::writeEventInfo()", "Couldn't copy EventInfo to output.");

    return EL::StatusCode::SUCCESS;
  }

} // namespace HG

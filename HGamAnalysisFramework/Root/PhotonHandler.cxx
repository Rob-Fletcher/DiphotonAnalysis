#include "HGamAnalysisFramework/PhotonHandler.h"

#include "PhotonVertexSelection/PhotonPointingTool.h"

#include "IsolationSelection/IsolationSelectionTool.h"
#include "ElectronPhotonShowerShapeFudgeTool/ElectronPhotonShowerShapeFudgeTool.h"
#include "ElectronPhotonSelectorTools/EGammaAmbiguityTool.h"

namespace HG {
  
  SG::AuxElement::Decorator<float> PhotonHandler::effSF("SF_IDeff");
  SG::AuxElement::Decorator<float> PhotonHandler::effSFunc("SF_IDeff_unc");
  SG::AuxElement::Decorator<float> PhotonHandler::scaleFactor("scaleFactor");
  SG::AuxElement::Decorator<float> PhotonHandler::Ecalib_ratio("Ecalib_ratio");
  SG::AuxElement::Decorator<float> PhotonHandler::r_SL1("r_SL1");
  SG::AuxElement::Decorator<float> PhotonHandler::z_SL1("z_SL1");
  SG::AuxElement::Decorator<float> PhotonHandler::etaS1("eta_s1");
  SG::AuxElement::Decorator<float> PhotonHandler::etaS2("eta_s2");
  SG::AuxElement::Decorator<float> PhotonHandler::cl_eta("cl_eta");
  SG::AuxElement::Decorator<float> PhotonHandler::cl_phi("cl_phi");
  SG::AuxElement::Decorator<float> PhotonHandler::relEreso("relEreso");
  SG::AuxElement::Decorator<char> PhotonHandler::isLoose("isLoose");
  SG::AuxElement::Decorator<char> PhotonHandler::isTight("isTight");
  SG::AuxElement::Decorator<char> PhotonHandler::isConv("isConv");
  SG::AuxElement::Decorator<char> PhotonHandler::passOQ("passOQ");

  //______________________________________________________________________________
  PhotonHandler::PhotonHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
  , m_photonCalibTool(nullptr)
  , m_pointTool(nullptr)
  , m_fudgeMC(nullptr)
  , m_fakeTool(nullptr)
  , m_photonSF(nullptr)
  { }

  //______________________________________________________________________________
  PhotonHandler::~PhotonHandler()
  {
    SafeDelete(m_photonCalibTool);
    SafeDelete(m_pointTool);
    SafeDelete(m_fudgeMC);
    SafeDelete(m_fakeTool);
    SafeDelete(m_photonSF);

    for (auto sel: m_photonSelectors) SafeDelete(sel.second);
    m_photonSelectors.clear();

    for (auto dec: m_pidAcc) SafeDelete(dec.second);
    m_pidAcc.clear();

    for (auto iso: m_isoTools) SafeDelete(iso.second);
    m_isoTools.clear();

    for (auto dec: m_isoAcc) SafeDelete(dec.second);
    m_isoAcc.clear();
  }

  //______________________________________________________________________________
  EL::StatusCode PhotonHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // General configs
    m_isAFII = config.getBool("IsAFII", false);

    // Photon energy scale calibration (data) and extra smearing correction (MC)
    m_photonCalibTool = new CP::EgammaCalibrationAndSmearingTool("EgammaCalibrationAndSmearingTool");

    CP_CHECK(m_name, m_photonCalibTool->setProperty("ESModel", config.getStr(m_name+".Calibration.ESModel").Data()));
    CP_CHECK(m_name, m_photonCalibTool->setProperty("decorrelationModel", config.getStr(m_name+".Calibration.decorrelationModel").Data()));
    if (m_isAFII)
      CP_CHECK(m_name, m_photonCalibTool->setProperty("useAFII", true));

    if (m_photonCalibTool->initialize().isFailure())
      fatal("Failed to initialize EgammaCalibrationAndSmearingTool");
    
    // Photon pointing tool
    m_pointTool = new CP::PhotonPointingTool("PhotonPointingTool");
    if (m_pointTool->initialize().isFailure()) {
      fatal("Couldn't initialize PhotonPointingTool, exiting.");
    }

    // photon selection
    m_doPidCut   = config.getBool(m_name+".Selection.ApplyPIDCut", true);
    m_pidCuts    = config.getStrV(m_name+".Selection.PID");
    if (m_pidCuts.size() < 1) fatal("You must specify at least one PID criteria for photons");

    // photon fudge tool
    if (m_isMC && not m_isAFII) {
      m_fudgeMC  = new ElectronPhotonShowerShapeFudgeTool("FudgeMCTool");
      m_fudgeSet = config.getNum(m_name+".Calibration.FFSet"  , 16);
      CP_CHECK(m_name, m_fudgeMC->setProperty("Preselection", m_fudgeSet));
      if (m_fudgeMC->initialize().isFailure()) {
        fatal("Failed to initialize FudgeMCTool");
      }
    }

    // Ambiguity tool (electron -> photon fakes)
    m_fakeTool = new EGammaAmbiguityTool("EGammaAmbiguityTool");
    if (m_fakeTool->initialize().isFailure())
      fatal("Failed to initialize EGammaAmbiguityTool");

    // isolation tools
    m_doIsoCut = config.getBool(m_name+".Selection.ApplyIsoCut", true);
    m_isoCuts  = config.getStrV(m_name+".Selection.IsoCriteria");
    if (m_doIsoCut && m_isoCuts.size() < 1)
      fatal("Isolation cut requested, but no working point supplied. Exiting!");

    // for each isolation cut
    for (size_t i = 0; i < m_isoCuts.size(); ++i) {
      HG::Iso::IsolationType iso = getIsoType(m_isoCuts[i]);
      m_isoAcc       [iso] = new SG::AuxElement::Accessor<char>(("isIso"+m_isoCuts[i]).Data());

      // first isolation in the list is the default one to apply
      if (i == 0) m_defaultIso = iso;
      
      m_isoTools[iso] = new CP::IsolationSelectionTool(m_isoCuts[i].Data());
      if (iso == HG::Iso::Cone20Higgs) {
        if (m_isoTools[iso]->initialize().isFailure())
          fatal("Failed to initialize IsolationSelectionTool with WP "+m_isoCuts[i]);

        std::vector<std::pair<xAOD::Iso::IsolationType, std::string> > myCuts;
        myCuts.push_back(std::make_pair<xAOD::Iso::IsolationType, std::string>(xAOD::Iso::ptcone20, "0.056*x"));
        myCuts.push_back(std::make_pair<xAOD::Iso::IsolationType, std::string>(xAOD::Iso::topoetcone20, "0.038*(x+30000)"));
        CP_CHECK(m_name, m_isoTools[iso]->addUserDefinedWP("Cone20Higgs", xAOD::Type::Photon, myCuts, "", CP::IsolationSelectionTool::Cut));
      } else {
        CP_CHECK(m_name, m_isoTools[iso]->setProperty("PhotonWP", m_isoCuts[i].Data()));
        if (m_isoTools[iso]->initialize().isFailure())
          fatal("Failed to initialize IsolationSelectionTool with WP "+m_isoCuts[i]);
      }

    }

    // loop over PID selections
    if (m_doPidCut && m_pidCuts.size() < 1)
      fatal("PID cut requested, but no working point supplied. Exiting!");
    for (size_t i = 0; i < m_pidCuts.size(); ++i) {
      unsigned int pid_mask = getPIDmask(m_pidCuts[i]);
      egammaPID::PID pid = getPID(m_pidCuts[i]);
      m_pidAcc[pid] = new SG::AuxElement::Accessor<char>(("is"+m_pidCuts[i]).Data());

      if (i == 0) m_defaultPid = pid;
      TString cfgFile(config.getStr(m_name+".Selection.ConfigFile."+m_pidCuts[i]));

      m_photonSelectors[pid] = new AsgPhotonIsEMSelector("PhotonSelector");
      CP_CHECK(m_name, m_photonSelectors[pid]->setProperty("isEMMask", pid_mask));
      CP_CHECK(m_name, m_photonSelectors[pid]->setProperty("ConfigFile",cfgFile.Data()));
      CP_CHECK(m_name, m_photonSelectors[pid]->setProperty("ForceConvertedPhotonPID", config.getBool(m_name+".Selection."+"ForceConvertedPhotonPID", false)));
      CP_CHECK(m_name, m_photonSelectors[pid]->setProperty("ForceNonConvertedPhotonPID", config.getBool(m_name+".Selection."+"ForceNonConvertedPhotonPID", false)));

      if (!m_photonSelectors[pid]->initialize().isSuccess()) {
        fatal("Failed to initialize PhotonTightIsEMSelector");
      }
    }

    // efficiency scale factor tool
    m_photonSF = new AsgPhotonEfficiencyCorrectionTool("AsgPhotonEfficiencyCorrectionTool");

    std::string file_unc = PathResolverFindCalibFile(config.getStr(m_name+".ScaleFactor.CorrectionFileNameConv").Data());
    if (m_isAFII)
      file_unc = PathResolverFindCalibFile(config.getStr(m_name+".ScaleFactor.CorrectionFileNameConvAFII").Data());

    std::string file_con = PathResolverFindCalibFile(config.getStr(m_name+".ScaleFactor.CorrectionFileNameUnconv").Data());
    if (m_isAFII)
      file_con = PathResolverFindCalibFile(config.getStr(m_name+".ScaleFactor.CorrectionFileNameUnconvAFII").Data());

    CP_CHECK(m_name, m_photonSF->setProperty("CorrectionFileNameConv", file_con));
    CP_CHECK(m_name, m_photonSF->setProperty("CorrectionFileNameUnconv", file_unc));
    // CP_CHECK(m_name, m_photonSF->setProperty("ForceDataType", m_isMC ? 1 : 0)); // FIXME need 3 for AF2
    CP_CHECK(m_name, m_photonSF->setProperty("ForceDataType", 1)); // FIXME need 3 for AF2, for data-driven FullSim = 1 (no data vs. MC)

    if (m_photonSF->initialize().isFailure()) {
      fatal("Failed to initialize AsgPhotonEfficiencyCorrectionTool");
    }

    // Read in configuration information
    m_containerName = config.getStr(m_name+".ContainerName"     , "Photons");

    m_etaCut      = config.getNum(m_name+".Selection.MaxAbsEta"  , 2.37);
    m_ptCut       = config.getNum(m_name+".Selection.PtPreCutGeV", 25.0)*GeV;

    m_crackReject = config.getBool(m_name+".Selection.ApplyCrackRejection", true);
    m_barrelMax   = config.getNum(m_name+".Selection.BarrelMaxAbsEta"  , 1.37);
    m_endcapMin   = config.getNum(m_name+".Selection.EndcapMinAbsEta"  , 1.52);

    m_doAuthor    = config.getBool(m_name+".Selection.ApplyAuthorCut", true);
    m_doQuality   = config.getBool(m_name+".Selection.ApplyQualityCut", true);

    m_doAmbCut   = config.getBool(m_name+".Selection.ApplyAmbiguityCut", false);

    m_correctIso  = config.getBool(m_name+".Selection.CorrectIsoVertex", false);

    // Temporary configs
    m_doCalib = config.getBool(m_name+".Calibration.DoCalibration", true);
    m_doFudge = config.getBool(m_name+".Calibration.DoFudgeFactor", true);
    m_doScale = config.getBool(m_name+".Calibration.DoScaleFactor", true);

    // If this is data, we don't want to fudge it
    if (m_isData) m_doFudge = false;
    
    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  void PhotonHandler::removeAuthorCut(xAOD::PhotonContainer &photons)
  {
    // Check if author cut should be applied
    if (!m_doAuthor) return;

    static SG::AuxElement::Accessor<uint16_t> author("author");

    for (auto ph = photons.begin(); ph != photons.end();) {
      // MxAODs don't save author variable, because cut is already applied
      if (!author.isAvailable(**ph)) {
        ++ph;
        continue;
      }

      // Remove if photon doesn't pass author cut
      if (passAuthorCut(*ph))
        ++ph;
      else                  
        ph = photons.erase(ph);
    }
  }

  //______________________________________________________________________________
  xAOD::PhotonContainer PhotonHandler::getCorrectedContainer()
  {
    // Get Shallow copy from TEvent
    bool calib = false;
    xAOD::PhotonContainer shallowContainer = getShallowContainer(calib);
    removeAuthorCut(shallowContainer);

    if (calib)
      return shallowContainer;

    // get the event info
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      fatal("Cannot access EventInfo");
    }

    // calibrate, correct, and smear photons
    if (m_doCalib) {
      for (auto photon: shallowContainer) {
        calibrateAndSmearPhoton(photon, eventInfo, m_photonCalibTool);
      }
    }

    // Apply PV correction
    HgammaHandler::correctContainerPV(shallowContainer);
    
    // Apply fudge factors, decorate various values
    for (auto photon: shallowContainer) {

      // Fudge photon shape variables
      if (m_doFudge) applyFudgeFactor(photon, eventInfo);


      // Efficiency scale factors
      if (m_doScale) applyScaleFactor(photon, eventInfo);

      decoratePID(*photon);
      decorateIso(*photon);
      decorateOQ(*photon);
      decorateAmbCut(*photon);

      isConv(*photon) = photon->vertex() == nullptr ? false : true;

      static SG::AuxElement::Decorator<unsigned int> dec_isEMTight("isEMTight");
      m_photonSelectors[egammaPID::IsEMTight]->accept(photon);
      dec_isEMTight(*photon) = m_photonSelectors[egammaPID::IsEMTight]->IsemValue();
    }

    // sort the photons
    shallowContainer.sort(comparePt);

    return shallowContainer;
  }

  //______________________________________________________________________________
  xAOD::PhotonContainer PhotonHandler::applyPreSelection(xAOD::PhotonContainer &container)
  {
    xAOD::PhotonContainer selected(SG::VIEW_ELEMENTS);
    for (auto photon: container) {

      // Check if the photon passes pre-selection cuts
      if (!passPreSelection(photon))
        continue;

      selected.push_back(photon);
    }
    return selected;
  }

  //______________________________________________________________________________
  xAOD::PhotonContainer PhotonHandler::applySelection(xAOD::PhotonContainer &container)
  {
    xAOD::PhotonContainer selected(SG::VIEW_ELEMENTS);
    for (auto photon: container) {

      // require photon away from bad calorimeter region
      if (m_doQuality && !passOQCut(photon))
        continue;

      // apply pT and eta selection cuts
      if (!passPtEtaCuts(photon))
        continue;

      // require ambiguity cut for electron fakes
      if (m_doAmbCut && !passAmbCut(photon))
        continue;
      
      // require ID
      if (m_doPidCut && !passPIDCut(photon))
        continue;

      // require Isolation
      if (m_doIsoCut && !passIsoCut(photon))
        continue;

      selected.push_back(photon);
    }
    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode PhotonHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool isAffected = false;
    for (auto var: sys) {
      if (m_photonCalibTool->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
      if (m_photonSF->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
    }

    if (isAffected) {
      CP_CHECK(m_name, m_photonCalibTool->applySystematicVariation(sys));
      CP_CHECK(m_name, m_photonSF       ->applySystematicVariation(sys));
      m_sysName = sys.name() == "" ? "" : "_"+sys.name();
    } else {
      CP_CHECK(m_name, m_photonCalibTool->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_photonSF       ->applySystematicVariation(CP::SystematicSet()));
      m_sysName = "";
    }
    
    return CP::SystematicCode::Ok;
  }

  //______________________________________________________________________________
  void PhotonHandler::calibrateAndSmearPhoton(xAOD::Photon *gam,
                                              const xAOD::EventInfo *evtInfo,
                                              CP::EgammaCalibrationAndSmearingTool *photonCalibTool) {
    // applies calibration to a photon

    // // Set random seed similarly in all code
    // photonCalibTool->setRandomSeed(evtInfo->eventNumber()*100+gam->index()); // done automatically
    
    // Calibrate the photon
    double E_before = gam->e();
    CP::CorrectionCode cc = photonCalibTool->applyCorrection(*gam);
    if (cc==CP::CorrectionCode::Error)
      Error("calibratedPhoton()","Error calibrating current photon");
    if (cc==CP::CorrectionCode::OutOfValidityRange)
      Warning("calibratedPhoton()","Current photon has no valid calibration due to out-of-range");
    
    // decorate the photon with the calibration factor
    // gam->auxdata< float >( "Ecalib_ratio" ) = gam->e()/E_before;
    Ecalib_ratio(*gam) = gam->e()/E_before;
    
    if (gam->caloCluster()==nullptr)
      fatal("PhotonHandler::calibrateAndSmearPhoton, photons has no associated caloCluster");

    etaS1(*gam)    = gam->caloCluster()->etaBE(1);
    etaS2(*gam)    = gam->caloCluster()->etaBE(2);
    cl_eta(*gam)   = gam->caloCluster()->eta();
    cl_phi(*gam)   = gam->caloCluster()->phi();
    relEreso(*gam) = photonCalibTool->resolution(gam->e(),cl_eta(*gam),PATCore::ParticleType::Photon);
  }
    
  //______________________________________________________________________________
  void PhotonHandler::correctPrimaryVertex(const xAOD::Vertex *vertex, xAOD::Photon &gam)
  {
    if (vertex == nullptr) {
      Warning("correctPrimaryVertex()", "Passed vertex = nullptr, skipping correction.");
      return;
    }

    // Correct pt/eta for primary vertex z
    if (m_pointTool->correctPrimaryVertex(gam, vertex->z()).isFailure())
      Warning("correctPrimaryVertex()", "Pointing tool unable to correct photon for PV");

    // float ptcone20 = 0.0, ptcone20_new = 0.0;
    // gam.isolationValue(ptcone20, xAOD::Iso::ptcone20);

    // Update isolation variables for primary vertex
    if (m_correctIso) {
      const std::set<const xAOD::TrackParticle*> tracksToExclude = xAOD::EgammaHelpers::getTrackParticles(&gam, true);
      m_trackIsoTool->decorateParticle(gam, m_isoT, m_corrList, vertex, &tracksToExclude);
    }

    // gam.isolationValue(ptcone20_new, xAOD::Iso::ptcone20);
    // Info("PhotonHandler::correctPrimaryVertex()", "old ptcone20 = %f, new ptcone20 = %f", ptcone20, ptcone20_new);
  }
    
  //______________________________________________________________________________
  void PhotonHandler::applyScaleFactor(xAOD::Photon *gam, const xAOD::EventInfo *evtInfo)
  {
    double _effSF = 1.0, _effSFunc = 1.0;
    if (!etaS2.isAvailable(*gam))
      fatal("PhotonHandler::correctPrimaryVertex etaS2 not availlable for photon");
    double aeta = std::abs(etaS2(*gam));
    if (evtInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        aeta < 2.47 && gam->pt() > 10.0*HG::GeV && gam->caloCluster()->e()/cosh(aeta) > 10.0*HG::GeV) {
      CP::CorrectionCode cc = m_photonSF->getEfficiencyScaleFactor     (*gam, _effSF   );
      if (cc == CP::CorrectionCode::Error)
        Error("applyScaleFactor()", "SF returned error");
      if (cc == CP::CorrectionCode::OutOfValidityRange)
        Warning("applyScaleFactor()", "Current photon has no valid SF due to out-of-range");

      cc = m_photonSF->getEfficiencyScaleFactorError(*gam, _effSFunc);
      if (cc == CP::CorrectionCode::Error)
        Error("applyScaleFactor()", "SF returned error");
      if (cc == CP::CorrectionCode::OutOfValidityRange)
        Warning("applyScaleFactor()", "Current photon has no valid SF due to out-of-range");
    }

    scaleFactor(*gam) = _effSF; 
    effSF(*gam) = _effSF; 
    effSFunc(*gam) = _effSFunc;
  }

  //______________________________________________________________________________
  void PhotonHandler::applyFudgeFactor(xAOD::Photon *gam, const xAOD::EventInfo *evtInfo)
  {
    // Yum... Fudge....

    // Decorate corrected isEMLoose bit and isTight before fudging
    static SG::AuxElement::Decorator<char> dec_isTight_nofudge("isTight_nofudge");
    dec_isTight_nofudge(*gam)   = m_photonSelectors[egammaPID::IsEMTight]->accept(gam);

    static SG::AuxElement::Decorator<unsigned int> dec_isEMTight_nofudge("isEMTight_nofudge");
    dec_isEMTight_nofudge(*gam) = m_photonSelectors[egammaPID::IsEMTight]->IsemValue();

    // Don't actually apply the fudging in AFII
    if (m_isAFII)
      return;

    // 3. fudge if MC
    int isMC = evtInfo->eventType(xAOD::EventInfo::IS_SIMULATION);
    if (isMC) {
      CP::CorrectionCode cc = m_fudgeMC->applyCorrection(*gam);
      if (cc==CP::CorrectionCode::Error)
        Error("calibrateAndSmearAndFudgePhoton()","Fudging returned error");
      if (cc==CP::CorrectionCode::OutOfValidityRange)
	Warning("calibrateAndSmearAndFudgePhoton()","Current photon has no valid fudging due to out-of-range");
    }

  }
  
  //______________________________________________________________________________
  bool PhotonHandler::passAuthorCut(const xAOD::Photon *gam)
  {
    uint16_t author = gam->author();
    if (author & xAOD::EgammaParameters::AuthorPhoton   ) return true;
    if (author & xAOD::EgammaParameters::AuthorAmbiguous) return true;
    return false;
  }

  //______________________________________________________________________________
  void PhotonHandler::decorateOQ(xAOD::Photon &gam)
  {
    if (gam.isGoodOQ(34214))
      passOQ(gam) = true;
    else
      passOQ(gam) = false;
  }

  /// Requires photon to pass preselection, defined by
  /// a) OQ, b) pT>25, |eta_s2| selection, c) Loose PID, and d) electron-photon ambiguity cut
  bool PhotonHandler::passPreSelection(const xAOD::Photon *gam)
  {
    if (!passOQCut(gam)) return false;
    if (!passPtEtaCuts(gam)) return false;
    if (!passPIDCut(gam,egammaPID::IsEMLoose)) return false;
    if (!passAmbCut(gam)) return false;

    return true;
  }
  
  //______________________________________________________________________________
  bool PhotonHandler::passOQCut(const xAOD::Photon *gam)
  {
    if (passOQ.isAvailable(*gam))
      return passOQ(*gam);

    return true;
  }

  //______________________________________________________________________________
  void PhotonHandler::decorateAmbCut(xAOD::Photon &gam)
  {
    static SG::AuxElement::Decorator<char> dec_passAmbCut("passAmbCut");
    dec_passAmbCut(gam) = true;

    if (m_doAmbCut)
      dec_passAmbCut(gam) = m_fakeTool->accept(gam);
  }

  //______________________________________________________________________________
  bool PhotonHandler::passAmbCut(const xAOD::Photon *gam)
  {
    static SG::AuxElement::ConstAccessor<char> acc_passAmbCut("passAmbCut");
    if (acc_passAmbCut.isAvailable(*gam))
      return acc_passAmbCut(*gam);

    return true;
  }

  //______________________________________________________________________________
  void PhotonHandler::decorateIso(xAOD::Photon &gam)
  {
    for (auto dec: m_isoAcc) {
      if (m_isoTools[dec.first]->accept(gam))
        (*dec.second)(gam) = true;
      else
        (*dec.second)(gam) = false;
    }
  }

  //______________________________________________________________________________
  bool PhotonHandler::passIsoCut(const xAOD::Photon *gam, HG::Iso::IsolationType iso)
  {
    /// applies PID cut specified in config file
    if (iso == HG::Iso::Undefined) {
      if (!m_isoAcc[m_defaultIso]->isAvailable(*gam))
        return true;
      return (*m_isoAcc[m_defaultIso])(*gam);
    }
    
    if (m_isoTools.find(iso) != m_isoTools.end()) {
      if (!m_isoAcc[iso]->isAvailable(*gam))
        return true;
      return (*m_isoAcc[iso])(*gam);
    }

    fatal("Isolation cut requested that wasn't specified in config file. Exiting.");
    return false;
  }

  //______________________________________________________________________________
  void PhotonHandler::decoratePID(xAOD::Photon &gam)
  {
    for (auto dec: m_pidAcc) {
      if (m_photonSelectors[dec.first]->accept(&gam))
        (*dec.second)(gam) = true;
      else
        (*dec.second)(gam) = false;
    }
  }

  //______________________________________________________________________________
  bool PhotonHandler::passPIDCut(const xAOD::Photon *gam, egammaPID::PID pid)
  {
    /// applies PID cut specified in config file
    if (pid == egammaPID::LastEgammaPID) {
      if (!m_pidAcc[m_defaultPid]->isAvailable(*gam))
        return true;
      return (*m_pidAcc[m_defaultPid])(*gam);
    }
    
    if (m_pidAcc.find(pid) != m_pidAcc.end()) {
      if (!m_pidAcc[pid]->isAvailable(*gam))
        return true;
      return (*m_pidAcc[pid])(*gam);
    }

    fatal("PID cut requested that wasn't specified in config file. Exiting.");
    return false;
  }

  //______________________________________________________________________________
  bool PhotonHandler::passPtEtaCuts(const xAOD::Photon *gam)
  {
    /// applies kinematic preselection cuts: not-in-crack + pT cut

    // eta cuts
    static SG::AuxElement::ConstAccessor<float> eta_s2("eta_s2");

    // if (!etaS2.isAvailable(*gam))
    if (!eta_s2.isAvailable(*gam))
      fatal("PhotonHandler::passPtEtaCuts etaS2 not availalbe!");
    double abs_eta_s2 = fabs(eta_s2(*gam));
    //fabs(gam->caloCluster()->etaBE(2));
    if (abs_eta_s2 > m_etaCut) return false;
    if (m_crackReject && (abs_eta_s2 > m_barrelMax && abs_eta_s2 < m_endcapMin)) return false;
  
    // pt cuts
    if (gam->pt() < m_ptCut) return false;
  
    return true;
  }
  
  //______________________________________________________________________________
  HG::Iso::IsolationType PhotonHandler::getIsoType(TString isoName) {
    if      (isoName == "Tight") return HG::Iso::Tight;
    else if (isoName == "Loose") return HG::Iso::Loose;
    else if (isoName == "Gradient") return HG::Iso::Gradient;
    else if (isoName == "Cone40CaloOnly") return HG::Iso::Cone40CaloOnly;
    else if (isoName == "Cone40") return HG::Iso::Cone40;
    else if (isoName == "Cone20") return HG::Iso::Cone20;
    else if (isoName == "Cone20Higgs") return HG::Iso::Cone20Higgs;
    else if (isoName == "UserDefined") return HG::Iso::UserDefined;
    else fatal("Isolation "+isoName+" read from: "+
               m_name+".Selection.IsolationCriteria is not Tight, Gradient, Loose, or UserDefined. Exiting.");
    return HG::Iso::Undefined;
  }

  //______________________________________________________________________________
  unsigned int PhotonHandler::getPIDmask(TString pidName) {
    if      (pidName=="Tight")  return egammaPID::PhotonTight;
    else if (pidName=="Medium") return egammaPID::PhotonMedium;
    else if (pidName=="Loose")  return egammaPID::PhotonLoose;
    else fatal("PhotonHandler::getPIDmask  Don't know about pid "+pidName);
    return 0;
  }
  
  //______________________________________________________________________________
  egammaPID::PID PhotonHandler::getPID(TString pidName) {
    if      (pidName=="Tight")  return egammaPID::IsEMTight;
    else if (pidName=="Medium") return egammaPID::IsEMMedium;
    else if (pidName=="Loose")  return egammaPID::IsEMLoose;
    else fatal("PhotonHandler::getPIDmask  Don't know about pid "+pidName);
    return egammaPID::LastEgammaPID;
  }
  
  
  
  //______________________________________________________________________________
  double PhotonHandler::PVz_corrected_eta(double eta_s1, double PVz)
  {
    /// Returns the photon eta corrected for a different primary vertex
    /// geometry: cosh(eta) = p/pT, sinh(eta) = pz/pT = z/rT
    ///   where rT is the transverse distance
    
    // Barrel: etaCorr = asinh( (z-PVz)/r ) = asinh( sinh(eta) - PVz/r )
    if ( fabs(eta_s1)<1.5 )
      return asinh( sinh(eta_s1) - PVz/BarrelR_s1(eta_s1) );
  
    // Endcap:  etaCorr = asinh( (z-PVz)/r ) = [ r = z/sinh(eta) ] = asinh( sinh(eta)*(1-PVz/z) )
    return asinh( sinh(eta_s1)*(1.0-PVz/EndcapZ_s1(eta_s1)) );
  }
  
  //______________________________________________________________________________
  double PhotonHandler::BarrelR_s1(double eta_s1)
  {
    /// barrel r coordinate [mm] in first sampling layer

    double aeta_s1=fabs(eta_s1);
    if (aeta_s1>1.5) fatal(Form("BarrelR_s1 eta=%.3f not in Barrel!",eta_s1));
    if (aeta_s1<0.8) return (1558.859292 - 4.990838*aeta_s1 - 21.144279*aeta_s1*aeta_s1);
    // 0.8<|eta|<1.5
    return 1522.775373 + 27.970192*aeta_s1 - 21.104108*aeta_s1*aeta_s1;
  }
  
  //______________________________________________________________________________
  double PhotonHandler::EndcapZ_s1(double eta_s1)
  {
    /// endcap z coordinate [mm] in first sampling layer

    double aeta_s1=fabs(eta_s1);
    double z = (aeta_s1<1.5) ? 12453.297448 - 5735.787116*aeta_s1 : 3790.671754;
    if (eta_s1<0) z*=-1;
    return z;
  }
  
  //______________________________________________________________________________
  double PhotonHandler::r_s1(double eta_s1)
  {
    /// r-coordinate in first sampling layer

    if (fabs(eta_s1)<1.5) return BarrelR_s1(eta_s1);
    return EndcapZ_s1(eta_s1)/sinh(eta_s1);
  }
  
  //______________________________________________________________________________
  double PhotonHandler::z_s1(double eta_s1)
  {
    /// z-coordinate in first sampling layer

    if (fabs(eta_s1)<1.5) return BarrelR_s1(eta_s1)*sinh(eta_s1);
    return EndcapZ_s1(eta_s1);
  }
  
  //______________________________________________________________________________
  const xAOD::Photon* PhotonHandler::findTrigMatchedPhoton(const xAOD::Photon *photon,
                                                    const xAOD::PhotonContainer* trigPhotons,
                                                    bool debug)
  {
    const xAOD::Photon *match=0;
    if (photon->caloCluster()==nullptr)
      fatal("PhotonHandler::findTrigMatchedPhoton calo cluster not availalbe - can't do trigger matching");
    const xAOD::CaloCluster &clust = *photon->caloCluster();
    double bestDpt = 999*GeV;
    for ( const xAOD::Photon *trigPhoton : *trigPhotons) {
      double dr  = clust.p4().DeltaR(trigPhoton->caloCluster()->p4());
      if (dr<0.05) {
        double dpT = fabs(clust.pt()-trigPhoton->caloCluster()->pt());
        if (debug&&bestDpt!=999*GeV) {
          printPhoton(photon,"Current photon match two trigger photons");
          printPhoton(match,"Matching trig-photon 1");
          printPhoton(trigPhoton,"Matching trig-photon 2");
        }
        if (dpT<bestDpt) { match=trigPhoton; bestDpt=dpT; }
      }
    }
    return match;
  }

  //______________________________________________________________________________
  void PhotonHandler::printPhoton(const xAOD::Photon *gam, TString comment)
  {
    // prints details about the photon

    printf("Photon %2zu  %s\n",
           gam->index(),comment.Data());
    
    // print the 4-vector
    printf("   (pT,eta,phi,m) = (%5.1f GeV,%6.3f,%6.3f,%4.1f GeV)\n",
           gam->pt()/GeV,gam->eta(),gam->phi(),gam->m()/GeV);
    
    // print (eta,phi) for first and second sampling layer
    //  const xAOD::CaloCluster* cluster  = gam->caloCluster();
    //  printf("   (eta,phi) s1: (%6.3f,%6.3f), s2: (%6.3f,%6.3f), cl: (%6.3f,%6.3f)\n",
    //         cluster->etaBE(1),cluster->phiBE(2),cluster->etaBE(2),cluster->phiBE(2),cluster->eta(),cluster->phi());
    
    // print some more information
    TString str;
    if ( Ecalib_ratio.isAvailable(*gam)) str+=Form("   calibFactor = %.3f", Ecalib_ratio(*gam));
    if ( r_SL1.isAvailable(*gam)       ) str+=Form("   r_s1 = %.0f mm", r_SL1(*gam));
    if ( z_SL1.isAvailable(*gam)       ) str+=Form("   z_s1 = %.0f mm", z_SL1(*gam));
    if (str.Sizeof()) printf("%s\n",str.Data());

    if (isLoose.isAvailable(*gam)&&isTight.isAvailable(*gam))
      printf("  isLoose: %d, isTight: %d\n",isLoose(*gam),isTight(*gam));
    
    float ptcone20=-99e3, etcone20=-99e3;
    bool has_iso  = gam->isolationValue(ptcone20,xAOD::Iso::ptcone20);
    has_iso      &= gam->isolationValue(etcone20,xAOD::Iso::topoetcone20);
    if (has_iso) {
      printf("  isolation (ptcone20,topoetcone20) = (%5.2f,%5.2f) GeV\n",
             ptcone20*HG::invGeV,etcone20*HG::invGeV);
    }
    
  }
}

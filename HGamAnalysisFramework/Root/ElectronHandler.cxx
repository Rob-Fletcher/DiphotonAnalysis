#include "HGamAnalysisFramework/ElectronHandler.h"

#ifdef __DC14__
#include "ElectronIsolationSelection/IsolationSelectionTool.h"
#else
#include "IsolationSelection/IsolationSelectionTool.h"
#endif

namespace HG {
  
  SG::AuxElement::Decorator<float> ElectronHandler::effIDSF("SF_IDeff");
  SG::AuxElement::Decorator<float> ElectronHandler::effRecoSF("SF_Recoeff");
  SG::AuxElement::Decorator<float> ElectronHandler::scaleFactor("scaleFactor");
  SG::AuxElement::Decorator<float> ElectronHandler::Ecalib_ratio("Ecalib_ratio");
  SG::AuxElement::Decorator<float> ElectronHandler::Ereso("Ereso");
  SG::AuxElement::Decorator<char>  ElectronHandler::passIPCut("passIPCut");
  SG::AuxElement::Decorator<char>  ElectronHandler::isTight("isTight");
  SG::AuxElement::Decorator<char>  ElectronHandler::isMedium("isMedium");
  SG::AuxElement::Decorator<char>  ElectronHandler::isLoose("isLoose");
  
  //______________________________________________________________________________
  ElectronHandler::ElectronHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
  , m_electronCalibTool(NULL)
  , m_electronIDSF(NULL)
  , m_electronRecoSF(NULL)
  { }

  //______________________________________________________________________________
  ElectronHandler::~ElectronHandler()
  {
    SafeDelete(m_electronCalibTool);
    SafeDelete(m_electronIDSF);
    SafeDelete(m_electronRecoSF);

    for (auto sel: m_electronSelectors) SafeDelete(sel.second);
    m_electronSelectors.clear();
    
    for (auto dec: m_electronSelDecorators) SafeDelete(dec.second);
    m_electronSelDecorators.clear();

    for (auto iso: m_isoTools) SafeDelete(iso.second);
    m_isoTools.clear();
    
    for (auto dec: m_isoDecorators) SafeDelete(dec.second);
    m_isoDecorators.clear();
  }

  //______________________________________________________________________________
  EL::StatusCode ElectronHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // egamma energy scale (data) and extra smearing correction (MC)
    m_electronCalibTool = new CP::EgammaCalibrationAndSmearingTool("ElectronCalibrationAndSmearingTool");

    // Let's loop over electron calibration properties and set them with the values
    // we read in from the config database (see  HGamAnalysisFramework/data/HgammaConfig.cfg)
#ifdef __DC14__    
    for (TString prop : {"ESModel", "ResolutionType"}) {
      CP_CHECK(m_name, m_electronCalibTool->setProperty(prop.Data(), config.getStr(m_name+".Calibration."+prop).Data()));
    }
#else
    for (TString prop : {"ESModel", "decorrelationModel"}) {
      CP_CHECK(m_name, m_electronCalibTool->setProperty(prop.Data(), config.getStr(m_name+".Calibration."+prop).Data()));
    }
#endif
    
    if (m_electronCalibTool->initialize().isFailure()) {
      fatal("Failed to initialize EgammaCalibrationAndSmearingTool");
    }

    //electron selection
    m_doPidCut   = config.getBool(m_name+".Selection.ApplyPIDCut", true);
    m_pidCuts    = config.getStrV(m_name+".Selection.PID");
    
    if (m_pidCuts.size() < 1) fatal("You must specify at least one PID criteria for electrons");

    // loop over PID selections
    if (m_doPidCut && m_pidCuts.size() < 1)
      fatal("Electron PID cut requested, but no working point supplied. Exiting!");
    
    for (size_t i = 0; i < m_pidCuts.size(); ++i) {
      TString pid = m_pidCuts[i];
      m_electronSelDecorators[pid] = new SG::AuxElement::Accessor<char>(("is"+m_pidCuts[i]).Data());
      
      if (i == 0) m_defaultPid = pid;
      
      TString cfgFile(config.getStr(m_name+".Selection.ConfigFile."+m_pidCuts[i]));

      m_electronSelectors[pid] = new AsgElectronLikelihoodTool("ElectronLikelihoodTool");

      CP_CHECK(m_name, m_electronSelectors[pid]->setProperty("primaryVertexContainer",config.getStr("PrimaryVertices.ContainerName","PrimaryVertices").Data()));
      CP_CHECK(m_name, m_electronSelectors[pid]->setProperty("ConfigFile",cfgFile.Data()));
      
      if (!m_electronSelectors[pid]->initialize().isSuccess()) {
        fatal(TString::Format("Failed to initialize %sElectronLikelihoodTool",pid.Data()));
      }
    }
    
    // isolation tools
    m_doIsoCut = config.getBool(m_name+".Selection.ApplyIsoCut", true);
    m_isoCuts  = config.getStrV(m_name+".Selection.IsoCriteria");
    if (m_doIsoCut && m_isoCuts.size() < 1)
      fatal("Isolation cut requested, but no working point supplied. Exiting!");
    for (size_t i = 0; i < m_isoCuts.size(); ++i) {
      HG::Iso::IsolationType iso = getIsoType(m_isoCuts[i]);
      m_isoDecorators[iso] = new SG::AuxElement::Accessor<char>(("isIso"+m_isoCuts[i]).Data());

      // first isolation in the list is the default one to apply
      if (i == 0) m_defaultIso = iso;

      m_isoTools[iso] = new CP::IsolationSelectionTool(m_isoCuts[i].Data());
#ifdef __DC14__
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("WorkingPoint", m_isoCuts[i].Data()));
#else
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("ElectronWP", m_isoCuts[i].Data()));
#endif

      if (m_isoTools[iso]->initialize().isFailure()) {
        fatal("Failed to initialize IsolationSelectionTool with WP"+m_isoCuts[i]);
      }
    }
    
    // efficiency scale factor tool
    m_electronIDSF = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyIDCorrectionTool");
    m_electronRecoSF = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyRecoCorrectionTool");
    
    std::string file_ID   = config.getStr(m_name+".ScaleFactor.IDCorrectionFileName").Data();
    std::string file_Reco = config.getStr(m_name+".ScaleFactor.RecoCorrectionFileName").Data();
    std::vector< std::string > correctionFilesID;
    std::vector< std::string > correctionFilesReco;
    correctionFilesID.push_back(file_ID);
    correctionFilesReco.push_back(file_Reco);
    
    CP_CHECK(m_name, m_electronIDSF->setProperty("CorrectionFileNameList", correctionFilesID));
    CP_CHECK(m_name, m_electronIDSF->setProperty("ForceDataType", m_isMC ? 1 : 0)); // need 3 for AFII
    CP_CHECK(m_name, m_electronRecoSF->setProperty("CorrectionFileNameList", correctionFilesReco));
    CP_CHECK(m_name, m_electronRecoSF->setProperty("ForceDataType", m_isMC ? 1 : 0)); // need 3 for AFII


    if (m_electronIDSF->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyIDCorrectionTool");
    }

    if (m_electronRecoSF->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyRecoIDCorrectionTool");
    }

    // Read in configuration information
    m_containerName = config.getStr(m_name+".ContainerName", "ElectronCollection");

    m_etaCut      = config.getNum(m_name+".Selection.MaxAbsEta"  , 2.47);
    m_ptCut       = config.getNum(m_name+".Selection.PtPreCutGeV", 25.0)*GeV;

    m_crackReject = config.getBool(m_name+".Selection.ApplyCrackRejection", true);
    m_barrelMax   = config.getNum(m_name+".Selection.BarrelMaxAbsEta"  , 1.37);
    m_endcapMin   = config.getNum(m_name+".Selection.EndcapMinAbsEta"  , 1.52);
    
    m_applyIPCuts   = config.getBool(m_name+".Selection.ApplyIPCuts", false);
    m_d0BySigd0Cut  = config.getNum(m_name+".Selection.d0BySigd0Max", 6.5);
    m_z0Cut         = config.getNum(m_name+".Selection.z0Max", 10.0);


    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::ElectronContainer ElectronHandler::getCorrectedContainer()
  {
    // get the event info
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      fatal("Cannot access EventInfo");
    }
      
    bool calib = false;
    xAOD::ElectronContainer shallowContainer = getShallowContainer(calib);
    if (calib) return shallowContainer;

    //PID on not calibrated objects, SF on calibrated objects, isolation on calibrated objects
    for (auto electron: shallowContainer) {
      calibrateAndSmearElectron(electron, eventInfo, m_electronCalibTool);
      decoratePID(*electron);
      decorateIso(*electron);
      decorateIPCut(*electron);
      applyScaleFactor(electron, eventInfo);
    }
    
    // sort the electrons
    shallowContainer.sort(comparePt);
    
    return shallowContainer;
  }

  //______________________________________________________________________________
  xAOD::ElectronContainer ElectronHandler::applySelection(xAOD::ElectronContainer &container)
  {
    xAOD::ElectronContainer selected(SG::VIEW_ELEMENTS);
    for (auto electron: container) {
      // pT and eta cuts
      if (!passPtEtaCuts(electron)) continue;

      // d0/z0 cuts
      if (m_applyIPCuts && !passIPCuts(electron)) continue;

      // PID LH identification
      if (m_doPidCut && !passPIDCut(electron)) continue;

      // isolation cuts
      if (m_doIsoCut && !passIsoCut(electron)) continue;

      selected.push_back(electron);
    }
    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode ElectronHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool isAffected = false;
    for (auto var: sys) {
      if (m_electronCalibTool->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
      if (m_electronIDSF->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
      if (m_electronRecoSF->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
    }

    if (isAffected) {
      CP_CHECK(m_name, m_electronCalibTool->applySystematicVariation(sys));
      CP_CHECK(m_name, m_electronIDSF     ->applySystematicVariation(sys));
      CP_CHECK(m_name, m_electronRecoSF   ->applySystematicVariation(sys));
      m_sysName = sys.name() == "" ? "" : "_"+sys.name();
    } else {
      CP_CHECK(m_name, m_electronCalibTool->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_electronIDSF     ->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_electronRecoSF   ->applySystematicVariation(CP::SystematicSet()));
      m_sysName = "";
    }
    
    return CP::SystematicCode::Ok;
  }

  //______________________________________________________________________________
  void ElectronHandler::calibrateAndSmearElectron(xAOD::Electron *ele,
                                                  const xAOD::EventInfo *evtInfo,
                                                  CP::EgammaCalibrationAndSmearingTool *electronCalibTool) {
    float cl_eta = 10.;
    const xAOD::CaloCluster* cluster = ele->caloCluster();
    if(cluster) cl_eta = cluster->eta();
    
    // Apply smearing 
    if(std::abs(cl_eta) < 2.47 && ele->pt() >= 20000.){
      electronCalibTool->setRandomSeed(evtInfo->eventNumber()*100+ele->index());
      
      // Calibrate the electron
      double E_before = ele->e();
      CP::CorrectionCode cc = electronCalibTool->applyCorrection(*ele);
      if (cc==CP::CorrectionCode::Error)
	Error("calibratedElectron()","Error calibrating current electron");
      if (cc==CP::CorrectionCode::OutOfValidityRange)
	Warning("calibratedElectron()","Current electron has no valid calibration due to out-of-range");
      
      Ecalib_ratio(*ele) = ele->e()/E_before;
#ifndef __DC14__
      Ereso(*ele) = electronCalibTool->getResolution(*ele);
#endif
    }
  }

  //______________________________________________________________________________
  void ElectronHandler::applyScaleFactor(xAOD::Electron *ele, const xAOD::EventInfo *evtInfo)
  {
    double _effIDSF = 1.0, _effRecoSF = 1.0;
    
    float cl_eta = 10.;
    const xAOD::CaloCluster* cluster = ele->caloCluster();
    if(cluster)  cl_eta = cluster->eta();
    
    if (m_isMC && std::abs(cl_eta) < 2.47 && ele->pt() >= 15000.){
      if(m_electronIDSF->getEfficiencyScaleFactor     (*ele, _effIDSF   ) == CP::CorrectionCode::Error)
        fatal("ElectronEfficiencyIDCorrection returned CP::CorrectionCode::Error");
      if(m_electronRecoSF->getEfficiencyScaleFactor   (*ele, _effRecoSF   ) == CP::CorrectionCode::Error)
	fatal("ElectronEfficiencyRecoCorrection returned CP::CorrectionCode::Error");
    }

    effIDSF(*ele) = _effIDSF;
    effRecoSF(*ele) = _effRecoSF;
    scaleFactor(*ele) = _effIDSF*_effRecoSF;
  }

  //______________________________________________________________________________
  bool ElectronHandler::passPtEtaCuts(const xAOD::Electron *ele) {

    double abs_eta = fabs(ele->eta());
    if (abs_eta > m_etaCut) return false;
    if (m_crackReject && (abs_eta > m_barrelMax && abs_eta < m_endcapMin)) return false;
    
    if (ele->pt() < m_ptCut) return false;
  
    return true;
  }

  //______________________________________________________________________________
  void ElectronHandler::decorateIso(xAOD::Electron &ele)
  {
    for (auto dec: m_isoDecorators) {
      if (m_isoTools[dec.first]->accept(ele))
        (*dec.second)(ele) = true;
      else
        (*dec.second)(ele) = false;
    }
  }

  //______________________________________________________________________________
  bool ElectronHandler::passIsoCut(const xAOD::Electron *ele, HG::Iso::IsolationType iso)
  {
    /// applies isolation cut specified in config file
    if (iso == HG::Iso::Undefined) {
      if (!m_isoDecorators[m_defaultIso]->isAvailable(*ele))
        return true;
      return (*m_isoDecorators[m_defaultIso])(*ele);
    }
    
    if (m_isoTools.find(iso) != m_isoTools.end()) {
      if (!m_isoDecorators[iso]->isAvailable(*ele))
        return true;
      return (*m_isoDecorators[iso])(*ele);
    }

    fatal("Isolation cut requested that wasn't specified in config file. Exiting.");
    return false;
  }

  //______________________________________________________________________________
  void ElectronHandler::decoratePID(xAOD::Electron &ele)
  {
    for (auto dec: m_electronSelDecorators) {
      if (m_electronSelectors[dec.first]->accept(&ele))
        (*dec.second)(ele) = true;
      else
        (*dec.second)(ele) = false;
    }
  }
  
  //______________________________________________________________________________
  bool ElectronHandler::passPIDCut(const xAOD::Electron *ele, TString pid)
  {

    // applies PID cut specified in config file
    if (pid == "Default") {
      if (!m_electronSelDecorators[m_defaultPid]->isAvailable(*ele))
        return true;
      return (*m_electronSelDecorators[m_defaultPid])(*ele);
    }
    
    if (m_electronSelDecorators.find(pid) != m_electronSelDecorators.end()) {
      if (!m_electronSelDecorators[pid]->isAvailable(*ele))
        return true;
      return (*m_electronSelDecorators[pid])(*ele);
    }

    fatal(TString::Format("%s: PID cut requested that wasn't specified in config file. Exiting.",m_name.Data()));
    return false;
  }

  //______________________________________________________________________________
  bool ElectronHandler::passIPCuts(const xAOD::Electron *ele) 
  {
    if(passIPCut.isAvailable(*ele) && !passIPCut(*ele)) return false;
    return true;
  }

  //______________________________________________________________________________
  void ElectronHandler::decorateIPCut(xAOD::Electron &ele)
  {
    passIPCut(ele) = false;
    const xAOD::TrackParticle *track = ele.trackParticle();
    if (track==nullptr) return;
    double d0sign = fabs(track->d0()) / std::sqrt(track->definingParametersCovMatrixVec()[0]);
    if (d0sign > m_d0BySigd0Cut) return;
    
    const xAOD::VertexContainer* vertexCont =0;
    if (m_event->retrieve(vertexCont,"PrimaryVertices").isFailure()) return;
    
    const xAOD::Vertex* pvx = xAOD::PVHelpers::getHardestVertex(vertexCont);
    if (pvx == nullptr) return;
    
    if (fabs(track->z0() + track->vz() - pvx->z()) > m_z0Cut) return;
    passIPCut(ele) = true;
  }
  
  //______________________________________________________________________________
  HG::Iso::IsolationType ElectronHandler::getIsoType(TString isoName) {
    if      (isoName == "Tight") return HG::Iso::Tight;
    else if (isoName == "Loose") return HG::Iso::Loose;
    else if (isoName == "Gradient") return HG::Iso::Gradient;
    else if (isoName == "Cone40CaloOnly") return HG::Iso::Cone40CaloOnly;
    else if (isoName == "Cone40") return HG::Iso::Cone40;
    else if (isoName == "Cone20") return HG::Iso::Cone20;
    else if (isoName == "UserDefined") return HG::Iso::UserDefined;
    else fatal("Isolation "+isoName+" read from: "+
               m_name+".Selection.IsolationCriteria is not Tight, Gradient, Loose, or UserDefined. Exiting.");
    return HG::Iso::Undefined;
  }
    
}

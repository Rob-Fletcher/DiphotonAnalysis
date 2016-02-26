#include "HGamAnalysisFramework/MuonHandler.h"

#ifdef __DC14__
#include "ElectronIsolationSelection/IsolationSelectionTool.h"
#else
#include "IsolationSelection/IsolationSelectionTool.h"
#endif

namespace HG {
  //______________________________________________________________________________
  SG::AuxElement::Decorator<float> MuonHandler::effSF("SF_eff");
  SG::AuxElement::Decorator<float> MuonHandler::scaleFactor("scaleFactor");
  SG::AuxElement::Decorator<char> MuonHandler::isAccepted("isAccepted");
  SG::AuxElement::Decorator<char> MuonHandler::passIPCut("passIPCut");

  //______________________________________________________________________________
  MuonHandler::MuonHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store) : HgammaHandler(name, event, store)
  { }

  //______________________________________________________________________________
  MuonHandler::~MuonHandler()
  {
      SafeDelete(m_muonEffScaleFactors);
      SafeDelete(m_muonTrigScaleFactors);
      SafeDelete(m_muonCalibTool);
      SafeDelete(m_muonSelectTool);

      for (auto iso: m_isoTools) SafeDelete(iso.second);
      m_isoTools.clear();

      for (auto dec: m_isoDecorators) SafeDelete(dec.second);
      m_isoDecorators.clear();
  }

  //______________________________________________________________________________
  EL::StatusCode MuonHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // Selecting muons
    m_muonSelectTool = new CP::MuonSelectionTool("MuonSelectionTool");
    m_pidCut = config.getStr(m_name+".Selection.PID", "Medium");
    if (m_pidCut == "Tight")
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::Tight)));
    }
    else if (m_pidCut == "Medium")
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::Medium)));
    }
    else if (m_pidCut == "Loose") 
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::Loose)));
    }
    else if (m_pidCut == "VeryLoose") 
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::VeryLoose)));
    }
    else
      fatal(TString::Format("Value: %s for key: %s.Selection.PID not vlid MuQuality. Exiting.", m_name.Data(), m_pidCut.Data()).Data());

    m_MaxEta = config.getNum(m_name+".Selection.MaxEta");
    CP_CHECK(m_name, m_muonSelectTool->setProperty("MaxEta", m_MaxEta));

    if (!m_muonSelectTool->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonSelection Tool. Exiting.");
    }
   

    // Efficiency scale factors
    m_muonEffScaleFactors = new CP::MuonEfficiencyScaleFactors("MuonEfficiencyScaleFactors");
    for (TString prop : {"WorkingPoint", "DataPeriod"}) 
    {
      CP_CHECK(m_name, m_muonEffScaleFactors->setProperty(prop.Data(), config.getStr(m_name+".Efficiency."+prop).Data()));
    }

    if (! m_muonEffScaleFactors->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonEfficiencyScaleFactors Tool. Exiting." );
    }


    // Trigger efficiency scale factors
    m_muonTrigScaleFactors = new CP::MuonTriggerScaleFactors("MuonTriggerScaleFactors");
#ifndef __DC14__
    CP_CHECK(m_name, m_muonTrigScaleFactors->setProperty("MuonQuality", m_pidCut.Data()));
#else
    CP_CHECK(m_name, m_muonTrigScaleFactors->setProperty("runNumber", config.getInt(m_name+".Efficiency.RunNumber", 900000) ));
#endif

    if (! m_muonTrigScaleFactors->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonTriggerScaleFactors Tool. Exiting." );
    }


    // Calibrate and smear
    m_muonCalibTool = new CP::MuonCalibrationAndSmearingTool( "MuonCalibrationAndSmearingTool" );
    for (TString prop : {"Year", "Algo", "SmearingType"}) 
    {
      CP_CHECK(m_name, m_muonCalibTool->setProperty(prop.Data(), config.getStr(m_name+".Calibration."+prop).Data()));
    }
    
    if (! m_muonCalibTool->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonCalibrationAndSmearingTool Tool. Exiting." );
    }

          
    // Isolation tools
    m_doIsoCut = config.getBool(m_name+".Selection.ApplyIsoCut", true);
    m_isoCuts  = config.getStrV(m_name+".Selection.IsoCriteria");
    if (m_doIsoCut && m_isoCuts.size() < 1)
      fatal("Isolation cut requested, but no working point supplied. Exiting!");
    for (size_t i = 0; i < m_isoCuts.size(); ++i) {
      HG::Iso::IsolationType iso = HG::Iso::Undefined;
      TString name = m_isoCuts[i];
      if (name == "Tight") 
      {
        iso = HG::Iso::Tight;
	m_isoDecorators[iso] = new SG::AuxElement::Accessor<char>("isIsoTight");
      } 
      else if (name == "Gradient") 
      {
        iso = HG::Iso::Gradient;
	m_isoDecorators[iso] = new SG::AuxElement::Accessor<char>("isIsoGradient");	
      }
      else if (name == "Loose") 
      {
        iso = HG::Iso::Loose;
	m_isoDecorators[iso] = new SG::AuxElement::Accessor<char>("isIsoLoose");
      }
      else if (name == "UserDefined") 
      {
        iso = HG::Iso::UserDefined;
	m_isoDecorators[iso] = new SG::AuxElement::Accessor<char>("isIsoUserDefined");
      } 
      else 
      {
        fatal(TString::Format("Value: %s for key: %s.Selection.IsolationCriteria not Tight, Gradient, Loose, or UserDefined. Exiting.", m_name.Data(), m_isoCuts[i].Data()).Data());
      }
      
      if (i == 0) m_defaultIso = iso;
      
      m_isoTools[iso] = new CP::IsolationSelectionTool(name.Data());
#ifdef __DC14__
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("WorkingPoint", name.Data()));
#else
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("MuonWP", name.Data()));
#endif
      
      if (m_isoTools[iso]->initialize().isFailure()) {
        fatal("Failed to initialize IsolationSelectionTool with WP: "+name);
      }
    }
    
    
    // Read in remaining configuration information
    m_containerName = config.getStr(m_name+".ContainerName", "Muons");
    m_ApplyPtCut    = config.getBool(m_name+".Selection.ApplyPtCut", true);
    m_PtCut         = config.getNum(m_name+".Selection.PtCutGeV",   10.0);

    m_ApplyIPCuts   = config.getBool(m_name+".Selection.ApplyIPCuts", false);
    m_d0Cut         = config.getNum(m_name+".Selection.d0Max",   1.0);
    m_z0Cut         = config.getNum(m_name+".Selection.z0Max",  10.0);

    
    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::MuonContainer MuonHandler::getCorrectedContainer()
  {
    // get the event info
    // const xAOD::EventInfo *eventInfo = 0;
    //if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) fatal("Cannot access EventInfo");

    bool calib = false;
    xAOD::MuonContainer shallowContainer = getShallowContainer(calib);
    if (calib) return shallowContainer;

    for (auto muon: shallowContainer)
    {
      if (fabs(muon->eta()) <= m_MaxEta) calibrateAndSmearMuon(muon, m_muonCalibTool);
      applyScaleFactor(muon);
      
      // Add selection decorations
      isAccepted(*muon) = m_muonSelectTool->accept(muon);
      decorateIPCut(muon);
      decorateIso(*muon);
    }
    

    // sort the muons
    shallowContainer.sort(comparePt);
    
    return shallowContainer;
  }

  //______________________________________________________________________________
  xAOD::MuonContainer MuonHandler::applySelection(xAOD::MuonContainer &container)
  {
    xAOD::MuonContainer selected(SG::VIEW_ELEMENTS);
    for (auto muon: container)
    {
      // Apply selection cuts
      if (!passSelection(muon)) continue;
      
      // require Isolation
      if (m_doIsoCut && !passIsoCut(muon)) continue;
      
      selected.push_back(muon);
    }
    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode MuonHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool canBeShifted = false;

    // Try shifting calibration tool
    if (m_muonCalibTool->applySystematicVariation(sys) == CP::SystematicCode::Ok)
      canBeShifted = true;
    else
      m_muonCalibTool->applySystematicVariation(CP::SystematicSet()).ignore();
    
    // If either was shifted, reflect that in container name
    if (canBeShifted)
      m_sysName = sys.name() == "" ? "" : TString::Format("_%s", sys.name().c_str());
    else
      m_sysName = "";
    
    return CP::SystematicCode::Ok;

  }

  
  //______________________________________________________________________________
  void MuonHandler::applyScaleFactor(xAOD::Muon *muon)
  {
    float EfficiencyScaleFactor = 1.0;
    float TriggerScaleFactor    = 1.0;
    
    CP::CorrectionCode cc = m_muonEffScaleFactors->applyEfficiencyScaleFactor( *muon );
    if (cc==CP::CorrectionCode::Error)
      Error("applyScaleFactor()","Error applying efficiency scale factor to current muon");
    
    if(m_muonEffScaleFactors->getEfficiencyScaleFactor(*muon, EfficiencyScaleFactor) == CP::CorrectionCode::Error)
      fatal("applyScaleFactor():  Error retrieving applied scale factor from muon");
    

    effSF(*muon) = EfficiencyScaleFactor;
    scaleFactor(*muon) = EfficiencyScaleFactor;
    //trigSF(*muon) = TriggerScaleFactor;
    //scaleFactor(*muon) = EfficiencyScaleFactor * TriggerScaleFactor;

  }
  
  
  //______________________________________________________________________________
  void MuonHandler::calibrateAndSmearMuon(xAOD::Muon *muon, CP::MuonCalibrationAndSmearingTool *muonCalibTool)
  {
    // applies calibration to a muon

    // Apply smearing (?)     <<<<<< Does not apply to muons??
    //photonCalibTool->setRandomSeed(evtInfo->eventNumber()*100+gam->index());
    
    // Calibrate the muon
    double E_before = muon->e();
    CP::CorrectionCode cc = muonCalibTool->applyCorrection( *muon );
    if (cc==CP::CorrectionCode::Error)
      Error("calibratedAndSmearMuon()","Error calibrating current muon");
    if (cc==CP::CorrectionCode::OutOfValidityRange)
      Warning("calibratedAndSmearMuon()","Current muon has no valid calibration due to out-of-range");
    
    // decorate the muon with the calibration factor
    muon->auxdata< float >( "Ecalib_ratio" ) = muon->e()/E_before;
  }
    
  
  //______________________________________________________________________________
  //bool MuonHandler::passPtEtaCuts(const xAOD::Muon *muon) {
    /// applies kinematic preselection cut
    // eta cuts
  //  if (fabs(muon->eta()) > m_etaCut) return false;
  
    // pt cuts
  //  if (muon->pt() < m_ptCut) return false;
  //  return true;
  //}
  
  //______________________________________________________________________________
  bool MuonHandler::passSelection(const xAOD::Muon *muon)
  {
    // Muon selector performs the following cuts:
    //if (!m_muonSelectTool->getQuality(muon) <= xAOD::Muon::Medium) return false;
    //if (!m_muonSelectTool->passedIDCuts(muon)) return false;
    //if (!m_muonSelectTool->passedHighPtCuts(muon)) return false;
    if (isAccepted.isAvailable(*muon) && !isAccepted(*muon)) return false;
    	  
    if (m_ApplyPtCut)
    {
      if (muon->pt()/HG::GeV < m_PtCut) return false;
    }

    if (m_ApplyIPCuts && passIPCut.isAvailable(*muon) && !passIPCut(*muon)) return false;

    return true;
  }
  
  //______________________________________________________________________________
  void MuonHandler::decorateIPCut(const xAOD::Muon *muon)
  {
    passIPCut(*muon) = false;

    //std::cout << "checking z0/d0... type = " << muon->muonType() << std::endl;
    //Info("checking z0/d0.. type=" +  muon->muonType());
    if (!(muon->muonType() == xAOD::Muon::MuonStandAlone))
    {
      if (fabs(muon->primaryTrackParticle()->d0()) > m_d0Cut) return;

      const xAOD::VertexContainer* vertexCont =0;
      if (m_event->retrieve(vertexCont,"PrimaryVertices").isFailure()) return;

      const xAOD::Vertex* pvx = xAOD::PVHelpers::getHardestVertex(vertexCont);
      if (pvx == nullptr) return;

      if (fabs(muon->primaryTrackParticle()->z0() + muon->primaryTrackParticle()->vz() - pvx->z()) > m_z0Cut) return;

      // Still here? It passes!
      passIPCut(*muon) = true;
    }
  }


  //______________________________________________________________________________
  bool MuonHandler::passIsoCut(const xAOD::Muon *muon, HG::Iso::IsolationType iso)
  {
    /// applies PID cut specified in config file
    if (iso == HG::Iso::Undefined) 
    {
      if (!m_isoDecorators[m_defaultIso]->isAvailable(*muon)) return true;
      return m_isoTools[m_defaultIso]->accept(*muon);
    }
        
    if (m_isoTools.find(iso) != m_isoTools.end()) 
    {
      if (!m_isoDecorators[iso]->isAvailable(*muon))
        return true;
      return (*m_isoDecorators[iso])(*muon);
    }

    fatal("Isolation cut requested that wasn't specified in config file. Exiting.");
    return false;
  }


  //______________________________________________________________________________
  void MuonHandler::decorateIso(xAOD::Muon &muon)
  {
    for (auto dec: m_isoDecorators) {
      if (m_isoTools[dec.first]->accept(muon))
        (*dec.second)(muon) = true;
      else
        (*dec.second)(muon) = false;
    }
  }
  
  
  //______________________________________________________________________________
  void MuonHandler::printMuon(const xAOD::Muon *muon, TString comment)
  {
    // prints details about the photon
    printf("Muon %2zu  %s\n", muon->index(), comment.Data());
    
    // print the 4-vector
    printf("   (pT,eta,phi,m) = (%5.1f GeV,%6.3f,%6.3f,%4.4f GeV)\n", muon->pt()/GeV, muon->eta(), muon->phi(), muon->m()/GeV);
    
    // print some more information
    TString str;
    if (muon->isAvailable<float>("Ecalib_ratio"))
      str+=Form("   calibFactor = %.3f", muon->auxdata<float>("Ecalib_ratio"));
    
    if (muon->isAvailable<float>("EfficiencyScaleFactor"))
      str+=Form("   scalefactor = %5.5f",muon->auxdata<float>("EfficiencyScaleFactor"));
    if (str.Sizeof()) printf("%s\n",str.Data());
  }

}

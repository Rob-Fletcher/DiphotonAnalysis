#include "HGamAnalysisFramework/JetHandler.h"

#include "HGamAnalysisFramework/HgammaUtils.h"

#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"

#include "PhotonVertexSelection/PhotonVertexHelpers.h"
#include "JetUncertainties/JetUncertaintiesTool.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"
#include "JetMomentTools/JetOriginCorrectionTool.h"

#include "JetResolution/JERTool.h"
#include "JetResolution/JERSmearingTool.h"

#include "JetCalibTools/JetCalibrationTool.h"
#include "JetRec/PseudoJetGetter.h"

namespace HG {

  SG::AuxElement::Accessor<std::vector<float> > JetHandler::JVF("JVF");
  SG::AuxElement::Accessor<float>  JetHandler::Jvt("Jvt");
  SG::AuxElement::Accessor<float>  JetHandler::DetectorEta("DetectorEta");
  SG::AuxElement::Decorator<float> JetHandler::Jvf("Jvf");
  SG::AuxElement::Decorator<char>  JetHandler::isClean("isClean");
  SG::AuxElement::Decorator<float> JetHandler::scaleFactor("scaleFactor");
  
  //______________________________________________________________________________
  JetHandler::JetHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
  , m_jetCalibTool(nullptr)
  , m_jetCleaning(nullptr)
  , m_jesProvider(nullptr)
  , m_jerTool(nullptr)
  , m_jerSmear(nullptr)
  , m_trackTool(nullptr)
  , m_jvtLikelihood(nullptr)
  , m_jvtTool(nullptr)
  , m_jetOriginTool(nullptr)
  { }

  //______________________________________________________________________________
  JetHandler::~JetHandler()
  {
    SafeDelete(m_jetCalibTool);
    SafeDelete(m_jetCleaning);
    SafeDelete(m_jesProvider);
    SafeDelete(m_jerTool);
    SafeDelete(m_jerSmear);
    SafeDelete(m_trackTool);
    SafeDelete(m_jvtLikelihood);
    SafeDelete(m_jvtTool);
    SafeDelete(m_jetOriginTool);

    for (auto name: m_bTagNames)
      for (auto tool: m_bTagEffTools[name])
        SafeDelete(tool);
  }

  //______________________________________________________________________________
  EL::StatusCode JetHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // General configs
    m_isAFII = config.getBool("IsAFII", false);

    // Read in configuration information
    m_containerName = config.getStr(m_name+".ContainerName");
    m_truthName     = config.getStr(m_name+".TruthContainerName", "AntiKt4TruthJets");

    m_rapidityCut 	= config.getNum(m_name+".Selection.MaxAbsRapidity" , 4.4);
    m_ptCut       	= config.getNum(m_name+".Selection.PtPreCutGeV"    , 25.0)*GeV;
    m_jvf               = config.getNum(m_name+".Selection.JVF"            , 0.25);
    m_jvt               = config.getNum(m_name+".Selection.JVT"            , 0.941);

    // B-tagging
    m_enableBTagging    = config.getBool(m_name+".EnableBTagging"          , false);
    m_bTagRapidityCut   = config.getNum(m_name+".BTagging.MaxAbsRapidity"  , 2.5);
    m_bTagNames     = config.getStrV(m_name+".BTagging.TaggerNames");
    TString bTagCDI = config.getStr(m_name+".BTagging.ScaleFactorFileName", "");

    for (auto name: m_bTagNames) {
      // Get efficiency and corresponding operating point for TaggerName
      m_bTagEffs[name] = config.getStrV(m_name+"."+name+".Efficiencies");
      m_bTagOPs [name] = config.getStrV(m_name+"."+name+".OperatingPoints");

      for (size_t i = 0; i < m_bTagEffs[name].size(); ++i) {
        TString op = m_bTagOPs[name].at(i);

        // Configure efficiency tool
        BTaggingEfficiencyTool *tempTool = new BTaggingEfficiencyTool("BTaggingEfficiencyTool");
        CP_CHECK(m_name, tempTool->setProperty("TaggerName"         , name           .Data()));
        CP_CHECK(m_name, tempTool->setProperty("OperatingPoint"     , op             .Data()));
        CP_CHECK(m_name, tempTool->setProperty("JetAuthor"          , m_containerName.Data()));
        CP_CHECK(m_name, tempTool->setProperty("ScaleFactorFileName", bTagCDI        .Data()));

        if (tempTool->initialize().isFailure())
          fatal("Couldn't initialize BTaggingEfficeincyTool, exiting.");

        m_bTagEffTools[name].push_back(tempTool);

        // Convert operating point to cut value
        op.ReplaceAll("_", ".");
        m_bTagCuts[name].push_back(op.Atof());
      }
    }

    // calibration tool
    TString jetAlgo   = m_containerName;
    jetAlgo.ReplaceAll("Jets", "");
    TString jetconfig = config.getStr(m_name+".Calibration.ConfigFile");
    if (m_isAFII)
      jetconfig = config.getStr(m_name+".CalibrationAFII.ConfigFile");
    TString calibSeq  = config.getStr(m_name+".Calibration.CalibSeq");
    if (m_isData) calibSeq += "_Insitu";

    m_correctVertex   = config.getBool(m_name+".Calibration.CorrectVertex");

    m_jetCalibTool = new JetCalibrationTool("JetCalibTool", jetAlgo, jetconfig, calibSeq, m_isData);

    if (m_jetCalibTool->initializeTool("JetCalibTool").isFailure()) {
      fatal("Failed to initialize JetCalibrationTool");
    }

    // smearing tool
    if (isMC()) {
      m_jerTool = new JERTool("JERTool");
      CP_CHECK(m_name, m_jerTool->setProperty("PlotFileName", config.getStr(m_name+".Resolution.PlotFileName").Data()));
      CP_CHECK(m_name, m_jerTool->setProperty("CollectionName", m_containerName.Data()));

      ToolHandle<IJERTool> jerHandle(m_jerTool->name());

      m_jerSmear = new JERSmearingTool("JERSmearingTool");
      CP_CHECK(m_name, m_jerSmear->setProperty("JERTool", jerHandle));
      CP_CHECK(m_name, m_jerSmear->setProperty("ApplyNominalSmearing", config.getBool(m_name+".Resolution.ApplyNominalSmearing")));
      CP_CHECK(m_name, m_jerSmear->setProperty("isMC", isMC()));
      CP_CHECK(m_name, m_jerSmear->setProperty("SystematicMode", config.getStr(m_name+".Resolution.SystematicMode").Data()));

      if (m_jerTool->initialize().isFailure())
        fatal("JERTool failed to initialize, exiting!");

      if (m_jerSmear->initialize().isFailure())
        fatal("JERSmearingTool failed to initialize, exiting!");
    }

    // cleaning tool
    m_jetCleaning = new JetCleaningTool("JetCleaning");
    CP_CHECK(m_name, m_jetCleaning->setProperty("CutLevel", config.getStr(m_name+".Selection.CutLevel").Data()));
    CP_CHECK(m_name, m_jetCleaning->setProperty("DoUgly", config.getBool(m_name+".Selection.DoUgly")));

    if (m_jetCleaning->initialize().isFailure()) {
      fatal("Failed to initialize JetCleaningTool");
    }

    // JES uncertainty provider
    m_jesProvider = new JetUncertaintiesTool("JESProvider");
    CP_CHECK(m_name, m_jesProvider->setProperty("JetDefinition", jetAlgo.Data()));
    CP_CHECK(m_name, m_jesProvider->setProperty("MCType"    , config.getStr(m_name+".Uncertainty.MCType"    ).Data()));
    CP_CHECK(m_name, m_jesProvider->setProperty("ConfigFile", config.getStr(m_name+".Uncertainty.ConfigFile").Data()));

    if (m_jesProvider->initialize().isFailure()) {
      fatal("Failed to initialize uncertainty tool");
    }

    // Track selection tool
    m_trackTool = new InDet::InDetTrackSelectionTool("TrackSelection");
    m_trackTool->setCutLevel(InDet::CutLevel::Loose);
    if (m_trackTool->initialize().isFailure())
      fatal("Failed to initialize InDetTrackSelectionTool");

    // JVT likelihood histogram
    TString jvtFile = "JetMomentTools/JVTlikelihood_20140805.root";
    TString jvtName = "JVTRootCore_kNN100trim_pt20to50_Likelihood";
    m_jvtLikelihood = (TH2F*)HG::getHistogramFromFile(jvtName, jvtFile);

    if (m_jvtLikelihood == nullptr)
      fatal("Failed to get JVT likelihood file, exiting.");

    // JVT re-scaling tool
    m_jvtTool = new JetVertexTaggerTool("JVTag");
    CP_CHECK(m_name, m_jvtTool->setProperty("JVTFileName", jvtFile.Data()));
    if (m_jvtTool->initialize().isFailure())
      fatal("Failed to initialize JetVertexTaggerTool, exiting.");

    // Jet origin correction tool
    m_jetOriginTool = new JetOriginCorrectionTool("JOT");

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::JetContainer JetHandler::getCorrectedContainer()
  {
    // Get Shallow copy from TEvent/TStore
    bool calib = false;
    xAOD::JetContainer shallowContainer = getShallowContainer(calib);
    if (calib) return shallowContainer;

    // Get the photon pointed vertex, if available
    if (m_correctVertex && m_store->contains<ConstDataVector<xAOD::VertexContainer> >("HGamVertices")) {

      ConstDataVector<xAOD::VertexContainer> *vertices = nullptr;
      if (m_store->retrieve(vertices, "HGamVertices").isFailure())
        Warning("JetHandler::recalculateJVT()", "TStore contains HGamVertices, but couldn't retrieve it.");

      const xAOD::Vertex *vertex = nullptr;
      if (vertices && vertices->size() > 0)
        vertex = (*vertices)[0];

      if (vertex && vertex->vertexType() != xAOD::VxType::PriVtx) {
        const xAOD::EventInfo *eventInfo = nullptr;
        if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
          HG::fatal("Cannot access EventInfo");

        static SG::AuxElement::Decorator<int> PVIndex("PVIndex");

        PVIndex(*eventInfo) = vertex->index();

        m_jetOriginTool->modify(shallowContainer);
      }

    }

    // calibrate and decorate jets
    for (auto jet : shallowContainer) {
      calibrateJet(jet);

      // Decorations necessary for MxAOD
      decorateJVF(jet);
      recalculateJVT(*jet);
      isClean(*jet) = m_jetCleaning->accept(*jet);
      if (m_enableBTagging)
        decorateBJet(jet);
      DetectorEta(*jet) = jet->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum").eta();
    }

    // sort the Jets
    shallowContainer.sort(comparePt);

    return shallowContainer;
  }

  //______________________________________________________________________________
  xAOD::JetContainer JetHandler::applySelection(xAOD::JetContainer &container)
  {
    xAOD::JetContainer selected(SG::VIEW_ELEMENTS);
    for (auto jet: container) {
      // apply pT and eta selection cuts
      if (!passPtEtaCuts(jet)) continue;

      // Apply cleaning cuts
      // If not available, it's an MxAOD where the cut is already applied
      if (isClean.isAvailable(*jet) && !isClean(*jet)) continue;

      // JVF cuts
      if (!passJVFCut(jet)) continue;
  
      // JVT cuts
      if (!passJVTCut(jet, m_jvt)) continue;

      scaleFactor(*jet) = 1.0;
  
      selected.push_back(jet);
    }
    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode JetHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool isAffected = false;
    for (auto var: sys) {
      if (m_jesProvider->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
      if (isMC() && m_jerSmear->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
      // if (isMC()) {
      //   for (auto name: m_bTagNames) {
      //     for (auto tool: m_bTagEffTools[name]) {
      //       if (tool->isAffectedBySystematic(var)) {
      //         isAffected = true;
      //         break;
      //       }
      //     }
      //   }
      // }

    }


    if (isAffected) {
      // This should mean that the jets can be shifted by this systematic
      CP_CHECK(m_name, m_jesProvider->applySystematicVariation(sys))
      if (isMC())
        CP_CHECK(m_name, m_jerSmear->applySystematicVariation(sys))

      // if (isMC()) {
      //   for (auto name: m_bTagNames) {
      //     for (auto tool: m_bTagEffTools[name]) {
      //       CP_CHECK(m_name, tool->applySystematicVariation(sys));
      //     }
      //   }
      // }

      m_sysName = sys.name() == "" ? "" : "_"+sys.name();
    } else {
      // Jets are not affected
      CP_CHECK(m_name, m_jesProvider->applySystematicVariation(CP::SystematicSet()));

      if (isMC())
        CP_CHECK(m_name, m_jerSmear->applySystematicVariation(CP::SystematicSet()));

      // if (isMC()) {
      //   for (auto name: m_bTagNames) {
      //     for (auto tool: m_bTagEffTools[name]) {
      //       CP_CHECK(m_name, tool->applySystematicVariation(CP::SystematicSet()));
      //     }
      //   }
      // }

      m_sysName = "";
    }

    return CP::SystematicCode::Ok;
  }

  //______________________________________________________________________________
  void JetHandler::calibrateJet(xAOD::Jet *jet)
  {
    // applies calibration to a jet
    float E_before = jet->e();
    m_jetCalibTool->applyCorrection(*jet).ignore();

    if (isMC())
      m_jerSmear->applyCorrection(*jet).ignore();

    // Uncertainty shifting
    if (jet->pt() > 20.0*HG::GeV && fabs(jet->rapidity()) < 4.4)
      m_jesProvider->applyCorrection(*jet).ignore();
      
    // decorate the photon with the calibration factor
    jet->auxdata< float >( "Ecalib_ratio" ) = jet->e()/E_before;
  }
    
  //______________________________________________________________________________
  bool JetHandler::passPtEtaCuts(const xAOD::Jet *jet)
  {
    /// applies kinematic preselection cuts: not-in-crack + pT cut

    // eta cuts
    if (fabs(jet->rapidity()) > m_rapidityCut) return false;
  
    // pt cuts
    if (jet->pt() < m_ptCut) return false;

    return true;
  }
  
  //______________________________________________________________________________
  void JetHandler::rescaleJVT(xAOD::Jet &jet)
  {
    m_jvtTool->updateJvt(jet);
  }

  //______________________________________________________________________________
  void JetHandler::recalculateJVT(xAOD::Jet &jet)
  {
    // Get vertex for calculation
    if (!m_store->contains<ConstDataVector<xAOD::VertexContainer> >("HGamVertices")) {
      rescaleJVT(jet);
      return;
    }

    ConstDataVector<xAOD::VertexContainer> *vertices = nullptr;
    if (m_store->retrieve(vertices, "HGamVertices").isFailure()) {
      Warning("JetHandler::recalculateJVT()", "TStore contains HGamVertices, but couldn't retrieve it.");
      rescaleJVT(jet);
      return;
    }

    const xAOD::Vertex *vertex = nullptr;
    if (vertices && vertices->size() > 0)
      vertex = (*vertices)[0];

    if (vertex == nullptr) {
      rescaleJVT(jet); 
      return;
    }

    if (vertex->vertexType() == xAOD::VxType::PriVtx) {
      rescaleJVT(jet); 
      return;
    }

    // // This gets the hardest vertex instead, for cross checks
    // const xAOD::VertexContainer *vertices = nullptr;
    // if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
    //   fatal("JetHandler::recalculateJVT() : Couldn't get PrimaryVertices, exiting.");

    // const xAOD::Vertex *vertex = xAOD::PVHelpers::getHardestVertex(vertices);
    // if (vertex == nullptr)
    //   fatal("Couldn't find hardest vertex? Exiting.");

    // Get track container
    const xAOD::TrackParticleContainer *tracks = nullptr;
    if (m_event->retrieve(tracks, "InDetTrackParticles").isFailure())
      fatal("JetHandler::recalculateJVT() : Couldn't get InDetTrackParticles, exiting");

    if (tracks == nullptr) {
      Warning("JetHandler::recalculateJVT()", "InDetTrackParticles retrieved as a nullptr?");
      rescaleJVT(jet);
      return;
    }

    int nPileupTracks = 0;
    for (auto track: *tracks) {
      if (track == nullptr)
        continue;

      if (m_trackTool->accept(*track, vertex)
          && track->vertex() 
          && track->vertex() != vertex
          && track->pt() < 30e3)
        nPileupTracks++;
    }
    if (nPileupTracks == 0)
      nPileupTracks = 1;

    std::vector<const xAOD::IParticle*> jetTracks;
    jet.getAssociatedObjects<xAOD::IParticle>(xAOD::JetAttribute::GhostTrack, jetTracks);

    double ptSum_all = 0.0, ptSum_pv = 0.0, ptSum_pu = 0.0;
    for (size_t i = 0; i < jetTracks.size(); ++i) {
      if (jetTracks[i] == nullptr)
        continue;

      const xAOD::TrackParticle *track = static_cast<const xAOD::TrackParticle*>(jetTracks[i]);

      bool accept = track->pt() > 500 && m_trackTool->accept(*track, vertex);

      if (accept) {
        ptSum_all += track->pt();

        if (track->vertex() == vertex ||
            (!track->vertex()
             && fabs((track->z0() + track->vz() - vertex->z())*sin(track->theta())) < 3.0))
          ptSum_pv += track->pt();

        if (track->vertex()
            && track->vertex() != vertex)
          ptSum_pu += track->pt();
      }
    }

    double myRpt     = ptSum_pv/jet.pt();
    double myCorrJVF = ptSum_pv + ptSum_pu > 0 ? ptSum_pv/(ptSum_pv + 100*ptSum_pu/nPileupTracks) : -1;
    double myJVT     = myCorrJVF >= 0 ? m_jvtLikelihood->Interpolate(myCorrJVF, std::min(myRpt, 1.0)) : -0.1;

    Jvt(jet) = myJVT;

  }

  //______________________________________________________________________________
  bool JetHandler::passJVTCut(const xAOD::Jet *jet, float jvtCut)
  {
    // Normal jet check
    if (jvtCut < 0)
      return true;

    if (jet->pt() < 50.0*HG::GeV &&
        fabs(jet->eta()) < 2.4   &&
        Jvt(*jet) < jvtCut) {
      return false;
    }

    return true;
  }

  //______________________________________________________________________________
  bool JetHandler::passJVFCut(const xAOD::Jet *jet, bool useBTagCut)
  {
    float jvf = useBTagCut ? m_bTagJvfCut : m_jvf;

    if (jet->pt() < 50.0*HG::GeV && fabs(jet->eta()) < 2.4 &&
        fabs(Jvf(*jet)) < jvf) {
      return false;
    }

    return true;
  }

  //______________________________________________________________________________
  void JetHandler::decorateJVF(xAOD::Jet *jet)
  {
    // JVF decoration
    const xAOD::VertexContainer* vertices = 0;
    if (!m_event->retrieve(vertices, "PrimaryVertices").isSuccess()) {
      fatal("Could not retrieve PrimaryVertices container.");
    }

    size_t pv = 0, npv = vertices->size();
    for (; pv < npv; ++pv) {
      if ((*vertices)[pv]->vertexType() == xAOD::VxType::VertexType::PriVtx) {
        Jvf(*jet) = JVF(*jet).at(pv);
        return;
      }
    }

    Jvf(*jet) = -1.0;
  }

  //______________________________________________________________________________
  void JetHandler::decorateBJet(xAOD::Jet *jet)
  {
    for (auto name: m_bTagNames) {
      for (auto eff: m_bTagEffs[name]) {
        jet->auxdata<char>((name+"_"+eff).Data()) = false;
        if (isMC()) {
          jet->auxdata<float>((name+"_"+eff+"_Eff").Data()) = 0.0;
          jet->auxdata<float>(("SF_"+name+"_"+eff).Data()) = 0.0;
        }
      }
    }
    
    if (fabs(jet->eta()) < m_bTagRapidityCut) {
      // Require JVT cut for MC15 tagging
      if(passJVTCut(jet, m_jvt)) {
        const xAOD::BTagging *bTagJet = jet->btagging();
        if (bTagJet == NULL)
          HG::fatal("Cannot access jet->btagging(). Missing in input file? If so, run with "+m_name+".EnableBTagging: NO");

        double disc = 0.0;
        float eff = 0.0;
        for (auto name: m_bTagNames) {
          for(size_t i = 0; i < m_bTagEffs[name].size(); i++) {

            bTagJet->MVx_discriminant(name.Data(), disc);
            if (disc > m_bTagCuts[name].at(i))
              jet->auxdata<char>((name+"_"+m_bTagEffs[name].at(i)).Data()) = true;

            if (isMC()) {
              CP::CorrectionCode result = m_bTagEffTools[name].at(i)->getMCEfficiency(*jet, eff);
              if (result == CP::CorrectionCode::Ok)
                jet->auxdata<float>((name+"_"+m_bTagEffs[name].at(i)+"_Eff").Data()) = eff;

              result = m_bTagEffTools[name].at(i)->getScaleFactor(*jet, eff);
              if (result == CP::CorrectionCode::Ok)
                jet->auxdata<float>(("SF_"+name+"_"+m_bTagEffs[name].at(i)).Data()) = eff;
            }

          }
        }
      }
    }

  }

  //______________________________________________________________________________
  void JetHandler::printJet(const xAOD::Jet *jet, TString comment)
  {
    // prints details about the photon

    printf("Jet %2zu  %s\n",
           jet->index(),comment.Data());
    
    // print the 4-vector
    printf("   (pT,eta,phi,m) = (%5.1f GeV,%6.3f,%6.3f,%4.1f GeV)\n",
           jet->pt()/GeV, jet->eta(), jet->phi(), jet->m()/GeV);
    
    // print some more information
    TString str;
    if (jet->isAvailable<float>("Ecalib_ratio"))
      str+=Form("   calibFactor = %.3f",jet->auxdata<float>("Ecalib_ratio"));

  }

}

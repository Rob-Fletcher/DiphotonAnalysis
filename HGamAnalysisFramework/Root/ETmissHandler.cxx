#include "HGamAnalysisFramework/ETmissHandler.h"

#include "HGamAnalysisFramework/HgammaUtils.h"

// MET EDM
#include "xAODMissingET/MissingET.h"
#include "xAODMissingET/MissingETContainer.h"

// MET Aux and association map
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETAssociationMap.h"

// METInterface includes
//#include "METInterface/IMETMaker.h"

// METMaker
//#include "METUtilities/METMaker.h"

#include "xAODBase/IParticleHelpers.h"

namespace HG {

  /*! \brief Class that
   *  \author Luis March
   */

  //______________________________________________________________________________
  ETmissHandler::ETmissHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
  , m_metMakerTool(nullptr)
  , m_metSysTool(nullptr)
  { }

  //______________________________________________________________________________
  ETmissHandler::~ETmissHandler()
  {
    SafeDelete(m_metMakerTool);
    SafeDelete(m_metSysTool);
  }

  //______________________________________________________________________________
  EL::StatusCode ETmissHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // Read in configuration information
    m_containerName = config.getStr(m_name+".ContainerName");
    m_truthName     = config.getStr(m_name+".TruthContainerName", "MET_Truth");

    m_assocMapName  = config.getStr(m_name+".METAssociactionMapName", "METAssoc_AntiKt4EMTopo");
    m_coreName      = config.getStr(m_name+".METCoreName", "METAssoc_AntiKt4EMTopo");

    m_metTypes      = config.getStrV(m_name+".METTypes");
    
    // METMaker tool
    m_metMakerTool = new met::METMaker("METMaker");   // Defining m_metMakerTool as re-builder of ETmiss tool
    if (m_metMakerTool->initialize().isFailure())
      fatal("Failed to initialize METMakerTool");

    // METSystematics tool
    m_metSysTool = new met::METSystematicsTool("METSystematicsTool");
    CP_CHECK(m_name, m_metSysTool->setProperty("ConfigJetTrkFile", "JetTrackSyst.config"                                                   ));
    CP_CHECK(m_name, m_metSysTool->setProperty("JetColl"         , config.getStr("JetHandler.JetContainerName", "AntiKt4EMTopoJets").Data()));
    if (m_metSysTool->initialize().isFailure())
      fatal("Failed to initialize METSystematicTool");

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer ETmissHandler::getCorrectedContainer()
  {
    xAOD::MissingETContainer dummyContainer;
    return dummyContainer;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer ETmissHandler::getCorrectedContainer(const xAOD::PhotonContainer   *photons  ,
                                                                const xAOD::JetContainer      *jets     ,
                                                                const xAOD::ElectronContainer *electrons,
                                                                const xAOD::MuonContainer     *muons    )
  {
    // Get Shallow copy from TEvent/TStore
    bool calib = false;
    // second argument false --> make empty raw contianer (not shallow copy of xAOD), since it's rebuilt below
    xAOD::MissingETContainer shallowContainer = getShallowContainer(calib, false);
    if (calib) return shallowContainer;

    // Retrieve the MET association map: Needed for METMaker tool
    const xAOD::MissingETAssociationMap *metMap = nullptr;
    if (m_event->retrieve(metMap, m_assocMapName).isFailure())
      fatal("Unable to retrieve MissingETAssociationMap from TEvent");

    // Retrieve the MET core container: Needed for METMaker tool
    const xAOD::MissingETContainer *coreMet  = nullptr;
    if (m_event->retrieve(coreMet, m_coreName).isFailure())
      fatal("Unable to retrieve coreMet from TEvent");

    // For MET rebuilding, need access to the actual container linked with the auxdata
    // which is already in TStore (HgammaHandler magic)
    TString shallowName = "Shallow" + m_containerName + m_sysName;
    xAOD::MissingETContainer *shallowMet  = nullptr;
    if (m_store->retrieve(shallowMet, shallowName.Data()).isFailure())
      fatal("Unable to retrieve ShallowMET from TEvent");

    // Reset the MET map before each building
    metMap->resetObjSelectionFlags();

    // Rebuild ETmiss RefGamma, RefEle, and Muons terms
    HG_CHECK(m_name, m_metMakerTool->rebuildMET("RefGamma", xAOD::Type::Photon  , shallowMet, photons  , metMap));
    HG_CHECK(m_name, m_metMakerTool->rebuildMET("RefEle"  , xAOD::Type::Electron, shallowMet, electrons, metMap));
    HG_CHECK(m_name, m_metMakerTool->rebuildMET("Muons"   , xAOD::Type::Muon    , shallowMet, muons    , metMap));

    // Jet and Soft Terms
    HG_CHECK(m_name, m_metMakerTool->rebuildJetMET("RefJet", "SoftClus", "PVSoftTrk", shallowMet, jets, coreMet, metMap, true));

    // Apply possible systematic uncertainty shifts
    xAOD::MissingET *softClusMet = (*shallowMet)["SoftClus"];
    if (softClusMet == nullptr)
      fatal("Couldn't retrieve SoftClus from shallowMet, exiting!");
    CC_CHECK(m_name, m_metSysTool->applyCorrection(*softClusMet));

    xAOD::MissingET *softTrkMet = (*shallowMet)["PVSoftTrk"];
    if (softTrkMet == nullptr)
      fatal("Couldn't retrieve PVSoftTrk from shallowMet, exiting!");
    CC_CHECK(m_name, m_metSysTool->applyCorrection(*softTrkMet));

    // Rebuild full MET
    HG_CHECK(m_name, m_metMakerTool->buildMETSum("CST", shallowMet, MissingETBase::Source::LCTopo));
    HG_CHECK(m_name, m_metMakerTool->buildMETSum("TST", shallowMet, MissingETBase::Source::Track ));

    // Return corrected MET
    xAOD::MissingETContainer corrected(shallowMet->begin(),
                                       shallowMet->end()  ,
                                       SG::VIEW_ELEMENTS  );
    return corrected;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer ETmissHandler::applySelection(xAOD::MissingETContainer &container)
  {
    xAOD::MissingETContainer selected(SG::VIEW_ELEMENTS);
    for (auto met: container) {
      // Limit MET to the types specified in config (TST, CST, ...)
      if (std::find(m_metTypes.begin(), m_metTypes.end(), met->name().c_str()) == m_metTypes.end())
        continue;
      
      selected.push_back(met);
    }

    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode ETmissHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool isAffected = false;
    for (auto var: sys) {
      if (m_metSysTool->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
    }

    if (isAffected) {
      CP_CHECK(m_name, m_metSysTool->applySystematicVariation(sys));
    } else {
      CP_CHECK(m_name, m_metSysTool->applySystematicVariation(CP::SystematicSet()));
    }

    // For MET, always make new container since shifting hard objects indirectly affects result
    m_sysName = sys.name() == "" ? "" : "_"+sys.name();

    return CP::SystematicCode::Ok;
  }

}

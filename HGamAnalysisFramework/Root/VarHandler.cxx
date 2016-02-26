#include "HGamAnalysisFramework/VarHandler.h"

#include "HGamAnalysisFramework/HgammaUtils.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODCore/AuxInfoBase.h"

namespace HG {

  //____________________________________________________________________________
  VarHandler* VarHandler::m_ptr = nullptr;

  //____________________________________________________________________________
  VarHandler::VarHandler()
  : m_sysName("")
  , m_MxAODName("HGam")
  , m_event(nullptr)
  , m_store(nullptr)
  { }

  //____________________________________________________________________________
  VarHandler::~VarHandler()
  { }

  //____________________________________________________________________________
  VarHandler* VarHandler::getInstance()
  {
    if (m_ptr == nullptr)
      m_ptr = new VarHandler();
    return m_ptr;
  }

  //____________________________________________________________________________
  CP::SystematicCode VarHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    m_MxAODName = "HGam";
    m_sysName = sys.name() == "" ? "" : "_" + sys.name();
    return CP::SystematicCode::Ok;
  }

  //____________________________________________________________________________
  const xAOD::EventInfo* VarHandler::getEventInfoFromEvent()
  {
    std::string name = getEventInfoName();

    const xAOD::EventInfo *eventInfo = nullptr;
    if (!m_event->contains<xAOD::EventInfo>(name))
      return nullptr;

    if (m_event->retrieve(eventInfo, name).isFailure())
      return nullptr;

    return eventInfo;
  }

  //____________________________________________________________________________
  xAOD::EventInfo* VarHandler::getEventInfoFromStore(bool createInfo)
  {
    std::string name = getEventInfoName();

    xAOD::EventInfo *eventInfo = nullptr;
    if (m_store->contains<xAOD::EventInfo>(name)) {
      if (m_store->retrieve(eventInfo, name).isFailure()) {
        return nullptr;
      }
      return eventInfo;
    }

    if (!createInfo) return nullptr;

    const xAOD::EventInfo *constInfo = nullptr;
    if (m_event->contains<xAOD::EventInfo>(name)) {
      if (m_event->retrieve(constInfo, name).isFailure())
        return nullptr;
    }

    eventInfo = new xAOD::EventInfo();
    xAOD::AuxInfoBase *eventInfoAux = new xAOD::AuxInfoBase();
    eventInfo->setStore(eventInfoAux);

    if (constInfo != nullptr)
      *eventInfo = *constInfo;

    if (!m_store->record(eventInfo, name))
      return nullptr;

    name += "Aux";
    if (!m_store->record(eventInfoAux, name))
      return nullptr;

    return eventInfo;
  }

  //____________________________________________________________________________
  const xAOD::EventInfo* VarHandler::getTruthEventInfoFromEvent()
  {
    // Save reco names
    std::string sysName   = m_sysName;
    std::string MxAODName = m_MxAODName;

    // Set true names
    m_sysName = "";
    m_MxAODName = "HGamTruth";
    
    // Retrieve eventInfo
    const xAOD::EventInfo *eventInfo =  getEventInfoFromEvent();

    // Reset reco names
    m_sysName = sysName;
    m_MxAODName = MxAODName;

    return eventInfo;
  }

  //____________________________________________________________________________
  xAOD::EventInfo* VarHandler::getTruthEventInfoFromStore(bool createInfo)
  {
    // Save reco names
    std::string sysName   = m_sysName;
    std::string MxAODName = m_MxAODName;

    // Set true names
    m_sysName = "";
    m_MxAODName = "HGamTruth";
    
    // Retrieve eventInfo
    xAOD::EventInfo *eventInfo =  getEventInfoFromStore(createInfo);

    // Reset reco names
    m_sysName = sysName;
    m_MxAODName = MxAODName;

    return eventInfo;
  }

  //____________________________________________________________________________
  std::string VarHandler::getEventInfoName() const
  {
    return m_MxAODName + "EventInfo" + m_sysName;
  }

  //____________________________________________________________________________
  const xAOD::IParticleContainer* VarHandler::getPhotons(bool truth) const
  {
    if (truth) return &m_truthPhotons;
    return &m_photons;
  }

  //____________________________________________________________________________
  const xAOD::IParticleContainer* VarHandler::getJets(bool truth) const
  {
    if (truth) return &m_truthJets;
    return &m_jets;
  }

  //____________________________________________________________________________
  const xAOD::IParticleContainer* VarHandler::getElectrons(bool truth) const
  {
    if (truth) return &m_truthElectrons;
    return &m_electrons;
  }

  //____________________________________________________________________________
  const xAOD::IParticleContainer* VarHandler::getMuons(bool truth) const
  {
    if (truth) return &m_truthMuons;
    return &m_muons;
  }

  //____________________________________________________________________________
  const xAOD::IParticleContainer* VarHandler::getHiggsBosons() const
  {
    return &m_higgsBosons;
  }

  //____________________________________________________________________________
  void VarHandler::setEventAndStore(xAOD::TEvent *event, xAOD::TStore *store)
  {
    m_event = event;
    m_store = store;
  }

  //____________________________________________________________________________
  void VarHandler::setContainers(const xAOD::IParticleContainer *photons  ,
                                 const xAOD::IParticleContainer *electrons,
                                 const xAOD::IParticleContainer *muons    ,
                                 const xAOD::IParticleContainer *jets     )
  {
    if (photons  ) m_photons   = *photons;
    if (electrons) m_electrons = *electrons;
    if (muons    ) m_muons     = *muons;
    if (jets     ) m_jets      = *jets;
  }

  //____________________________________________________________________________
  void VarHandler::setTruthContainers(const xAOD::IParticleContainer *photons  ,
                                      const xAOD::IParticleContainer *electrons,
                                      const xAOD::IParticleContainer *muons    ,
                                      const xAOD::IParticleContainer *jets     )
  {
    if (photons  ) m_truthPhotons   = *photons;
    if (electrons) m_truthElectrons = *electrons;
    if (muons    ) m_truthMuons     = *muons;
    if (jets     ) m_truthJets      = *jets;
  }

  //____________________________________________________________________________
  void VarHandler::setHiggsBosons(const xAOD::IParticleContainer *higgs)
  {
    if (higgs) m_higgsBosons = *higgs;
  }

  //____________________________________________________________________________
  void VarHandler::clearContainers()
  {
    // Clear reco containers
    m_photons.clear();
    m_electrons.clear();
    m_muons.clear();
    m_jets.clear();

    // Clear truth containers
    m_truthPhotons.clear();
    m_truthElectrons.clear();
    m_truthMuons.clear();
    m_truthJets.clear();
  }

  //____________________________________________________________________________
  EL::StatusCode VarHandler::write()
  {
    xAOD::EventInfo *eventInfo = getEventInfoFromStore();
    if (eventInfo == nullptr)
      return EL::StatusCode::FAILURE;

    xAOD::EventInfo *copy = new xAOD::EventInfo();
    xAOD::AuxInfoBase *copyAux = new xAOD::AuxInfoBase();
    copy->setStore(copyAux);
    *copy = *eventInfo;

    std::string name = getEventInfoName();
    if (!m_event->record(copy, name))
      return EL::StatusCode::FAILURE;

    name += "Aux.";
    if (!m_event->record(copyAux, name))
      return EL::StatusCode::FAILURE;

    return EL::StatusCode::SUCCESS;
  }

  //____________________________________________________________________________
  EL::StatusCode VarHandler::writeTruth()
  {
    // Save reco names
    std::string sysName   = m_sysName;
    std::string MxAODName = m_MxAODName;

    // Set true names
    m_sysName = "";
    m_MxAODName = "HGamTruth";
    
    // Retrieve eventInfo
    EL::StatusCode code = write();

    // Reset reco names
    m_sysName = sysName;
    m_MxAODName = MxAODName;

    return code;
  }

}

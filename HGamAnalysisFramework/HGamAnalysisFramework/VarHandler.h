#ifndef HGamAnalysisFramework_VarHandler_H
#define HGamAnalysisFramework_VarHandler_H

#include "EventLoop/StatusCode.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"

#include "PATInterfaces/SystematicRegistry.h"

#include <string>

namespace xAOD {
  class TEvent;
  class TStore;
  class AuxInfoBase;
}

namespace HG {

  class VarHandler {
  private:
    static VarHandler        *m_ptr;

    std::string               m_sysName;
    std::string               m_MxAODName;

    xAOD::TEvent             *m_event;
    xAOD::TStore             *m_store;

    xAOD::IParticleContainer  m_photons;
    xAOD::IParticleContainer  m_electrons;
    xAOD::IParticleContainer  m_muons;
    xAOD::IParticleContainer  m_jets;

    xAOD::IParticleContainer  m_truthPhotons;
    xAOD::IParticleContainer  m_truthElectrons;
    xAOD::IParticleContainer  m_truthMuons;
    xAOD::IParticleContainer  m_truthJets;
    xAOD::IParticleContainer  m_higgsBosons;


  public:
    /// Get instance of singleton class
    static VarHandler* getInstance();

    /// Get MxAOD EventInfo name for event
    std::string        getEventInfoName() const;

    /// Get pointer to collection
    const xAOD::IParticleContainer* getPhotons    (bool truth = false) const;
    const xAOD::IParticleContainer* getJets       (bool truth = false) const;
    const xAOD::IParticleContainer* getElectrons  (bool truth = false) const;
    const xAOD::IParticleContainer* getMuons      (bool truth = false) const;
    const xAOD::IParticleContainer* getHiggsBosons() const;

    /// Get MxAOD EventInfo for event
    xAOD::EventInfo*       getEventInfoFromStore(bool createInfo = true);
    const xAOD::EventInfo* getEventInfoFromEvent();

    /// Get MxAOD TruthEventInfo for event
    xAOD::EventInfo*       getTruthEventInfoFromStore(bool createInfo = true);
    const xAOD::EventInfo* getTruthEventInfoFromEvent();

    /// Set TEvent and TStore
    void               setEventAndStore(xAOD::TEvent *event, xAOD::TStore *store);

    /// Set current systematic variation
    CP::SystematicCode applySystematicVariation(const CP::SystematicSet &sys);

    /// Set object containers
    void setContainers     (const xAOD::IParticleContainer *photons   = nullptr,
                            const xAOD::IParticleContainer *electrons = nullptr,
                            const xAOD::IParticleContainer *muons     = nullptr,
                            const xAOD::IParticleContainer *jets      = nullptr);

    /// Set truth object containers
    void setTruthContainers(const xAOD::IParticleContainer *photons   = nullptr,
                            const xAOD::IParticleContainer *electrons = nullptr,
                            const xAOD::IParticleContainer *muons     = nullptr,
                            const xAOD::IParticleContainer *jets      = nullptr);
    void setHiggsBosons    (const xAOD::IParticleContainer *higgs);

    /// Reset containers to null pointers to avoid carry-over from previous event
    void               clearContainers();

    /// Write MxAOD EventInfo to output
    EL::StatusCode     write();

    /// Write MxAOD TruthEventInfo to output
    EL::StatusCode     writeTruth();



  private:
    /// Default constructor
    VarHandler();

    /// Default detructor
    ~VarHandler();



  }; // class VarHandler


  template <class T>
  class VarBase {
  protected:
    T                                m_default;
    std::string                      m_name;
    SG::AuxElement::Decorator<T>     decVal;
    SG::AuxElement::Accessor<T>      accVal;



  protected:



  protected:
    /// Calculate variable of interest, should be defined by inherited class
    virtual T calculateValue(bool truth = false);

    /// Get variable of interest, first by checking TEvent, then calculateVarlue()
    T getValue(bool truth = false);

    /// Check if variable exists in TEvent. If so return true, and set value by reference
    bool checkInEvent(T &value, bool truth = false);
    bool checkInStore(T &value, bool truth = false);



  public:
    /// Default constructor
    VarBase(const std::string &name);

    /// Default constructor
    virtual ~VarBase();

    /// Add variable to the current EventInfo
    /// @truth: if true, calculates value using truth containers
    void addToStore(bool truth);

    /// Set value manually
    void setValue(const T &value);
    void setTruthValue(const T &value);

    /// Check if value exists in TEvent or TStore
    bool exists();

    /// Get variable of interest
    T operator()();

    /// Get truth variable of interest
    T truth();



  }; // class VarBase

}

#include "HGamAnalysisFramework/VarHandler.hpp"

#endif // HGamAnalysisFramework_VarHandler_H

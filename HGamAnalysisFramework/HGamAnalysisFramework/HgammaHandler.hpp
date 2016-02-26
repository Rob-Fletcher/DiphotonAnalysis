#ifndef HGamAnalysisFramework_ObjectHandler_HPP
#define HGamAnalysisFramework_ObjectHandler_HPP

#include "xAODBase/IParticleHelpers.h"

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"

namespace HG {
  
  template <class partType, class partContainer, class auxContainer>
  HgammaHandler<partType, partContainer, auxContainer>::HgammaHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : m_name(name)
  , m_event(event)
  , m_store(store)
  , m_trackIsoTool(nullptr)
  , m_sysName("nominal")
  , m_MxAODname("HGam")
  , m_isMxAOD(false)
  , m_isData(false)
  , m_isMC(false)
  , m_isVtxCorr(false)
  { }

  template <class partType, class partContainer, class auxContainer>
  HgammaHandler<partType, partContainer, auxContainer>::~HgammaHandler()
  {
    SafeDelete(m_trackIsoTool);
  }

  template <class partType, class partContainer, class auxContainer>
  EL::StatusCode
  HgammaHandler<partType, partContainer, auxContainer>::initialize(Config &config)
  {
    // Check whether it's simulation or data
    m_isMC   = isMC();
    m_isData = !m_isMC;

    // Check if it's an MxAOD being run over
    m_isMxAOD = config.getBool("IsMxAOD", false);

    // Track selection tool used by isolation correction tool
    InDet::InDetTrackSelectionTool *trackSel = new InDet::InDetTrackSelectionTool("TrackSel");
    CP_CHECK(m_name, trackSel->setProperty("CutLevel"     , "Loose"));
    CP_CHECK(m_name, trackSel->setProperty("minPt"        , 1000.0 ));
    CP_CHECK(m_name, trackSel->setProperty("maxZ0SinTheta", 3.0    ));

    if (trackSel->initialize().isFailure())
      fatal("Cannot initialize track selection tool in HgammaHandler, exiting");

    ToolHandle<InDet::IInDetTrackSelectionTool> th(trackSel);

    // Track isolation correction tool
    m_trackIsoTool = new xAOD::TrackIsolationTool("trackIsoTool");
    CP_CHECK(m_name, m_trackIsoTool->setProperty("TrackSelectionTool", th));

    if (m_trackIsoTool->initialize().isFailure())
      fatal("Cannot initialize track isolation tool in HgammaHandler, exiting");


    // Variables needed by isolation correction tool
    m_isoT = {xAOD::Iso::ptcone40, xAOD::Iso::ptcone30, xAOD::Iso::ptcone20};
    m_corrList.trackbitset.set(static_cast<unsigned int>(xAOD::Iso::coreTrackPtr));

    return EL::StatusCode::SUCCESS;
  }

  template <class partType, class partContainer, class auxContainer>
  void
  HgammaHandler<partType, partContainer, auxContainer>::correctPrimaryVertex(const xAOD::Vertex *,
                                                                             partType &)
  {
    // Should be implemented by inheriting class for something to happen
    // By default, nothing happens (and shouldn't crash or complain)
  }

  template <class partType, class partContainer, class auxContainer>
  void
  HgammaHandler<partType, partContainer, auxContainer>::correctContainerPV(partContainer &cont)
  {
    // Check if correction is already applied, then if the updated hard-scatter vertex exists
    if (m_isVtxCorr) return;
    if (!m_store->contains<ConstDataVector<xAOD::VertexContainer> >("HGamVertices"))
      return;

    // Retrieve updated hard-scatter vertex from TStore
    ConstDataVector<xAOD::VertexContainer> *vertices = nullptr;
    if (m_store->retrieve(vertices, "HGamVertices").isFailure()) {
      fatal("Cannot access HGamVertices container");
    }

    if (vertices->size() == 0) return;
    const xAOD::Vertex *vertex = (*vertices)[0];

    // Correct objects in collection
    for (auto part: cont)
      correctPrimaryVertex(vertex, *part);

    return;
  }

  template <class partType, class partContainer, class auxContainer>
  partContainer
  HgammaHandler<partType, partContainer, auxContainer>::getShallowContainer(bool &calibStatus, bool makeShallowCopy)
  {
    calibStatus = true;

    TString rawName = m_containerName + m_sysName;
    TString aodName = m_MxAODname + m_containerName + m_sysName;
    TString name = "";

    // 1. Check for MxAOD collection (assumed calibrated) in TStore, stored by a previous
    // call to this Handler instance
    name = "Shallow" + aodName;
    if (m_isMxAOD && m_store->contains<partContainer>(name.Data())) {
      return getStoreContainer(name.Data());
    }

    // 2. Check for MxAOD collection (assumed calibrated) in TEvent (not yet accessed)
    // and make a shallow copy for TStore if it exists
    name = aodName;
    if (m_isMxAOD && m_event->contains<partContainer>(name.Data())) {
      return getEventContainer(name.Data());
    }

    // 3. Check for already calibrated "raw" collection in TStore, from a previous
    // call to this Handler instance
    name = "Shallow" + rawName;
    if (m_store->contains<partContainer>(name.Data())) {
      partContainer preContainer = getStoreContainer(name.Data());
      correctContainerPV(preContainer);
      return preContainer;
    }

    calibStatus = false;

    // 4. If all else fails, get the uncalibrated collection from TEvent and make a shallow
    // copy for TStore (calibStatus flag is now false)
    rawName.ReplaceAll(m_sysName, "");
    return getEventContainer(rawName.Data(), m_sysName.Data(), makeShallowCopy);
  }

  template <class partType, class partContainer, class auxContainer>
  partContainer
  HgammaHandler<partType, partContainer, auxContainer>::getStoreContainer(std::string name)
  {
    // Get shallow copy from TStore
    partContainer *storeContainer = NULL;
    if (m_store->retrieve(storeContainer, name).isFailure()) {
      fatal("Cannot access container");
    }

    // Make the container from the TStore'd shallow copy
    partContainer container(storeContainer->begin(),
                            storeContainer->end(),
                            SG::VIEW_ELEMENTS);
    return container;
  }

  template <class partType, class partContainer, class auxContainer>
  partContainer
  HgammaHandler<partType, partContainer, auxContainer>::getContainer(std::string name)
  {
    return getEventContainer(name, "", true);
  }

  template <class partType, class partContainer, class auxContainer>
  partContainer
  HgammaHandler<partType, partContainer, auxContainer>::getEventContainer(std::string name, std::string post, bool makeShallowCopy)
  {
    // Get raw container from TEvent
    const partContainer *rawContainer = NULL;
    if (!m_event->contains<partContainer>(name))
      fatal("In getEventContainer, cannot access \""+name+"\" in the input file");

    if (m_event->retrieve(rawContainer, name).isFailure())
      fatal("Cannot access container");

    // Make a shallow copy of the container
    std::pair<partContainer*, xAOD::ShallowAuxContainer*> shallowCopy;
    auxContainer *fullStore = nullptr;
    if (makeShallowCopy) {
      shallowCopy = xAOD::shallowCopyContainer(*rawContainer);
    } else {
      // shallowCopy = std::make_pair<partContainer*, xAOD::ShallowAuxContainer*>(new partContainer(), nullptr);
      shallowCopy.first = new partContainer();
      fullStore         = new auxContainer();
      shallowCopy.first->setStore(fullStore);
    }

#ifndef __HGamMET__
    static SG::AuxElement::ConstAccessor<ElementLink<xAOD::IParticleContainer> > acc("originalObjectLink");
    if (rawContainer->size() > 0 && !acc.isAvailable(*rawContainer->at(0))) {
      if (!xAOD::setOriginalObjectLink(*rawContainer, *shallowCopy.first)) {
        fatal("Failed to add original object links to shallow copy");
      }
    }
#endif

    TString shallowName = TString::Format("Shallow%s%s", name.c_str(), post.c_str());

    // Add shallow copy to TStore, so it handles memory cleanup
    if (m_store->record(shallowCopy.first, shallowName.Data()).isFailure()) {
      fatal("Cannot store container in TStore");
    }
    shallowName += "Aux";
    if (makeShallowCopy && m_store->record(shallowCopy.second, shallowName.Data()).isFailure())
      fatal("Cannot store aux container in TStore");
    if (!makeShallowCopy && m_store->record(fullStore, shallowName.Data()).isFailure())
      fatal("Cannot store aux container in TStore");
    
    // Make the container from the TStore'd shallow copy
    partContainer container(shallowCopy.first->begin(),
                            shallowCopy.first->end(),
                            SG::VIEW_ELEMENTS);

    return container;
  }

  template <class partType, class partContainer, class auxContainer>
  partContainer
  HgammaHandler<partType, partContainer, auxContainer>::getShallowContainer(std::string name)
  {
    // DEPRECATED!!!

    // By default, use name specified in config file
    if (name == "") name = m_containerName.Data();

    // get raw container
    const partContainer *rawContainer = NULL;
    if (m_event->retrieve(rawContainer, name).isFailure()) {
      fatal("Cannot access container");
    }

    // make a shallow copy of the container
    std::pair<partContainer*, xAOD::ShallowAuxContainer*> shallowCopy = xAOD::shallowCopyContainer(*rawContainer);

    // add shallow copy to TStore, so it handles memory cleanup
    TString shallowName = "Shallow" + name + m_sysName;
    if (m_store->record(shallowCopy.first, shallowName.Data()).isFailure())
      fatal("Cannot store container in TStore");

    shallowName += "Aux";
    if (m_store->record(shallowCopy.second, shallowName.Data()).isFailure())
      fatal("Cannot store aux container in TStore");
    
    // make the container
    partContainer container(shallowCopy.first->begin(),
                            shallowCopy.first->end(),
                            SG::VIEW_ELEMENTS);

    return container;
  }

  template <class partType, class partContainer, class auxContainer>
  EL::StatusCode
  HgammaHandler<partType, partContainer, auxContainer>::writeTruthContainer(partContainer &container, TString name)
  {
    if (name == "") {
      name = m_MxAODname + m_truthName;
    }
    return writeContainer(container, name);
  }

  template <class partType, class partContainer, class auxContainer>
  EL::StatusCode
  HgammaHandler<partType, partContainer, auxContainer>::writeContainer(partContainer &container, TString name)
  {
    // Check for name of collection to write
    if (name == "") {
      name = m_MxAODname + m_containerName + m_sysName;
    }

    // Create the new container and its auxilliary store
    partContainer          *output    = new partContainer();
    xAOD::AuxContainerBase *outputAux = new xAOD::AuxContainerBase();
    output->setStore(outputAux);

    for (auto part: container) {
      // Copy to output container
      partType *save = new partType();
      output->push_back(save);
      *save = *part;
    }

    if (m_event->record(output, name.Data()).isFailure())
      return EL::StatusCode::FAILURE;

    name += "Aux.";
    if (m_event->record(outputAux, name.Data()).isFailure())
      return EL::StatusCode::FAILURE;

    return EL::StatusCode::SUCCESS;
  }

  template <class partType, class partContainer, class auxContainer>
  bool
  HgammaHandler<partType, partContainer, auxContainer>::isMC()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      fatal("Cannot access EventInfo");
    }
    return eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION);
  }
  
  template <class partType, class partContainer, class auxContainer>
  bool
  HgammaHandler<partType, partContainer, auxContainer>::comparePt(const partType *a, const partType *b)
  {
     return a->pt() > b->pt();
  }
  
}



#endif // HGamAnalysisFramework_ObjectHandler_HPP

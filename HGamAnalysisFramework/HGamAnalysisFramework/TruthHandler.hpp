#ifndef HGamAnalysisFramework_TruthHandler_HPP
#define HGamAnalysisFramework_TruthHandler_HPP

#include "xAODBase/IParticleHelpers.h"

namespace HG {

  //____________________________________________________________________________
  template <class T>
  void
  TruthHandler::setP4(const xAOD::IParticle *p, T &copy)
  { copy.setP4(p->pt(), p->eta(), p->phi(), p->m()); }

  //____________________________________________________________________________
  template <>
  void
  TruthHandler::setP4<xAOD::Muon>(const xAOD::IParticle *p, xAOD::Muon &copy);

  //____________________________________________________________________________
  template <>
  void
  TruthHandler::setP4<xAOD::TruthParticle>(const xAOD::IParticle *p, xAOD::TruthParticle &copy);

  //____________________________________________________________________________
  template <class T>
  void
  TruthHandler::setOriginalObjectLink(const T *orig, T *copy)
  {
    xAOD::setOriginalObjectLink(*orig, *copy);
  }

  //____________________________________________________________________________
  template <>
  void
  TruthHandler::setOriginalObjectLink<xAOD::MissingET>(const xAOD::MissingET *orig, xAOD::MissingET *copy);

  //____________________________________________________________________________
  template <class T>
  void
  TruthHandler::setOriginalObjectLink(const DataVector<T> *orig, DataVector<T> *copy)
  {
    xAOD::setOriginalObjectLink(*orig, *copy);
  }

  //____________________________________________________________________________
  template <>
  void
  TruthHandler::setOriginalObjectLink<xAOD::MissingET>(const DataVector<xAOD::MissingET> *orig, DataVector<xAOD::MissingET> *copy);

  //____________________________________________________________________________
  template <class T>
  DataVector<T>
  TruthHandler::getDeepCopy(const DataVector<T> *parts, std::string name)
  {
    // Create the new container and its auxilliary store
    DataVector<T>          *copy    = new DataVector<T>();
    xAOD::AuxContainerBase *copyAux = new xAOD::AuxContainerBase();
    copy->setStore(copyAux);

    for (auto part: *parts) {
      // Copy to output container
      T *element = new T();
      copy->push_back(element);
      *element = *part;
      setOriginalObjectLink(part, element);
    }

    name = "Deep" + name;
    if (m_store->record(copy, name).isFailure())
      fatal("Cannot store deep copy of Truth to TStore, exiting.");

    name = name + "Aux";
    if (m_store->record(copyAux, name).isFailure())
      fatal("Cannot store deep copy of TruthAux to TStore, exiting.");

    // Make a container for returning
    DataVector<T> cont(copy->begin(),
                       copy->end(),
                       SG::VIEW_ELEMENTS);

    return cont;
  }

  //____________________________________________________________________________
  template <class T>
  DataVector<T>
  TruthHandler::getShallowCopy(const DataVector<T> *parts, std::string name)
  {
    // Make a shallow copy
    std::pair<DataVector<T>*, xAOD::ShallowAuxContainer*> shallowCopy = xAOD::shallowCopyContainer(*parts);
    setOriginalObjectLink(parts, shallowCopy.first); 

    // Add shallowcopy to TStore
    name = "Shallow" + name;
    if (m_store->record(shallowCopy.first, name).isFailure())
      fatal("Cannot store shallow copy of Truth to TStore, exiting.");

    name = name + "Aux";
    if (m_store->record(shallowCopy.second, name).isFailure())
      fatal("Cannot store shallow copy of TruthAux to TStore, exiting.");

    // Make a container for returning
    DataVector<T> cont(shallowCopy.first->begin(),
                       shallowCopy.first->end(),
                       SG::VIEW_ELEMENTS);

    return cont;
  }

  //____________________________________________________________________________
  template <class T>
  DataVector<T>
  TruthHandler::getContainer(const TruthContainer *parts, std::string name)
  {
    DataVector<T> *cont = new DataVector<T>();
    xAOD::AuxContainerBase *contAux = new xAOD::AuxContainerBase();
    cont->setStore(contAux);

    for (auto part: *parts) {
      T *truth = new T();
      cont->push_back(truth);
      setP4<T>(part, *truth);
    }

    if (m_store->record(cont, name).isFailure())
      HG::fatal("Couldn't store Truth container in TStore, exiting.");

    name += "Aux";
    if (m_store->record(contAux, name).isFailure())
      HG::fatal("Couldn't store TruthAux container in TStore, exiting.");

    DataVector<T> container(cont->begin(), cont->end(), SG::VIEW_ELEMENTS);
    return container;
  }

  //____________________________________________________________________________
  template <class T>
  bool
  TruthHandler::checkEventAndStore(DataVector<T> &cont, std::string name)
  {
    // First check in TStore
    if (m_store->contains<DataVector<T> >(name)) {
      // Get shallow copy from TStore
      DataVector<T> *storeContainer = nullptr;
      if (m_store->retrieve(storeContainer, name).isFailure())
        fatal("Cannot access Truth container in TStore, exiting.");

      // Make a container for returning
      cont = DataVector<T>(storeContainer->begin(),
                           storeContainer->end(),
                           SG::VIEW_ELEMENTS);
      return true;
    }

    // Now check in TEvent
    if (m_event->contains<DataVector<T> >(name)) {
      // Get container from TEvent
      const DataVector<T> *eventContainer = nullptr;
      if (m_event->retrieve(eventContainer, name).isFailure())
        fatal("Cannot access Truth container in TEvent, exiting.");

      // Make a shallow copy
      cont = getShallowCopy(eventContainer, name);
      return true;
    }

    // Doesn't exist, return false
    return false;
  }

  //______________________________________________________________________________
  template <class T>
  EL::StatusCode
  TruthHandler::writeContainer(DataVector<T> &container, std::string name)
  {
    // Create the new container and its auxilliary store
    DataVector<T>          *output    = new DataVector<T>();
    xAOD::AuxContainerBase *outputAux = new xAOD::AuxContainerBase();
    output->setStore(outputAux);

    for (auto part: container) {
      // Copy to output container
      T *save = new T();
      output->push_back(save);
      *save = *part;
    }

    if (m_event->record(output, name).isFailure())
      return EL::StatusCode::FAILURE;

    name += "Aux.";
    if (m_event->record(outputAux, name).isFailure())
      return EL::StatusCode::FAILURE;

    return EL::StatusCode::SUCCESS;
  }

}

#endif // HGamAnalysisFramework_TruthHandler_H

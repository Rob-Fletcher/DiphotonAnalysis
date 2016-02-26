#ifndef HGamAnalysisFramework_EventHandler_HPP
#define HGamAnalysisFramework_EventHandler_HPP

#include "HGamAnalysisFramework/VarHandler.h"

namespace HG {
  
  //______________________________________________________________________________
  template <typename T>
  T EventHandler::getTruthVar(const char *name)
  {
    xAOD::EventInfo *eventInfo = HG::VarHandler::getInstance()->getTruthEventInfoFromStore();
    if (eventInfo != nullptr) {
      return eventInfo->auxdata<T>(name);
    }

    const xAOD::EventInfo *constEventInfo = HG::VarHandler::getInstance()->getTruthEventInfoFromEvent();
    if (constEventInfo == nullptr) {
      HG::fatal("EventHandler::getVar() cannot access EventInfo");
    }

    return constEventInfo->auxdata<T>(name);
  }

  //______________________________________________________________________________
  template <typename T>
  T EventHandler::getVar(const char *name)
  {
    xAOD::EventInfo *eventInfo = HG::VarHandler::getInstance()->getEventInfoFromStore();
    if (eventInfo != nullptr) {
      return eventInfo->auxdata<T>(name);
    }

    const xAOD::EventInfo *constEventInfo = HG::VarHandler::getInstance()->getEventInfoFromEvent();
    if (constEventInfo == nullptr) {
      HG::fatal("EventHandler::getVar() cannot access EventInfo");
    }

    return constEventInfo->auxdata<T>(name);
  }

  //______________________________________________________________________________
  template <typename T>
  void EventHandler::storeVariable(const char *name, T value)
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      fatal("Cannot access EventInfo");
    }
    eventInfo->auxdecor<T>(name) = value;
  }

  //______________________________________________________________________________
  template <typename T>
  void EventHandler::storeTruthVar(const char *name, T value)
  {
    xAOD::EventInfo *eventInfo = HG::VarHandler::getInstance()->getTruthEventInfoFromStore();
    if (eventInfo == nullptr) {
      HG::fatal("EventHandler cannot access EventInfo");
    }
    eventInfo->auxdecor<T>(name) = value;
  }

  //______________________________________________________________________________
  template <typename T>
  void EventHandler::storeVar(const char *name, T value)
  {
    xAOD::EventInfo *eventInfo = HG::VarHandler::getInstance()->getEventInfoFromStore();
    if (eventInfo == nullptr) {
      HG::fatal("EventHandler cannot access EventInfo");
    }
    eventInfo->auxdecor<T>(name) = value;
  }

  //______________________________________________________________________________
  template <typename T>
  void EventHandler::storeVar(SG::AuxElement::Decorator<T> &decor, T value)
  {
    xAOD::EventInfo *eventInfo = HG::VarHandler::getInstance()->getEventInfoFromStore();
    if (eventInfo == nullptr) {
      HG::fatal("EventHandler cannot access EventInfo");
    }
    decor(*eventInfo) = value;
  }

}

#endif // HGamAnalysisFramework_EventHandler_HPP

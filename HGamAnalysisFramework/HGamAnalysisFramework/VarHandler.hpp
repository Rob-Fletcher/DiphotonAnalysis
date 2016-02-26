#ifndef HGamAnalysisFramework_VarHandler_HPP
#define HGamAnalysisFramework_VarHandler_HPP

#include "xAODEventInfo/EventInfo.h"

namespace HG {

  template <class T>
  VarBase<T>::VarBase(const std::string &name)
  : m_name(name)
  , decVal(name)
  , accVal(name)
  { }

  template <class T>
  VarBase<T>::~VarBase()
  { }

  template <class T>
  void
  VarBase<T>::addToStore(bool truth)
  {
    // Decorate to the TStore
    xAOD::EventInfo *eventInfo = nullptr;
    if (truth) eventInfo = VarHandler::getInstance()->getTruthEventInfoFromStore();
    else       eventInfo = VarHandler::getInstance()->getEventInfoFromStore();

    // First check TEvent, if not there then calculate it
    T value;
    if (checkInEvent(value)) {
      decVal(*eventInfo) = value;
    } else {
      decVal(*eventInfo) = this->calculateValue(truth);
    }

  }

  template <class T>
  T
  VarBase<T>::calculateValue(bool truth)
  { return m_default; }

  template <class T>
  T
  VarBase<T>::getValue(bool truth)
  {
    // Check in TStore and TEvent first.
    // Here TStore is assumed to have the more up to date value
    T value;
    if (checkInStore(value, truth)) return value;
    if (checkInEvent(value, truth)) return value;

    // Otherwise, calculate the value, add to TStore, then return
    xAOD::EventInfo *eventInfo = nullptr;
    if (truth) eventInfo = VarHandler::getInstance()->getTruthEventInfoFromStore();
    else       eventInfo = VarHandler::getInstance()->getEventInfoFromStore();
    decVal(*eventInfo) = this->calculateValue(truth);

    return accVal(*eventInfo);
  }

  template <class T>
  T
  VarBase<T>::operator()()
  { return this->getValue(); }

  template <class T>
  T
  VarBase<T>::truth()
  { return this->getValue(true); }

  template <class T>
  bool
  VarBase<T>::checkInEvent(T &value, bool truth)
  {
    const xAOD::EventInfo *eventInfo = nullptr;
    if (truth) eventInfo = VarHandler::getInstance()->getTruthEventInfoFromEvent();
    else       eventInfo = VarHandler::getInstance()->getEventInfoFromEvent();

    if (eventInfo == nullptr)
      return false;

    if (accVal.isAvailable(*eventInfo)) {
      value = accVal(*eventInfo);
      return true;
    }

    return false;
  }

  template <class T>
  bool
  VarBase<T>::checkInStore(T &value, bool truth)
  {
    xAOD::EventInfo *eventInfo = nullptr;
    if (truth) eventInfo = VarHandler::getInstance()->getTruthEventInfoFromStore(false);
    else       eventInfo = VarHandler::getInstance()->getEventInfoFromStore(false);

    if (eventInfo == nullptr)
      return false;

    if (accVal.isAvailable(*eventInfo)) {
      value = accVal(*eventInfo);
      return true;
    }

    return false;
  }

  template <class T>
  void
  VarBase<T>::setValue(const T &value)
  {
    // Decorate to the TStore
    xAOD::EventInfo *eventInfo = VarHandler::getInstance()->getEventInfoFromStore();
    if (eventInfo == nullptr) return;

    decVal(*eventInfo) = value;
  }

  template <class T>
  void
  VarBase<T>::setTruthValue(const T &value)
  {
    // Decorate to the TStore
    xAOD::EventInfo *eventInfo = VarHandler::getInstance()->getTruthEventInfoFromStore();
    if (eventInfo == nullptr) return;

    decVal(*eventInfo) = value;
  }

  template <class T>
  bool
  VarBase<T>::exists()
  {
    T value;
    if (checkInStore(value)) return true;
    if (checkInEvent(value)) return true;

    return false;
  }

}

#endif // HGamAnalysisFramework_VarHandler_HPP

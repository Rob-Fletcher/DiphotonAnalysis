#include "HGamAnalysisFramework/Config.h"
#include "HGamAnalysisFramework/HgammaUtils.h"

#include "TSystem.h"
#include "THashList.h"

//! Constructor

ClassImp(HG::Config)

namespace HG {

  Config::Config()
  : m_env("env")
  {
    // Must have no pointer initialization, for CINT
    m_env.IgnoreDuplicates(true);
  }
  
  Config::Config(const Config &config)
  : Config()
  {
    config.m_env.Copy(m_env);
    copyTable(config.m_env);
  }
  
  Config::Config(TString fileName)
  : Config()
  {
    addFile(fileName);
  }

  /*
  Config::Config(TEnv *env) : Config() {
    m_env->Copy(env);
  }
   */
  
  void Config::ensureDefined(TString key) {
    if (!isDefined(key)) fatal("HG::Config no value found for "+key);
  }
  
  bool Config::isDefined(TString key) {
    return m_env.Defined(key);
  }
  
  TString Config::getStr(TString key) {
    ensureDefined(key);
    return gSystem->ExpandPathName(m_env.GetValue(key,""));
  }
  
  TString Config::getStr(TString key, TString dflt) {
    return m_env.GetValue(key, dflt);
  }
  
  int Config::getInt(TString key) {
    ensureDefined(key);
    return m_env.GetValue(key,-99);
  }
  
  int Config::getInt(TString key, int dflt) {
    return m_env.GetValue(key,dflt);
  }
  
  bool Config::getBool(TString key, bool dflt) {
    return m_env.GetValue(key,dflt);
  }
  
  bool Config::getBool(TString key) {
    ensureDefined(key); return getBool(key,false);
  }
  
  double Config::getNum(TString key) {
    ensureDefined(key);
    return m_env.GetValue(key,-99.0);
  }
  
  double Config::getNum(TString key, double dflt) {
    return m_env.GetValue(key, dflt);
  }
  
  StrV Config::getStrV(TString key) {
    ensureDefined(key);
    return vectorize(m_env.GetValue(key,"")," \t");
  }
  
  StrV Config::getStrV(TString key, StrV dflt) {
    if (isDefined(key)) {
      return vectorize(m_env.GetValue(key,"")," \t");
    }
    return dflt;
  }

  NumV Config::getNumV(TString key) {
    ensureDefined(key);
    return vectorizeNum(m_env.GetValue(key,"")," \t");
  }
  
  NumV Config::getNumV(TString key, NumV dflt) {
    if (isDefined(key)) {
      return vectorizeNum(m_env.GetValue(key,"")," \t");
    }
    return dflt;
  }

  void Config::printDB() {
    TIter next(m_env.GetTable());
    while (TEnvRec *er = (TEnvRec*) next())
      printf("  %-60s%s\n",Form("%s:",er->GetName()),er->GetValue());
  }
  
  void Config::addFile(TString fileName) {
    TString path(fileName);
    if (!fileExist(path))
      path = PathResolverFindCalibFile(fileName.Data());
    if (path == "")
      fatal("Cannot find settings file "+fileName+"\n  also searched in "+path);

    // settings read in by files should not overwrite values set by setValue()
    TEnv env;
    int status = env.ReadFile(path.Data(),EEnvLevel(0));
    if (status!=0) fatal("Cannot read settings file "+fileName);
    TIter next(env.GetTable());
    while (TEnvRec *er = (TEnvRec*) next())
      if (!isDefined(er->GetName())) setValue(er->GetName(),er->GetValue());
  }
  
}

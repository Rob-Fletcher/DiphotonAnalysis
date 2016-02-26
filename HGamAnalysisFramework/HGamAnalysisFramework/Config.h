#ifndef HGamAnalysisFramework_Config_H
#define HGamAnalysisFramework_Config_H

/********************
 * Config: 
 *   class to read settings from text files
 *   relies on root's TEnv
 *
 * Usage:
 *   Config settings("Hgamma.config");
 *   TString gamContainerName = settings.getStr("PhotonContainer");
 *   TString elContainerName  = settings.getStr("ElectronContainer");
 *   vector<TString> systShifts = settings.getStrV("Systematics");
 *
 */

#include <vector>
#include "TString.h"
#include "TEnv.h"
#include "THashList.h"

// \brief Hgamma namespace
namespace HG {

  //! \brief typedef for a vector of strings (to save some typing)
  typedef std::vector<TString> StrV;
  typedef std::vector<double>  NumV;
  
  /*! \brief Class that handles reading input from a configuration text file
   *        based on root's TEnv format, i.e. key-value pairs
   *  \author Dag Gillberg
   */
  class Config {
  private:
    TEnv m_env; // TEnv objects holding the settings database

  public:
    
    //! \brief Config constructor
    //! \param fileName name of text file with user-specified settings
    Config(TString fileName);
    Config(const Config &config);
    Config();
    virtual ~Config() { }

    Config &operator=(const Config &rhs); // assignment operator
    
    //! \brief Access a string value from the config database. Exception thrown if no entry exist.
    TString getStr(TString key);
    
    //! \brief Access a string value from the config database. Default value used if no entry exist.
    TString getStr(TString key, TString dflt);

    //! \brief Access a vector of strings from the config database
    StrV  getStrV(TString key);
    StrV  getStrV(TString key, StrV dflt);

    //! \brief Access an integer from the config database. Exception thrown if no entry exist.
    int getInt(TString key, int dflt);

    //! \brief Access an integer from the config database. Default value used if no entry exist.
    int getInt(TString key);

    //! \brief Access a boolean from the config database. Exception thrown if no entry exist.
    bool getBool(TString key);

    //! \brief Access a boolean from the config database. Default value used if no entry exist
    bool getBool(TString key, bool dflt);

    //! \brief Access a real number from the config database
    double getNum(TString key);
    double getNum(TString key, double dflt);

    //! \brief Access a vector of integers from the config database
    // std::vector<int>    getIntV(TString key);
    
    //! \brief Access a vector of doubles from the config database
    NumV getNumV(TString key);
    NumV getNumV(TString key, NumV dflt);
    
    //! \brief returns true if the key is defined
    bool isDefined(TString key);
    
    //! \brief Add more user specified settings to the
    void addFile(TString fileName);

    //! \brief Set value
    inline void setValue(TString key, TString value) { m_env.SetValue(key,value); };
    
    //! \brief accessor to the TEnv database
    inline const TEnv* getDB() { return &m_env; }
    
    //! \brief prints the TEnv database to screen
    void printDB();
    
    
  private:
    
    // //! \brief method to abort program with error message
    //    void fatal(TString msg);

    //! \brief ensures that there is a value in the database assocated with key
    //!        if not, abort with error message
    void ensureDefined(TString key);
    inline void copyTable(const TEnv &env);

    ClassDef(Config, 1);
  };

  inline
  Config& Config::operator=(const Config &rhs)
  {
     rhs.m_env.Copy(m_env);
     copyTable(rhs.m_env);
     return *this;
  }

  inline
  void Config::copyTable(const TEnv &env)
  {
    m_env.GetTable()->Delete();
    THashList* hl=0;
    if ( (hl = env.GetTable()) ){
      TIter next(hl);
      TEnvRec* env_rec = 0;
      while ( (env_rec = (TEnvRec* )next.Next()) ){
        m_env.SetValue(env_rec->GetName(),env_rec->GetValue());
      }
    }
  }

}

#endif

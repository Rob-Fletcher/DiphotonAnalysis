#include "HGamAnalysisFramework/HgammaUtils.h"
#include <TObjString.h>
#include <TObjArray.h>

#include "TFile.h"
#include "TH1.h"

// See header file for documentation

namespace HG {

  void fatal(TString msg) {
    printf("\nFATAL\n  %s\n\n",msg.Data());
    abort();
  }
  
  StrV vectorize(TString str, TString sep) {
    StrV result;
    TObjArray *strings = str.Tokenize(sep.Data());
    if (strings->GetEntries()==0) { delete strings; return result; }
    TIter istr(strings);
    while (TObjString* os=(TObjString*)istr()) {
      // the number sign and everything after is treated as a comment
      if (os->GetString()[0]=='#') break;
      result.push_back(os->GetString());
    }
    delete strings;
    return result;
  }
  
  // convert a text line containing a list of numbers to a vector<double>
  NumV vectorizeNum(TString str, TString sep) {
    NumV result; StrV vecS = vectorize(str,sep);
    for (uint i=0;i<vecS.size();++i)
      result.push_back(atof(vecS[i]));
    return result;
  }
  
  // checks if a given file or directory exist
  bool fileExist(TString fn) {
    return !(gSystem->AccessPathName(fn.Data()));
  }

  
  TString fourVecAsText(const TLorentzVector &p4) {
    return TString::Format("(pT,y,phi,m) = (%6.1f GeV,%6.3f,%6.3f,%5.1f GeV )",
                           p4.Pt()*invGeV,p4.Rapidity(),p4.Phi(),p4.M()*invGeV);
  }

  TString fourVecAsText(const xAOD::IParticle *p) {
    return fourVecAsText(p->p4());
  }

  TH1* getHistogramFromFile(TString hname, TString fname)
  {
    fname = PathResolverFindCalibFile(fname.Data());
    TFile *file = TFile::Open(fname.Data(), "READ");

    if (file == nullptr) {
      std::cout << "HgammaUtils::getHistogramFromFile() : Couldn't open file "
                << fname.Data() << ", returning nullptr." << std::endl;
      return nullptr;
    }

    TH1 *temp = (TH1*)file->Get(hname.Data());

    if (temp == nullptr) {
      std::cout << "HgammaUtils::getHistogramFromFile() : Couldn't find histogram "
                << hname.Data() << " in file "
                << fname.Data() << ", returning nullptr." << std::endl;
      return nullptr;
    }

    bool status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    hname = "cloned_" + hname;
    TH1 *hist = (TH1*)temp->Clone(hname.Data());
    SafeDelete(file);
    TH1::AddDirectory(status);

    return hist;
  }
  
  /*
  const char* fourVecAsCharStar(const xAOD::IParticle *p) {
    return fourVecAsText(p).Data();
  }
  */
}

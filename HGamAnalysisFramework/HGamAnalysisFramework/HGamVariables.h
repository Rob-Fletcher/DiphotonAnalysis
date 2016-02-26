#ifndef HGamAnalysisFramework_HGamVars_H
#define HGamAnalysisFramework_HGamVars_H

#include "HGamAnalysisFramework/VarHandler.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

namespace HG {

  //____________________________________________________________________________
  class pT_h1 : public VarBase<float> {
  public:
    pT_h1() : VarBase("pT_h1") { m_default = -99; }
    ~pT_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
        return m_default;
      const xAOD::IParticleContainer *higgs = HG::VarHandler::getInstance()->getHiggsBosons();
      if (higgs->size() < 1)
        return m_default;
      return (*higgs)[0]->pt();
    }
  };

  //____________________________________________________________________________
  class pT_h2 : public VarBase<float> {
  public:
    pT_h2() : VarBase("pT_h2") { m_default = -99; }
    ~pT_h2() { }

    float calculateValue(bool truth)
    {
      if (not truth)
        return m_default;
      const xAOD::IParticleContainer *higgs = HG::VarHandler::getInstance()->getHiggsBosons();
      if (higgs->size() < 2)
        return m_default;
      return (*higgs)[1]->pt();
    }
  };

  //____________________________________________________________________________
  class y_h1 : public VarBase<float> {
  public:
    y_h1() : VarBase("y_h1") { m_default = -99; }
    ~y_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
        return m_default;
      const xAOD::IParticleContainer *higgs = HG::VarHandler::getInstance()->getHiggsBosons();
      if (higgs->size() < 1)
        return m_default;
      return (*higgs)[0]->rapidity();
    }
  };

  //____________________________________________________________________________
  class y_h2 : public VarBase<float> {
  public:
    y_h2() : VarBase("y_h2") { m_default = -99; }
    ~y_h2() { }

    float calculateValue(bool truth)
    {
      if (not truth)
        return m_default;
      const xAOD::IParticleContainer *higgs = HG::VarHandler::getInstance()->getHiggsBosons();
      if (higgs->size() < 2)
        return m_default;
      return (*higgs)[1]->rapidity();
    }
  };

  //____________________________________________________________________________
  class m_h1 : public VarBase<float> {
  public:
    m_h1() : VarBase("m_h1") { m_default = -99; }
    ~m_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
        return m_default;
      const xAOD::IParticleContainer *higgs = HG::VarHandler::getInstance()->getHiggsBosons();
      if (higgs->size() < 1)
        return m_default;
      return (*higgs)[0]->m();
    }
  };

  //____________________________________________________________________________
  class m_h2 : public VarBase<float> {
  public:
    m_h2() : VarBase("m_h2") { m_default = -99; }
    ~m_h2() { }

    float calculateValue(bool truth)
    {
      if (not truth)
        return m_default;
      const xAOD::IParticleContainer *higgs = HG::VarHandler::getInstance()->getHiggsBosons();
      if (higgs->size() < 2)
        return m_default;
      return (*higgs)[1]->m();
    }
  };

  //____________________________________________________________________________
  class yAbs_yy : public VarBase<float> {
  public:
    yAbs_yy() : VarBase("yAbs_yy") { m_default = -99; }
    ~yAbs_yy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      return fabs(((*photons)[0]->p4() + (*photons)[1]->p4()).Rapidity());
    }
  };

  //____________________________________________________________________________
  class pTt_yy : public VarBase<float> {
  public:
    pTt_yy() : VarBase("pTt_yy") { m_default = -99; }
    ~pTt_yy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      TLorentzVector g1 = (*photons)[0]->p4(), g2 = (*photons)[1]->p4();
      return fabs(g1.Px()*g2.Py() - g2.Px()*g1.Py())/(g1-g2).Pt()*2.0;
    }
  };

  //____________________________________________________________________________
  class m_yy : public VarBase<float> {
  public:
    m_yy() : VarBase("m_yy") { m_default = -99; }
    ~m_yy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      return ((*photons)[0]->p4() + (*photons)[1]->p4()).M();
    }
  };

  //____________________________________________________________________________
  class pT_y1 : public VarBase<float> {
  public:
    pT_y1() : VarBase("pT_y1") { m_default = -99; }
    ~pT_y1() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 1)
        return m_default;
      return (*photons)[0]->pt();
    }
  };

  //____________________________________________________________________________
  class pT_y2 : public VarBase<float> {
  public:
    pT_y2() : VarBase("pT_y2") { m_default = -99; }
    ~pT_y2() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      return (*photons)[1]->pt();
    }
  };

  //____________________________________________________________________________
  class pT_yy : public VarBase<float> {
  public:
    pT_yy() : VarBase("pT_yy") { m_default = -99; }
    ~pT_yy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      return ((*photons)[0]->p4() + (*photons)[1]->p4()).Pt();
    }
  };

  //____________________________________________________________________________
  class cosTS_yy : public VarBase<float> {
  public:
    cosTS_yy() : VarBase("cosTS_yy") { m_default = -99; }
    ~cosTS_yy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (gams->size() < 2)
        return m_default;
      const TLorentzVector &y1 = (*gams)[0]->p4(), &y2 = (*gams)[1]->p4();
      return std::abs( ( (y1.E()+y1.Pz())*(y2.E()-y2.Pz()) - (y1.E()-y1.Pz())*(y2.E()+y2.Pz()) )
                      / ((y1+y2).M()*sqrt(pow((y1+y2).M(),2)+pow((y1+y2).Pt(),2)) ) );
    }
  };

  //____________________________________________________________________________
  class phiStar_yy : public VarBase<float> {
  public:
    phiStar_yy() : VarBase("phiStar_yy") { m_default = -99; }
    ~phiStar_yy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      const TLorentzVector &g1 = (*photons)[0]->p4(), &g2 = (*photons)[1]->p4();
      return tan((TMath::Pi() - fabs(g1.DeltaPhi(g2)))/2.0) *
             sqrt(1.0 - pow(tanh((g1.Eta() - g2.Eta())/2.0), 2.0)); // FIXME?
    }
  };

  //____________________________________________________________________________
  class Dphi_y_y : public VarBase<float> {
  public:
    Dphi_y_y() : VarBase("Dphi_y_y") { m_default = -99; }
    ~Dphi_y_y() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      TLorentzVector g1 = (*photons)[0]->p4(), g2 = (*photons)[1]->p4();
      return fabs(g1.DeltaPhi(g2));
    }
  };

  //____________________________________________________________________________
  class Dy_y_y : public VarBase<float> {
  public:
    Dy_y_y() : VarBase("Dy_y_y") { m_default = -99; }
    ~Dy_y_y() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      if (photons->size() < 2)
        return m_default;
      return fabs((*photons)[0]->rapidity() - (*photons)[1]->rapidity());
    }
  };

  //____________________________________________________________________________
  class Njets : public VarBase<float> {
  public:
    Njets() : VarBase("Njets") { m_default = -99; }
    ~Njets() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      return jets->size() > 3 ? 3 : jets->size();
    }
  };

  //____________________________________________________________________________
  class pT_j1 : public VarBase<float> {
  public:
    pT_j1() : VarBase("pT_j1") { m_default = -99; }
    ~pT_j1() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 1)
        return m_default;
      return (*jets)[0]->pt();
    }
  };

  //____________________________________________________________________________
  class yAbs_j1 : public VarBase<float> {
  public:
    yAbs_j1() : VarBase("yAbs_j1") { m_default = -99; }
    ~yAbs_j1() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 1)
        return m_default;
      return fabs((*jets)[0]->rapidity());
    }
  };

  //____________________________________________________________________________
  class pT_yyj : public VarBase<float> {
  public:
    pT_yyj() : VarBase("pT_yyj") { m_default = -99; }
    ~pT_yyj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (photons->size() < 2)
        return m_default;
      if (jets->size() < 1)
        return m_default;
      return ((*photons)[0]->p4() + (*photons)[1]->p4() + (*jets)[0]->p4()).Pt();
    }
  };

  //____________________________________________________________________________
  class m_yyj : public VarBase<float> {
  public:
    m_yyj() : VarBase("m_yyj") { m_default = -99; }
    ~m_yyj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *photons = HG::VarHandler::getInstance()->getPhotons(truth);
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (photons->size() < 2)
        return m_default;
      if (jets->size() < 1)
        return m_default;
      return ((*photons)[0]->p4() + (*photons)[1]->p4() + (*jets)[0]->p4()).M();
    }
  };

  //____________________________________________________________________________
  class pT_j2 : public VarBase<float> {
  public:
    pT_j2() : VarBase("pT_j2") { m_default = -99; }
    ~pT_j2() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 2)
        return m_default;
      return (*jets)[1]->pt();
    }
  };

  //____________________________________________________________________________
  class yAbs_j2 : public VarBase<float> {
  public:
    yAbs_j2() : VarBase("yAbs_j2") { m_default = -99; }
    ~yAbs_j2() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 2)
        return m_default;
      return fabs((*jets)[1]->rapidity());
    }
  };

  //____________________________________________________________________________
  class m_jj : public VarBase<float> {
  public:
    m_jj() : VarBase("m_jj") { m_default = -99; }
    ~m_jj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 2)
        return m_default;
      return ((*jets)[0]->p4() + (*jets)[1]->p4()).M();
    }
  };

  //____________________________________________________________________________
  class Dy_j_j : public VarBase<float> {
  public:
    Dy_j_j() : VarBase("Dy_j_j") { m_default = -99; }
    ~Dy_j_j() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 2)
        return m_default;
      return fabs((*jets)[0]->rapidity() - (*jets)[1]->rapidity());
    }
  };

  //____________________________________________________________________________
  class Dy_yy_jj : public VarBase<float> {
  public:
    Dy_yy_jj() : VarBase("Dy_yy_jj") { m_default = -99; }
    ~Dy_yy_jj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *ps = HG::VarHandler::getInstance()->getPhotons(truth);
      if (ps->size() < 2)
        return m_default;
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);
      if (js->size() < 2)
        return m_default;
      return fabs(((*ps)[0]->p4() + (*ps)[1]->p4()).Rapidity() - ((*js)[0]->p4() + (*js)[1]->p4()).Rapidity());
    }
  };

  //____________________________________________________________________________
  class Dphi_j_j : public VarBase<float> {
  public:
    Dphi_j_j() : VarBase("Dphi_j_j") { m_default = -99; }
    ~Dphi_j_j() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (jets->size() < 2)
        return m_default;
      return (*jets)[0]->p4().DeltaPhi((*jets)[1]->p4());
    }
  };

  //____________________________________________________________________________
  class Dphi_yy_jj : public VarBase<float> {
  public:
    Dphi_yy_jj() : VarBase("Dphi_yy_jj") { m_default = -99; }
    ~Dphi_yy_jj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *ps = HG::VarHandler::getInstance()->getPhotons(truth);
      if (ps->size() < 2)
        return m_default;
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);
      if (js->size() < 2)
        return m_default;
      return fabs(((*ps)[0]->p4() + (*ps)[1]->p4()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
    }
  };

  //____________________________________________________________________________
  class m_yyjj : public VarBase<float> {
  public:
    m_yyjj() : VarBase("m_yyjj") { m_default = -99; }
    ~m_yyjj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *ps = HG::VarHandler::getInstance()->getPhotons(truth);
      if (ps->size() < 2)
        return m_default;
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);
      if (js->size() < 2)
        return m_default;
      return ((*ps)[0]->p4() + (*ps)[1]->p4() + (*js)[0]->p4() + (*js)[1]->p4()).M();
    }
  };

  //____________________________________________________________________________
  class pT_yyjj : public VarBase<float> {
  public:
    pT_yyjj() : VarBase("pT_yyjj") { m_default = -99; }
    ~pT_yyjj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *ps = HG::VarHandler::getInstance()->getPhotons(truth);
      if (ps->size() < 2)
        return m_default;
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);
      if (js->size() < 2)
        return m_default;
      return ((*ps)[0]->p4() + (*ps)[1]->p4() + (*js)[0]->p4() + (*js)[1]->p4()).Pt();
    }
  };

  //____________________________________________________________________________
  class Zepp : public VarBase<float> {
  public:
    Zepp() : VarBase("Zepp") { m_default = -99; }
    ~Zepp() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *ps = HG::VarHandler::getInstance()->getPhotons(truth);
      if (ps->size() < 2)
        return m_default;
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);
      if (js->size() < 2)
        return m_default;
      return ((*ps)[0]->p4() + (*ps)[1]->p4()).Eta() - ((*js)[0]->eta() + (*js)[1]->eta())/2.0;
    }
  };

  //____________________________________________________________________________
  class isPassedBasic : public VarBase<char> {
  public:
    isPassedBasic() : VarBase("isPassedBasic") { m_default = false; }
    ~isPassedBasic() { }
  };

  class isDalitzEvent : public VarBase<char> {
  public:
    isDalitzEvent() : VarBase("isDalitz") { m_default = false; }
    ~isDalitzEvent() { }
  };

  
  //____________________________________________________________________________
  class isPassed : public VarBase<char> {
  public:
    isPassed() : VarBase("isPassed") { m_default = false; }
    ~isPassed() { }
  };

  //____________________________________________________________________________
  class isPassedEventClean : public VarBase<char> {
  public:
    isPassedEventClean() : VarBase("isPassedEventClean") { m_default = false; }
    ~isPassedEventClean() { }
  };

  //____________________________________________________________________________
  class isFiducial : public VarBase<char> {
  public:
    isFiducial() : VarBase("isFiducial") { m_default = false; }
    ~isFiducial() { }
  };

  //____________________________________________________________________________
  class isFiducialKinOnly : public VarBase<char> {
  public:
    isFiducialKinOnly() : VarBase("isFiducialKinOnly") { m_default = false; }
    ~isFiducialKinOnly() { }
  };

  //____________________________________________________________________________
  class cutFlow : public VarBase<int> {
  public:
    cutFlow() : VarBase("cutFlow") { m_default = -99; }
    ~cutFlow() { }
  };

  //____________________________________________________________________________
  class weightInitial : public VarBase<float> {
  public:
    weightInitial() : VarBase("weightInitial") { m_default = 1.0; }
    ~weightInitial() { }
  };

  //____________________________________________________________________________
  class weight : public VarBase<float> {
  public:
    weight() : VarBase("weight") { m_default = 1.0; }
    ~weight() { }
  };

  //____________________________________________________________________________
  class weightCategory : public VarBase<float> {
  public:
    weightCategory() : VarBase("weightCategory") { m_default = 1.0; }
    ~weightCategory() { }
  };

  //____________________________________________________________________________
  class category : public VarBase<int> {
  public:
    category() : VarBase("category") { m_default = -99; }
    ~category() { }
  };

  //____________________________________________________________________________
  class numberOfPrimaryVertices : public VarBase<int> {
  public:
    numberOfPrimaryVertices() : VarBase("numberOfPrimaryVertices") { m_default = -99; }
    ~numberOfPrimaryVertices() { }
  };

  //____________________________________________________________________________
  class selectedVertexZ : public VarBase<float> {
  public:
    selectedVertexZ() : VarBase("selectedVertexZ") { m_default = -999; }
    ~selectedVertexZ() { }
  };

  //____________________________________________________________________________
  class hardestVertexZ : public VarBase<float> {
  public:
    hardestVertexZ() : VarBase("hardestVertexZ") { m_default = -999; }
    ~hardestVertexZ() { }
  };

  //____________________________________________________________________________
  class eventShapeDensity : public VarBase<float> {
  public:
    eventShapeDensity() : VarBase("eventShapeDensity") { m_default = -99; }
    ~eventShapeDensity() { }
  };

  //____________________________________________________________________________
  class mu : public VarBase<float> {
  public:
    mu() : VarBase("mu") { m_default = -99; }
    ~mu() { }
  };

  //____________________________________________________________________________
  class truthVertexZ : public VarBase<float> {
  public:
    truthVertexZ() : VarBase("truthVertexZ") { m_default = -999; }
    ~truthVertexZ() { }
  };

}

namespace var {
  extern HG::pT_h1 pT_h1;
  extern HG::pT_h2 pT_h2;
  extern HG::y_h1 y_h1;
  extern HG::y_h2 y_h2;
  extern HG::m_h1 m_h1;
  extern HG::m_h2 m_h2;
  extern HG::yAbs_yy yAbs_yy;
  extern HG::pTt_yy pTt_yy;
  extern HG::m_yy m_yy;
  extern HG::pT_yy pT_yy;
  extern HG::pT_y1 pT_y1;
  extern HG::pT_y2 pT_y2;
  extern HG::cosTS_yy cosTS_yy;
  extern HG::phiStar_yy phiStar_yy;
  extern HG::Dphi_y_y Dphi_y_y;
  extern HG::Dy_y_y Dy_y_y;
  extern HG::Njets Njets;
  extern HG::pT_j1 pT_j1;
  extern HG::pT_j2 pT_j2;
  extern HG::yAbs_j1 yAbs_j1;
  extern HG::yAbs_j2 yAbs_j2;
  extern HG::pT_yyj pT_yyj;
  extern HG::m_yyj m_yyj;
  extern HG::m_jj m_jj;
  extern HG::Dy_j_j Dy_j_j;
  extern HG::Dy_yy_jj Dy_yy_jj;
  extern HG::Dphi_j_j Dphi_j_j;
  extern HG::Dphi_yy_jj Dphi_yy_jj;
  extern HG::m_yyjj m_yyjj;
  extern HG::pT_yyjj pT_yyjj;
  extern HG::Zepp Zepp;

  extern HG::isPassedBasic isPassedBasic;
  extern HG::isDalitzEvent isDalitzEvent;
  extern HG::isPassed isPassed;
  extern HG::isPassedEventClean isPassedEventClean;
  extern HG::isFiducial isFiducial;
  extern HG::isFiducialKinOnly isFiducialKinOnly;
  extern HG::cutFlow cutFlow;
  extern HG::weightInitial weightInitial;
  extern HG::weight weight;
  extern HG::weightCategory weightCategory;
  extern HG::category category;
  extern HG::numberOfPrimaryVertices numberOfPrimaryVertices;
  extern HG::selectedVertexZ selectedVertexZ;
  extern HG::hardestVertexZ hardestVertexZ;
  extern HG::truthVertexZ truthVertexZ;
  extern HG::eventShapeDensity eventShapeDensity;
  extern HG::mu mu;
}

#endif // HGamAnalysisFramework_HGamVars_H

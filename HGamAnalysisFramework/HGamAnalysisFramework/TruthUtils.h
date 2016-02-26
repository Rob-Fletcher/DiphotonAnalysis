#ifndef HGamAnalysisFramework_TruthUtils
#define HGamAnalysisFramework_TruthUtils

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "TruthUtils/PIDUtils.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "AthContainers/ConstDataVector.h"
#include <TLorentzVector.h>

#include "HGamAnalysisFramework/HgammaUtils.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "AthContainers/DataVector.h"

//! \brief Hgamma namespace
namespace HG {
  
  typedef ConstDataVector<xAOD::TruthParticleContainer> TruthPtcls;
  typedef ConstDataVector<xAOD::JetContainer> TruthJets;
  
  //! \name Helper methods to deal with truth particles.
  //!       These method should all follow the truth harmonization recommendation:
  //!       https://cds.cern.ch/record/1700541/
  //! \author Stephen Menary
  //! \author Michaela Queitsch-Maitland
  //! \author Dag Gillberg
  //@{

  /// print details about the truth particle to the screen
  void printTruthPtcl(const xAOD::TruthParticle *ptcl, TString comment="", 
		      int NchildDepth=0, int NparentDepth=0, int currentDepth=0);
  
  //! /brief returns true if this is a stable "generator "particle according to the ATALS defintion.
  //!        status shoudl be =1 (stable) and barcode < 200k (non-GEANT partcle)
  bool isStable(const xAOD::TruthParticle *ptcl);
  
  //! \brief rejects events that decay according to H -> y y* -> y f fbar
  //!        which is turned off by default in Pythia8
  //!        (is turned off by: )
  bool isDalitz(const xAOD::TruthParticleContainer *truthPtcls);
  
  //! /brief returns the Et within DR of the passed particle
  float getTruthIsolation(const xAOD::TruthParticle *ptcl, const xAOD::TruthParticleContainer *truthPtcls, double dr = 0.4, std::vector<int> ignorePdgId = {});

  //! /brief returns the sum of four-momenta of all stable particles
  TLorentzVector getStableParticle4VectorSum(const xAOD::TruthParticleContainer * truthPtcls);

  //! /brief returns false if the particle originates from a hadron
  bool notFromHadron(const xAOD::TruthParticle *ptcl);

  //! /brief returns true if the particle is a final Higgs boson
  bool isFinalHiggs(const xAOD::TruthParticle *part);

  //! /brief returns pointer to the first final Higgs boson found in the event
  TruthPtcls getFinalHiggsBosons(const xAOD::TruthParticleContainer *truthParticles);

  //! /brief returns true if the particle originates from a Higgs boson (recursively)
  bool isFromHiggs(const xAOD::TruthParticle *ptcl);
  
  //! /brief returns true if the particle originates from a Bhadron (recursively)
  bool isFromBhadron(const xAOD::TruthParticle *ptcl);
  
  //! /brief whether the paritlce is a stable truth photon, not decaying from a hadron
  bool isGoodTruthPhoton(const xAOD::TruthParticle *ptcl);

  //! /brief whether the paritlce is a stable truth electron, not decaying from a hadron
  bool isGoodTruthElectron(const xAOD::TruthParticle *ptcl);
  
  //! /brief whether the paritlce is a stable truth electron, not decaying from a hadron
  bool isGoodTruthMuon(const xAOD::TruthParticle *ptcl);
  
  //! /brief returns a collection of good truth photons
  //  xAOD::TruthParticleContainer getGoodTruthPhotons(const xAOD::TruthParticleContainer * truthPtcls );
  ::std::vector<const xAOD::TruthParticle*> getGoodTruthPhotonsOld(const xAOD::TruthParticleContainer *truthPtcls);

  //! /brief returns a collection of good truth photons
  TruthPtcls getGoodTruthPhotons( const xAOD::TruthParticleContainer * truthPtcls );

  //! /brief returns a collection of good truth photons
  TruthPtcls getHiggsDecayProducts( const xAOD::TruthParticleContainer * truthPtcls );
  
  //! /brief returns all stable electrons that do not originate from hadrons
  TruthPtcls getGoodTruthElectrons( const xAOD::TruthParticleContainer * truthPtcls );

  //! /brief returns all stable electrons that do not originate from hadrons
  TruthPtcls getGoodTruthMuons( const xAOD::TruthParticleContainer * truthPtcls );

  //! /brief follows ancestry and returns all stalbe particles pointing back to ptcls
  TruthPtcls getStableDecayProducts( const xAOD::TruthParticle *ptcl );
  
  //! /brief returns all stable particles that either are hadrons or originate from the decay of a hadron
  TruthPtcls getHadronsAndTheirDecay( const xAOD::TruthParticleContainer * truthPtcls );

  //! /brief returns a list of Bhadrons.
  //!        If pTcut is specified (positive), then only Bhadrons above this value are selected
  TruthPtcls getBHadrons( const xAOD::TruthParticleContainer *truthPtcls, double pTcut = -1.0 );
  
  //! /brief returns a list of Dhadrons.
  //!        If pTcut is specified (positive), then only Dhadrons above this value are selected
  TruthPtcls getDHadrons( const xAOD::TruthParticleContainer *truthPtcls, double pTcut = -1.0 );

  //! /brief returns the photons from the Higgs
  TruthPtcls getPhotonsFromHiggs( const xAOD::TruthParticleContainer *truthPtcls );

  //! /brief returns list of truth muons originating from a Bhadron
  TruthPtcls getMuonsFromBs( const xAOD::TruthParticleContainer *truthPtcls );

  // Full reconstruction of all TruthParticles
  struct TruthParticleStruct {
    // all stable particles are classifed as, photons, muons, electrons or "hadrons"
    // hadrons includes taus and all particles produced in hadronic decays
    // i.e. photons from pi0 are in the "hadrons" collection
    TruthPtcls photons;
    TruthPtcls muons;     // TO-DO: apply dressing
    TruthPtcls electrons; // TO-DO: apply dressing
    TruthPtcls hadrons;
    
    // The above containters makes "a complete set"
    // the below conntainers contain additional information
    TruthPtcls photonsFromHiggs;
    TruthPtcls HiggsDecay;
    TruthPtcls Bhadrons;
    TruthPtcls Dhadrons;
    TruthPtcls muonsFromBs;
    
    // all truth jets - after photon, electron overlap removal
    TruthJets jets;

    // all truth jets split into b-, c- or light jets
    // this is based on hadrons
    TruthJets bJets;
    TruthJets cJets;
    TruthJets lightJets;

    // light jets further split into light-quark, gluon, and unknown-flavour jets
    // WARNING: this splitting is "unphysical" - it relies on off-shell partons
    TruthJets LQjets;
    TruthJets gluonJets;
    TruthJets UFjets;
  };

  TruthParticleStruct identifyTruthParticles(xAOD::TEvent *event,
                                             double jet_pTcut = 10*HG::GeV);

  TruthParticleStruct identifyTruthParticles(const xAOD::TruthParticleContainer *truthPtcls,
                                             const xAOD::JetContainer *truthJets,
                                             double jet_pTcut = 10*HG::GeV);
  void printTruthParticles( const TruthParticleStruct &truthPtcls );

  void removeTruthOverlap(DataVector<xAOD::IParticle> &photons  ,
                          DataVector<xAOD::IParticle> &electrons,
                          DataVector<xAOD::IParticle> &muons    ,
                          DataVector<xAOD::IParticle> &jets     ,
                          double jet_pTcut = 10*HG::GeV         );

  template<class T, class contT>
  contT
  getContainer(const TruthPtcls *particles, TString name,  xAOD::TStore *store);
  
  //@}
  
}

namespace HG {

  template<class T, class contT>
  contT
  getContainer(const TruthPtcls *particles, TString name, xAOD::TStore *store)
  {
    contT *truths = new contT();
    xAOD::AuxContainerBase *truthsAux = new xAOD::AuxContainerBase();
    truths->setStore(truthsAux);

    for (auto part: *particles) {
      T *truth = new T();
      truths->push_back(truth);
      truth->setP4(part->pt(), part->eta(), part->phi(), part->m());
    }

    if (!store->record(truths, name.Data()))
      std::cout << "Couldn't store " + name << std::endl;

    name += "Aux";
    if (!store->record(truthsAux, name.Data()))
      std::cout << "Couldn't store " + name << std::endl;

    contT container(truths->begin(), truths->end(), SG::VIEW_ELEMENTS);
    return container;
  }

}
#endif // HGamAnalysisFramework_TruthUtils

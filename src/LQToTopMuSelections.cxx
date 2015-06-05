#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;


DijetSelection::DijetSelection(float dphi_min_, float third_frac_max_): dphi_min(dphi_min_), third_frac_max(third_frac_max_){}
bool DijetSelection::passes(const Event & event){
    assert(event.jets); // if this fails, it probably means jets are not read in
    if(event.jets->size() < 2) return false;
    const auto & jet0 = event.jets->at(0);
    const auto & jet1 = event.jets->at(1);
    auto dphi = deltaPhi(jet0, jet1);
    if(dphi < dphi_min) return false;
    if(event.jets->size() == 2) return true;
    const auto & jet2 = event.jets->at(2);
    auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
    return third_jet_frac < third_frac_max;
}

HtSelection::HtSelection(double ht_min_, double ht_max_):ht_min(ht_min_), ht_max(ht_max_){}
bool HtSelection::passes(const Event & event){
  auto met = event.met->pt();

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;

  for(const auto & jet : *event.jets){
    ht += jet.pt();
    }
  //Bedeutung der for-Schleife
  /*const auto jets = event.jets;
  for(unsigned int i=0; i<jets.size();i++){
    Jet jet=jets[i];
    ht +=jet.pt();
    }*/
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
  }

  ht = ht_lep + ht_jets + met;

  bool pass = false;
  pass = ht > ht_min && (ht_max < 0 || ht < ht_max);
  return pass;
}

InvMass2MuVeto::InvMass2MuVeto(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2MuVeto::passes(const Event & event){

  bool pass = true;
  int Nmuons = event.muons->size();
  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j > i){
	M_mumu = (muons[i] + muons[j]).M();
	if(M_mumu > m_min && M_mumu < m_max){pass = false;}
      }
    }
  }
  return pass;
}

PtLeadingMuonSelection::PtLeadingMuonSelection(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool PtLeadingMuonSelection::passes(const Event & event){

  bool pass = true;
  double pt_leadingmu = event.muons->at(0).pt();
  pass = pt_leadingmu >= pt_min && (pt_leadingmu <= pt_max || pt_max < 0);
  return pass;
}

Pt2ndMuonSelection::Pt2ndMuonSelection(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool Pt2ndMuonSelection::passes(const Event & event){

  bool pass = true;
  double pt_2ndmu = event.muons->at(1).pt();
  pass = pt_2ndmu >= pt_min && (pt_2ndmu <= pt_max || pt_max < 0);
  return pass;
}

PtLeadingJetSelection::PtLeadingJetSelection(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool PtLeadingJetSelection::passes(const Event & event){

  bool pass = true;
  double pt_leadingjet = event.jets->at(0).pt();
  pass = pt_leadingjet >= pt_min && (pt_leadingjet <= pt_max || pt_max < 0);
  return pass;
}

Pt2ndJetSelection::Pt2ndJetSelection(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool Pt2ndJetSelection::passes(const Event & event){

  bool pass = true;
  double pt_2ndjet = event.jets->at(1).pt();
  pass = pt_2ndjet >= pt_min && (pt_2ndjet <= pt_max || pt_max < 0);
  return pass;
}

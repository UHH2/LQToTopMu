#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

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
    ht_jets += jet.pt();
    }
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  /*for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
    }*/

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
	if(M_mumu > m_min && M_mumu < m_max){
	  pass = false;
	}
      }
    }
  }
  return pass;
}

InvMass2MuVetoInverted::InvMass2MuVetoInverted(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2MuVetoInverted::passes(const Event & event){

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
	if(M_mumu < m_min || M_mumu > m_max){
	  pass = false;
	}
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

PtRelMuJetSelection::PtRelMuJetSelection(double ptrel_min_, double ptrel_max_) : ptrel_min(ptrel_min_), ptrel_max(ptrel_max_){}
bool PtRelMuJetSelection::passes(const Event & event){
  
  bool pass = true;
  double ptrel = 0;
  for(const auto & muon : *event.muons){
    auto nextjet = nextJet(muon, *event.jets);
    ptrel = pTrel(muon, nextjet);
    if(ptrel >= ptrel_min && (ptrel <= ptrel_max || ptrel_max < 0)){
      pass = true;
    }
    else{
      pass = false;
      return pass; 
    }
  }
  return pass;
}

METSelection::METSelection(double met_min_, double met_max_) : met_min(met_min_), met_max(met_max_){}
bool METSelection::passes(const Event & event){

  auto met = event.met->pt();
  bool pass = true;

  pass = met >= met_min && (met <= met_max || met_max < 0);
  return pass;

}

GenLvlZMuMuSelection::GenLvlZMuMuSelection(){}
bool GenLvlZMuMuSelection::passes(const Event & event){

  bool pass = true;
  for (const auto & genpart : *event.genparticles){
    if (genpart.pdgId() == 23){
      if(abs(genpart.daughter(event.genparticles, 1)->pdgId()) == 13){
	pass = true;
      }
      else{
	pass = false;
      }
    }
  }
  return pass;
}

PtRelMu1JetSelection::PtRelMu1JetSelection(double ptrel_min_, double ptrel_max_) : ptrel_min(ptrel_min_), ptrel_max(ptrel_max_){}
bool PtRelMu1JetSelection::passes(const Event & event){
  
  bool pass = true;
  auto muon = event.muons->at(0);
  auto nextjet = nextJet(muon, *event.jets);
  double ptrel = pTrel(muon, nextjet);

  pass = ptrel > ptrel_min && (ptrel < ptrel_max || ptrel_max < 0);
  return pass;
}

HTJetsSelection::HTJetsSelection(double ht_min_, double ht_max_) : ht_min(ht_min_), ht_max(ht_max_){}
bool HTJetsSelection::passes(const Event & event){
  
  bool pass = true;
  double ht_jets = 0.0;
  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
    }

  pass = ht_jets > ht_min && (ht_jets < ht_max || ht_max < 0);
  return pass;
}

HTLeptSelection::HTLeptSelection(double ht_min_, double ht_max_) : ht_min(ht_min_), ht_max(ht_max_){}
bool HTLeptSelection::passes(const Event & event){
  
  bool pass = true;
  double ht_lept = 0.0;
  for(const auto & ele : *event.electrons){
    ht_lept += ele.pt();
    }
  for(const auto & mu : *event.muons){
    ht_lept += mu.pt();
    }

  pass = ht_lept > ht_min && (ht_lept < ht_max || ht_max < 0);
  return pass;
}

EtaLeadingJetSelection::EtaLeadingJetSelection(double eta_max_) : eta_max(eta_max_){}
bool EtaLeadingJetSelection::passes(const Event & event){
  
  bool pass = true;
  double etajet1 = fabs(event.jets->at(0).eta());

  pass = etajet1 < eta_max || eta_max < 0;
  return pass;
}

InvMassMuEleVeto::InvMassMuEleVeto(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassMuEleVeto::passes(const Event & event){

  bool pass = true;
  int Nmuons = event.muons->size();
  int Nelectrons = event.electrons->size();
  double M_muele;
  LorentzVector muons[Nmuons];
  LorentzVector electrons[Nelectrons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nelectrons; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nelectrons; j++){
      M_muele = (muons[i] + electrons[j]).M();
      if(M_muele > m_min && M_muele < m_max){
	pass = false;
      } 
    }
  }
  return pass;
}



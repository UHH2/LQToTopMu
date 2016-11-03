#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;




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




#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

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



InvMass2MuVeto_Inverted::InvMass2MuVeto_Inverted(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2MuVeto_Inverted::passes(const Event & event){

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
	if(M_mumu <= m_min || M_mumu >= m_max){
	  pass = false;
	}
      }
    }
  }
  return pass;
}



InvMass2EleVeto_Inverted::InvMass2EleVeto_Inverted(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2EleVeto_Inverted::passes(const Event & event){

  int Nelectrons = event.electrons->size();
  double M_eleele;
  LorentzVector electrons[Nelectrons];
  for(int i=0; i<Nelectrons; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nelectrons; i++){
    for(int j=0; j<Nelectrons; j++){
      if(j > i){
	M_eleele = (electrons[i] + electrons[j]).M();
	/*
	  if(M_eleele <= m_min || M_eleele >= m_max){
	  pass = false;
	  }
	*/
	if(M_eleele > m_min && M_eleele < m_max){
	  return true;
	}
      }
    }
  }
  return false;
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



METSelection::METSelection(double MET_min_, double MET_max_):MET_min(MET_min_), MET_max(MET_max_){}
bool METSelection::passes(const Event & event){

  bool pass = true;
  float met = event.met->pt();

  if(met < MET_min || met > MET_max){
    pass = false;
  }
  return pass;
}



dRSelection::dRSelection(double dR_min_, double dR_max_):dR_min(dR_min_), dR_max(dR_max_){}
bool dRSelection::passes(const Event & event) {

  bool pass = false;
  std::vector<Jet>* jets = event.jets;

  int Neles = event.electrons->size();
  double M_eleele;
  LorentzVector electrons[Neles];
  double ele_eta[Neles];
  double ele_phi[Neles];

  for(int g=0; g<Neles; g++){
    electrons[g] = event.electrons->at(g).v4();
    ele_eta[g] = event.electrons->at(g).eta();
    ele_phi[g] = event.electrons->at(g).phi();
  }

  for(unsigned int i=0; i<event.jets->size(); i++){
    double jet_eta = jets->at(i).eta();
    double jet_phi = jets->at(i).phi();
    for(int j=0; j<Neles; j++){
      for(int k=0; k<Neles; k++){
	if(k > j){
	  M_eleele = (electrons[j] + electrons[k]).M();
	  if(71 < M_eleele && M_eleele < 111){
	    for(int l=0; l<Neles; l++){
	      if(electrons[l] != electrons[j] && electrons[l] != electrons[k]){
		double dR_jetelefakeele;
		double deltaetafakeele = jet_eta - ele_eta[l];
		double deltaphifakeele = fabs(jet_phi - ele_phi[l]);
		if(deltaphifakeele > M_PI) deltaphifakeele = 2* M_PI - deltaphifakeele;
		dR_jetelefakeele = sqrt(deltaetafakeele * deltaetafakeele + deltaphifakeele * deltaphifakeele);
		pass = dR_jetelefakeele <= dR_max;
		break;
	      }
	    }
	  }
	}
	if(pass) break;      
      }
      if(pass) break;    
    }
    if(pass) break;  
  }

  return pass;
}

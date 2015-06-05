#pragma once

#include "UHH2/core/include/Particle.h"
#include <map>

/**
 *  @short container class to store the results of the top quark reconstruction
 * 
 * The top quarks reconstruction only applied to semileptonic ttbar events. A
 * LQReconstructionHypothesis then consists of finding the lepton, and assigning the jets of the event
 * to either the leptonically decaying top or the hadronically decaying top. In addition
 * to accessing these information (i.e., which jets are assigned to which side, etc.), each
 * hypothesis can have a number of associated *discriminators*. A discriminator is identified
 * by name and is a floating point value meant to measure how 'good' the hypothesis is according to some criterion;
 * see LQReconstructionHypothesisDiscriminators.h for different criteria and to fill the discriminators.
 */
class LQReconstructionHypothesis {
public:
  LorentzVector toplep_v4() const{return m_toplep_v4;}
  LorentzVector tophad_v4() const{return m_tophad_v4;} 
  LorentzVector neutrino_v4() const{return m_neutrino_v4;} 
  Particle electron() const{return m_electron;}
  LorentzVector electron_v4() const{return m_electron_v4;}
  LorentzVector mu_lep_v4() const{return m_mu_lep_v4;}
  LorentzVector mu_had_v4() const{return m_mu_had_v4;}
  Particle mu_lep() const{return m_mu_lep;}
  Particle mu_had() const{return m_mu_had;}

  std::vector<Particle> tophad_jets() const{return m_tophad_jets;}
  std::vector<Particle> toplep_jets() const{return m_toplep_jets;}
  LorentzVector top_v4() const{ return m_electron.charge() > 0 ? m_toplep_v4 : m_tophad_v4;}
  LorentzVector antitop_v4() const{ return m_electron.charge() < 0 ? m_toplep_v4 : m_tophad_v4;}
  LorentzVector wlep_v4() const{ return m_neutrino_v4+m_electron.v4();}
  LorentzVector blep_v4() const{return m_blep_v4;}

  /// get the discriminator value for this hypothesis; thows a runtime_error if it does not exist.
  float discriminator(const std::string & l) const {
      auto it = m_discriminators.find(l);
      if(it == m_discriminators.end()){
          throw std::runtime_error("LQReconstructionHypothesis::discriminator: discriminator with label '" + l + "' not set");
      }
      return it->second;
  }
  
  /// test if a discriminator value with a certian label has already been added
  bool has_discriminator(const std::string & label) const {
      return m_discriminators.find(label) != m_discriminators.end();
  }
  
  void set_blep_v4(LorentzVector v4){m_blep_v4=v4;}
  void set_toplep_v4(LorentzVector v4){m_toplep_v4=v4;}
  void set_tophad_v4(LorentzVector v4){m_tophad_v4=v4;} 
  void set_neutrino_v4(LorentzVector v4){m_neutrino_v4=v4;}
  void add_toplep_jet(const Particle & j){m_toplep_jets.push_back(j);}
  void add_tophad_jet(const Particle & j){m_tophad_jets.push_back(j);}
  void set_electron(const Particle & l){m_electron = l;}
  void set_mu_lep(const Particle & l){m_mu_lep = l;}
  void set_mu_had(const Particle & l){m_mu_had = l;}

  void set_electron_v4(LorentzVector v4){m_electron_v4=v4;}
  void set_mu_lep_v4(LorentzVector v4){m_mu_lep_v4=v4;}
  void set_mu_had_v4(LorentzVector v4){m_mu_had_v4=v4;}

  void set_discriminator(const std::string & label, float discr){
      m_discriminators[label] = discr;
  }
  
private:
  LorentzVector m_blep_v4;
  LorentzVector m_toplep_v4;
  LorentzVector m_tophad_v4;
  LorentzVector m_neutrino_v4;

  LorentzVector m_electron_v4;
  LorentzVector m_mu_lep_v4;
  LorentzVector m_mu_had_v4;

  std::vector<Particle> m_tophad_jets;
  std::vector<Particle> m_toplep_jets;
  Particle m_electron;
  Particle m_mu_lep;
  Particle m_mu_had;

  std::map<std::string, float> m_discriminators;
};


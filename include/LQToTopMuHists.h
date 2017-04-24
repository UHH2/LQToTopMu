#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesis.h"
#include "UHH2/LQToTopMu/include/LQGen.h"


namespace uhh2examples {

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class LQToTopMuHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    LQToTopMuHists(uhh2::Context & ctx, const std::string & dirname, const bool & _isSR);

    virtual void fill(const uhh2::Event & ev) override;

  protected:
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_hadr_hyps;

    uhh2::Event::Handle<std::vector<LorentzVector>> h_jet_ele_ele;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_jet_ele_fele;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_jet_fele_fele;

    uhh2::Event::Handle<std::vector<LorentzVector>> h_ele;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_jet;

    uhh2::Event::Handle<std::vector<LorentzVector>> h_fele_nogen;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fjet_nogen;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fele_nogen_inZ;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fjet_nogen_inZ;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fele_nogen_notZ;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fjet_nogen_notZ;

    uhh2::Event::Handle<std::vector<LorentzVector>> h_ele_close;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fele_nogen_close;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_ele_notclose;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fele_nogen_notclose;

    uhh2::Event::Handle<std::vector<LorentzVector>> h_jet_close;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fjet_nogen_close;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_jet_notclose;
    uhh2::Event::Handle<std::vector<LorentzVector>> h_fjet_nogen_notclose;

    uhh2::Event::Handle<std::vector<double>> h_M_fele_ele;
    uhh2::Event::Handle<std::vector<double>> h_M_ele_ele_wofele;

    uhh2::Event::Handle<std::vector<double>> h_weights;

    uhh2::Event::Handle<bool> h_is_pt_weight;


    std::string m_discriminator_name;
    bool is_mc, isSR, is_apply_SF, cross;


    virtual ~LQToTopMuHists();
  };

}

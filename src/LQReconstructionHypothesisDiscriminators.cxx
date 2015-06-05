#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/core/include/Utils.h"

#include <set>

using namespace uhh2;
using namespace std;

namespace {
    
// invariant mass of a lorentzVector, but save for timelike / spacelike vectors
float inv_mass(const LorentzVector & p4){
    if(p4.isTimelike()){
            return p4.mass();
    }
    else{
        return -sqrt(-p4.mass2());
    }
}

}


const LQReconstructionHypothesis * get_best_hypothesis(const std::vector<LQReconstructionHypothesis> & hyps, const std::string & label){
    const LQReconstructionHypothesis * best = nullptr;
    float current_best_disc = numeric_limits<float>::infinity();
    for(const auto & hyp : hyps){
        if(!hyp.has_discriminator(label)) continue;
        auto disc = hyp.discriminator(label);
        if(disc < current_best_disc){
            best = &hyp;
            current_best_disc = disc;
        }
    }
    if(std::isfinite(current_best_disc)){
        return best;
    }
    else{
        return nullptr;
    }
}

LQChi2Discriminator::LQChi2Discriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>(rechyps_name);
}


bool LQChi2Discriminator::process(uhh2::Event & event){
    auto & hyps = event.get(h_hyps);
    const double mass_thad = 181;
    const double mass_thad_sigma = 15;
    const double mass_tlep = 174;
    const double mass_tlep_sigma = 18;
    for(auto & hyp: hyps){
        double mass_thad_rec = inv_mass(hyp.tophad_v4());
        double mass_tlep_rec = inv_mass(hyp.toplep_v4());
        double chi2_thad = pow((mass_thad_rec - mass_thad) / mass_thad_sigma, 2);
        double chi2_tlep = pow((mass_tlep_rec - mass_tlep) / mass_tlep_sigma, 2);
        hyp.set_discriminator(config.discriminator_label, chi2_thad + chi2_tlep);
        hyp.set_discriminator(config.discriminator_label + "_tlep", chi2_tlep);
        hyp.set_discriminator(config.discriminator_label + "_thad", chi2_thad);
  }
  return true;
}


LQTopDRMCDiscriminator::LQTopDRMCDiscriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>(rechyps_name);
    h_ttbargen = ctx.get_handle<TTbarGen>(config.ttbargen_name);
}


bool LQTopDRMCDiscriminator::process(uhh2::Event & event){
    auto & hyps = event.get(h_hyps);
    const auto & ttbargen = event.get(h_ttbargen);
    for(auto & hyp: hyps){
        auto deltar_sum = deltaR(ttbargen.Top(), hyp.top_v4()) + deltaR(ttbargen.Antitop(), hyp.antitop_v4());
        hyp.set_discriminator(config.discriminator_label, deltar_sum);
    }
    return true;
}



LQCorrectMatchDiscriminator::LQCorrectMatchDiscriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>(rechyps_name);
    h_ttbargen = ctx.get_handle<TTbarGen>(config.ttbargen_name);
    h_LQLQbargen = ctx.get_handle<LQGen>(config.LQLQbargen_name);
}

namespace {

// match particle p to one of the jets (Delta R < 0.3); return the deltaR
// of the match.
template<typename T> // T should inherit from Particle
float match_dr(const Particle & p, const std::vector<T> & jets, int& index){
  float mindr = infinity;
  index = -1;
  for(unsigned int i=0; i<jets.size(); ++i){
    float dR = deltaR(p, jets.at(i));
    if( dR <0.3 && dR<mindr) {
      mindr=dR;
      index=i;
    }
  }
  return mindr;
}


}

bool LQCorrectMatchDiscriminator::process(uhh2::Event & event){
    auto & hyps = event.get(h_hyps);
    const auto & ttbargen = event.get(h_ttbargen);
    const auto & LQLQbargen = event.get(h_LQLQbargen);
    auto dec = ttbargen.DecayChannel();
    if(dec != TTbarGen::e_ehad){ // if not 1x Electron, 1x hadronic
        for(auto & hyp: hyps){
            hyp.set_discriminator(config.discriminator_label, infinity);
        }
        return true;
    }
    
    // note that it is allowed that two partons from the hadronic ttbar decay match the same jet.
    for(auto & hyp: hyps){
        auto hadr_jets = hyp.tophad_jets();
        auto lept_jets = hyp.toplep_jets();
	auto hadr_mu = hyp.mu_had();
	auto lept_mu = hyp.mu_lep();
	//auto ele = hyp.electron();
        
        if(lept_jets.size() != 1){
            hyp.set_discriminator(config.discriminator_label, infinity);
            continue;
        }
        if(hadr_jets.size() > 3){ // < 3 is allowed ...
            hyp.set_discriminator(config.discriminator_label, infinity);
            continue;
        }

        //index lists of jets that can be matched to partons
        std::set<int> matched_hadr_jets;

        // match b jets
        int index_l, index_h;
        float correct_dr = match_dr(ttbargen.BLep(), lept_jets, index_l) + match_dr(ttbargen.BHad(), hadr_jets, index_h);
        if(index_h >= 0) matched_hadr_jets.insert(index_h);
        //match quarks from W decays
        correct_dr += match_dr(ttbargen.Q1(), hadr_jets, index_h);
        if(index_h >= 0) matched_hadr_jets.insert(index_h);
        correct_dr += match_dr(ttbargen.Q2(), hadr_jets, index_h);
        if(index_h >= 0) matched_hadr_jets.insert(index_h);
        
        // if not all jets of the hadronic side of the reconstruction could be matched: infinite
        // value:
        if(matched_hadr_jets.size() != hadr_jets.size()){
            hyp.set_discriminator(config.discriminator_label, infinity);
            continue;
        }

	//match muons: always 1 lept. & 1 had. muon
	float dR_mu_lept1 = deltaR(LQLQbargen.muLQ(), lept_mu);
	float dR_mu_lept2 = deltaR(LQLQbargen.muAntiLQ(), lept_mu);
	float dR_mu_hadr1 = deltaR(LQLQbargen.muLQ(), hadr_mu);
	float dR_mu_hadr2 = deltaR(LQLQbargen.muAntiLQ(), hadr_mu);
	//float dR_min_lept = infinity;
	//float dR_min_hadr = infinity;
	float dR_min_ges = infinity;

	//calculate dR for all possible matches (R=0.3): exactly 1 muon is assigned to each LQ
	//lept1 : muLQ - hypothesis
	if(dR_mu_lept1 <= 0.3 && dR_mu_hadr2 <= 0.3){
	  if(dR_mu_lept1 + dR_mu_hadr2 < dR_min_ges){
	    dR_min_ges = dR_mu_lept1 + dR_mu_hadr2;
	    //dR_min_lept = dR_mu_lept1;
	    //dR_min_hadr = dR_mu_hadr2;
	  }
	}
	//lept2 : muLQ - hypothesis
	if(dR_mu_lept2 <= 0.3 && dR_mu_hadr1 <= 0.3){
	  if(dR_mu_lept2 + dR_mu_hadr1 < dR_min_ges){
	    dR_min_ges = dR_mu_lept2 + dR_mu_hadr1;
	    //dR_min_lept = dR_mu_lept2;
	    //dR_min_hadr = dR_mu_hadr1; 
	  }
	}

	//kick out events, where the reconstructed muons are too close to each other or a jet
	int dummie_index;
	if(match_dr(lept_mu, hadr_jets, dummie_index) <= 0.3) dR_min_ges = infinity;
	if(match_dr(lept_mu, lept_jets, dummie_index) <= 0.3) dR_min_ges = infinity;
	if(match_dr(hadr_mu, hadr_jets, dummie_index) <= 0.3) dR_min_ges = infinity;
	if(match_dr(hadr_mu, lept_jets, dummie_index) <= 0.3) dR_min_ges = infinity;
	if(deltaR(lept_mu, hadr_mu) <= 0.3) dR_min_ges = infinity;

	//add minimum of dR. which muon was assigned to which LQ is irrelevant.
	correct_dr += dR_min_ges;

        //add deltaR between reconstructed and true neutrino
        correct_dr += deltaR(ttbargen.Neutrino(), hyp.neutrino_v4());
        hyp.set_discriminator(config.discriminator_label, correct_dr);
    }
    return true;
}


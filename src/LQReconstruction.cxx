#include "UHH2/LQToTopMu/include/LQReconstruction.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include <cassert>

using namespace uhh2;
using namespace std;

LQPrimaryLepton::LQPrimaryLepton(Context & ctx) {
    h_primlep = ctx.get_handle<FlavorParticle>("LQPrimaryLepton");
}

bool LQPrimaryLepton::process(uhh2::Event & event) {
  assert(/*event.muons || */event.electrons);
    double ptmax = -infinity;
    FlavorParticle primlep;
    if(event.electrons) {
        for(const auto & ele : *event.electrons) {
            if(ele.pt() > ptmax) {
                ptmax = ele.pt();
                primlep = ele;
            }
        }
    }
    /*if(event.muons) {
      for(const auto & mu : *event.muons) {
      if(mu.pt() > ptmax) {
      ptmax = mu.pt();
      primlep = mu;
      }
      }
      }*/
    event.set(h_primlep, std::move(primlep));
    return true;
}

LQPrimaryLepton::~LQPrimaryLepton() {}

HighMassLQReconstruction::HighMassLQReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label): m_neutrinofunction(neutrinofunction) {
    h_recohyps = ctx.declare_event_output<vector<LQReconstructionHypothesis>>(label);
    h_primlep = ctx.get_handle<FlavorParticle>("LQPrimaryLepton");
}

HighMassLQReconstruction::~HighMassLQReconstruction() {}

bool HighMassLQReconstruction::process(uhh2::Event & event) {
    assert(event.jets);
    assert(event.met);
    //find primary charged lepton
    const Particle & electron = event.get(h_primlep); //always an electron, as only electrons are considered possible primary leptons
    LorentzVector electron_v4 = electron.v4();
    std::vector<LQReconstructionHypothesis> recoHyps;
    //reconstruct neutrino
    std::vector<LorentzVector> neutrinos = m_neutrinofunction( electron.v4(), event.met->v4());
    unsigned int n_muons = event.muons->size();
    unsigned int n_jets = event.jets->size();
    if(n_jets>10) n_jets=10; //avoid crashes in events with many jets
    // idea: loop over 3^Njet possibilities and write the current loop
    // index j in the 3-base system. The Njets digits represent whether
    // to assign each jet to the hadronic side (0), leptonic side (1),
    // or none of them (2).
    const unsigned int max_j = pow(3, n_jets); 

    //loop over neutrino solutions and jet assignments to fill hyotheses
    for(const auto & neutrino_p4 : neutrinos) {
      const LorentzVector wlep_v4 = electron.v4() + neutrino_p4;
      for (unsigned int j=0; j < max_j; j++) {
	LorentzVector tophad_v4;
	LorentzVector toplep_v4 = wlep_v4;
	int hadjets=0;
	int lepjets=0;
	int num = j;
	LQReconstructionHypothesis hyp;
	hyp.set_electron(electron);
	hyp.set_electron_v4(electron_v4);
	hyp.set_neutrino_v4(neutrino_p4);
	for (unsigned int k=0; k<n_jets; k++) {
	  if(num%3==0) {
	    tophad_v4 = tophad_v4 + event.jets->at(k).v4();
	    hyp.add_tophad_jet(event.jets->at(k));
	    hadjets++;
	  }
	  
	  if(num%3==1) {
	    toplep_v4 = toplep_v4 + event.jets->at(k).v4();
	    hyp.add_toplep_jet(event.jets->at(k));
	    lepjets++;
	  }
	  //in case num%3==2 do not take this jet at all
	  //shift the trigits of num to the right:
	  num /= 3;
	}
	
	//search jet with highest pt assigned to leptonic top
	int blep_idx(-1);
	float maxpt(-1.);
	for(unsigned int i=0; i<hyp.toplep_jets().size(); ++i){
	  if(maxpt< hyp.toplep_jets().at(i).pt()){
	    maxpt = hyp.toplep_jets().at(i).pt();
	    blep_idx = i;
	  }
	}
	if(blep_idx != -1) hyp.set_blep_v4(hyp.toplep_jets().at(blep_idx).v4());
	
	//fill only hypotheses with at least one jet assigned to each top quark
	if(hadjets>0 && lepjets>0) {
	  int max_i = pow(3,n_muons); // analogous to jet combinations
	  for(int i=0; i<max_i; i++){ // for each jet comb loop over all possible muon combs
	    LorentzVector mu1_v4;
	    LorentzVector mu2_v4;
	    int hadmu=0;
	    int lepmu=0;
	    int num = i;
	    for(unsigned int k=0; k<n_muons; k++){
	      if(num%3==0){
		mu1_v4 = event.muons->at(k).v4();
		hyp.set_mu_had(event.muons->at(k)); // had is only a hypothesis not having considered the charge!!!
		hadmu++;
	      }
	      if(num%3==1){
		mu2_v4 = event.muons->at(k).v4();
		hyp.set_mu_lep(event.muons->at(k)); // lep is only a hypothesis not having considered the charge!!!
		lepmu++;
	      }
	      //for num%3==2 do nothing
	      num /= 3;
	    }
	    
	    Particle mu_2 = hyp.mu_lep();
	    if(hadmu==1 && lepmu==1){ //require exactly 1 muon assigned to each top
	      if(mu_2.charge() != electron.charge()){ //electron and leptonic mu must have opposite charges
		hyp.set_mu_had_v4(mu1_v4);
		hyp.set_mu_lep_v4(mu2_v4); // mu2 really is the leptonic one.
		hyp.set_tophad_v4(tophad_v4);
		hyp.set_toplep_v4(toplep_v4);
		recoHyps.emplace_back(move(hyp));
	      } // charge comparison
	      else{
		hyp.set_mu_had_v4(mu2_v4);
		hyp.set_mu_lep_v4(mu1_v4); // mu1 really is the leptonic one, the original hypothesis was wrong.
		hyp.set_tophad_v4(tophad_v4);
		hyp.set_toplep_v4(toplep_v4);
		recoHyps.emplace_back(move(hyp));	
	      } // charge comparison_2
	    } // 1 muon per top
	  } // muon combs for-loop
	} // if at least 1 jet is assigned to each top quark
      } // 3^n_jets jet combinations * n_muon muon combinations
    } // neutrinos
    event.set(h_recohyps, move(recoHyps));
    return true;
}




HighMassHadronicLQReconstruction::HighMassHadronicLQReconstruction(Context & ctx, const string & label) {
    h_hadr_recohyps = ctx.declare_event_output<vector<LQReconstructionHypothesis>>(label);
}

HighMassHadronicLQReconstruction::~HighMassHadronicLQReconstruction() {}

bool HighMassHadronicLQReconstruction::process(uhh2::Event & event) {
    assert(event.jets);
    std::vector<LQReconstructionHypothesis> recoHyps;

    unsigned int n_jets = event.jets->size();
    if(n_jets>12) n_jets=12; //avoid crashes in events with many jets
    // idea: loop over 3^Njet possibilities and write the current loop
    // index j in the 3-base system. The Njets digits represent whether
    // to assign each jet to the hadronic side (0), leptonic side (1),
    // or none of them (2).
    const unsigned int max_j = pow(3, n_jets);

    //loop over neutrino solutions and jet assignments to fill hyotheses
    for (unsigned int j=0; j < max_j; j++) {
      LorentzVector tophad1_v4;
      LorentzVector tophad2_v4;
      int hadjets1=0;
      int hadjets2=0;
      int num = j;
      LQReconstructionHypothesis hyp;
      for (unsigned int k=0; k<n_jets; k++) {
	if(num%3==0) {
	  tophad1_v4 = tophad1_v4 + event.jets->at(k).v4();
	  hyp.add_tophad1_jet(event.jets->at(k));
	  hadjets1++;
	}

	if(num%3==1) {
	  tophad2_v4 = tophad2_v4 + event.jets->at(k).v4();
	  hyp.add_tophad2_jet(event.jets->at(k));
	  hadjets2++;
	}
	//in case num%3==2 do not take this jet at all
	//shift the trigits of num to the right:
	num /= 3;
      }

      //fill only hypotheses with at least one jet assigned to each top quark
      if(hadjets1>0 && hadjets2>0) {
	LorentzVector mu1_v4;
	LorentzVector mu2_v4;

	for(int i=0; i<2; i++){
	hyp.set_mu_had1(event.muons->at(i));
	mu1_v4 = event.muons->at(i).v4();
	hyp.set_mu_had2(event.muons->at(1-i));
	mu2_v4 = event.muons->at(1-i).v4();
	
	hyp.set_mu_had1_v4(mu1_v4);
	hyp.set_mu_had2_v4(mu2_v4);
	hyp.set_tophad1_v4(tophad1_v4);
	hyp.set_tophad2_v4(tophad2_v4);
	recoHyps.emplace_back(move(hyp));
	} // 2 possible muon combinations	
      } // if at least 1 jet is assigned to each top quark
    } // 3^n_jets jet combinations
    event.set(h_hadr_recohyps, move(recoHyps));
    return true;
}





















LQTopTagReconstruction::LQTopTagReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label, TopJetId tjetid, float minDR_tj_j):
  m_neutrinofunction(neutrinofunction), topjetID_(tjetid), minDR_topjet_jet_(minDR_tj_j) {

  h_recohyps = ctx.declare_event_output<vector<LQReconstructionHypothesis>>(label);
  h_primlep = ctx.get_handle<FlavorParticle>("LQPrimaryLepton");
}

bool LQTopTagReconstruction::process(uhh2::Event & event) {

  assert(event.jets && event.topjets);
  assert(event.met);

  std::vector<LQReconstructionHypothesis> recoHyps;

  const Particle& lepton = event.get(h_primlep);
  std::vector<LorentzVector> neutrinos = m_neutrinofunction(lepton.v4(), event.met->v4());

  for(const auto& tj : *event.topjets){

    if(!topjetID_(tj, event)) continue;

    // jet candidates for leptonic-top (not overlapping with top-tagged jet)
    std::vector<const Jet*> tlep_jets;
    tlep_jets.reserve(event.jets->size());
    for(const auto & jet : *event.jets)
      if(deltaR(tj, jet) > minDR_topjet_jet_) tlep_jets.push_back(&jet);

    const unsigned int jet_combs = pow(2, tlep_jets.size());

    for(const auto& neutrino_p4 : neutrinos){

      for(unsigned int i=1; i<jet_combs; ++i){

        LQReconstructionHypothesis hyp;
        hyp.set_electron(lepton);
        hyp.set_neutrino_v4(neutrino_p4);

        LorentzVector tophad_v4(tj.v4());
        hyp.add_tophad_jet(tj);

        LorentzVector toplep_v4(lepton.v4() + neutrino_p4);

        for(unsigned int j=0; j<tlep_jets.size(); ++j){
          // index for jet assignment to top leg (0=none, 1=leptonic-top)
          int jet_topidx = int(i/(pow(2,j))) % 2;

          if(jet_topidx == 1){
            toplep_v4 += tlep_jets.at(j)->v4();
            hyp.add_toplep_jet(*tlep_jets.at(j));
          }
        }

        // b-jet of leptonic top (pt-leading)
        int blep_idx(-1);
        float maxpt(-1.);
        for(unsigned int i=0; i<hyp.toplep_jets().size(); ++i){
          if(maxpt< hyp.toplep_jets().at(i).pt()){
            maxpt = hyp.toplep_jets().at(i).pt();
            blep_idx = i;
          }
        }
        if(blep_idx != -1) hyp.set_blep_v4(hyp.toplep_jets().at(blep_idx).v4());

        if(hyp.tophad_jets().size() && hyp.toplep_jets().size()){
          hyp.set_tophad_v4(tophad_v4);
          hyp.set_toplep_v4(toplep_v4);
          recoHyps.emplace_back(std::move(hyp));
        }
      }
    }
  }

  event.set(h_recohyps, std::move(recoHyps));
  return true;
}

std::vector<LorentzVector> LQNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met) {
    TVector3 lepton_pT = toVector(lepton);
    lepton_pT.SetZ(0);
    TVector3 neutrino_pT = toVector(met);
    neutrino_pT.SetZ(0);
    constexpr float mass_w = 80.399f;
    float mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
    float A = - (lepton_pT * lepton_pT);
    float B = mu * lepton.pz();
    float C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
    float discriminant = B * B - A * C;
    std::vector<LorentzVector> solutions;
    if (0 >= discriminant) {
        // Take only real part of the solution for pz:
        LorentzVectorXYZE solution (met.Px(),met.Py(),-B / A,0);
        solution.SetE(solution.P());
        solutions.emplace_back(toPtEtaPhi(solution));
    }
    else {
        discriminant = sqrt(discriminant);
        LorentzVectorXYZE solution (met.Px(),met.Py(),(-B - discriminant) / A,0);
        solution.SetE(solution.P());
        solutions.emplace_back(toPtEtaPhi(solution));

        LorentzVectorXYZE solution2 (met.Px(),met.Py(),(-B + discriminant) / A,0);
        solution2.SetE(solution2.P());
        solutions.emplace_back(toPtEtaPhi(solution2));
    }
    return solutions;
}

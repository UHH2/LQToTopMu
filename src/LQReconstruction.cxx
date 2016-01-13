#include "UHH2/LQToTopMu/include/LQReconstruction.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include <cassert>

using namespace uhh2;
using namespace std;

float inv_mass(const LorentzVector & p4){
    if(p4.isTimelike()){
            return p4.mass();
    }
    else{
        return -sqrt(-p4.mass2());
    }
}

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
    int n_final_hyps = 0;
    /*double chi2 = 9999999999;
    double chi2min = chi2;
    const double mass_thad = 181;
    const double mass_thad_sigma = 15;
    const double mass_tlep = 174;
    const double mass_tlep_sigma = 18;
    const double mass_LQ_diff = -12; // from Histo-sum: -12, from M500: -13, from M1300: -12
    const double mass_LQ_diff_sigma = 46; // from Histo-sum: 46, from M500: 63, from M1300: 152
    */

    //find primary charged lepton
    const Particle & electron = event.get(h_primlep); //always an electron, as only electrons are considered possible primary leptons
    LorentzVector electron_v4 = electron.v4();
    std::vector<LQReconstructionHypothesis> recoHyps;

    //reconstruct neutrino
    std::vector<LorentzVector> neutrinos = m_neutrinofunction( electron.v4(), event.met->v4());

    unsigned int n_muons = event.muons->size();
    unsigned int n_jets = event.jets->size();
    if(n_jets>7) n_jets=7; //avoid crashes in events with many jets
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
	    Particle mu_1 = hyp.mu_had();
	    if(hadmu==1 && lepmu==1 && (mu_1.charge()!=mu_2.charge())){ //require exactly 1 muon assigned to each top, opposite charges
	      if(mu_2.charge() != electron.charge()){ //electron and leptonic mu must have opposite charges
		hyp.set_mu_had_v4(mu1_v4);
		hyp.set_mu_lep_v4(mu2_v4); // mu2 really is the leptonic one.
		hyp.set_mu_had(hyp.mu_had());
		hyp.set_mu_lep(hyp.mu_lep());
		hyp.set_tophad_v4(tophad_v4);
		hyp.set_toplep_v4(toplep_v4);
	      } // charge comparison
	      else{
		hyp.set_mu_had_v4(mu2_v4);
		hyp.set_mu_lep_v4(mu1_v4); // mu1 really is the leptonic one, the original hypothesis was wrong.
		hyp.set_mu_had(hyp.mu_lep());
		hyp.set_mu_lep(hyp.mu_had());
		hyp.set_tophad_v4(tophad_v4);
		hyp.set_toplep_v4(toplep_v4);
	      } // charge comparison_2
	      
	      // orginally: save all hypotheses and calculate chi2 in a second step -> now calculate chi2 on the fly and only save the best hypothesis
	      // makes get_best_hypothesis obsolete for this reconstruction
	      //recoHyps.emplace_back(move(hyp));

	      //calculate chi2 for this hypothesis
	      /*double mass_thad_rec = inv_mass(hyp.tophad_v4());
	      double mass_tlep_rec = inv_mass(hyp.toplep_v4());
	      double mass_LQ_had_rec = inv_mass(hyp.LQhad_v4()); // added
	      double mass_LQ_lep_rec = inv_mass(hyp.LQlep_v4()); // added
	      double chi2_thad = pow((mass_thad_rec - mass_thad) / mass_thad_sigma, 2);
	      double chi2_tlep = pow((mass_tlep_rec - mass_tlep) / mass_tlep_sigma, 2);
	      double chi2_MLQdiff = pow(((mass_LQ_had_rec - mass_LQ_lep_rec) - mass_LQ_diff) / mass_LQ_diff_sigma, 2); // added
	      chi2 = chi2_thad + chi2_tlep + chi2_MLQdiff;
	      if(chi2<chi2min){
		chi2min = chi2;
		if(recoHyps.size() == 1){
		  //cout << "recoHyps-size = 1, popping back the element" << endl;
		  recoHyps.pop_back();
		  //cout << "emplacing back new best hypothesis" << endl;
		  recoHyps.emplace_back(move(hyp));
		}
		else if(recoHyps.size() == 0){*/
		  //cout << "recoHyps empty, emplacing back best (=first) hypothesis" << endl;
	      recoHyps.emplace_back(move(hyp));
	      /*cout << "in hyp no. " << n_final_hyps << ": mtoplep: " << hyp.toplep_v4().M() << " with " << hyp.toplep_jets().size() << " leptonic jets." << endl;
		  cout << "same hyp: lep_v4 timelike?: " << hyp.toplep_v4().isTimelike() << endl;
		  cout << "ele 4 vec:     " << hyp.electron_v4() << "N electrons: " << event.electrons->size() << endl;
		  cout << "nu  4 vec:     " << hyp.neutrino_v4() << endl;
		  cout << "w   4 vec:     " << hyp.wlep_v4() << endl;
		  cout << "top 4 vec:     " << hyp.toplep_v4() << endl;
		  cout << "mulep vec:     " << hyp.mu_lep_v4() << endl;
		  cout << "LQlep vec:     " << hyp.LQlep_v4() << endl;*/
		  n_final_hyps++;
		  /*}
		  else throw runtime_error("size of recoHyps neither 0 nor 1");
		  }*/

	    } // 1 muon per top
	  } // muon combs for-loop
	} // if at least 1 jet is assigned to each top quark
      } // 3^n_jets jet combinations * n_muon muon combinations
    } // neutrinos

    //cout << "# of considered reco hyps: " << n_final_hyps << endl;
    event.set(h_recohyps, move(recoHyps));
    return true;
}




/*HighMassHadronicLQReconstruction::HighMassHadronicLQReconstruction(Context & ctx, const string & label) {
    h_hadr_recohyps = ctx.declare_event_output<vector<LQReconstructionHypothesis>>(label);
}

HighMassHadronicLQReconstruction::~HighMassHadronicLQReconstruction() {}

bool HighMassHadronicLQReconstruction::process(uhh2::Event & event) {
    assert(event.jets);
    double chi2 = 9999999999;
    double chi2min = chi2;
    const double mass_thad1 = 181;
    const double mass_thad1_sigma = 15;
    const double mass_LQ_diff = -12; // from Histo-sum: -12, from M500: -13, from M1300: -12
    const double mass_LQ_diff_sigma = 46; // from Histo-sum: 46, from M500: 63, from M1300: 152
    //const double mass_W = 80.3;
    //const double mass_W_sigma = 2.1;
    const double PTLQLQ = 0;
    const double PTLQLQ_sigma = 150;


    std::vector<LQReconstructionHypothesis> recoHyps;

    unsigned int n_jets = event.jets->size();
    if(n_jets>7) n_jets=7; //avoid crashes in events with many jets
    // idea: loop over 3^Njet possibilities and write the current loop
    // index j in the 3-base system. The Njets digits represent whether
    // to assign each jet to the hadronic side (0), leptonic side (1),
    // or none of them (2).
    const unsigned int max_j = pow(3, n_jets);

    //loop over neutrino solutions and jet assignments to fill hyotheses
    for (unsigned int j=0; j < max_j; j++) {
      LorentzVector tophad1_v4;
      LorentzVector tophad2_v4;
      //LorentzVector w1_v4;
      //LorentzVector w2_v4;
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

	if((unsigned int)hadjets1 != hyp.tophad1_jets().size() || (unsigned int)hadjets2 != hyp.tophad2_jets().size() ) throw runtime_error("0mismatch between hadjets and real size()");

	LorentzVector mu1_v4;
	LorentzVector mu2_v4;
 
	if(event.muons->size() < 2) return false;
	for(int i=0; i<2; i++){
	  hyp.set_mu_had1(event.muons->at(i));
	  mu1_v4 = event.muons->at(i).v4();
	  hyp.set_mu_had2(event.muons->at(1-i));
	  mu2_v4 = event.muons->at(1-i).v4();
	
	  hyp.set_mu_had1_v4(mu1_v4);
	  hyp.set_mu_had2_v4(mu2_v4);
	  hyp.set_tophad1_v4(tophad1_v4);
	  hyp.set_tophad2_v4(tophad2_v4);*/

	  //find W-mass inside of each top hypothesis

	  //tophad1
	  /*double chi2_w1 = 99999999999;
	  double chi2min_w1 = chi2_w1;
	  int n_w1jets = 0;
	  int combs_max = pow(2,hyp.tophad1_jets().size()-1);
	  for(int n=0; n<combs_max; n++){
	    n_w1jets = 0;
	    int iter = n;
	    for(unsigned int m=0; m<hyp.tophad1_jets().size()-1; m++){
	      if(iter%2==0){
		w1_v4 = w1_v4 + hyp.tophad1_jets().at(m).v4();
		n_w1jets++;
	      }

	      iter/=2;
	    }

	    //get best w-hypothesis
	    if(n_w1jets > 0){
	      chi2_w1 = pow((w1_v4.M()-mass_W) / mass_W_sigma,2);
	      if(chi2_w1 < chi2min_w1){
		chi2min_w1 = chi2_w1;
		hyp.set_whad1_v4(w1_v4);
	      }
	    }
	  }*/

	  //tophad2
	  /*double chi2_w2 = 99999999999;
	  double chi2min_w2 = chi2_w2;
	  int n_w2jets = 0;
	  combs_max = pow(2,hyp.tophad2_jets().size()-1);
	  for(int n=0; n<combs_max; n++){
	    n_w2jets = 0;
	    int iter = n;
	    for(unsigned int m=0; m<hyp.tophad2_jets().size()-1; m++){
	      if(iter%2==0){
		w2_v4 = w2_v4 + hyp.tophad2_jets().at(m).v4();
		n_w2jets++;
	      }

	      iter/=2;
	    }

	    //get best w-hypothesis
	    if(n_w2jets > 0){
	      chi2_w2 = pow((w2_v4.M()-mass_W) / mass_W_sigma,2);
	      if(chi2_w2 < chi2min_w2){
		chi2min_w2 = chi2_w2;
		hyp.set_whad2_v4(w2_v4);
	      }
	    }
	  }*/


	  //original way, cf above
	  //recoHyps.emplace_back(move(hyp));
	
	  //calculate chi2 for this hypothesis
	
	  /*double mass_thad1_rec = inv_mass(hyp.tophad1_v4());
	  double mass_thad2_rec = inv_mass(hyp.tophad2_v4());
	  double mass_LQ_had1_rec = inv_mass(hyp.LQhad1_v4()); // added
	  double mass_LQ_had2_rec = inv_mass(hyp.LQhad2_v4()); // added
	  double PTLQLQ_rec = (hyp.LQhad1_v4()+hyp.LQhad2_v4()).Pt();
	  //double mass_whad1_rec = inv_mass(hyp.whad1_v4());
	  //double mass_whad2_rec = inv_mass(hyp.whad2_v4());
	  double chi2_thad1 = pow((mass_thad1_rec - mass_thad1) / mass_thad1_sigma, 2);
	  double chi2_thad2 = pow((mass_thad2_rec - mass_thad1) / mass_thad1_sigma, 2);
	  double chi2_MLQdiff = pow(((mass_LQ_had1_rec - mass_LQ_had2_rec) - mass_LQ_diff) / mass_LQ_diff_sigma, 2); // added
	  double chi2_PTLQLQ = pow((PTLQLQ_rec - PTLQLQ) / PTLQLQ_sigma,2);
	  //double chi2_whad1 = pow((mass_whad1_rec-mass_W) / mass_W_sigma,2);
	  //double chi2_whad2 = pow((mass_whad2_rec-mass_W) / mass_W_sigma,2);
	  chi2 = chi2_thad1 + chi2_thad2 + chi2_MLQdiff + chi2_PTLQLQ + chi2_whad1 + chi2_whad2;


	  if(chi2<chi2min){
	    chi2min = chi2;
	    if(recoHyps.size() == 1){
	      //cout << "recoHyps-size = 1, popping back the element" << endl;
	      recoHyps.pop_back();
	      //cout << "emplacing back new best hypothesis" << endl;
	      recoHyps.push_back(hyp);
	    }
	    else if(recoHyps.size() == 0){
	      //cout << "recoHyps empty, emplacing back best (=first) hypothesis" << endl;
	      recoHyps.push_back(hyp);
	    }
	    else throw runtime_error("size of recoHyps neither 0 nor 1");
	  }
	} // 2 possible muon combinations	
      } // if at least 1 jet is assigned to each top quark
    } // 3^n_jets jet combinations

    event.set(h_hadr_recohyps, move(recoHyps));
    return true;
}*/





















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

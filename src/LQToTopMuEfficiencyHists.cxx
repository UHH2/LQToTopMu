#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

#include "TH1F.h"
#include "TH1D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuEfficiencyHists::LQToTopMuEfficiencyHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1D>("N_ele", "N_{e}", 10, -0.5, 9.5); 
  book<TH1D>("pt_ele", "p_{T}^{e} [GeV/c]", 100, 0, 1000);
  book<TH1D>("eta_ele", "#eta^{e}", 40, -2.5, 2.5);
  book<TH1D>("reliso_ele", "ele rel. Iso", 40, 0, 0.5);
  book<TH1D>("Eff_ele_int", "Electrons    -1: unmatched   +1: matched", 3, -1.5, 1.5);
  book<TH1D>("Eff_ele_pt", "e efficiency(p_{T}^{e})", 100, 0, 1000);
  book<TH1D>("Eff_ele_ht", "e efficiency(H_{T}^{event})", 80, 0, 7000);

  book<TH1D>("N_mu", "N_{#mu}", 10, -0.5, 9.5); 
  book<TH1D>("pt_mu", "p_{T}^{#mu} [GeV/c]", 100, 0, 1000);
  book<TH1D>("eta_mu", "#eta^{#mu}", 40, -2.5, 2.5);
  book<TH1D>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1D>("Eff_mu_int", "Muons    -1: unmatched   +1: matched", 3, -1.5, 1.5);
  book<TH1D>("Eff_mu_pt", "#mu efficiency(p_{T}^{#mu})", 100, 0, 1000);
  book<TH1D>("Eff_mu_ht", "#mu efficiency(H_{T}^{event})", 80, 0, 7000);


  // general
  book<TH1D>("H_T", "H_{T}", 80, 0, 7000);
  book<TH1D>("H_T_1Mu", "H_{T}, N_{#mu} #geq 1", 80, 0, 7000);
  book<TH1D>("H_T_1Ele", "H_{T}, N_{ele} #geq 1", 80, 0, 7000);
}


void LQToTopMuEfficiencyHists::fill(const Event & event){
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      hist("pt_mu")->Fill(thismu.pt(), weight);
      hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }
int Nele = event.electrons->size();
  hist("N_ele")->Fill(Nele, weight);
 for (const Electron & thisele : *event.electrons){
      hist("pt_ele")->Fill(thisele.pt(), weight);
      hist("eta_ele")->Fill(thisele.eta(), weight);
      hist("reliso_ele")->Fill(thisele.relIso(), weight);
  }
  //HT
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
  hist("H_T")->Fill(ht, weight);

  if(Nele>0){
  // electron efficiency
  // integrated
  int N_gen_ele = 0;
  for(const auto & gp : *event.genparticles){
    int id = abs(gp.pdgId());
    if(id == 11){//gen electrons
      N_gen_ele++;
    }
  }

  int e_matched = 0;
  int e_unmatched = 0;

  if(N_gen_ele>0){
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11){//gen electrons
	//loop over reco electrons   
	bool thise_matched = false;
	for(const auto & thisele: *event.electrons){
	  if(deltaR(gp, thisele) <= 0.1 && thise_matched == false){
	    e_matched ++;
	    thise_matched = true;
	    hist("Eff_ele_int")->Fill(1,weight);
	    hist("Eff_ele_pt")->Fill(thisele.pt(),weight);
	  }
	}
	if(thise_matched == false){
	  e_unmatched++;
	  hist("Eff_ele_int")->Fill(-1,weight);
	}
      }
    }
  
    if(e_matched+e_unmatched != N_gen_ele) throw runtime_error("In Efficiency calculation: N_matched + N_unmatched != N_gen_ele");

    //cout << "e_matched: " << e_matched << endl;
    //cout << "e_unmatched: " << e_unmatched << endl;
    double e_tot_eff = e_matched / (e_matched+e_unmatched);
    //cout << "e_tot_eff: " << e_tot_eff << endl;
    double fillweight = e_tot_eff * weight;
    hist("Eff_ele_ht")->Fill(ht,fillweight);
  }

  hist("H_T_1Ele")->Fill(ht,weight);
  } // Nele > 0


  if(Nmuons>0){
  // muon efficiency
  // integrated
  int N_gen_mu = 0;
  for(const auto & gp : *event.genparticles){
    int id = abs(gp.pdgId());
    if(id == 13){//gen electrons
      N_gen_mu++;
    }
  }

  int mu_matched = 0;
  int mu_unmatched = 0;

  if(N_gen_mu>0){
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 13){//gen electrons
	//loop over reco electrons   
	bool thismu_matched = false;
	for(const auto & thismu: *event.muons){
	  if(deltaR(gp, thismu) <= 0.1 && thismu_matched == false){
	    mu_matched ++;
	    thismu_matched = true;
	    hist("Eff_mu_int")->Fill(1,weight);
	    hist("Eff_mu_pt")->Fill(thismu.pt(),weight);
	  }
	}
	if(thismu_matched == false){
	  mu_unmatched++;
	  hist("Eff_mu_int")->Fill(-1,weight);
	}
      }
    }
  
    if(mu_matched+mu_unmatched != N_gen_mu) throw runtime_error("In Efficiency calculation: N_matched + N_unmatched != N_gen_mu");
    double mu_tot_eff = mu_matched / (mu_matched+mu_unmatched);
    hist("Eff_mu_ht")->Fill(ht,mu_tot_eff*weight);
  }

  hist("H_T_1Mu")->Fill(ht,weight);
  } // Nmu > 0




} //Methode


LQToTopMuEfficiencyHists::~LQToTopMuEfficiencyHists(){}

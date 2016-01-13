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
  book<TH1D>("h_N_gen_ele", "N_{e} generated", 10, -0.5, 9.5); 
  book<TH1D>("pt_1ele", "p_{T}^{e} [GeV/c] (sideband)", 100, 0, 1000);
  book<TH1D>("pt_1ele_fakerate_side", "p_{T}^{e} [GeV/c] (Sideband)", 100, 0, 1000);
  book<TH1D>("Eff_elemu_int", "Electrons    -1: unmatched   +1: matched", 3, -1.5, 1.5);
  book<TH1D>("Eff_1ele_pt", "e efficiency(p_{T}^{e}) (Siddeband)", 100, 0, 1000);
  book<TH1D>("Eff_elemu_ht", "e efficiency(H_{T}^{event})", 80, 0, 7000);
  book<TH1D>("Eff_elemu_int_fakerate", "Electrons    -1: unmatched   +1: matched, for fakes", 3, -1.5, 1.5);
  book<TH1D>("Eff_ele_pt_fakerate_side", "e efficiency(p_{T}^{e}), for fakes (Sideband)", 100, 0, 1000);
  book<TH1D>("Eff_elemu_ht_fakerate", "e efficiency(H_{T}^{event}), for fakes", 80, 0, 7000);

  book<TH1D>("N_mu", "N_{#mu}", 10, -0.5, 9.5); 
  book<TH1D>("h_N_gen_mu", "N_{#mu} generated", 10, -0.5, 9.5); 
  book<TH1D>("pt_1mu", "p_{T}^{#mu} [GeV/c] (Sideband)", 100, 0, 1000);
  book<TH1D>("pt_2mu", "p_{T}^{#mu} [GeV/c] (Signal region)", 100, 0, 1000);
  book<TH1D>("pt_1mu_fakerate_side", "p_{T}^{#mu} [GeV/c] (Sideband)", 100, 0, 1000);
  book<TH1D>("pt_1mu_fakerate_sig", "p_{T}^{#mu} [GeV/c] (Signal region)", 100, 0, 1000);
  book<TH1D>("Eff_mumu_int", "Muons    -1: unmatched   +1: matched", 3, -1.5, 1.5);
  book<TH1D>("Eff_1mu_pt", "#mu efficiency(p_{T}^{#mu}) (Sideband)", 100, 0, 1000);
  book<TH1D>("Eff_2mu_pt", "#mu efficiency(p_{T}^{#mu}) (Signal region)", 100, 0, 1000);
  book<TH1D>("Eff_mumu_ht", "#mu efficiency(H_{T}^{event})", 80, 0, 7000);
  book<TH1D>("Eff_mumu_int_fakerate", "Muons    -1: unmatched   +1: matched, for fakes", 3, -1.5, 1.5);
  book<TH1D>("Eff_mu_pt_fakerate_side", "#mu efficiency(p_{T}^{#mu}), for fakes (Sideband)", 100, 0, 1000);
  book<TH1D>("Eff_mu_pt_fakerate_sig", "#mu efficiency(p_{T}^{#mu}), for fakes (Signal region)", 100, 0, 1000);
  book<TH1D>("Eff_mumu_ht_fakerate", "#mu efficiency(H_{T}^{event}), for fakes", 80, 0, 7000);


  // general
  book<TH1D>("H_T", "H_{T}", 80, 0, 7000);
  //book<TH1D>("H_T_1Mu", "H_{T}, N_{#mu} #geq 1", 80, 0, 7000);
  book<TH1D>("H_T_EleMu", "H_{T}, N_{#mu} = 1 and  N_{ele} = 1 (for sideband)", 80, 0, 7000);
  book<TH1D>("H_T_MuMu", "H_{T}, N_{#mu} = 2 (no events with > 2) (for signal region)", 80, 0, 7000);
  //book<TH1D>("H_T_1Ele", "H_{T}, N_{ele} #geq 1", 80, 0, 7000);
  book<TH1D>("H_T_MuMu_fakerate", "H_{T}, N_{#mu} #geq 1, for fakes", 80, 0, 7000);
  book<TH1D>("H_T_EleMu_fakerate", "H_{T}, N_{ele} #geq 1, for fakes", 80, 0, 7000);

  is_mc = ctx.get("dataset_type") == "MC";
}


void LQToTopMuEfficiencyHists::fill(const Event & event){
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);

  int Nele = event.electrons->size();
  hist("N_ele")->Fill(Nele, weight);

  int N_gen_ele = 0, N_gen_mu = 0;
  if(is_mc){
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11) N_gen_ele++;
      if(id == 13) N_gen_mu++;
    }
    hist("h_N_gen_ele")->Fill(N_gen_ele,weight);
    hist("h_N_gen_mu")->Fill(N_gen_mu,weight);
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



 if(is_mc) {


   /*

  //Reco - Effizienz Elektronen

  bool reco_ele = false;
  if(N_gen_ele==1){
    //in wie wahrscheinlich wird in den events, in denen ein elektron vorhanden ist, auch eins rekonstruiert?
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11){
	for(const auto & thisele: *event.electrons){
	  if(deltaR(gp,thisele) < 0.1){ //wurde rekonstruiert
	    reco_ele = true;
	  }
	}
	hist("pt_1ele")->Fill(gp.pt(),weight); //wahres pt aller geneles
	if(reco_ele) hist("Eff_ele_pt")->Fill(gp.pt(),weight); //wahres pt der geneles, die auch rekonstruiert werden
      }
    }
    hist("H_T_1Ele")->Fill(ht,weight); //ht, wenn N_gen_ele == 1
    if(reco_ele) hist("Eff_ele_ht")->Fill(ht,weight); //ht, wenn es auch rekonstruiert wurde  -->ht_reco / ht = eff(ht)
    if(reco_ele) hist("Eff_ele_int")->Fill(1,weight); // 1, wenn rekonstruktion erfolgreich
    else         hist("Eff_ele_int")->Fill(-1,weight);//-1, sonst --> normiert ergibt das gesamteffizienz
  }


  // Reco - Effizienz Muonen

  bool reco_mu = false;
  if(N_gen_mu==1){
    //in wie wahrscheinlich wird in den events, in denen ein muon vorhanden ist, auch eins rekonstruiert?
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 13){
	for(const auto & thismu: *event.muons){
	  if(deltaR(gp,thismu) < 0.1){ //wurde rekonstruiert
	    reco_mu = true;
	  }
	}
	hist("pt_1mu")->Fill(gp.pt(),weight); //wahres pt aller genmus
	if(reco_mu) hist("Eff_mu_pt")->Fill(gp.pt(),weight); //wahres pt der genmus, die auch rekonstruiert werden
      }
    }
    hist("H_T_1Mu")->Fill(ht,weight); //ht, wenn N_gen_mu == 1
    if(reco_mu) hist("Eff_mu_ht")->Fill(ht,weight); //ht, wenn es auch rekonstruiert wurde  -->ht_reco / ht = eff(ht)
    if(reco_mu) hist("Eff_mu_int")->Fill(1,weight); // 1, wenn rekonstruktion erfolgreich	       
    else        hist("Eff_mu_int")->Fill(-1,weight);//-1, sonst --> normiert ergibt das gesamteffizienz
  }




  // Fakerate Elektronen
  bool not_fake_ele = false;
  if(Nele == 1){ //1 Ele wurde rekonstruiert
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11){ //es ist auch eines vorhanden gewesen
	if(deltaR(event.electrons->at(0),gp) < 0.1){ // koennen sogar gematcht werden --> richtige reco
	  not_fake_ele = true;
	}
      }
    }
    if(not_fake_ele) hist("Eff_ele_pt_fakerate")->Fill(event.electrons->at(0).pt(),weight); //alle richtig rekonstruierten elektronen
    hist("pt_1ele_fakerate")->Fill(event.electrons->at(0).pt(),weight); //alle vorhandenen elektronen

    if(not_fake_ele) hist("Eff_ele_ht_fakerate")->Fill(ht,weight); // s.o.
    hist("H_T_1Ele_fakerate")->Fill(ht,weight); // s.o.

    if(not_fake_ele) hist("Eff_ele_int_fakerate")->Fill(1,weight);  // 1, wenn rekonstruktion erfolgreich	       
    else             hist("Eff_ele_int_fakerate")->Fill(-1,weight); //-1, sonst --> normiert ergibt das gesamteffizienz
  }


  // Fakerate Muonen
  bool not_fake_mu = false;
  if(Nmuons == 1){ //1 Mu wurde rekonstruiert
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 13){ //es ist auch eines vorhanden gewesen
	if(deltaR(event.muons->at(0),gp) < 0.1){ // koennen sogar gematcht werden --> richtige reco
	  not_fake_mu = true;
	}
      }
    }
    if(not_fake_mu) hist("Eff_mu_pt_fakerate")->Fill(event.muons->at(0).pt(),weight); //alle richtig rekonstruierten muonen
    hist("pt_1mu_fakerate")->Fill(event.muons->at(0).pt(),weight); //alle vorhandenen muonen

    if(not_fake_mu) hist("Eff_mu_ht_fakerate")->Fill(ht,weight); // s.o.
    hist("H_T_1Mu_fakerate")->Fill(ht,weight); // s.o.

    if(not_fake_mu) hist("Eff_mu_int_fakerate")->Fill(1,weight);  // 1, wenn rekonstruktion erfolgreich	       
    else            hist("Eff_mu_int_fakerate")->Fill(-1,weight); //-1, sonst --> normiert ergibt das gesamteffizienz
  }

*/





   //Effizienz fuer Sideband
  if(N_gen_ele == 1 && N_gen_mu == 1){
    bool reco_ele = false;
    bool reco_mu = false;
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11){ //gp = ele --> suche nach match in reco-eles 
	for(const auto & thisele: *event.electrons){
	  if(deltaR(gp,thisele) < 0.1){ //wurde rekonstruiert
	    reco_ele = true;
	  }
	}
	hist("pt_1ele")->Fill(gp.pt(),weight); //wahres pt aller geneles
	if(reco_ele) hist("Eff_1ele_pt")->Fill(gp.pt(),weight); //wahres pt der geneles, die auch rekonstruiert werden
      }
      else if(id == 13){
	for(const auto & thismu: *event.muons){
	  if(deltaR(gp,thismu) < 0.1){ //wurde rekonstruiert
	    reco_mu = true;
	  }
	}
	hist("pt_1mu")->Fill(gp.pt(),weight); //wahres pt aller genmus
	if(reco_mu) hist("Eff_1mu_pt")->Fill(gp.pt(),weight); //wahres pt der genmus, die auch rekonstruiert werden
      }
    }
    hist("H_T_EleMu")->Fill(ht,weight); //ht, wenn N_gen_ele == 1 && N_gen_mu == 1
    if(reco_ele && reco_mu) hist("Eff_elemu_ht")->Fill(ht,weight); //ht, wenn sie auch rekonstruiert wurde  -->ht_reco / ht = eff(ht)
    if(reco_ele && reco_mu) hist("Eff_elemu_int")->Fill(1,weight); // 1, wenn rekonstruktion erfolgreich
    else                    hist("Eff_elemu_int")->Fill(-1,weight);//-1, sonst --> normiert ergibt das gesamteffizienz

  }

  //Effizienz fuer Signal Region
  if(N_gen_mu == 2){
    int idx = 0;
    int idx_genp = 0;
    bool reco_mu1 = false;
    bool reco_mu2 = false;
    for(const auto & gp : *event.genparticles){
      idx = 0;
      int id = abs(gp.pdgId());
      if(id == 13){
	idx_genp++;
	//cout << "for gp " << idx_genp << ":" << endl;
	for(const auto & thismu: *event.muons){
	  idx ++;
	  if(deltaR(gp,thismu) < 0.1 && idx_genp == 1){ //wurde rekonstruiert
	    reco_mu1 = true;
	    //cout << "idx1: " << idx << endl;
	  }
	  if(deltaR(gp,thismu) < 0.1 && idx_genp == 2){
	    reco_mu2 = true;
	    //cout << "idx2: " << idx << endl;
	  }
	}
	hist("pt_2mu")->Fill(gp.pt(),weight); //wahres pt aller genmus
	if(reco_mu1 && reco_mu2) hist("Eff_2mu_pt")->Fill(gp.pt(),weight); //wahres pt der genmus, von denen beide auch rekonstruiert werden
      }
    }
    //cout << "1: " << reco_mu1 << ", 2: " << reco_mu2 << endl << endl;
    hist("H_T_MuMu")->Fill(ht,weight); //ht, wenn N_gen_mu == 2
    if(reco_mu1 && reco_mu2) hist("Eff_mumu_ht")->Fill(ht,weight); //ht, wenn sie auch rekonstruiert wurden  -->ht_reco / ht = eff(ht)
    if(reco_mu1 && reco_mu2) hist("Eff_mumu_int")->Fill(1,weight); // 1, wenn rekonstruktion erfolgreich
    else                     hist("Eff_mumu_int")->Fill(-1,weight);//-1, sonst --> normiert ergibt das gesamteffizienz

  }


  //Fakerate fuer Sideband
  if(Nele == 1 && Nmuons == 1){
    bool not_fake_ele = false;
    bool not_fake_mu = false;
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11){ //es ist auch ein ele vorhanden gewesen
	if(deltaR(event.electrons->at(0),gp) < 0.1){ // koennen sogar gematcht werden --> richtige reco
	  not_fake_ele = true;
	}
      }
      if(id == 13){ //es ist auch ein mu vorhanden gewesen
	if(deltaR(event.muons->at(0),gp) < 0.1){ // koennen sogar gematcht werden --> richtige reco
	  not_fake_mu = true;
	}
      }
    }
    if(not_fake_ele) hist("Eff_ele_pt_fakerate_side")->Fill(event.electrons->at(0).pt(),weight); //alle richtig rekonstruierten elektronen
                     hist("pt_1ele_fakerate_side")->Fill(event.electrons->at(0).pt(),weight);    //alle vorhandenen elektronen
    if(not_fake_mu)  hist("Eff_mu_pt_fakerate_side")->Fill(event.muons->at(0).pt(),weight);      //alle richtig rekonstruierten muonen
                     hist("pt_1mu_fakerate_side")->Fill(event.muons->at(0).pt(),weight);         //alle vorhandenen muonen

    if(not_fake_mu && not_fake_ele) hist("Eff_elemu_ht_fakerate")->Fill(ht,weight); 
                                    hist("H_T_EleMu_fakerate")->Fill(ht,weight); 

    if(not_fake_mu && not_fake_ele) hist("Eff_elemu_int_fakerate")->Fill(1,weight);  // 1, wenn rekonstruktion erfolgreich	       
    else                            hist("Eff_elemu_int_fakerate")->Fill(-1,weight); //-1, sonst --> normiert ergibt das gesamteffizienz

  }


  //Fakerate fuer Signal Region
  if(Nmuons == 2){
    bool not_fake_mu1 = false;
    bool not_fake_mu2 = false;
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 13){ //es ist auch ein mu vorhanden gewesen
	if(deltaR(event.muons->at(0),gp) < 0.1){ // koennen sogar gematcht werden --> richtige reco
	  not_fake_mu1 = true;
	}
	else if(deltaR(event.muons->at(1),gp) < 0.1){ //moeglich, dass ein gp das erste mu matcht und ein zweites beide mu's matcht und durch das 'else' hier nicht betrachtet wird...aber sehr unwahrscheinlich. dann gilt das zweite faelschlicherweise trotzdem als unmatched...aber sehr unwahrscheinlich
	  not_fake_mu2 = true;
	}
      }
    }
    if(not_fake_mu1) hist("Eff_mu_pt_fakerate_sig")->Fill(event.muons->at(0).pt(),weight); //alle richtig rekonstruierten elektronen
                     hist("pt_1mu_fakerate_sig")->Fill(event.muons->at(0).pt(),weight);    //alle vorhandenen elektronen
    if(not_fake_mu2) hist("Eff_mu_pt_fakerate_sig")->Fill(event.muons->at(1).pt(),weight);      //alle richtig rekonstruierten muonen
                     hist("pt_1mu_fakerate_sig")->Fill(event.muons->at(1).pt(),weight);         //alle vorhandenen muonen

    if(not_fake_mu1 && not_fake_mu2) hist("Eff_mumu_ht_fakerate")->Fill(ht,weight); 
                                     hist("H_T_MuMu_fakerate")->Fill(ht,weight); 

    if(not_fake_mu1 && not_fake_mu2) hist("Eff_mumu_int_fakerate")->Fill(1,weight);  // 1, wenn rekonstruktion erfolgreich	       
    else                             hist("Eff_mumu_int_fakerate")->Fill(-1,weight); //-1, sonst --> normiert ergibt das gesamteffizienz
  }







 } //is_mc
} //Methode


LQToTopMuEfficiencyHists::~LQToTopMuEfficiencyHists(){}

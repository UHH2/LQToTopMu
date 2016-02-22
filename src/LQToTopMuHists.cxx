#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuHists::LQToTopMuHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);  
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);
  book<TH1F>("pt_jets", "p_{T}^{jets} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3} [GeV]", 50, 20, 1500);
  book<TH1F>("N_bJets_loose", "#N_{Bjets}^{loose}", 10, 0, 10);
  book<TH1F>("N_bJets_med", "#N_{Bjets}^{medium}", 10, 0, 10);
  book<TH1F>("N_bJets_tight", "#N_{Bjets}^{tight}", 10, 0, 10);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV]", 50, 0, 1500);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1F>("M_mumu", "M_(#mu#mu) [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("Pt_mu_sum", "#Sum p_{T}^{#mu} [GeV]", 50, 0, 7000);
  double bins_HTlept_low[6] = {0, 300, 600, 900, 1200, 7000};
  book<TH1F>("Pt_mu_sum_rebin", "#Sum p_{T}^{#mu} [GeV]", 5, bins_HTlept_low);
  book<TH1F>("Pt_lept1", "leading lepton p_{T [GeV]}", 75, 0, 1500);
  book<TH1F>("Pt_lept2", "subleading lepton p_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_lept12", "leading lepton p_{T} + subleading lepton p_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_lept12_rebin", "leading lepton p_{T} + subleading lepton p_{T} [GeV]", 5,bins_HTlept_low);
  book<TH1F>("Pt_mu1", "p_{T}^{leading #mu} [GeV]", 75, 0, 1500);
  double bins_pt_low[26] = {0,30,60,90,120,150,180,210,240,270,300,350,400,450,500,550,600,650,700,750,800,900,1000,1100,1300,1500};
  book<TH1F>("Pt_mu1_rebin", "P_{T}^{leading #mu} [GeV]", 25, bins_pt_low);
  book<TH1F>("Pt_mu1_NoEle", "p_{T}^{leading #mu}, no Ele [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_mu1_NoEle_rebin", "P_{T}^{leading #mu}, no Ele [GeV]", 25, bins_pt_low);

  book<TH1F>("pt_ele", "p_{T}^{ele} [GeV] [GeV]", 50, 0, 1500);


  // general
  book<TH1F>("N_pv", "N_{PV}", 50, 0, 50);
  book<TH1F>("E_Tmiss", "missing E_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("E_Tmiss_0Ele2Mu", "missing E_{T} [GeV] for N_{e}=0, N_{#mu}=2", 75, 0,1500);
  book<TH1F>("H_T", "H_{T} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_from350", "H_{T} [GeV]", 24, 0, 4200);
  book<TH1F>("H_T_from350_rebin", "H_{T} [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_from350_rebin2", "H_{T} [GeV]", 12, 0, 4200);
  book<TH1F>("Parton_H_T", "H_{T} [GeV] on parton level", 80,0,7000);
  double bins_low_1Ele[12] = {0,350,500,700,900,1100,1300,1500,1750,2000,2500,7000};
  double bins_low_NoEle[23] = {0,200,350,500,650,800,950,1100,1250,1400,1550,1700,1850,2000,2150,2300,2450,2600,2750,2900,3050,3200,7000};
  double bins_low_NoEle2[11] = {0,350,500,650,800,950,1100,1250,1450,1750,2050};
  book<TH1F>("H_T_rebin", "H_{T} [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_rebin2", "H_{T} [GeV]", 100, 0, 7000);
  book<TH1F>("H_T_rebin3", "H_{T} [GeV]", 10,bins_low_NoEle2);
  book<TH1F>("H_T_jets", "H_{T}^{jets} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_lept", "H_{T}^{leptons} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_jets_rebin", "H_{T}^{jets} rebinned [GeV]", 5, bins_HTlept_low);
  book<TH1F>("H_T_lept_rebin", "H_{T}^{leptons} rebinned [GeV]", 5, bins_HTlept_low);
  book<TH1F>("H_T_comb_NoEle", "H_{T}, no Ele [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_NoEle_from350", "H_{T}, no Ele [GeV]", 24, 0, 4200);
  book<TH1F>("H_T_comb_NoEle_from350_rebin", "H_{T}, no Ele [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_comb_NoEle_from350_rebin2", "H_{T}, no Ele [GeV]", 12, 0, 4200);
  book<TH1F>("H_T_comb_NoEle_rebin", "H_{T}, no Ele [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_comb_NoEle_rebin2", "H_{T}, no Ele [GeV]", 10, bins_low_NoEle2);
  book<TH1F>("Integral_NoEle", "BinContent = sum(eventweights), NoEle", 1, 0.5, 1.5);
  book<TH1F>("H_T_comb_1Ele", "H_{T}, N_{Ele} #geq 1 [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_1Ele_from350", "H_{T}, N_{Ele} #geq 1 [GeV]", 24, 0, 4200);
  book<TH1F>("H_T_comb_1Ele_from350_rebin", "H_{T}, N_{Ele} #geq 1 [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_comb_1Ele_from350_rebin2", "H_{T}, N_{Ele} #geq 1 [GeV]", 12, 0, 4200);
  book<TH1F>("H_T_comb_1Ele_rebin", "H_{T}, N_{Ele} #geq 1 [GeV]", 11, bins_low_1Ele);
  book<TH1F>("H_T_comb_1Ele_rebin2", "H_{T}, N_{Ele} #geq 1, same binning as for N_{Ele} = 0 [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_comb_1Ele_rebin3", "H_{T}, N_{Ele} #geq 1, same binning as for N_{Ele} = 0 [GeV]", 10, bins_low_NoEle2);
  book<TH1F>("Integral_1Ele", "BinContent = sum(eventweights), 1Ele", 1, 0.5, 1.5);
  book<TH1F>("M_LQ_comb", "M_{LQ,mean} [GeV]", 60, 0, 3000);
  double bins_mlq_low[17] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
  book<TH1F>("M_LQ_comb_rebin", "M_{LQ,mean} [GeV]", 16, bins_mlq_low);
  double bins_mlq_low2[6] = {0,200,400,600,800,1000};
  book<TH1F>("M_LQ_comb_rebin2", "M_{LQ,mean} [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV]", 50, -500, 500);
  book<TH1F>("M_LQ_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep})/M_{LQ,mean} [GeV]", 50, -0.5, 0.5);
  book<TH1F>("M_LQLQ", "M_{LQLQ} [GeV]", 100, 0, 5000);
  double bins_mlqlq_low[12] = {0,400,500,600,700,800,900,1000,1250,1500,2000,5000};
  book<TH1F>("M_LQLQ_rebin", "M_{LQLQ} [GeV]", 11, bins_mlqlq_low);

  book<TH1F>("dR_toplep_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_tophad_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_toplep_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_muX", "1: t^{had} closer to #mu^{had}, -1: closer to #mu^{lep}", 3,-1.5, 1.5);
  //book<TH1F>("dummy","dummy",50,0,5);

 
  //electron fakes
  book<TH1F>("ele_type", "0 real ele, 1 fake ele", 2,-0.5,1.5);

  //event weights: sum
  book<TH1F>("sum_event_weights", "BinContent = sum(eventweights)", 1, 0.5, 1.5);

 

  //For MLQ reconstruction
  h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
  //h_hadr_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassHadronicLQReconstruction");
  m_discriminator_name ="Chi2";
  //m_discriminator_name ="CorrectMatch";

  is_mc = ctx.get("dataset_type") == "MC";

}


void LQToTopMuHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.

  double weight = event.weight;
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);

  for(unsigned int i=0; i<event.jets->size(); i++){
    hist("pt_jets")->Fill(jets->at(i).pt(),weight);
  }
  
  if(Njets>=1){
    hist("eta_jet1")->Fill(jets->at(0).eta(), weight);
  }
  if(Njets>=2){
    hist("eta_jet2")->Fill(jets->at(1).eta(), weight);
  }
  if(Njets>=3){
    hist("eta_jet3")->Fill(jets->at(2).eta(), weight);
  }
  if(Njets>=4){
    hist("eta_jet4")->Fill(jets->at(3).eta(), weight);
  }
  if(Njets>=1){
    hist("pt_jet1")->Fill(jets->at(0).pt(), weight);
  }
  if(Njets>=2){
    hist("pt_jet2")->Fill(jets->at(1).pt(), weight);
  }
  if(Njets>=3){
    hist("pt_jet3")->Fill(jets->at(2).pt(), weight);
  }

  //# b-jets
  std::vector<Jet> bjets_loose, bjets_med, bjets_tight;
  CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  CSVBTag Btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  CSVBTag Btag_tight = CSVBTag(CSVBTag::WP_TIGHT);


  for (unsigned int i =0; i<jets->size(); ++i) {
    if(Btag_loose(jets->at(i),event)) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_loose.push_back(jets->at(i));
    }
    if(Btag_medium(jets->at(i),event)) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_med.push_back(jets->at(i));
    }
    if(Btag_tight(jets->at(i),event)) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_tight.push_back(jets->at(i));
    }
  }

  int NbJets_loose = bjets_loose.size();
  hist("N_bJets_loose")->Fill(NbJets_loose,weight);
  int NbJets_med = bjets_med.size();
  hist("N_bJets_med")->Fill(NbJets_med,weight);
  int NbJets_tight = bjets_tight.size();
  hist("N_bJets_tight")->Fill(NbJets_tight,weight);


  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  double sum_mu_pt = 0;
  for (const Muon & thismu : *event.muons){
    hist("pt_mu")->Fill(thismu.pt(), weight);
    hist("eta_mu")->Fill(thismu.eta(), weight);
    hist("reliso_mu")->Fill(thismu.relIso(), weight);
    sum_mu_pt += thismu.pt();
  }
  hist("Pt_mu_sum")->Fill(sum_mu_pt,weight);
  hist("Pt_mu_sum_rebin")->Fill(sum_mu_pt,weight);
  hist("Pt_mu1")->Fill(event.muons->at(0).pt(), weight);
  hist("Pt_mu1_rebin")->Fill(event.muons->at(0).pt(), weight);


  double pt_lept1 = 0, pt_lept2 = 0;
  if(event.electrons->size() > 1){
    if(event.muons->at(0).pt()>event.electrons->at(0).pt()){
      pt_lept1 = event.muons->at(0).pt();
      if(event.muons->size() > 1){
	if(event.muons->at(1).pt()>event.electrons->at(0).pt()) pt_lept2 = event.muons->at(1).pt();
	else pt_lept2 = event.electrons->at(0).pt();
      }
      else pt_lept2 = event.electrons->at(0).pt();
    }
    else{
      pt_lept1 = event.electrons->at(0).pt();
      if(event.electrons->at(1).pt()>=event.muons->at(0).pt()) pt_lept2 = event.electrons->at(1).pt();
      else pt_lept2 = event.muons->at(0).pt();
    }
  }
 if(event.electrons->size() == 1){
    if(event.muons->at(0).pt()>event.electrons->at(0).pt()){
      pt_lept1 = event.muons->at(0).pt();
      if(event.muons->size() > 1){
	if(event.muons->at(1).pt()>event.electrons->at(0).pt()) pt_lept2 = event.muons->at(1).pt();
	else pt_lept2 = event.electrons->at(0).pt();
      }
      else pt_lept2 = event.electrons->at(0).pt();
    }
    else{
      pt_lept1 = event.electrons->at(0).pt();
      pt_lept2 = event.muons->at(0).pt();
    }
  }
  if(event.electrons->size() == 0){
   pt_lept1 = event.muons->at(0).pt();
   if(event.muons->size() > 1) pt_lept2 = event.muons->at(1).pt();
  }



  hist("Pt_lept1")->Fill(pt_lept1,weight);
  if(event.electrons->size()+event.muons->size() > 1)hist("Pt_lept2")->Fill(pt_lept2,weight);
  hist("Pt_lept12")->Fill(pt_lept1+pt_lept2,weight);
  hist("Pt_lept12_rebin")->Fill(pt_lept1+pt_lept2,weight);

  for (const Electron & thisele : *event.electrons){
    hist("pt_ele")->Fill(thisele.pt(), weight);

  }
  
  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

  //HT
  auto met = event.met->pt();
  hist("E_Tmiss")->Fill(met, weight);
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
  hist("H_T_jets")->Fill(ht_jets,weight);
  hist("H_T_lept")->Fill(ht_lep,weight);
  hist("H_T_jets_rebin")->Fill(ht_jets,weight);
  hist("H_T_lept_rebin")->Fill(ht_lep,weight);
  hist("H_T")->Fill(ht, weight);
  hist("H_T_from350")->Fill(ht, weight);
  hist("H_T_from350_rebin")->Fill(ht, weight);
  hist("H_T_from350_rebin2")->Fill(ht, weight);
  hist("H_T_rebin")->Fill(ht, weight);
  hist("H_T_rebin2")->Fill(ht,weight);
  if(ht <= 2000) hist("H_T_rebin3")->Fill(ht,weight);
  else hist("H_T_rebin3")->Fill(2000,weight);

  //partonlvl HT:
  if(is_mc){
    double partonHT = 0;
    constexpr const int invalid_daughter = (unsigned short)(-1);
    for(const auto & gp : *event.genparticles){
      if(gp.daughter1() != invalid_daughter || gp.daughter2() != invalid_daughter) continue;
      // if we are here, it means we have a final state particle.
      // Add to HT in cas it is a parton (quark -- including b but not top as tops are never final state particles -- or gluon -- or ele/mu -- or its respective neutrino).
      // Note that the exact HT definition depends on the madgraph configuration, but this
      // should cover the most common case.
      int id = abs(gp.pdgId());
      if((id >= 1 && id <= 5) || (id == 21) || (id>=11 && id <= 14)){
	partonHT += gp.pt();
      }
    }
    hist("Parton_H_T")->Fill(partonHT, weight);
  }

  // M_mumu Invariant Mass
  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j > i){
	M_mumu = (muons[i] + muons[j]).M();
	hist("M_mumu")->Fill(M_mumu, weight);	
      }
    }
  }



  //Fill HT, if Nele = 0, else
  //reconstruct MLQ and fill MLQmean
  int Nele = event.electrons->size();
  if(Nele == 0){
    hist("H_T_comb_NoEle")->Fill(ht, weight);
    hist("H_T_comb_NoEle_from350")->Fill(ht, weight);
    hist("H_T_comb_NoEle_from350_rebin")->Fill(ht, weight);
    hist("H_T_comb_NoEle_from350_rebin2")->Fill(ht, weight);
    hist("H_T_comb_NoEle_rebin")->Fill(ht, weight);
    if(ht <= 2000) hist("H_T_comb_NoEle_rebin2")->Fill(ht, weight);
    else hist("H_T_comb_NoEle_rebin2")->Fill(2000., weight);
    hist("Pt_mu1_NoEle")->Fill(event.muons->at(0).pt(), weight);
    hist("Pt_mu1_NoEle_rebin")->Fill(event.muons->at(0).pt(), weight);
    hist("Integral_NoEle")->Fill(1,weight);
  }
  
  //check for at least 1 muon pair with opposite charge
  bool charge_opposite = false;
  for(unsigned int i=0; i<event.muons->size(); i++){
    for(unsigned int j=0; j<event.muons->size(); j++){
      if(j>i){
	if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
	  charge_opposite = true;
	}
      }
    }
  }
  //if(charge_opposite) cout << "opposite charges detected" << endl;
  //else cout << "NO opposite charges detected!!" << endl;

  if(Nele >= 1 && event.muons->size() >= 2 && charge_opposite){   
    std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps); 
    const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );

    double mLQlep_rec = 0;
    double mLQhad_rec = 0;
    double mLQmed_rec = 0;
    double mLQdiff = 0;
    double mLQdiff_rel = 0;
    double mLQLQ = 0;

    if( (hyp->LQlep_v4()).isTimelike() ) {mLQlep_rec = (hyp->LQlep_v4()).M();}
    else {mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());}
    if( (hyp->LQhad_v4()).isTimelike() ) {mLQhad_rec = (hyp->LQhad_v4()).M();}
    else {mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());}
    

    mLQdiff = mLQhad_rec - mLQlep_rec;
    mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;
    hist("M_LQ_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin")->Fill(mLQmed_rec, weight);
    if(mLQmed_rec < 900)   hist("M_LQ_comb_rebin2")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_comb_rebin2")->Fill(900., weight);
    hist("M_LQ_diff")->Fill(mLQdiff, weight);
    mLQdiff_rel = mLQdiff / mLQmed_rec;
    hist("M_LQ_diff_rel")->Fill(mLQdiff_rel,weight);
    mLQLQ = (hyp->LQlep_v4()+hyp->LQhad_v4()).M();
    hist("M_LQLQ")->Fill(mLQLQ,weight);
    hist("M_LQLQ_rebin")->Fill(mLQLQ,weight);
    
    //deltaR between tops and associated muons
    double dR_lep = deltaR(hyp->toplep_v4(),hyp->mu_lep_v4());
    double dR_had = deltaR(hyp->tophad_v4(),hyp->mu_had_v4());
    double dR_lephad = deltaR(hyp->toplep_v4(),hyp->mu_had_v4());
    double dR_hadlep = deltaR(hyp->tophad_v4(),hyp->mu_lep_v4());

    hist("dR_toplep_mulep")->Fill(dR_lep,weight);
    hist("dR_tophad_muhad")->Fill(dR_had,weight);
    hist("dR_toplep_muhad")->Fill(dR_lephad,weight);
    hist("dR_tophad_mulep")->Fill(dR_hadlep,weight);
    if(dR_had<dR_hadlep)  hist("dR_tophad_muX")->Fill(1,weight);
    else                  hist("dR_tophad_muX")->Fill(-1,weight);

  }

  if(Nele >= 1){
    hist("H_T_comb_1Ele")->Fill(ht, weight);
    hist("H_T_comb_1Ele_from350")->Fill(ht, weight);
    hist("H_T_comb_1Ele_from350_rebin")->Fill(ht, weight);
    hist("H_T_comb_1Ele_from350_rebin2")->Fill(ht, weight);
    hist("H_T_comb_1Ele_rebin")->Fill(ht, weight);
    hist("H_T_comb_1Ele_rebin2")->Fill(ht, weight);
    if(ht <= 2000) hist("H_T_comb_1Ele_rebin3")->Fill(ht, weight);
    else hist("H_T_comb_1Ele_rebin3")->Fill(2000, weight);
    hist("Integral_1Ele")->Fill(1,weight);
  }

 
 



    hist("sum_event_weights")->Fill(1, weight);

} //Methode



LQToTopMuHists::~LQToTopMuHists(){}

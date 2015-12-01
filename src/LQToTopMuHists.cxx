#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
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
  book<TH1F>("pt_jet1", "#p_{T}^{jet 1} [GeV/c]", 100, 0, 3000);
  book<TH1F>("pt_jet2", "#p_{T}^{jet 2} [GeV/c]", 100, 0, 3000);
  book<TH1F>("pt_jet3", "#p_{T}^{jet 3} [GeV/c]", 100, 0, 3000);
  book<TH1F>("N_bJets_loose", "#N_{Bjets}^{loose}", 10, 0, 10);
  book<TH1F>("N_bJets_med", "#N_{Bjets}^{medium}", 10, 0, 10);
  book<TH1F>("N_bJets_tight", "#N_{Bjets}^{tight}", 10, 0, 10);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 50, 0, 1500);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1F>("M_mumu", "M_(#mu#mu) [GeV/c^{2}]",50 , 0, 1000);
  book<TH1F>("Pt_mu_sum", "#Sum p_{T}^{#mu} [GeV/c]", 50, 0, 7000);
  double bins_HTlept_low[6] = {0, 300, 600, 900, 1200, 7000};
  book<TH1F>("Pt_mu_sum_rebin", "#Sum p_{T}^{#mu} [GeV/c]", 5, bins_HTlept_low);
  book<TH1F>("Pt_lept1", "leading lepton p_{T [GeV/c]}", 75, 0, 1500);
  book<TH1F>("Pt_lept2", "subleading lepton p_{T} [GeV/c]", 75, 0, 1500);
  book<TH1F>("Pt_lept12", "leading lepton p_{T} + subleading lepton p_{T} [GeV/c]", 75, 0, 1500);
  book<TH1F>("Pt_lept12_rebin", "leading lepton p_{T} + subleading lepton p_{T} [GeV/c]", 5,bins_HTlept_low);
  book<TH1F>("Pt_mu1", "p_{T}^{leading #mu} [GeV/c]", 75, 0, 1500);
  double bins_pt_low[26] = {0,30,60,90,120,150,180,210,240,270,300,350,400,450,500,550,600,650,700,750,800,900,1000,1100,1300,1500};
  book<TH1F>("Pt_mu1_rebin", "P_{T}^{leading #mu} [GeV/c]", 25, bins_pt_low);
  book<TH1F>("Pt_mu1_NoEle", "p_{T}^{leading #mu}, no Ele [GeV/c]", 75, 0, 1500);
  book<TH1F>("Pt_mu1_NoEle_rebin", "P_{T}^{leading #mu}, no Ele [GeV/c]", 25, bins_pt_low);

  book<TH1F>("pt_ele", "p_{T}^{ele} [GeV/c] [GeV/c]", 50, 0, 1500);


  // general
  book<TH1F>("N_pv", "N_{PV}", 50, 0, 50);
  book<TH1F>("E_Tmiss", "missing E_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("E_Tmiss_0Ele2Mu", "missing E_{T} [GeV] for N_{e}=0, N_{#mu}=2", 75, 0,1500);
  book<TH1F>("H_T", "H_{T} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_from350", "H_{T} [GeV] (from 350)", 40, 0,7000);
  book<TH1F>("H_T_from350_rebin", "H_{T} [GeV] (from 350)", 80, 0,7000);
  book<TH1F>("Parton_H_T", "H_{T} [GeV] (from 350) on parton level", 80,0,7000);
  double bins_low_1Ele[12] = {0,350,500,700,900,1100,1300,1500,1750,2000,2500,7000};
  double bins_low_NoEle[23] = {0,200,350,500,650,800,950,1100,1250,1400,1550,1700,1850,2000,2150,2300,2450,2600,2750,2900,3050,3200,7000};
  book<TH1F>("H_T_rebin", "H_{T} [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_rebin2", "H_{T} [GeV]", 100, 0, 7000);
  book<TH1F>("H_T_jets", "H_{T}^{jets} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_lept", "H_{T}^{leptons} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_jets_rebin", "H_{T}^{jets} rebinned [GeV]", 5, bins_HTlept_low);
  book<TH1F>("H_T_lept_rebin", "H_{T}^{leptons} rebinned [GeV]", 5, bins_HTlept_low);
  book<TH1F>("H_T_comb_NoEle", "H_{T}, no Ele [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_NoEle_from350", "H_{T}, no Ele [GeV] (from 350)", 40, 0, 7000);
  book<TH1F>("H_T_comb_NoEle_from350_rebin", "H_{T}, no Ele [GeV] (from 350)", 80, 0, 7000);
  book<TH1F>("H_T_comb_NoEle_rebin", "H_{T}, no Ele [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_comb_1Ele", "H_{T}, N_{Ele} #geq 1 [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_1Ele_from350", "H_{T}, N_{Ele} #geq 1 [GeV] (from 350)", 40, 0, 7000);
  book<TH1F>("H_T_comb_1Ele_from350_rebin", "H_{T}, N_{Ele} #geq 1 [GeV] (from 350)", 80, 0, 7000);
  book<TH1F>("H_T_comb_1Ele_rebin", "H_{T}, N_{Ele} #geq 1 [GeV]", 11, bins_low_1Ele);
  book<TH1F>("H_T_comb_1Ele_rebin2", "H_{T}, N_{Ele} #geq 1, same binning as for N_{Ele} = 0 [GeV]", 22, bins_low_NoEle);
  book<TH1F>("M_LQ_comb", "M_{LQ,mean} [GeV/c^{2}]", 60, 0, 3000);
  double bins_mlq_low[17] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
  book<TH1F>("M_LQ_comb_rebin", "M_{LQ,mean} [GeV/c^{2}]", 16, bins_mlq_low);
  double bins_mlq_low2[6] = {100,200,300,500,800,2000};
  book<TH1F>("M_LQ_comb_rebin2", "M_{LQ,mean} [GeV/c^{2}]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV/c^{2}]", 50, -500, 500);
  book<TH1F>("M_LQ_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep})/M_{LQ,mean} [GeV/c]", 50, -0.5, 0.5);
  book<TH1F>("M_LQLQ", "M_{LQLQ} [GeV/c^{2}]", 100, 0, 5000);
  double bins_mlqlq_low[11] = {400,500,600,700,800,900,1000,1250,1500,2000,5000};
  book<TH1F>("M_LQLQ_rebin", "M_{LQLQ} [GeV/c^{2}]", 10, bins_mlqlq_low);

  book<TH1F>("dR_toplep_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_tophad_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_toplep_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_muX", "1: t^{had} closer to #mu^{had}, -1: closer to #mu^{lep}", 3,-1.5, 1.5);
  //book<TH1F>("dummy","dummy",50,0,5);

  //substructure related
  book<TH1F>("M_jet", "M_{Jet} [GeV/c^{2}]", 100, 0, 2000);
  book<TH1F>("N_subjets", "N_{Subjets} in a Topjet", 11, -0.5, 10.5);
  book<TH1F>("min_mDisubjet", "Min(m_{ij}) [GeV/c^{2}]", 50, 0, 1000);
  book<TH1F>("N_TopTags", "Number of CMSTopTags",16 ,-0.5, 15.5 );

  //electron fakes
  book<TH1F>("ele_type", "0 real ele, 1 fake ele", 2,-0.5,1.5);

  //event weights: sum
  book<TH1F>("sum_event_weights", "BinContent = sum(eventweights)", 1, 0.5, 1.5);

  /*  book <TH1F>("M_t_had", "M_{t,had} total", 50, 0, 500);
  book <TH1F>("M_t_had1", "M_{t1,had}", 50, 0, 500);
  book <TH1F>("M_t_had2", "M_{t2,had}", 50, 0, 500);

  book <TH1F>("M_LQ_had", "M_{LQ,had} total", 60, 0, 3000);
  book <TH1F>("M_LQ_had1", "M_{LQ1,had}", 60, 0, 3000);
  book <TH1F>("M_LQ_had2", "M_{LQ2,had}", 60, 0, 3000);
  book <TH1F>("M_LQ_had_mean", "M_{LQ,had}^{mean}", 60, 0, 3000);

  book <TH1F>("Chi2Had", "#chi^{2}_{had}", 100, 0, 500);
  book <TH1F>("Chi2Had_top1", "#chi^{2}_{had} top1", 150, 0, 300);
  book <TH1F>("Chi2Had_top2", "#chi^{2}_{had} top2", 150, 0, 300);
  book <TH1F>("Chi2Had_MLQdiff", "#chi^{2}_{had} M_{LQ}^{diff}", 150, 0, 300);
  book <TH1F>("Chi2Had_PTLQLQ", "#chi^{2}_{had} p_{T}^{LQLQ}", 150, 0, 300);
  book <TH1F>("Chi2Had_w1", "#chi^{2}_{had} w1", 100, 0, 500);
  book <TH1F>("Chi2Had_w2", "#chi^{2}_{had} w2", 100, 0, 500);
  book <TH1F>("Chi2Had_Top", "#chi^{2}_{had} from top", 200,0,1000);
  book <TH1F>("Chi2Had_LQ", "#chi^{2}_{had} from LQ", 150,0,300);
  book <TH2D>("Chi2Had_top1_vs_NJetsTop1","#chi^{2}_{had} top1 vs N_{Jets}^{top1}",150,0,300,10,-0.5,9.5);
  book <TH1F>("PTLQhadLQhad", "p_{T}^{LQ#bar{LQ}}", 100, 0, 1500);
  book <TH2D>("Chi2vsMLQ","#chi^{2} vs M_{LQ,had}^{mean}",150,0,300, 60,0,3000);
  book <TH2D>("PTvsMLQ","p_{T}^{LQ#bar{LQ}} vs M_{LQ,had}^{mean}",100,0,1500, 60,0,3000);
  book <TH2D>("Chi2vsPT","#chi^{2} vs p_{T}^{LQ#bar{LQ}}",150,0,300, 100,0,1500); */

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
  //auto ptmu2 = event.muons->at(1).pt();
  //double corrfactor_metreweight = (7.76-0.169*met+0.00128*met*met-0.00000352*met*met*met+0.00000000323*met*met*met*met) * 4.6087; // MET reweight
  //double corrfactor_ptmu2reweight = 0.08502+0.001746*ptmu2-0.000002253*ptmu2*ptmu2;
  double weight = event.weight;
  //if(met >= 100 && met <= 580){weight = event.weight * corrfactor_metreweight;}
  //weight = event.weight * corrfactor_ptmu2reweight;
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);
  
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
  for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.605) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_loose.push_back(jets->at(i));
    }
    if(jets->at(i).btag_combinedSecondaryVertex()>0.890) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_med.push_back(jets->at(i));
    }
    if(jets->at(i).btag_combinedSecondaryVertex()>0.970) { //loose: >0.605, medium: >0.890, tight: >0.970
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
  hist("H_T_rebin")->Fill(ht, weight);
  hist("H_T_rebin2")->Fill(ht,weight);

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
    hist("H_T_comb_NoEle_rebin")->Fill(ht, weight);
    hist("Pt_mu1_NoEle")->Fill(event.muons->at(0).pt(), weight);
    hist("Pt_mu1_NoEle_rebin")->Fill(event.muons->at(0).pt(), weight);
  }
  
  if(Nele >= 1 && event.muons->size() >= 2){   
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
    //mLQmed_rec = (1/hyp->discriminator("Chi2_tlep")*mLQhad_rec + 1/hyp->discriminator("Chi2_thad")*mLQlep_rec)/(2*(1/hyp->discriminator("Chi2_tlep")+1/hyp->discriminator("Chi2_thad"))/2);
    mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;
    hist("M_LQ_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin2")->Fill(mLQmed_rec, weight);
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
    hist("H_T_comb_1Ele_rebin")->Fill(ht, weight);
    hist("H_T_comb_1Ele_rebin2")->Fill(ht, weight);
  }

  //all-hadronic LQ reco
  /*if(Nele == 0 && event.muons->size() == 2){

    hist("E_Tmiss_0Ele2Mu")->Fill(met,weight);

    std::vector<LQReconstructionHypothesis> hadr_hyps = event.get(h_hadr_hyps); 
    const LQReconstructionHypothesis* hadr_hyp = get_best_hypothesis( hadr_hyps, "Chi2Hadronic" );

    double mTopHad1 = hadr_hyp->tophad1_v4().M();
    double mTopHad2 = hadr_hyp->tophad2_v4().M();
    double mLQHad1 = hadr_hyp->LQhad1_v4().M();
    double mLQHad2 = hadr_hyp->LQhad2_v4().M();
    double mLQHadMean = (mLQHad1 + mLQHad2)/2;
    double PTLQLQ = (hadr_hyp->LQhad1_v4()+ hadr_hyp->LQhad2_v4()).Pt();

    if(hadr_hyp->discriminator("Chi2Hadronic_thad1") <= 100 && hadr_hyp->discriminator("Chi2Hadronic_thad2") <= 100){
      hist("M_t_had1")->Fill(mTopHad1,weight);
      hist("M_t_had2")->Fill(mTopHad2,weight);
      hist("M_t_had")->Fill(mTopHad1,weight);
      hist("M_t_had")->Fill(mTopHad2,weight);

      hist("M_LQ_had1")->Fill(mLQHad1,weight);
      hist("M_LQ_had2")->Fill(mLQHad2,weight);
      hist("M_LQ_had")->Fill(mLQHad1,weight);
      hist("M_LQ_had")->Fill(mLQHad2,weight);
      hist("M_LQ_had_mean")->Fill(mLQHadMean,weight);
      hist("Chi2Had")->Fill(hadr_hyp->discriminator("Chi2Hadronic"), weight);
      hist("Chi2Had_top1")->Fill(hadr_hyp->discriminator("Chi2Hadronic_thad1"),weight);
      hist("Chi2Had_top2")->Fill(hadr_hyp->discriminator("Chi2Hadronic_thad2"),weight);
      hist("Chi2Had_MLQdiff")->Fill(hadr_hyp->discriminator("Chi2Hadronic_MLQdiff"),weight);
      hist("Chi2Had_PTLQLQ")->Fill(hadr_hyp->discriminator("Chi2Hadronic_PTLQLQ"),weight);
      //hist("Chi2Had_w1")->Fill(hadr_hyp->discriminator("Chi2Hadronic_whad1"),weight);
      //hist("Chi2Had_w2")->Fill(hadr_hyp->discriminator("Chi2Hadronic_whad2"),weight);
      hist("Chi2Had_Top")->Fill(hadr_hyp->discriminator("Chi2Hadronic_thad1")+hadr_hyp->discriminator("Chi2Hadronic_thad2"),weight);
      hist("Chi2Had_LQ")->Fill(hadr_hyp->discriminator("Chi2Hadronic_PTLQLQ")+hadr_hyp->discriminator("Chi2Hadronic_MLQdiff"),weight);
      hist("PTLQhadLQhad")->Fill(PTLQLQ,weight);
      ((TH2D*)hist("Chi2Had_top1_vs_NJetsTop1"))->Fill(hadr_hyp->discriminator("Chi2Hadronic_thad1"),hadr_hyp->tophad1_jets().size(),weight);
      ((TH2D*)hist("Chi2vsMLQ"))->Fill(hadr_hyp->discriminator("Chi2Hadronic"),mLQHadMean,weight);
      ((TH2D*)hist("PTvsMLQ"))->Fill(PTLQLQ,mLQHadMean,weight);
      ((TH2D*)hist("Chi2vsPT"))->Fill(hadr_hyp->discriminator("Chi2Hadronic"),PTLQLQ,weight);
    }

      if(hadr_hyp->tophad1_jets().size() == 0 || hadr_hyp->tophad2_jets().size() == 0) throw runtime_error("0 jets assinged to at least one top quark");
      
      
      if((100 < hadr_hyp->discriminator("Chi2Hadronic_thad1") ) || (100 < hadr_hyp->discriminator("Chi2Hadronic_thad2") ) ){
	cout << "Chi2_top1:    " << hadr_hyp->discriminator("Chi2Hadronic_thad1") << endl;
	cout << "Chi2_top2:    " << hadr_hyp->discriminator("Chi2Hadronic_thad2") << endl;
	cout << "Chi2_w1:      " << hadr_hyp->discriminator("Chi2Hadronic_whad1") << endl;
	cout << "Chi2_w2:      " << hadr_hyp->discriminator("Chi2Hadronic_whad2") << endl;
	cout << "Chi2_MLQdiff: " << hadr_hyp->discriminator("Chi2Hadronic_MLQdiff") << endl;
	cout << "Chi2_PTLQLQ:  " << hadr_hyp->discriminator("Chi2Hadronic_PTLQLQ") << endl << endl;

	cout << "Chi2_TOP:     " << hadr_hyp->discriminator("Chi2Hadronic_thad1")+hadr_hyp->discriminator("Chi2Hadronic_thad2")+hadr_hyp->discriminator("Chi2Hadronic_whad1")+hadr_hyp->discriminator("Chi2Hadronic_whad2") << endl;
	cout << "Chi2_LQ:      " << hadr_hyp->discriminator("Chi2Hadronic_PTLQLQ")+hadr_hyp->discriminator("Chi2Hadronic_MLQdiff") << endl << endl;

	cout << "Chi2_TOTAL:   " << hadr_hyp->discriminator("Chi2Hadronic") << endl;
	cout << "# Jets in t1: " << hadr_hyp->tophad1_jets().size() << endl;
	cout << "# Jets in t2: " << hadr_hyp->tophad2_jets().size() << endl;
	cout << "MTopHad1:     " << mTopHad1 << endl;
	cout << "MTopHad2:     " << mTopHad2 << endl << endl << endl << endl;
      }
      

      //cout << "TopMass1: " << mTopHad1 << endl;
      //cout << "TopMass2: " << mTopHad2 << endl;
      //}

  }*/
 
  //CMSTopTags

double mDiminLower = 50., mjetLower = 140., mjetUpper = 250.;
 //std::vector<TopJet>* topjets = event.topjets;
 //std::vector<TopJet> taggedtopjets;
 int N_toptaggedjets = 0;
 bool CMSTopTag = true;

 double m_disubjet_min = 0.;

 for(const auto & topjet : *event.topjets){

   auto subjets = topjet.subjets();
   
   if(subjets.size() < 2) m_disubjet_min = 0.0;
   
   // only need to sort if subjets there are more than 3 subjets, as
   // otherwise, we use all 3 anyway.
   if(subjets.size() > 3) sort_by_pt(subjets);
   
   double m01 = 0;

   if(subjets.size() >= 3){
     auto sum01 = subjets[0].v4()+subjets[1].v4();
     if(sum01.isTimelike())  m01 = sum01.M();
   
     if(subjets.size() < 3) m_disubjet_min = m01;
   
     double m02 = 0;
     auto sum02 = subjets[0].v4()+subjets[2].v4();
     if( sum02.isTimelike() ) m02 = sum02.M();
   
     double m12 = 0;
     auto sum12 = subjets[1].v4()+subjets[2].v4();
     if( sum12.isTimelike() )  m12 = sum12.M();
     m_disubjet_min = std::min(m01,std::min(m02, m12));
   }


   hist("min_mDisubjet")->Fill(m_disubjet_min, weight);
   if(m_disubjet_min < mDiminLower) CMSTopTag = false;

   auto mjet = topjet.v4().M();
   hist("M_jet")->Fill(mjet, weight);
   if(mjet < mjetLower) CMSTopTag = false;
   if(mjet > mjetUpper) CMSTopTag = false;

   hist("N_subjets")->Fill(subjets.size(), weight);
   if(subjets.size() < 3) CMSTopTag = false;

   //if (CMSTopTag) taggedtopjets.push_back(topjet); 
   if (CMSTopTag) N_toptaggedjets++; 
 }


 //int N_toptaggedjets = taggedtopjets.size();
 hist("N_TopTags")->Fill(N_toptaggedjets, weight);

 //fake electrons
 /* std::vector<double> VecdR;
 double dR = 1000;
  for(const auto & ele : *event.electrons){
    dR = 1000;
    for(auto genp : *event.genparticles){
      if(abs(genp.pdgId())==11){
	double tmp = deltaR(ele,genp);
	if(tmp<dR){
	  dR = tmp;
	}
      }
    }
    VecdR.push_back(dR);
  }

  for (unsigned int i=0; i<VecdR.size(); i++){
    if(VecdR.at(i) <= 0.3){
      hist("ele_type")->Fill(0);
    }
    else{
      hist("ele_type")->Fill(1);
    }
  }*/



    hist("sum_event_weights")->Fill(1, weight);

} //Methode



LQToTopMuHists::~LQToTopMuHists(){}

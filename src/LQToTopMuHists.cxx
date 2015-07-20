#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

#include "TH1F.h"
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
  book<TH1F>("pt_jet1", "#p_{T}^{jet 1}", 100, 0, 3000);
  book<TH1F>("pt_jet2", "#p_{T}^{jet 2}", 100, 0, 3000);
  book<TH1F>("pt_jet3", "#p_{T}^{jet 3}", 100, 0, 3000);
  book<TH1F>("N_bJets_loose", "#N_{Bjets}^{loose}", 10, 0, 10);
  book<TH1F>("N_bJets_med", "#N_{Bjets}^{medium}", 10, 0, 10);
  book<TH1F>("N_bJets_tight", "#N_{Bjets}^{tight}", 10, 0, 10);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 300, 0, 1500);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1F>("M_mumu", "M_(#mu#mu)",50 , 0, 1000);
  book<TH1F>("Pt_mu1", "p_{T}^{leading #mu}", 75, 0, 1500);
  double bins_pt_low[26] = {0,30,60,90,120,150,180,210,240,270,300,350,400,450,500,550,600,650,700,750,800,900,1000,1100,1300,1500};
  book<TH1F>("Pt_mu1_rebin", "P_{T}^{leading #mu}", 25, bins_pt_low);
  book<TH1F>("Pt_mu1_NoEle", "p_{T}^{leading #mu}, no Ele", 75, 0, 1500);
  book<TH1F>("Pt_mu1_NoEle_rebin", "P_{T}^{leading #mu}, no Ele", 25, bins_pt_low);


  // general
  book<TH1F>("N_pv", "N_{PV}", 50, 0, 50);
  book<TH1F>("E_Tmiss", "missing E_{T}", 75, 0, 1500);
  book<TH1F>("H_T", "H_{T}", 50, 0, 7000);
  double bins_low[12] = {0,350,500,700,900,1100,1300,1500,1750,2000,2500,7000};
  book<TH1F>("H_T_rebin", "H_{T}", 11, bins_low);
  book<TH1F>("H_T_jets", "H_{T}^{jets}", 50, 0, 7000);
  book<TH1F>("H_T_lept", "H_{T}^{leptons}", 50, 0, 7000);
  double bins_HTlept_low[6] = {0, 300, 600, 900, 1200, 7000};
  book<TH1F>("H_T_jets_rebin", "H_{T}^{jets} rebinned", 5, bins_HTlept_low);
  book<TH1F>("H_T_lept_rebin", "H_{T}^{leptons} rebinned", 5, bins_HTlept_low);

  book<TH1F>("H_T_comb_NoEle", "H_{T}, no Ele", 50, 0, 7000);
  book<TH1F>("H_T_comb_NoEle_rebin", "H_{T}, no Ele", 11, bins_low);
  book<TH1F>("H_T_comb_1Ele", "H_{T}, N_{Ele} #geq 1", 50, 0, 7000);
  book<TH1F>("H_T_comb_1Ele_rebin", "H_{T}, N_{Ele} #geq 1", 11, bins_low);
  book<TH1F>("M_LQ_comb", "M_{LQ,mean}", 40, 0, 2000);
  double bins_mlq_low[17] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
  book<TH1F>("M_LQ_comb_rebin", "M_{LQ,mean}", 16, bins_mlq_low);
  double bins_mlq_low2[6] = {100,200,300,500,800,2000};
  book<TH1F>("M_LQ_comb_rebin2", "M_{LQ,mean}", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_diff", "M_{LQ}^{had} - M_{LQ}^{lep}", 50, -500, 500);

  //substructure related
  book<TH1F>("M_jet", "M_{Jet}", 100, 0, 2000);
  book<TH1F>("N_subjets", "N_{Subjets} in a Topjet", 11, -0.5, 10.5);
  book<TH1F>("min_mDisubjet", "Min(m_{ij})", 50, 0, 1000);
  book<TH1F>("N_TopTags", "Number of CMSTopTags",16 ,-0.5, 15.5 );

  //book <TH1F>("M_t_had", "M_{t,had}", 50, 0, 500);
  //book <TH1F>("M_t_lep", "M_{t,lep}", 70, 0, 700);
  //book <TH1F>("M_ttbar", "M_{t#bar{t}}", 100, 0, 5000);

  //For MLQ reconstruction
  h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
  m_discriminator_name ="Chi2";
  //m_discriminator_name ="CorrectMatch";

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
    if(jets->at(i).btag_combinedSecondaryVertex()>0.423) { //loose: >0.423, medium: >0.814, tight: >0.914
      bjets_loose.push_back(jets->at(i));
    }
    if(jets->at(i).btag_combinedSecondaryVertex()>0.814) { //loose: >0.423, medium: >0.814, tight: >0.914
      bjets_med.push_back(jets->at(i));
    }
    if(jets->at(i).btag_combinedSecondaryVertex()>0.914) { //loose: >0.423, medium: >0.814, tight: >0.914
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
  for (const Muon & thismu : *event.muons){
    hist("pt_mu")->Fill(thismu.pt(), weight);
    hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }
  hist("Pt_mu1")->Fill(event.muons->at(0).pt(), weight);
  hist("Pt_mu1_rebin")->Fill(event.muons->at(0).pt(), weight);
  
  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

  //HT
  auto met = event.met->pt();
  hist("E_Tmiss")->Fill(met, weight);
  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    //ht += jet.pt();
    ht_jets += jet.pt();
  }

  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  /*for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
    }*/
  ht = ht_lep + ht_jets + met;
  hist("H_T_jets")->Fill(ht_jets,weight);
  hist("H_T_lept")->Fill(ht_lep,weight);
  hist("H_T_jets_rebin")->Fill(ht_jets,weight);
  hist("H_T_lept_rebin")->Fill(ht_lep,weight);
  hist("H_T")->Fill(ht, weight);
  hist("H_T_rebin")->Fill(ht, weight);

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

  //HT / MLQ Mix

  //Fill HT, if Nele = 0, else
  //reconstruct MLQ and fill MLQmean
  int Nele = event.electrons->size();
  if(Nele == 0){
    hist("H_T_comb_NoEle")->Fill(ht, weight);
    hist("H_T_comb_NoEle_rebin")->Fill(ht, weight);
    hist("Pt_mu1_NoEle")->Fill(event.muons->at(0).pt(), weight);
    hist("Pt_mu1_NoEle_rebin")->Fill(event.muons->at(0).pt(), weight);
  }
  
if(Nele >= 1){   
  std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps); 
  const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
 
    
    //Combine Top and Muon (Electron and Muon Charge have to be opposite)
    double mLQlep_rec = 0;
    double mLQhad_rec = 0;
    double mLQmed_rec = 0;
    double mLQdiff = 0;



  if( (hyp->LQlep_v4()).isTimelike() ) {mLQlep_rec = (hyp->LQlep_v4()).M();}
  else {mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());}
  if( (hyp->LQhad_v4()).isTimelike() ) {mLQhad_rec = (hyp->LQhad_v4()).M();}
  else {mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());}
    

    
   mLQmed_rec = (mLQhad_rec + mLQlep_rec)/2;
    hist("M_LQ_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin2")->Fill(mLQmed_rec, weight);
    hist("H_T_comb_1Ele")->Fill(ht, weight);
    hist("H_T_comb_1Ele_rebin")->Fill(ht, weight);
    hist("M_LQ_diff")->Fill(mLQdiff, weight);

 }
 
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

} //Methode


LQToTopMuHists::~LQToTopMuHists(){}

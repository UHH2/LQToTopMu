#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <stdexcept>


#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuHists::LQToTopMuHists(Context & ctx, const string & dirname, const bool & _isSR): Hists(ctx, dirname) {
  // book all histograms here
  // jets

  double bins_fele_pt[3] = {20,100,1500};
  double bins_pt_1[6]    = {20,80,140,200,260,800};
  double bins_pt_2[4]    = {20,100,200,800};
  double bins_eta[15]    = {-2.5, -2.0, -1.653, -1.305, -1.044, -0.783, -0.4, 0, 0.4, 0.783, 1.044, 1.305, 1.653, 2.0, 2.5};


  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);
  book<TH1F>("eta_jets", "#eta^{jets}", 25, -2.5, 2.5);
  book<TH1F>("eta_jets_rebin", "#eta^{jets}", 14, bins_eta);
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);
  book<TH1F>("pt_jets", "p_{T}^{jets} [GeV]", 2, bins_fele_pt);
  book<TH1F>("pt_jets_rebin", "p_{T}^{jets} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jets_rebin_2", "p_{T}^{jets} [GeV]", 5, bins_pt_1);
  book<TH1F>("pt_jets_rebin_3", "p_{T}^{jets} [GeV]", 3, bins_pt_2);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3} [GeV]", 50, 20, 1500);
  book<TH1F>("N_bJets_loose", "N_{Bjets}^{loose}", 11, -0.5, 10.5);
  book<TH1F>("N_bJets_med", "N_{Bjets}^{medium}", 11, -0.5, 10.5);
  book<TH1F>("N_bJets_tight", "N_{Bjets}^{tight}", 11, -0.5, 10.5);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_mu_zoom", "p_{T}^{#mu} [GeV]", 34, 0, 1020);
  book<TH1F>("eta_mu", "#eta^{#mu}", 25, -2.5, 2.5);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  book<TH1F>("M_mumu", "M_(#mu#mu) [GeV^{2}]",40 , 51, 131);
  book<TH1F>("M_mumu_full", "M_(#mu#mu) [GeV^{2}]",20 , 0, 1000);
  book<TH1F>("M_mumu_rebin1", "M_(#mu#mu) [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("M_mumu_rebin2", "M_(#mu#mu) [GeV^{2}]",100 , 0, 1000);

  book<TH1F>("M_eleele", "M_(#ele#ele) [GeV^{2}]",40 , 51, 131);
  book<TH1F>("M_eleele_full", "M_(#ele#ele) [GeV^{2}]",20 , 0, 1000);
  book<TH1F>("M_eleele_rebin1", "M_(#ele#ele) [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("M_eleele_rebin2", "M_(#ele#ele) [GeV^{2}]",100 , 0, 1000);

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

  book<TH1F>("pt_ele", "p_{T}^{ele} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_ele_zoom", "p_{T}^{ele} [GeV]", 34, 0, 1020);


  // general
  book<TH1F>("N_pv", "number of primary vertices", 51, -0.50, 50.5);
  book<TH1F>("N_pv_zoom", "number of primary vertices", 31, -0.50, 30.5);
  book<TH1F>("E_Tmiss", "missing E_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("E_Tmiss_0Ele2Mu", "missing E_{T} [GeV] for N_{e}=0, N_{#mu}=2", 75, 0,1500);
  book<TH1F>("H_T", "H_{T} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_from350", "H_{T} [GeV]", 24, 0, 4200);
  double bins_from350[21] = {0,175,350,525,700,875,1050,1225,1400,1575,1750,1925,2100,2275,2450,2625,2800,2975,3325,3675,4200}; //same binning as _from350 up to 2975, then two double-size and one triple-size bin
  book<TH1F>("H_T_from350_all_filled", "H_{T} [GeV]", 20, bins_from350);
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
  book<TH1F>("H_T_lept_zoom", "H_{T}^{leptons} [GeV]", 40, 0, 4000);
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
  book<TH1F>("M_LQ_comb_rebin3", "M_{LQ,mean} [GeV]", 30, 0, 3000);
  book<TH1F>("M_LQ_comb_rebin4", "M_{LQ,mean} [GeV]", 20, 0, 3000);
  double bins_mlq_low[17] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
  book<TH1F>("M_LQ_comb_rebin", "M_{LQ,mean} [GeV]", 16, bins_mlq_low);
  double bins_mlq_low2[6] = {0,200,400,600,800,1000};
  book<TH1F>("M_LQ_comb_rebin2", "M_{LQ,mean} [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV]", 50, -500, 500);
  book<TH1F>("M_LQ_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep})/M_{LQ,mean} [GeV]", 50, -0.5, 0.5);
  book<TH1F>("M_LQLQ", "M_{LQLQ} [GeV]", 100, 0, 5000);
  double bins_mlqlq_low[12] = {0,400,500,600,700,800,900,1000,1250,1500,2000,5000};
  book<TH1F>("M_LQLQ_rebin", "M_{LQLQ} [GeV]", 11, bins_mlqlq_low);
  book<TH1F>("N_jets_had","number of jets in the hadronic top hypothesis",6,-0.5,5.5);
  book<TH1F>("N_jets_lep","number of jets in the leptonic top hypothesis",6,-0.5,5.5);

  book<TH1F>("dR_toplep_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_tophad_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_toplep_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_muX", "1: t^{had} closer to #mu^{had}, -1: closer to #mu^{lep}", 3,-1.5, 1.5);
  book<TH1F>("dR_jet_ele1","#Delta R(Jet,Ele1)",240,0,3);
  book<TH1F>("dR_jet_ele2","#Delta R(Jet,Ele2)",240,0,3);
  book<TH1F>("dR_jet_ele3","#Delta R(Jet,Ele3)",240,0,3);
  book<TH1F>("dR_jet_fakeele","#Delta R(Jet,fake Ele)",240,0,3);

  book<TH1F>("M_fele_ele", "M_{fele, ele} 1 fele with every ele [GeV^{2}]",40 , 51, 131);
  book<TH1F>("M_fele_ele_full", "M_{fele, ele} 1 fele with every ele [GeV^{2}]",20 , 0, 1000);
  book<TH1F>("M_fele_ele_rebin1", "M_{fele, ele} 1 fele with every ele [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("M_fele_ele_rebin2", "M_{fele, ele} 1 fele with every ele [GeV^{2}]",100 , 0, 1000);

  book<TH1F>("M_no_fele", "M_{ele, ele} w/o fele [GeV^{2}]",40 , 51, 131);
  book<TH1F>("M_no_fele_full", "M_{ele, ele} w/o fele [GeV^{2}]",20 , 0, 1000);
  book<TH1F>("M_no_fele_rebin1", "M_{ele, ele} w/o fele [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("M_no_fele_rebin2", "M_{ele, ele} w/o fele [GeV^{2}]",100 , 0, 1000);
  //book<TH1F>("dummy","dummy",50,0,5);

 
  //electron fakes

  book<TH1F>("ele_type", "0 real ele, 1 fake ele", 2,-0.5,1.5);

  
  book<TH1F>("pt_fele_nogen_1", "p_{T}^{fake ele} no gen 1 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_1", "#eta^{fake ele} no gen 1", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_1", "p_{T}^{fake jet} no gen 1 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_1", "#eta^{fake jet} no gen 1", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_1", "#Delta R(jet,fake ele1) no gen 1", 250,0,3);

  book<TH1F>("pt_fele_nogen_inZ_1", "p_{T}^{fake ele} no gen but in Z 1 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_inZ_1", "#eta^{fake ele} no gen but in Z 1", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_inZ_1", "p_{T}^{fake jet} no gen but in Z 1 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_inZ_1", "#eta^{fake jet} no gen but in Z 1", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_inZ_1", "#Delta R(jet,fake ele1) no gen but in Z 1", 250,0,3);

  book<TH1F>("pt_fele_nogen_notZ_1", "p_{T}^{fake ele} no gen and not in Z 1 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_notZ_1", "#eta^{fake ele} no gen and not in Z 1", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_notZ_1", "p_{T}^{fake jet} no gen and not in Z 1 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_notZ_1", "#eta^{fake jet} no gen and not in Z 1", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_notZ_1", "#Delta R(jet,fake ele) no gen and not in Z 1", 250,0,3);

  book<TH1F>("fele_not_1", "+1 no gen and not in Z, -1 no gen but in Z 1", 3,-1.5,1.5);



  book<TH1F>("pt_fele_nogen_2", "p_{T}^{fake ele} no gen 2 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_2", "#eta^{fake ele} no gen 2", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_2", "p_{T}^{fake jet} no gen 2 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_2", "#eta^{fake jet} no gen 2", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_2", "#Delta R(jet,fake ele2) no gen 2", 250,0,3);

  book<TH1F>("pt_fele_nogen_inZ_2", "p_{T}^{fake ele} no gen but in Z 2 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_inZ_2", "#eta^{fake ele} no gen but in Z 2", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_inZ_2", "p_{T}^{fake jet} no gen but in Z 2 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_inZ_2", "#eta^{fake jet} no gen but in Z 2", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_inZ_2", "#Delta R(jet,fake ele) no gen but in Z 2", 250,0,3);

  book<TH1F>("pt_fele_nogen_notZ_2", "p_{T}^{fake ele} no gen and not in Z 2 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_notZ_2", "#eta^{fake ele} no gen and not in Z 2", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_notZ_2", "p_{T}^{fake jet} no gen and not in Z 2 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_notZ_2", "#eta^{fake jet} no gen and not in Z 2", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_notZ_2", "#Delta R(jet,fake ele) no gen and not in Z 2", 250,0,3);

  book<TH1F>("fele_not_2", "+1 no gen and not in Z, -1 no gen but in Z 2", 3,-1.5,1.5);



  book<TH1F>("pt_fele_nogen_3", "p_{T}^{fake ele} no gen 3 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_3", "#eta^{fake ele} no gen 3", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_3", "p_{T}^{fake jet} no gen 3 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_3", "#eta^{fake jet} no gen 3", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_3", "#Delta R(jet,fake ele) no gen 3", 250,0,3);

  book<TH1F>("pt_fele_nogen_inZ_3", "p_{T}^{fake ele} no gen but in Z 3 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_inZ_3", "#eta^{fake ele} no gen but in Z 3", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_inZ_3", "p_{T}^{fake jet} no gen but in Z 3 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_inZ_3", "#eta^{fake jet} no gen but in Z 3", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_inZ_3", "#Delta R(jet,fake ele) no gen but in Z 3", 250,0,3);

  book<TH1F>("pt_fele_nogen_notZ_3", "p_{T}^{fake ele} no gen and not in Z 3 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fele_nogen_notZ_3", "#eta^{fake ele} no gen and not in Z 3", 25,-2.5,2.5);
  book<TH1F>("pt_fjet_nogen_notZ_3", "p_{T}^{fake jet} no gen and not in Z 3 [GeV]", 2, bins_fele_pt);
  book<TH1F>("eta_fjet_nogen_notZ_3", "#eta^{fake jet} no gen and not in Z 3", 25,-2.5,2.5);
  book<TH1F>("dR_jetfele_nogen_notZ_3", "#Delta R(jet,fake ele) no gen and not in Z 1", 250,0,3);

  book<TH1F>("fele_not_3", "+1 no gen and not in Z, -1 no gen but in Z 3", 3,-1.5,1.5);

  book<TH1F>("hypo", "-2: Gen, best pair, -1: Gen, not best pair, +1: No Gen, best pair, +2: No Gen, not best pair", 5,-2.5,2.5);
  book<TH1F>("hypo2", "-2: Gen, best pair, -1: Gen, not best pair, +1: No Gen, best pair, +2: No Gen, not best pair", 5,-2.5,2.5);
  book<TH1F>("hypo3", "-2: Gen, best pair, -1: Gen, not best pair, +1: No Gen, best pair, +2: No Gen, not best pair", 5,-2.5,2.5);

  book<TH1F>("pt_fjet_nogen_notclose", "p_{T}^{faking jet} no gen and not in best pair [GeV]", 2, bins_fele_pt);
  book<TH1F>("pt_fjet_nogen_notclose_rebin", "p_{T}^{faking jet} no gen and not in best pair [GeV]", 50, 20, 1500);
  book<TH1F>("pt_fjet_nogen_notclose_rebin_2", "p_{T}^{faking jet} no gen and not in best pair [GeV]", 5, bins_pt_1);
  book<TH1F>("pt_fjet_nogen_notclose_rebin_3", "p_{T}^{faking jet} no gen and not in best pair [GeV]", 3, bins_pt_2);
  book<TH1F>("eta_fjet_nogen_notclose", "#eta^{faking jet} no gen and not in best pair", 25,-2.5,2.5);
  book<TH1F>("eta_fjet_nogen_notclose_rebin", "#eta^{faking jet} no gen and not in best pair", 14,bins_eta);

  book<TH1F>("pt_jet_notclose", "p_{T}^{jet} gen but not in best pair [GeV]", 2, bins_fele_pt);
  book<TH1F>("pt_jet_notclose_rebin", "p_{T}^{jet} gen but not in best pair [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet_notclose_rebin_2", "p_{T}^{jet} gen but not in best pair [GeV]", 5, bins_pt_1);
  book<TH1F>("pt_jet_notclose_rebin_3", "p_{T}^{jet} gen but not in best pair [GeV]", 3, bins_pt_2);
  book<TH1F>("eta_jet_notclose", "#eta^{jet} gen but not in best pair", 25,-2.5,2.5);
  book<TH1F>("eta_jet_notclose_rebin", "#eta^{jet} gen but not in best pair", 14,bins_eta);

  book<TH1F>("pt_jet_data", "p_{T}^{jet} not in best pair [GeV]", 2, bins_fele_pt);
  book<TH1F>("pt_jet_data_rebin", "p_{T}^{jet} not in best pair [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet_data_rebin_2", "p_{T}^{jet} not in best pair [GeV]", 5, bins_pt_1);
  book<TH1F>("pt_jet_data_rebin_3", "p_{T}^{jet} not in best pair [GeV]", 3, bins_pt_2);
  book<TH1F>("eta_jet_data", "#eta^{jet} not in best pair", 25,-2.5,2.5);
  book<TH1F>("eta_jet_data_rebin", "#eta^{jet} not in best pair", 14,bins_eta);

  book<TH1F>("2_ele", "1: 2 real ele, 2: 1 real & 1 fake ele , 3: 2 fake ele", 3,0.5,3.5);

  book<TH1F>("cos_phi", "cos(#Delta #phi)", 20, -1, 1);
  book<TH1F>("cos_phi_rebin", "cos(#Delta #phi)", 40, -1, 1);
  book<TH1F>("cos_phi_rebin2", "cos(#Delta #phi)", 80, -1, 1);

  book<TH1F>("m_t", "M_{T}", 10, 0, 400);
  book<TH1F>("m_t_rebin", "M_{T}", 40, 0, 400);
  book<TH1F>("m_t_rebin2", "M_{T}", 100, 0, 400);

  //event weights: sum
  book<TH1F>("sum_event_weights", "BinContent = sum(eventweights)", 1, 0.5, 1.5);

 

  //For MLQ reconstruction
  h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
  //h_hadr_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassHadronicLQReconstruction");
  m_discriminator_name ="Chi2";
  //m_discriminator_name ="CorrectMatch";


  h_jet_ele_ele =ctx.get_handle<vector<LorentzVector>>("h_jet_ele_ele");
  h_jet_ele_fele =ctx.get_handle<vector<LorentzVector>>("h_jet_ele_fele");
  h_jet_fele_fele =ctx.get_handle<vector<LorentzVector>>("h_jet_fele_fele");

  h_ele =ctx.get_handle<vector<LorentzVector>>("h_ele");
  h_jet =ctx.get_handle<vector<LorentzVector>>("h_jet");

  h_fele_nogen =ctx.get_handle<vector<LorentzVector>>("h_fele_nogen");
  h_fjet_nogen =ctx.get_handle<vector<LorentzVector>>("h_fjet_nogen");
  h_fele_nogen_inZ =ctx.get_handle<vector<LorentzVector>>("h_fele_nogen_inZ");
  h_fjet_nogen_inZ =ctx.get_handle<vector<LorentzVector>>("h_fjet_nogen_inZ");
  h_fele_nogen_notZ =ctx.get_handle<vector<LorentzVector>>("h_fele_nogen_notZ");
  h_fjet_nogen_notZ =ctx.get_handle<vector<LorentzVector>>("h_fjet_nogen_notZ");

  h_ele_close = ctx.get_handle<vector<LorentzVector>>("h_ele_close");
  h_fele_nogen_close = ctx.get_handle<vector<LorentzVector>>("h_fele_nogen_close");
  h_ele_notclose = ctx.get_handle<vector<LorentzVector>>("h_ele_notclose");
  h_fele_nogen_notclose = ctx.get_handle<vector<LorentzVector>>("h_fele_nogen_notclose");

  h_jet_close = ctx.get_handle<vector<LorentzVector>>("h_jet_close");
  h_fjet_nogen_close = ctx.get_handle<vector<LorentzVector>>("h_fjet_nogen_close");
  h_jet_notclose = ctx.get_handle<vector<LorentzVector>>("h_jet_notclose");
  h_fjet_nogen_notclose = ctx.get_handle<vector<LorentzVector>>("h_fjet_nogen_notclose");

  h_M_fele_ele =ctx.get_handle<vector<double>>("h_M_fele_ele");
  h_M_ele_ele_wofele =ctx.get_handle<vector<double>>("h_M_ele_ele_wofele");

  h_weights =ctx.get_handle<vector<double>>("h_weights");

  h_is_pt_weight =ctx.get_handle<bool>("h_is_pt_scale");

  is_mc = ctx.get("dataset_type") == "MC";  
  is_apply_SF = ctx.get("apply_SF") == "true";
  cross = ctx.get("cross_check") == "true";
 
  isSR = _isSR;
}

void LQToTopMuHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  double alt_weight = event.weight;

  vector<double> new_weights;
  if(event.is_valid(h_weights)) new_weights = event.get(h_weights);

  bool is_pt_weight=false;
  if(event.is_valid(h_is_pt_weight)) is_pt_weight = event.get(h_is_pt_weight);
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();

  hist("N_jets")->Fill(Njets, weight);

  for(unsigned int i=0; i<jets->size(); i++){
    if(is_mc && is_apply_SF) alt_weight = new_weights.at(3);
    
    if(is_mc && is_apply_SF && is_pt_weight && !cross){
       if(jets->at(i).pt() >= 20 && jets->at(i).pt() < 100) alt_weight = new_weights.at(3) * new_weights.at(0);
       if(jets->at(i).pt() >= 100 && jets->at(i).pt() < 200) alt_weight = new_weights.at(3) * new_weights.at(1);
       if(jets->at(i).pt() >= 200) alt_weight = new_weights.at(3) * new_weights.at(2);
    }
    hist("pt_jets")->Fill(jets->at(i).pt(),alt_weight);
    hist("pt_jets_rebin")->Fill(jets->at(i).pt(),alt_weight);
    hist("pt_jets_rebin_2")->Fill(jets->at(i).pt(),alt_weight);
    hist("pt_jets_rebin_3")->Fill(jets->at(i).pt(),alt_weight);
    hist("eta_jets")->Fill(jets->at(i).eta(),alt_weight);
    hist("eta_jets_rebin")->Fill(jets->at(i).eta(),alt_weight);
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


  // b-jets
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

  
  if(isSR){
    int Nmuons = event.muons->size();
    hist("N_mu")->Fill(Nmuons, weight);
    double sum_mu_pt = 0;
    for (const Muon & thismu : *event.muons){
      hist("pt_mu")->Fill(thismu.pt(), weight);
      hist("pt_mu_zoom")->Fill(thismu.pt(), weight);
      hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
      sum_mu_pt += thismu.pt();
    }
    hist("Pt_mu_sum")->Fill(sum_mu_pt,weight);
    hist("Pt_mu_sum_rebin")->Fill(sum_mu_pt,weight);
    hist("Pt_mu1")->Fill(event.muons->at(0).pt(), weight);
    hist("Pt_mu1_rebin")->Fill(event.muons->at(0).pt(), weight);
  }

  double pt_lept1 = 0, pt_lept2 = 0;
  if(isSR){
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
  }

  hist("Pt_lept1")->Fill(pt_lept1,weight);
  if(event.electrons->size()+event.muons->size() > 1)hist("Pt_lept2")->Fill(pt_lept2,weight);
  hist("Pt_lept12")->Fill(pt_lept1+pt_lept2,weight);
  hist("Pt_lept12_rebin")->Fill(pt_lept1+pt_lept2,weight);

  for (const Electron & thisele : *event.electrons){
    hist("pt_ele")->Fill(thisele.pt(), weight);
    hist("pt_ele_zoom")->Fill(thisele.pt(), weight);

  }
  
  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
  hist("N_pv_zoom")->Fill(Npvs, weight);

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
  
  if(isSR){
    for(const auto & muon : *event.muons){
      ht_lep += muon.pt();
    }
  }

  ht = ht_lep + ht_jets + met;

  hist("H_T_jets")->Fill(ht_jets,weight);
  hist("H_T_lept")->Fill(ht_lep,weight);
  hist("H_T_lept_zoom")->Fill(ht_lep,weight);
  hist("H_T_jets_rebin")->Fill(ht_jets,weight);
  hist("H_T_lept_rebin")->Fill(ht_lep,weight);
  hist("H_T")->Fill(ht, weight);
  hist("H_T_from350")->Fill(ht, weight);
  if(ht <= 4000) hist("H_T_from350_all_filled")->Fill(ht, weight);
  else hist("H_T_from350_all_filled")->Fill(4000, weight);
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
      // Add to HT in case it is a parton (quark -- including b but not top as tops are never final state particles -- or gluon -- or ele/mu -- or its respective neutrino).
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
  if(isSR){
    int Nmuons = event.muons->size();
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
	  hist("M_mumu_full")->Fill(M_mumu, weight);	
	  hist("M_mumu_rebin1")->Fill(M_mumu, weight);	
	  hist("M_mumu_rebin2")->Fill(M_mumu, weight);	
	}
      }
    }
  }

  // M_eleele Invariant Mass
  if(!isSR){
    int Neles = event.electrons->size();
    double M_eleele;
    LorentzVector electrons[Neles];
    for(int i=0; i<Neles; i++){
      electrons[i] = event.electrons->at(i).v4();
    }
    for(int i=0; i<Neles; i++){
      for(int j=0; j<Neles; j++){
	if(j > i){
	  M_eleele = (electrons[i] + electrons[j]).M();
	  hist("M_eleele")->Fill(M_eleele, weight);
	  hist("M_eleele_full")->Fill(M_eleele, weight);
	  hist("M_eleele_rebin1")->Fill(M_eleele, weight);
	  hist("M_eleele_rebin2")->Fill(M_eleele, weight);	
	}
      }
    }
  }

  vector<double> M_fele_ele, M_ele_ele_wofele;

  if(event.is_valid(h_M_fele_ele)) M_fele_ele = event.get(h_M_fele_ele);
  if(event.is_valid(h_M_ele_ele_wofele)) M_ele_ele_wofele = event.get(h_M_ele_ele_wofele);

  for(unsigned int i=0; i<M_fele_ele.size(); i++){
    hist("M_fele_ele")->Fill(M_fele_ele.at(i), weight);
    hist("M_fele_ele_full")->Fill(M_fele_ele.at(i), weight);
    hist("M_fele_ele_rebin1")->Fill(M_fele_ele.at(i), weight);
    hist("M_fele_ele_rebin2")->Fill(M_fele_ele.at(i), weight);
  }

  for(unsigned int i=0; i<M_ele_ele_wofele.size(); i++){
    hist("M_no_fele")->Fill(M_ele_ele_wofele.at(i), weight);
    hist("M_no_fele_full")->Fill(M_ele_ele_wofele.at(i), weight);
    hist("M_no_fele_rebin1")->Fill(M_ele_ele_wofele.at(i), weight);
    hist("M_no_fele_rebin2")->Fill(M_ele_ele_wofele.at(i), weight);
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
    if(isSR){
      hist("Pt_mu1_NoEle")->Fill(event.muons->at(0).pt(), weight);
      hist("Pt_mu1_NoEle_rebin")->Fill(event.muons->at(0).pt(), weight);
    }
    hist("Integral_NoEle")->Fill(1,weight);
  }
  
  //check for at least 1 muon pair with opposite charge
  
  bool charge_opposite = false;
  if(isSR){
    for(unsigned int i=0; i<event.muons->size(); i++){
      for(unsigned int j=0; j<event.muons->size(); j++){
	if(j>i){
	  if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
	    charge_opposite = true;
	  }
	}
      }
    }
  }
  
  //dR of Jet and ele
  for(unsigned int i=0; i<event.jets->size(); i++){
    int Neles = event.electrons->size();
    double M_eleele;
    LorentzVector electrons[Neles];
    for(int g=0; g<Neles; g++){
      electrons[g] = event.electrons->at(g).v4();
    }
    for(int j=0; j<Neles; j++){
      for(int k=0; k<Neles; k++){
	if(k > j){
	  M_eleele = (electrons[j] + electrons[k]).M();

	  double dR_jetele1;
	  dR_jetele1 = deltaR(jets->at(i), electrons[0]);
	  hist("dR_jet_ele1")->Fill(dR_jetele1,weight);

	  double dR_jetele2;
	  dR_jetele2 = deltaR(jets->at(i), electrons[1]);
	  hist("dR_jet_ele2")->Fill(dR_jetele2,weight);

	  double dR_jetele3;
	  dR_jetele3 = deltaR(jets->at(i), electrons[2]);
	  hist("dR_jet_ele3")->Fill(dR_jetele3,weight);

	  if(71 < M_eleele && M_eleele < 111){
	    for(int l=0; l<Neles; l++){
	      if(electrons[l] != electrons[j] && electrons[l] != electrons[k]){
		double dR_jetelefakeele;
		dR_jetelefakeele = deltaR(jets->at(i), electrons[l]);
		hist("dR_jet_fakeele")->Fill(dR_jetelefakeele,weight);
	      }
	    }
	  }
	}
      }
    }
  }


  //Hists for fake ele
  LorentzVector zero;

  vector<LorentzVector> fele_nogen, fjet_nogen, fele_nogen_inZ, fjet_nogen_inZ, fele_nogen_notZ, fjet_nogen_notZ;
  fele_nogen = fjet_nogen = fele_nogen_inZ = fjet_nogen_inZ = fele_nogen_notZ = fjet_nogen_notZ  = {{0,0,0,0}};

  int idx=0;

  if(event.is_valid(h_fele_nogen)) fele_nogen = event.get(h_fele_nogen);
  if(event.is_valid(h_fjet_nogen)) fjet_nogen = event.get(h_fjet_nogen);
  if(event.is_valid(h_fele_nogen_inZ)) fele_nogen_inZ = event.get(h_fele_nogen_inZ);
  if(event.is_valid(h_fjet_nogen_inZ)) fjet_nogen_inZ = event.get(h_fjet_nogen_inZ);
  if(event.is_valid(h_fele_nogen_notZ)) fele_nogen_notZ = event.get(h_fele_nogen_notZ);
  if(event.is_valid(h_fjet_nogen_notZ)) fjet_nogen_notZ = event.get(h_fjet_nogen_notZ);

  if(is_mc){
    //hist for nogen + check of the hypothesis
    for(unsigned int i=0; i<fele_nogen.size(); i++){
      switch(idx){
      case 0:
	if(fele_nogen.at(i) != zero && fjet_nogen.at(i) != zero){
	  hist("pt_fele_nogen_1")->Fill(fele_nogen.at(i).Pt(), weight);
	  hist("eta_fele_nogen_1")->Fill(fele_nogen.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_1")->Fill(fjet_nogen.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_1")->Fill(fjet_nogen.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_1")->Fill(deltaR(fjet_nogen.at(i), fele_nogen.at(i)), weight);
	}
	idx++;
	break;

      case 1:
	if(fele_nogen.at(i) != zero && fjet_nogen.at(i) != zero){
	  hist("pt_fele_nogen_2")->Fill(fele_nogen.at(i).Pt(), weight);
	  hist("eta_fele_nogen_2")->Fill(fele_nogen.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_2")->Fill(fjet_nogen.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_2")->Fill(fjet_nogen.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_2")->Fill(deltaR(fjet_nogen.at(i), fele_nogen.at(i)), weight);
	}
	idx++;
	break;

      case 2:
	if(fele_nogen.at(i) != zero && fjet_nogen.at(i) != zero){
	  hist("pt_fele_nogen_3")->Fill(fele_nogen.at(i).Pt(), weight);
	  hist("eta_fele_nogen_3")->Fill(fele_nogen.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_3")->Fill(fjet_nogen.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_3")->Fill(fjet_nogen.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_3")->Fill(deltaR(fjet_nogen.at(i), fele_nogen.at(i)), weight);
	}
	idx++;	  
	break;
      case 3:
	  throw std::runtime_error ("A 4th Fakeelektron was found! Please extend your Histogramms!");
	idx++;
	break;
      }
    }
      

    //hists for no ge-ele and in Z
    idx=0;
    for(unsigned int i=0; i<fele_nogen_inZ.size(); i++){
      switch(idx){
      case 0:
	if(fele_nogen_inZ.at(i) != zero && fjet_nogen_inZ.at(i) != zero){
	  hist("pt_fele_nogen_inZ_1")->Fill(fele_nogen_inZ.at(i).Pt(), weight);
	  hist("eta_fele_nogen_inZ_1")->Fill(fele_nogen_inZ.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_inZ_1")->Fill(fjet_nogen_inZ.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_inZ_1")->Fill(fjet_nogen_inZ.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_inZ_1")->Fill(deltaR(fjet_nogen_inZ.at(i), fele_nogen_inZ.at(i)), weight);
	}
	idx++;
	break;

      case 1:
	if(fele_nogen_inZ.at(i) != zero && fjet_nogen_inZ.at(i) != zero){
	  hist("pt_fele_nogen_inZ_2")->Fill(fele_nogen_inZ.at(i).Pt(), weight);
	  hist("eta_fele_nogen_inZ_2")->Fill(fele_nogen_inZ.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_inZ_2")->Fill(fjet_nogen_inZ.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_inZ_2")->Fill(fjet_nogen_inZ.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_inZ_2")->Fill(deltaR(fjet_nogen_inZ.at(i), fele_nogen_inZ.at(i)), weight);
	}
	idx++;
	break;

      case 2:
	if(fele_nogen_inZ.at(i) != zero && fjet_nogen_inZ.at(i) != zero){
	  hist("pt_fele_nogen_inZ_3")->Fill(fele_nogen_inZ.at(i).Pt(), weight);
	  hist("eta_fele_nogen_inZ_3")->Fill(fele_nogen_inZ.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_inZ_3")->Fill(fjet_nogen_inZ.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_inZ_3")->Fill(fjet_nogen_inZ.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_inZ_3")->Fill(deltaR(fjet_nogen_inZ.at(i), fele_nogen_inZ.at(i)), weight);
	}
	idx++;	  
	break;
      case 3:
	  throw std::runtime_error ("A 4th Fakeelektron was found! Please extend your Histogramms!");
	idx++;
	break;
      }
    }

    //hists for no ge-ele and not in Z
    idx=0;
    for(unsigned int i=0; i<fele_nogen_notZ.size(); i++){
      switch(idx){
      case 0:
	if(fele_nogen_notZ.at(i) != zero && fjet_nogen_notZ.at(i) != zero){
	  hist("pt_fele_nogen_notZ_1")->Fill(fele_nogen_notZ.at(i).Pt(), weight);
	  hist("eta_fele_nogen_notZ_1")->Fill(fele_nogen_notZ.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_notZ_1")->Fill(fjet_nogen_notZ.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_notZ_1")->Fill(fjet_nogen_notZ.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_notZ_1")->Fill(deltaR(fjet_nogen_notZ.at(i), fele_nogen_notZ.at(i)), weight);
	}
	idx++;
	break;

      case 1:
	if(fele_nogen_notZ.at(i) != zero && fjet_nogen_notZ.at(i) != zero){
	  hist("pt_fele_nogen_notZ_2")->Fill(fele_nogen_notZ.at(i).Pt(), weight);
	  hist("eta_fele_nogen_notZ_2")->Fill(fele_nogen_notZ.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_notZ_2")->Fill(fjet_nogen_notZ.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_notZ_2")->Fill(fjet_nogen_notZ.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_notZ_2")->Fill(deltaR(fjet_nogen_notZ.at(i), fele_nogen_notZ.at(i)), weight);
	}
	idx++;
	break;

      case 2:
	if(fele_nogen_notZ.at(i) != zero && fjet_nogen_notZ.at(i) != zero){
	  hist("pt_fele_nogen_notZ_3")->Fill(fele_nogen_notZ.at(i).Pt(), weight);
	  hist("eta_fele_nogen_notZ_3")->Fill(fele_nogen_notZ.at(i).Eta(), weight);
	  hist("pt_fjet_nogen_notZ_3")->Fill(fjet_nogen_notZ.at(i).Pt(), weight);
	  hist("eta_fjet_nogen_notZ_3")->Fill(fjet_nogen_notZ.at(i).Eta(), weight);
	  //cout << "In LQToTopMuHists: Filling dR between fake ele and fake jet: " << deltaR(fjet_Z, fele_Z) << endl;
	  hist("dR_jetfele_nogen_notZ_3")->Fill(deltaR(fjet_nogen_notZ.at(i), fele_nogen_notZ.at(i)), weight);
	}
	idx++;	  
	break;
      case 3:
	  throw std::runtime_error ("A 4th Fakeelektron was found! Please extend your Histogramms!");
	idx++;
	break;
      }
    }

    //Check hypothesis
    idx=0;
    for(unsigned int i=0; i<fele_nogen_notZ.size(); i++){
      switch(idx){
      case 0:
	if(fele_nogen_notZ.at(i)!= zero && fjet_nogen_notZ.at(i) != zero){
	  hist("fele_not_1")->Fill(1,weight);
	}
	if(fele_nogen_inZ.at(i)!= zero && fjet_nogen_inZ.at(i) != zero){
	  hist("fele_not_1")->Fill(-1,weight);
	}
	idx++;
	break;

      case 1:
	if(fele_nogen_notZ.at(i)!= zero && fjet_nogen_notZ.at(i) != zero){
	  hist("fele_not_2")->Fill(1,weight);
	}
	if(fele_nogen_inZ.at(i)!= zero && fjet_nogen_inZ.at(i) != zero){
	  hist("fele_not_2")->Fill(-1,weight);
	}
	idx++;
	break;

      case 2:
	if(fele_nogen_notZ.at(i)!= zero && fjet_nogen_notZ.at(i) != zero){
	  hist("fele_not_3")->Fill(1,weight);
	}
	if(fele_nogen_inZ.at(i)!= zero && fjet_nogen_inZ.at(i) != zero){
	  hist("fele_not_3")->Fill(-1,weight);
	}
	idx++;	  
	break;
      case 3:
	  throw std::runtime_error ("A 4th Fakeelektron was found! Please extend your Histogramms!");
	idx++;
	break;
      }
    }


    vector<LorentzVector> ele_close, fele_nogen_close, ele_notclose, fele_nogen_notclose;
    vector<LorentzVector> jet_close, fjet_nogen_close, jet_notclose, fjet_nogen_notclose;
    ele_close = fele_nogen_close = ele_notclose = fele_nogen_notclose = {{0,0,0,0}};
    jet_close = fjet_nogen_close = jet_notclose = fjet_nogen_notclose = {{0,0,0,0}};

    if(event.is_valid(h_ele_close)) ele_close = event.get(h_ele_close);
    if(event.is_valid(h_fele_nogen_close)) fele_nogen_close = event.get(h_fele_nogen_close);
    if(event.is_valid(h_ele_notclose)) ele_notclose = event.get(h_ele_notclose);
    if(event.is_valid(h_fele_nogen_notclose)) fele_nogen_notclose = event.get(h_fele_nogen_notclose);

    if(event.is_valid(h_jet_close)) jet_close = event.get(h_jet_close);
    if(event.is_valid(h_fjet_nogen_close)) fjet_nogen_close = event.get(h_fjet_nogen_close);
    if(event.is_valid(h_jet_notclose)) jet_notclose = event.get(h_jet_notclose);
    if(event.is_valid(h_fjet_nogen_notclose)) fjet_nogen_notclose = event.get(h_fjet_nogen_notclose);

    for(unsigned int i=0; i<ele_close.size(); i++){
      if(ele_close.at(i)!= zero){
	hist("hypo")->Fill(-2,weight);
      }
      if(fele_nogen_close.at(i)!= zero){
	hist("hypo")->Fill(1,weight);
      }
      if(ele_notclose.at(i)!= zero){
	hist("hypo")->Fill(-1,weight);
	hist("hypo2")->Fill(-1,weight);
	//match each ele which is not in the best combi with a jet within R=0.4 (this jet is responsible for faking that electron) 
	vector<Jet> faking_jets_notclose;
	double dr_min_notclose = 999999;
	int idx_matching_jet_notclose = -1;
	for(unsigned int j=0; j<event.jets->size(); j++){
	  double dr_tmp_notclose = deltaR(event.jets->at(j), ele_notclose.at(i));
	  if(dr_tmp_notclose < dr_min_notclose){
	    dr_min_notclose = dr_tmp_notclose;
	    idx_matching_jet_notclose = j;
	  }
	}
	//save jets corresponding to fake eles
	if(dr_min_notclose < 0.4){
	  faking_jets_notclose.push_back(event.jets->at(idx_matching_jet_notclose));
	  // cout << "dR between ele with genparticle and jet: " << dr_min_notclose << endl;
	  // cout << "Elektron w Genparticle and not close: " << ele_notclose.at(i) << endl;
	  // cout << endl;
	  hist("hypo3")->Fill(-1,weight);
	}
      }

      if(fele_nogen_notclose.at(i)!= zero){
	hist("hypo")->Fill(2,weight);
	hist("hypo2")->Fill(2,weight);

	//match each ele which is not in the best comb with a jet within R=0.4 (this jet is responsible for faking that electron) 
	vector<Jet> faking_jets_nogen_notclose;
	double dr_min_nogen_notclose = 999999;
	int idx_matching_jet_nogen_notclose = -1;
	for(unsigned int j=0; j<event.jets->size(); j++){
	  double dr_tmp_nogen_notclose = deltaR(event.jets->at(j), fele_nogen_notclose.at(i));
	  if(dr_tmp_nogen_notclose < dr_min_nogen_notclose){
	    dr_min_nogen_notclose = dr_tmp_nogen_notclose;
	    idx_matching_jet_nogen_notclose = j;
	  }
	}
	//save jets corresponding to fake eles
	if(dr_min_nogen_notclose < 0.4){
	  faking_jets_nogen_notclose.push_back(event.jets->at(idx_matching_jet_nogen_notclose));	
	  // cout << "dR between fake ele and jets: " << dr_min_nogen_notclose << endl;
	  // cout << "Elektron w/o Genparticle and not close: " << fele_nogen_notclose.at(i) << endl;
	  // cout << endl;
	  hist("hypo3")->Fill(2,weight);
	}
      }
    }


    for(unsigned int i=0; i<fele_nogen_notclose.size(); i++){
      if(is_mc && is_apply_SF) alt_weight = new_weights.at(3);

      if(fjet_nogen_notclose.at(i) != zero){
	if(is_mc && is_apply_SF && is_pt_weight){
	  if(fjet_nogen_notclose.at(i).Pt() >= 20 && fjet_nogen_notclose.at(i).Pt() < 100) alt_weight = new_weights.at(3) * new_weights.at(0);
	  if(fjet_nogen_notclose.at(i).Pt() >= 100 && fjet_nogen_notclose.at(i).Pt() < 200) alt_weight = new_weights.at(3) * new_weights.at(1);
	  if(fjet_nogen_notclose.at(i).Pt() >= 200) alt_weight = new_weights.at(3) * new_weights.at(2);
	}
	hist("pt_fjet_nogen_notclose")->Fill(fjet_nogen_notclose.at(i).Pt(), alt_weight);
	hist("pt_fjet_nogen_notclose_rebin")->Fill(fjet_nogen_notclose.at(i).Pt(), alt_weight);
	hist("pt_fjet_nogen_notclose_rebin_2")->Fill(fjet_nogen_notclose.at(i).Pt(), alt_weight);
	hist("pt_fjet_nogen_notclose_rebin_3")->Fill(fjet_nogen_notclose.at(i).Pt(), alt_weight);


	hist("eta_fjet_nogen_notclose")->Fill(fjet_nogen_notclose.at(i).Eta(), alt_weight);
	hist("eta_fjet_nogen_notclose_rebin")->Fill(fjet_nogen_notclose.at(i).Eta(), alt_weight);
      }
      if(jet_notclose.at(i) != zero){
	if(is_mc && is_apply_SF && is_pt_weight && !cross){
	  if(jet_notclose.at(i).Pt() >= 20 && jet_notclose.at(i).Pt() < 100) alt_weight = new_weights.at(3) * new_weights.at(0);
	  if(jet_notclose.at(i).Pt() >= 100 && jet_notclose.at(i).Pt() < 200) alt_weight = new_weights.at(3) * new_weights.at(1);
	  if(jet_notclose.at(i).Pt() >= 200) alt_weight = new_weights.at(3) * new_weights.at(2);
	}
	hist("pt_jet_notclose")->Fill(jet_notclose.at(i).Pt(), alt_weight);
	hist("pt_jet_notclose_rebin")->Fill(jet_notclose.at(i).Pt(), alt_weight);
	hist("pt_jet_notclose_rebin_2")->Fill(jet_notclose.at(i).Pt(), alt_weight);
	hist("pt_jet_notclose_rebin_3")->Fill(jet_notclose.at(i).Pt(), alt_weight);
	hist("eta_jet_notclose")->Fill(jet_notclose.at(i).Eta(), alt_weight);
	hist("eta_jet_notclose_rebin")->Fill(jet_notclose.at(i).Eta(), alt_weight);
      }
    }
  }


  vector<LorentzVector> jet, ele;
  jet = ele  = {{0,0,0,0}};

  if(event.is_valid(h_ele)) ele = event.get(h_ele);
  if(event.is_valid(h_jet)) jet = event.get(h_jet);

  for(unsigned int i=0; i<ele.size(); i++){
    if(is_mc && is_apply_SF) alt_weight = new_weights.at(3);
    if(ele.at(i) != zero){
      if(is_mc && is_apply_SF && is_pt_weight && !cross){
	if(jet.at(i).Pt() >= 20 && jet.at(i).Pt() < 100) alt_weight = new_weights.at(3) * new_weights.at(0);
	if(jet.at(i).Pt() >= 100 && jet.at(i).Pt() < 200) alt_weight = new_weights.at(3) * new_weights.at(1);
	if(jet.at(i).Pt() >= 200) alt_weight = new_weights.at(3) * new_weights.at(2);
      }
      hist("pt_jet_data")->Fill(jet.at(i).Pt(), alt_weight);
      hist("pt_jet_data_rebin")->Fill(jet.at(i).Pt(), alt_weight);
      hist("pt_jet_data_rebin_2")->Fill(jet.at(i).Pt(), alt_weight);
      hist("pt_jet_data_rebin_3")->Fill(jet.at(i).Pt(), alt_weight);
      hist("eta_jet_data")->Fill(jet.at(i).Eta(), alt_weight);
      hist("eta_jet_data_rebin")->Fill(jet.at(i).Eta(), alt_weight);
    }
  }

  vector<LorentzVector> jet_ele_ele, jet_ele_fele, jet_fele_fele;
  jet_ele_ele = jet_ele_fele = jet_fele_fele  = {{0,0,0,0}};

  if(event.is_valid(h_jet_ele_ele)) jet_ele_ele = event.get(h_jet_ele_ele);
  if(event.is_valid(h_jet_ele_fele)) jet_ele_fele = event.get(h_jet_ele_fele);
  if(event.is_valid(h_jet_fele_fele)) jet_fele_fele = event.get(h_jet_fele_fele);

  for(unsigned int i=0; i<jet_ele_ele.size(); i++){
    if(jet_ele_ele.at(i) != zero){
      hist("2_ele")->Fill(1, weight);
    }
    if(jet_ele_fele.at(i) != zero){
      hist("2_ele")->Fill(2, weight);
    }
    if(jet_fele_fele.at(i) != zero){
      hist("2_ele")->Fill(3, weight);
    }
  }


  double dphi, c_dphi, m_t;

  for(unsigned int i=0; i<ele.size(); i++){
    if(ele.at(i) != zero){
      for(unsigned int j=0; j<event.electrons->size(); j++){
	if(ele.at(i) == event.electrons->at(j).v4()){

	  dphi = deltaPhi(event.electrons->at(j),*event.met);
	  c_dphi = cos(dphi);
	  m_t = sqrt(2 * event.electrons->at(j).pt() * event.met->pt() * (1-c_dphi));

	  hist("cos_phi")->Fill(c_dphi, weight);
	  hist("cos_phi_rebin")->Fill(c_dphi, weight);
	  hist("cos_phi_rebin2")->Fill(c_dphi, weight);

	  hist("m_t")->Fill(m_t, weight);
	  hist("m_t_rebin")->Fill(m_t, weight);
	  hist("m_t_rebin2")->Fill(m_t, weight);
	}
      }
    }
  }




  //if(charge_opposite) cout << "opposite charges detected" << endl;
  //else cout << "NO opposite charges detected!!" << endl;

  hist("sum_event_weights")->Fill(1, weight);

  if(isSR) {
    if(!(Nele >= 1)) return;
  }
  else{
    if(!(Nele == 1)) return;
  }
  if((!isSR && Nele >= 1) || (isSR && event.muons->size() >= 2 && charge_opposite)) {   
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
    hist("M_LQ_comb_rebin3")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin4")->Fill(mLQmed_rec, weight);
    if(mLQmed_rec < 900)   hist("M_LQ_comb_rebin2")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_comb_rebin2")->Fill(900., weight);
    hist("M_LQ_diff")->Fill(mLQdiff, weight);
    mLQdiff_rel = mLQdiff / mLQmed_rec;
    hist("M_LQ_diff_rel")->Fill(mLQdiff_rel,weight);
    mLQLQ = (hyp->LQlep_v4()+hyp->LQhad_v4()).M();
    hist("M_LQLQ")->Fill(mLQLQ,weight);
    hist("M_LQLQ_rebin")->Fill(mLQLQ,weight);

    double n_jets_had = hyp->tophad_jets().size();
    double n_jets_lep = hyp->toplep_jets().size();
    hist("N_jets_had")->Fill(n_jets_had,weight);
    hist("N_jets_lep")->Fill(n_jets_lep,weight);
    
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

  }


} //Methode



LQToTopMuHists::~LQToTopMuHists(){}

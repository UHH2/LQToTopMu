#include "TFile.h"
#include "TTree.h"
#include <TH1F.h>

using namespace std;

void Fakerate_Check_HT350()
{
  //Fakerate of pt

  //setup canvas

  TCanvas *c1 = new TCanvas("c1","pt of Fake Jet / pt of Jets [MC & DATA] Rebin_3_newPre",200,10,800,600);
  TCanvas *c2 = new TCanvas("c2","Fakerate (MC)/ Fakerate (DATA) [pt] Rebin_3_newPre",200,10,800,600);
  // TCanvas *c3 = new TCanvas("c3","pt of Fake Jet / pt of Jets [MC & DATA] Check",200,10,800,600);
  // TCanvas *c4 = new TCanvas("c4","Fakerate (MC)/ Fakerate (DATA) [pt] Check",200,10,800,600);
  TCanvas *c5 = new TCanvas("c5","pt of Fake Jet / pt of Jets [MC] SCALE_newPre",200,10,800,600);
  TCanvas *c6 = new TCanvas("c6","pt of Fake Jet / pt of Jets [DATA] SCALE_newPre",200,10,800,600);
  TCanvas *c7 = new TCanvas("c7","Fakerate (DATA)/ Fakerate (MC) [pt] SCALE_newPre",200,10,800,600);
  TCanvas *sf = new TCanvas("sf","Scalefactor",200,10,800,600);

  // c1->cd();

  //load DATA
  TFile *DY = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.MC.DYJets.root", "READ");
  TFile *Diboson = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.MC.Diboson.root", "READ");
  TFile *QCD = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.MC.QCD.root", "READ");
  TFile *SingleTop = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.MC.SingleTop.root", "READ");
  TFile *TTbar = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.MC.TTbar.root", "READ");
  TFile *WJets = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.MC.WJets.root", "READ");

  TFile *DATA = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_check_HT350/uhh2.AnalysisModuleRunner.DATA.DATA.root", "READ");

  c1->cd();
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Rebin 3  
  //create hists
  // get pt hists of the jets which are the closest to the real electrons which are not part of the best pair, this needs to be substracted from DATA and MC
  TH1D* DY_real_ele_rebin3 = (TH1D*)DY->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* Diboson_real_ele_rebin3 = (TH1D*)Diboson->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* QCD_real_ele_rebin3 = (TH1D*)QCD->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* SingleTop_real_ele_rebin3 = (TH1D*)SingleTop->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* TTbar_real_ele_rebin3 = (TH1D*)TTbar->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* WJets_real_ele_rebin3 = (TH1D*)WJets->Get("FinalSelection/pt_jet_notclose_rebin_3");

  TH1D* DY_real_ele_check = (TH1D*) DY_real_ele_rebin3->Clone("DY_real_ele_check");
  TH1D* Diboson_real_ele_check = (TH1D*) Diboson_real_ele_rebin3->Clone("Diboson_real_ele_check");
  TH1D* QCD_real_ele_check = (TH1D*) QCD_real_ele_rebin3->Clone("QCD_real_ele_check");
  TH1D* SingleTop_real_ele_check = (TH1D*) SingleTop_real_ele_rebin3->Clone("SingleTop_real_ele_check");
  TH1D* TTbar_real_ele_check = (TH1D*) TTbar_real_ele_rebin3->Clone("TTbar_real_ele_check");
  TH1D* WJets_real_ele_check = (TH1D*) WJets_real_ele_rebin3->Clone("WJets_real_ele_check");

  TH1D* DY_real_ele_scale_up = (TH1D*) DY_real_ele_rebin3->Clone("DY_real_ele_scale_up");
  TH1D* Diboson_real_ele_scale_up = (TH1D*) Diboson_real_ele_rebin3->Clone("Diboson_real_ele_scale_up");
  TH1D* QCD_real_ele_scale_up = (TH1D*) QCD_real_ele_rebin3->Clone("QCD_real_ele_scale_up");
  TH1D* SingleTop_real_ele_scale_up = (TH1D*) SingleTop_real_ele_rebin3->Clone("SingleTop_real_ele_scale_up");
  TH1D* TTbar_real_ele_scale_up = (TH1D*) TTbar_real_ele_rebin3->Clone("TTbar_real_ele_scale_up");
  TH1D* WJets_real_ele_scale_up = (TH1D*) WJets_real_ele_rebin3->Clone("WJets_real_ele_scale_up");

  TH1D* DY_real_ele_scale_down = (TH1D*) DY_real_ele_rebin3->Clone("DY_real_ele_scale_down");
  TH1D* Diboson_real_ele_scale_down = (TH1D*) Diboson_real_ele_rebin3->Clone("Diboson_real_ele_scale_down");
  TH1D* QCD_real_ele_scale_down = (TH1D*) QCD_real_ele_rebin3->Clone("QCD_real_ele_scale_down");
  TH1D* SingleTop_real_ele_scale_down = (TH1D*) SingleTop_real_ele_rebin3->Clone("SingleTop_real_ele_scale_down");
  TH1D* TTbar_real_ele_scale_down = (TH1D*) TTbar_real_ele_rebin3->Clone("TTbar_real_ele_scale_down");
  TH1D* WJets_real_ele_scale_down = (TH1D*) WJets_real_ele_rebin3->Clone("WJets_real_ele_scale_down");

  // get pt hists of all jets for Nele >= 3
  TH1D* h_DY_fele_pt_rebin3 = (TH1D*)DY->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_Diboson_fele_pt_rebin3 = (TH1D*)Diboson->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_QCD_fele_pt_rebin3 = (TH1D*)QCD->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_SingleTop_fele_pt_rebin3 = (TH1D*)SingleTop->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_TTbar_fele_pt_rebin3 = (TH1D*)TTbar->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_WJets_fele_pt_rebin3 = (TH1D*)WJets->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_DATA_fele_pt_rebin3 = (TH1D*)DATA->Get("FinalSelection/pt_jet_data_rebin_3");

  TH1D* h_DY_fele_pt_check = (TH1D*) h_DY_fele_pt_rebin3->Clone("h_DY_fele_pt_check");
  TH1D* h_Diboson_fele_pt_check = (TH1D*) h_Diboson_fele_pt_rebin3->Clone("h_Diboson_fele_pt_check");
  TH1D* h_QCD_fele_pt_check = (TH1D*) h_QCD_fele_pt_rebin3->Clone("h_QCD_fele_pt_check");
  TH1D* h_SingleTop_fele_pt_check = (TH1D*) h_SingleTop_fele_pt_rebin3->Clone("h_SingleTop_fele_pt_check");
  TH1D* h_TTbar_fele_pt_check = (TH1D*) h_TTbar_fele_pt_rebin3->Clone("h_TTbar_fele_pt_check");
  TH1D* h_WJets_fele_pt_check = (TH1D*) h_WJets_fele_pt_rebin3->Clone("h_WJets_fele_pt_check");
  TH1D* h_DATA_fele_pt_check = (TH1D*) h_DATA_fele_pt_rebin3->Clone("h_DATA_fele_pt_check");

  TH1D* h_DY_fele_pt_scale_up = (TH1D*) h_DY_fele_pt_rebin3->Clone("h_DY_fele_pt_scale_up");
  TH1D* h_Diboson_fele_pt_scale_up = (TH1D*) h_Diboson_fele_pt_rebin3->Clone("h_Diboson_fele_pt_scale_up");
  TH1D* h_QCD_fele_pt_scale_up = (TH1D*) h_QCD_fele_pt_rebin3->Clone("h_QCD_fele_pt_scale_up");
  TH1D* h_SingleTop_fele_pt_scale_up = (TH1D*) h_SingleTop_fele_pt_rebin3->Clone("h_SingleTop_fele_pt_scale_up");
  TH1D* h_TTbar_fele_pt_scale_up = (TH1D*) h_TTbar_fele_pt_rebin3->Clone("h_TTbar_fele_pt_scale_up");
  TH1D* h_WJets_fele_pt_scale_up = (TH1D*) h_WJets_fele_pt_rebin3->Clone("h_WJets_fele_pt_scale_up");
  TH1D* h_DATA_fele_pt_scale_up = (TH1D*) h_DATA_fele_pt_rebin3->Clone("h_DATA_fele_pt_scale_up");

  TH1D* h_DY_fele_pt_scale_down = (TH1D*) h_DY_fele_pt_rebin3->Clone("h_DY_fele_pt_scale_down");
  TH1D* h_Diboson_fele_pt_scale_down = (TH1D*) h_Diboson_fele_pt_rebin3->Clone("h_Diboson_fele_pt_scale_down");
  TH1D* h_QCD_fele_pt_scale_down = (TH1D*) h_QCD_fele_pt_rebin3->Clone("h_QCD_fele_pt_scale_down");
  TH1D* h_SingleTop_fele_pt_scale_down = (TH1D*) h_SingleTop_fele_pt_rebin3->Clone("h_SingleTop_fele_pt_scale_down");
  TH1D* h_TTbar_fele_pt_scale_down = (TH1D*) h_TTbar_fele_pt_rebin3->Clone("h_TTbar_fele_pt_scale_down");
  TH1D* h_WJets_fele_pt_scale_down = (TH1D*) h_WJets_fele_pt_rebin3->Clone("h_WJets_fele_pt_scale_down");
  TH1D* h_DATA_fele_pt_scale_down = (TH1D*) h_DATA_fele_pt_rebin3->Clone("h_DATA_fele_pt_scale_down");

  // get pt hists of all jets for Nele >= 2
  TH1D* h_DY_jet_pt_rebin3 = (TH1D*)DY->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_Diboson_jet_pt_rebin3 = (TH1D*)Diboson->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_QCD_jet_pt_rebin3 = (TH1D*)QCD->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_SingleTop_jet_pt_rebin3 = (TH1D*)SingleTop->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_TTbar_jet_pt_rebin3 = (TH1D*)TTbar->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_WJets_jet_pt_rebin3 = (TH1D*)WJets->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_DATA_jet_pt_rebin3 = (TH1D*)DATA->Get("0bJetLoose/pt_jets_rebin_3");

  TH1D* h_DY_jet_pt_check = (TH1D*) h_DY_jet_pt_rebin3->Clone("h_DY_jet_pt_check");
  TH1D* h_Diboson_jet_pt_check = (TH1D*) h_Diboson_jet_pt_rebin3->Clone("h_Diboson_jet_pt_check");
  TH1D* h_QCD_jet_pt_check = (TH1D*) h_QCD_jet_pt_rebin3->Clone("h_QCD_jet_pt_check");
  TH1D* h_SingleTop_jet_pt_check = (TH1D*) h_SingleTop_jet_pt_rebin3->Clone("h_SingleTop_jet_pt_check");
  TH1D* h_TTbar_jet_pt_check = (TH1D*) h_TTbar_jet_pt_rebin3->Clone("h_TTbar_jet_pt_check");
  TH1D* h_WJets_jet_pt_check = (TH1D*) h_WJets_jet_pt_rebin3->Clone("h_WJets_jet_pt_check");
  TH1D* h_DATA_jet_pt_check = (TH1D*) h_DATA_jet_pt_rebin3->Clone("h_DATA_jet_pt_check");

  TH1D* h_DY_jet_pt_scale_up = (TH1D*) h_DY_jet_pt_rebin3->Clone("h_DY_jet_pt_scale_up");
  TH1D* h_Diboson_jet_pt_scale_up = (TH1D*) h_Diboson_jet_pt_rebin3->Clone("h_Diboson_jet_pt_scale_up");
  TH1D* h_QCD_jet_pt_scale_up = (TH1D*) h_QCD_jet_pt_rebin3->Clone("h_QCD_jet_pt_scale_up");
  TH1D* h_SingleTop_jet_pt_scale_up = (TH1D*) h_SingleTop_jet_pt_rebin3->Clone("h_SingleTop_jet_pt_scale_up");
  TH1D* h_TTbar_jet_pt_scale_up = (TH1D*) h_TTbar_jet_pt_rebin3->Clone("h_TTbar_jet_pt_scale_up");
  TH1D* h_WJets_jet_pt_scale_up = (TH1D*) h_WJets_jet_pt_rebin3->Clone("h_WJets_jet_pt_scale_up");
  TH1D* h_DATA_jet_pt_scale_up = (TH1D*) h_DATA_jet_pt_rebin3->Clone("h_DATA_jet_pt_scale_up");

  TH1D* h_DY_jet_pt_scale_down = (TH1D*) h_DY_jet_pt_rebin3->Clone("h_DY_jet_pt_scale_down");
  TH1D* h_Diboson_jet_pt_scale_down = (TH1D*) h_Diboson_jet_pt_rebin3->Clone("h_Diboson_jet_pt_scale_down");
  TH1D* h_QCD_jet_pt_scale_down = (TH1D*) h_QCD_jet_pt_rebin3->Clone("h_QCD_jet_pt_scale_down");
  TH1D* h_SingleTop_jet_pt_scale_down = (TH1D*) h_SingleTop_jet_pt_rebin3->Clone("h_SingleTop_jet_pt_scale_down");
  TH1D* h_TTbar_jet_pt_scale_down = (TH1D*) h_TTbar_jet_pt_rebin3->Clone("h_TTbar_jet_pt_scale_down");
  TH1D* h_WJets_jet_pt_scale_down = (TH1D*) h_WJets_jet_pt_rebin3->Clone("h_WJets_jet_pt_scale_down");
  TH1D* h_DATA_jet_pt_scale_down = (TH1D*) h_DATA_jet_pt_rebin3->Clone("h_DATA_jet_pt_scale_down");

  //for adding hists
  TH1D* h_sum_MC_fele_pt_rebin3 = (TH1D*)DY->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_sum_MC_jet_pt_rebin3 = (TH1D*)DY->Get("0bJetLoose/pt_jets_rebin_3");

  TH1D* h_sum_MC_fele_pt_check = (TH1D*) h_sum_MC_fele_pt_rebin3->Clone("h_sum_MC_fele_pt_check");
  TH1D* h_sum_MC_jet_pt_check = (TH1D*) h_sum_MC_jet_pt_rebin3->Clone("h_sum_MC_jet_pt_check");

  TH1D* h_sum_MC_fele_pt_scale_up = (TH1D*) h_sum_MC_fele_pt_rebin3->Clone("h_sum_MC_fele_pt_scale_up");
  TH1D* h_sum_MC_jet_pt_scale_up = (TH1D*) h_sum_MC_jet_pt_rebin3->Clone("h_sum_MC_jet_pt_scale_up");

  TH1D* h_sum_MC_fele_pt_scale_down = (TH1D*) h_sum_MC_fele_pt_rebin3->Clone("h_sum_MC_fele_pt_scale_down");
  TH1D* h_sum_MC_jet_pt_scale_down = (TH1D*) h_sum_MC_jet_pt_rebin3->Clone("h_sum_MC_jet_pt_scale_down");


  gStyle->SetOptStat(0);

  //Add pt hists for the jets of the electrons w/o a gen particle (MC) hists for Nele >= 3
  h_sum_MC_fele_pt_rebin3->Add(h_QCD_fele_pt_rebin3);
  h_sum_MC_fele_pt_rebin3->Add(h_Diboson_fele_pt_rebin3);
  h_sum_MC_fele_pt_rebin3->Add(h_SingleTop_fele_pt_rebin3);
  h_sum_MC_fele_pt_rebin3->Add(h_TTbar_fele_pt_rebin3);
  h_sum_MC_fele_pt_rebin3->Add(h_WJets_fele_pt_rebin3);

  //Add all MC hists for Nele >= 2
  h_sum_MC_jet_pt_rebin3->Add(h_QCD_jet_pt_rebin3);
  h_sum_MC_jet_pt_rebin3->Add(h_Diboson_jet_pt_rebin3);
  h_sum_MC_jet_pt_rebin3->Add(h_SingleTop_jet_pt_rebin3);
  h_sum_MC_jet_pt_rebin3->Add(h_TTbar_jet_pt_rebin3);
  h_sum_MC_jet_pt_rebin3->Add(h_WJets_jet_pt_rebin3);

  //Substract the pt of jets which belong to real electrons, Nele >= 3
  cout << endl;
  cout << "h_DATA before subtracting Diboson: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;
  h_DATA_fele_pt_rebin3->Add(Diboson_real_ele_rebin3, -1);
  cout << "h_DATA after subtracting Diboson: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;
  h_DATA_fele_pt_rebin3->Add(DY_real_ele_rebin3, -1);
  cout << "h_DATA after subtracting DY: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;
  h_DATA_fele_pt_rebin3->Add(QCD_real_ele_rebin3, -1);
  cout << "h_DATA after subtracting QCD: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;
  h_DATA_fele_pt_rebin3->Add(SingleTop_real_ele_rebin3, -1);
  cout << "h_DATA after subtracting SingleTop: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;
  h_DATA_fele_pt_rebin3->Add(TTbar_real_ele_rebin3, -1);
  cout << "h_DATA after subtracting TTbar: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;
  h_DATA_fele_pt_rebin3->Add(WJets_real_ele_rebin3, -1);
  cout << "h_DATA after subtracting WJets: Bin 1: " << h_DATA_fele_pt_rebin3->GetBinContent(1) << " Bin 2: " << h_DATA_fele_pt_rebin3->GetBinContent(2) << " Bin 3: " << h_DATA_fele_pt_rebin3->GetBinContent(3) << endl;

  cout << endl;
  //divide fake ele pt and jet pt
  TGraphAsymmErrors* ratio_MC = new TGraphAsymmErrors(h_sum_MC_fele_pt_rebin3, h_sum_MC_jet_pt_rebin3);
  h_sum_MC_fele_pt_rebin3->Divide(h_sum_MC_jet_pt_rebin3);
  TH1D* dummy_MC = (TH1D*) h_sum_MC_fele_pt_rebin3->Clone("dummy_MC");
  ratio_MC->GetXaxis()->SetRangeUser(20,800);
  ratio_MC->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_MC->GetYaxis()->SetRangeUser(0.000001,0.0014);
  ratio_MC->GetYaxis()->SetTitle("#epsilon");
  ratio_MC->SetMarkerColor(kBlack);
  ratio_MC->SetLineColor(kBlack);
  dummy_MC->SetMarkerColor(kBlack);
  dummy_MC->SetLineColor(kBlack);
  ratio_MC->SetMarkerStyle(8);
  dummy_MC->SetMarkerStyle(8);
  ratio_MC->SetMarkerSize(0.6);
  dummy_MC->SetMarkerSize(0.6);
  ratio_MC->SetLineWidth(2);
  dummy_MC->SetLineWidth(2);
  ratio_MC->SetTitle("p_{t} fakerate [MC & DATA]");
  ratio_MC->Draw("ap");

  //h_DATA_fele_pt_rebin3->Divide(h_DATA_jet_pt_rebin3);
  TGraphAsymmErrors* ratio_DATA = new TGraphAsymmErrors(h_DATA_fele_pt_rebin3, h_DATA_jet_pt_rebin3);
  h_DATA_fele_pt_rebin3->Divide(h_DATA_jet_pt_rebin3);
  TH1D* dummy_DATA = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA");
  ratio_DATA->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA->GetYaxis()->SetRangeUser(0.000001,0.0014);
  ratio_DATA->SetMarkerColor(kGreen+1);
  ratio_DATA->SetLineColor(kGreen+1);
  dummy_DATA->SetMarkerColor(kGreen+1);
  dummy_DATA->SetLineColor(kGreen+1);
  ratio_DATA->SetMarkerStyle(25);
  dummy_DATA->SetMarkerStyle(25);
  ratio_DATA->SetMarkerSize(0.6);
  dummy_DATA->SetMarkerSize(0.6);
  ratio_DATA->SetLineWidth(2);
  dummy_DATA->SetLineWidth(2);
  ratio_DATA->Draw("e1psame");

  leg3 = new TLegend(0.639098,0.810105,0.803258,0.874564);
  leg3->AddEntry(dummy_MC,"MC","lep");
  leg3->AddEntry(dummy_DATA,"DATA","lep");
  leg3->Draw();
  
  c2->cd();
  //h_DATA_fele_pt_rebin3->Divide(h_sum_MC_fele_pt_rebin3);
  

  TGraphAsymmErrors* ratio_DATA_MC = new TGraphAsymmErrors(dummy_DATA, dummy_MC, "pois");
  TH1D* dummy_DATA_MC = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_MC");
  cout << endl;
  cout << "Used option 'pois' for the DATA/MC plot, can you use it?" << endl;
  cout << endl;
  ratio_DATA_MC->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC->GetYaxis()->SetRangeUser(0,3.5);
  ratio_DATA_MC->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_MC->SetMarkerColor(kBlue+1);
  ratio_DATA_MC->SetLineColor(kBlue+1);
  dummy_DATA_MC->SetMarkerColor(kBlue+1);
  dummy_DATA_MC->SetLineColor(kBlue+1);
  ratio_DATA_MC->SetMarkerStyle(8);
  dummy_DATA_MC->SetMarkerStyle(8);
  ratio_DATA_MC->SetMarkerSize(0.6);
  dummy_DATA_MC->SetMarkerSize(0.6);
  ratio_DATA_MC->SetLineWidth(2);
  dummy_DATA_MC->SetLineWidth(2);
  ratio_DATA_MC->SetTitle("p_{t} fakerate [DATA/MC]");
  ratio_DATA_MC->Draw("ap");

  leg4 = new TLegend(0.714286,0.824042,0.884712,0.881533);
  leg4->AddEntry(dummy_DATA_MC,"DATA/MC","lep");
  leg4->Draw();
  
  //c3->cd();
  
  ratio_DATA_MC_nice = (TGraphErrors*)ratio_DATA_MC->Clone("ratio_DATA_MC_nice");

  ratio_DATA_nice = (TGraphErrors*)ratio_DATA->Clone("ratio_DATA_nice");
  ratio_MC_nice = (TGraphErrors*)ratio_MC->Clone("ratio_MC_nice");

  // lets save the graphs
  // c1->SaveAs("fakerate_data_and_mc.eps");
  // c2->SaveAs("fakerate_data_mc.eps");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* ***********  Nice Ratio plots from fits ********** */
  gStyle->SetOptTitle(0);

  TH1D* dummy_MC_nice = (TH1D*) h_sum_MC_fele_pt_rebin3->Clone("dummy_MC_nice");
  TH1D* dummy_DATA_nice = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_nice");

  dummy_MC_nice->SetLineWidth(2);
  dummy_DATA_nice->SetLineWidth(2);

  dummy_MC_nice->SetMarkerStyle(25);
  dummy_DATA_nice->SetMarkerStyle(8);

  dummy_MC_nice->SetMarkerSize(0.6);
  dummy_DATA_nice->SetMarkerSize(0.6);

  dummy_MC_nice->SetMarkerColor(kGreen+1);
  dummy_MC_nice->SetLineColor(kGreen+1);

  dummy_DATA_nice->SetMarkerColor(kBlack);
  dummy_DATA_nice->SetLineColor(kBlack);


  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 400;
  CanHeight = 400;

  TCanvas* cfitsigold = new TCanvas("cfitsigold", "Fake rate", CanWidth, CanHeight);


  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //
  Float_t y1, y2, y3;                           //  y3 +-------------+
  y3 = 0.99;                                    //     |             |
  y2 = y3-yplot;                                //     | pad1        |
  y1 = y2-yratio;                               //  y2 |-------------|
  Float_t x1, x2;                               //     | rp1         |
  x1 = 0.01;                                    //  y1 +-------------+
  x2 = 0.99;                                    //     x1            x2
                                                //
                                                // No Pad 2!


  TPad* m_rp1_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  TPad* m_rp1 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  TPad* m_pad1 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);

  TPad* m_rp2_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  TPad* m_rp2 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  TPad* m_pad2 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);


  m_pad1->SetTopMargin(0.05); m_pad1->SetBottomMargin(0.16);  m_pad1->SetLeftMargin(0.19); m_pad1->SetRightMargin(0.05);
  m_pad2->SetTopMargin(0.05); m_pad2->SetBottomMargin(0.16);  m_pad2->SetLeftMargin(0.19); m_pad2->SetRightMargin(0.05);

  m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.0); m_rp1_top->SetLeftMargin(0.19); m_rp1_top->SetRightMargin(0.05);
  m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.0); m_rp2_top->SetLeftMargin(0.19); m_rp2_top->SetRightMargin(0.05);
  m_rp1->SetTopMargin(0.0); m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.19); m_rp1->SetRightMargin(0.05);
  m_rp2->SetTopMargin(0.0); m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.19); m_rp2->SetRightMargin(0.05);

  // m_rp1_top->SetFillColor(kYellow);
  // m_rp2_top->SetFillColor(kOrange);
  // m_rp1->SetFillColor(kGray);
  // m_rp2->SetFillColor(kGray);

  m_pad1->Draw();
  m_rp1_top->Draw();
  m_rp1->Draw();

  m_rp1_top->cd();


  TPaveText *textbox9 = new TPaveText(0.6,0.77,0.75,0.87,"NDC");
  textbox9->SetFillColor(0);
  textbox9->SetLineColor(0);
  TText *line9a = textbox9->AddText("Control region");
  line9a->SetTextColor(1);
  line9a->SetTextAlign(12);//12
  line9a->SetTextFont(43);
  line9a->SetTextSizePixels(20);
  textbox9->SetBorderSize(1);
  textbox9->Draw();



  gPad->SetTicks(1,1);
  gPad->RedrawAxis();
  m_rp1->cd();

  //calculate ratios
  gPad->SetTicks(1,1);
  // ratio_DATA_MC->SetTitle("ratio of fake rates");
  ratio_DATA_MC_nice->Draw("ap");
  ratio_DATA_MC_nice->GetYaxis()->SetTitle("DATA/MC");
  ratio_DATA_MC_nice->GetYaxis()->SetTitleSize(0.07);
  ratio_DATA_MC_nice->GetYaxis()->SetLabelSize(0.08);
  ratio_DATA_MC_nice->GetYaxis()->SetTickLength(0.015);
  ratio_DATA_MC_nice->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC_nice->GetXaxis()->SetTitleSize(0.1);
  ratio_DATA_MC_nice->GetXaxis()->SetLabelSize(0.1);
  ratio_DATA_MC_nice->GetXaxis()->SetTickLength(0.055);
  ratio_DATA_MC_nice->SetMarkerColor(kBlack);
  ratio_DATA_MC_nice->SetLineColor(kBlack);
  ratio_DATA_MC_nice->GetYaxis()->SetTitleOffset(0.9);
  ratio_DATA_MC_nice->GetYaxis()->CenterTitle(true);
  ratio_DATA_MC_nice->GetYaxis()->SetRangeUser(0,3.4);

  TLine *line = new TLine(20,1,800,1);
  line->SetLineColor(14);
  line->SetLineStyle(2);
  line->Draw();
                                        
  blub = (TGraphErrors*)ratio_DATA_MC->Clone("blub");

  sf->cd();
  double x_SF,sf_1,sf_2,sf_3;
  ratio_DATA_MC->GetPoint(0,x_SF,sf_1);
  cout << "Scalefactor for Bin 1: " << sf_1<< endl;
  ratio_DATA_MC->GetPoint(1,x_SF,sf_2);
  cout << "Scalefactor for Bin 2: " << sf_2<< endl;
  ratio_DATA_MC->GetPoint(2,x_SF,sf_3);
  cout << "Scalefactor for Bin 3: " << sf_3<< endl;


  const Int_t n = 3;
  Double_t x[n]   = {60, 150, 500};
  Double_t y[n]   = {sf_1, sf_2, sf_3};
  Double_t exl[n] = {0, 0, 0};
  Double_t eyl[n] = {0, 0, 0};
  Double_t exh[n] = {0, 0, 0};
  Double_t eyh[n] = {0, 0, 0};
  gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(1.3);
  gr->GetYaxis()->SetTickLength(0.015);
  gr->GetXaxis()->SetTickLength(0.055);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
  gr->SetLineWidth(2);
  gr->GetYaxis()->SetTitleOffset(0.65);
  gr->GetYaxis()->SetRangeUser(0,3.4);
  gr->Draw("e1p same");


  gStyle->SetEndErrorSize(5);
  gPad->SetTicks(1,1);
  blub->GetYaxis()->SetTitle("DATA/MC");
  //blub->GetYaxis()->SetTitleSize(0.07);                                                                                                                                                                   
  //blub->GetYaxis()->SetLabelSize(0.055);                                                                                                                                                                  
  blub->GetYaxis()->SetTickLength(0.015);
  blub->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  //blub->GetXaxis()->SetTitleSize(0.1);                                                                                                                                                                    
  //blub->GetXaxis()->SetLabelSize(0.07);                                                                                                                                                                   
  blub->GetXaxis()->SetTickLength(0.055);
  blub->SetMarkerColor(kBlack);
  blub->SetLineColor(kBlack);
  blub->GetYaxis()->SetTitleOffset(1);
  blub->GetYaxis()->SetRangeUser(0,3.4);
  blub->Draw("ap");
  gr->Draw("e1p same");
  line->Draw();

  //sf->SaveAs("Fake_Check.eps");

  m_rp1_top->cd();

  ratio_DATA_nice->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_nice->GetYaxis()->SetLabelSize(0.04);
  ratio_MC_nice->GetYaxis()->SetLabelSize(0.04);
  ratio_DATA_nice->GetYaxis()->SetTitleSize(0.07);
  ratio_DATA_nice->GetYaxis()->SetTitleOffset(0.9);
  ratio_DATA_nice->Draw("ap");
  ratio_MC_nice->Draw("e1p same");

  leg_nice = new TLegend(0.746429,0.727759,0.90998,0.888294);
  leg_nice->AddEntry(dummy_MC_nice,"MC","lep");
  leg_nice->AddEntry(dummy_DATA_nice,"DATA","lep");
  leg_nice->SetLineWidth(0);
  leg_nice->Draw();



  cfitsigold->SaveAs("check.eps");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // scale Diboson up and down (50%)
  gStyle->SetOptTitle(0);
  c5->cd();

  gStyle->SetOptTitle(0);

  //scale the Diboson plots
  Diboson_real_ele_scale_up->Scale(1.1118);
  h_Diboson_fele_pt_scale_up->Scale(1.1118);
  h_Diboson_jet_pt_scale_up->Scale(1.1118);

  Diboson_real_ele_scale_down->Scale(0.8882);
  h_Diboson_fele_pt_scale_down->Scale(0.8882);
  h_Diboson_jet_pt_scale_down->Scale(0.8882);


  //Add pt hists for the jets of the electrons w/o a gen particle (MC) hists for Nele >= 3
  h_sum_MC_fele_pt_scale_up->Add(h_QCD_fele_pt_scale_up);
  h_sum_MC_fele_pt_scale_up->Add(h_Diboson_fele_pt_scale_up);
  h_sum_MC_fele_pt_scale_up->Add(h_SingleTop_fele_pt_scale_up);
  h_sum_MC_fele_pt_scale_up->Add(h_TTbar_fele_pt_scale_up);
  h_sum_MC_fele_pt_scale_up->Add(h_WJets_fele_pt_scale_up);

  //Add all MC hists for Nele >= 2
  h_sum_MC_jet_pt_scale_up->Add(h_QCD_jet_pt_scale_up);
  h_sum_MC_jet_pt_scale_up->Add(h_Diboson_jet_pt_scale_up);
  h_sum_MC_jet_pt_scale_up->Add(h_SingleTop_jet_pt_scale_up);
  h_sum_MC_jet_pt_scale_up->Add(h_TTbar_jet_pt_scale_up);
  h_sum_MC_jet_pt_scale_up->Add(h_WJets_jet_pt_scale_up);

  //Substract the pt of jets which belong to real electrons, Nele >= 3
  h_DATA_fele_pt_scale_up->Add(Diboson_real_ele_scale_up, -1);
  h_DATA_fele_pt_scale_up->Add(DY_real_ele_scale_up, -1);
  h_DATA_fele_pt_scale_up->Add(QCD_real_ele_scale_up, -1);
  h_DATA_fele_pt_scale_up->Add(SingleTop_real_ele_scale_up, -1);
  h_DATA_fele_pt_scale_up->Add(TTbar_real_ele_scale_up, -1);
  h_DATA_fele_pt_scale_up->Add(WJets_real_ele_scale_up, -1);



  //Add pt hists for the jets of the electrons w/o a gen particle (MC) hists for Nele >= 3
  h_sum_MC_fele_pt_scale_down->Add(h_QCD_fele_pt_scale_down);
  h_sum_MC_fele_pt_scale_down->Add(h_Diboson_fele_pt_scale_down);
  h_sum_MC_fele_pt_scale_down->Add(h_SingleTop_fele_pt_scale_down);
  h_sum_MC_fele_pt_scale_down->Add(h_TTbar_fele_pt_scale_down);
  h_sum_MC_fele_pt_scale_down->Add(h_WJets_fele_pt_scale_down);

  //Add all MC hists for Nele >= 2
  h_sum_MC_jet_pt_scale_down->Add(h_QCD_jet_pt_scale_down);
  h_sum_MC_jet_pt_scale_down->Add(h_Diboson_jet_pt_scale_down);
  h_sum_MC_jet_pt_scale_down->Add(h_SingleTop_jet_pt_scale_down);
  h_sum_MC_jet_pt_scale_down->Add(h_TTbar_jet_pt_scale_down);
  h_sum_MC_jet_pt_scale_down->Add(h_WJets_jet_pt_scale_down);

  //Substract the pt of jets which belong to real electrons, Nele >= 3
  h_DATA_fele_pt_scale_down->Add(Diboson_real_ele_scale_down, -1);
  h_DATA_fele_pt_scale_down->Add(DY_real_ele_scale_down, -1);
  h_DATA_fele_pt_scale_down->Add(QCD_real_ele_scale_down, -1);
  h_DATA_fele_pt_scale_down->Add(SingleTop_real_ele_scale_down, -1);
  h_DATA_fele_pt_scale_down->Add(TTbar_real_ele_scale_down, -1);
  h_DATA_fele_pt_scale_down->Add(WJets_real_ele_scale_down, -1);

  //divide fake ele pt and jet pt
  TGraphAsymmErrors* ratio_MC_scale_up = new TGraphAsymmErrors(h_sum_MC_fele_pt_scale_up, h_sum_MC_jet_pt_scale_up);
  TGraphAsymmErrors* ratio_MC_scale_down = new TGraphAsymmErrors(h_sum_MC_fele_pt_scale_down, h_sum_MC_jet_pt_scale_down);
  h_sum_MC_fele_pt_scale_up->Divide(h_sum_MC_jet_pt_scale_up);
  h_sum_MC_fele_pt_scale_down->Divide(h_sum_MC_jet_pt_scale_down);
  TH1D* dummy_MC_scale_up = (TH1D*) h_sum_MC_fele_pt_scale_up->Clone("dummy_MC_scale_up");
  TH1D* dummy_MC_scale_down = (TH1D*) h_sum_MC_fele_pt_scale_down->Clone("dummy_MC_scale_down");
  ratio_MC_scale_up->GetXaxis()->SetRangeUser(20,800);
  ratio_MC_scale_up->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_MC_scale_up->GetYaxis()->SetRangeUser(0,0.0014);
  ratio_MC_scale_up->GetYaxis()->SetTitle("#epsilon");
  ratio_MC_scale_up->SetMarkerColor(kBlue+1);
  ratio_MC_scale_up->SetLineColor(kBlue+1);
  dummy_MC_scale_up->SetMarkerColor(kBlue+1);
  dummy_MC_scale_up->SetLineColor(kBlue+1);
  ratio_MC_scale_up->SetMarkerStyle(8);
  dummy_MC_scale_up->SetMarkerStyle(8);
  ratio_MC_scale_up->SetMarkerSize(1.3);
  dummy_MC_scale_up->SetMarkerSize(1.3);
  ratio_MC_scale_up->SetLineWidth(2);
  dummy_MC_scale_up->SetLineWidth(2);
  ratio_MC_scale_up->Draw("ap");

  ratio_MC_scale_down->GetXaxis()->SetRangeUser(20,800);
  ratio_MC_scale_down->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_MC_scale_down->GetYaxis()->SetRangeUser(0,0.0014);
  ratio_MC_scale_down->GetYaxis()->SetTitle("#epsilon");
  ratio_MC_scale_down->SetMarkerColor(kGreen+1);
  ratio_MC_scale_down->SetLineColor(kGreen+1);
  dummy_MC_scale_down->SetMarkerColor(kGreen+1);
  dummy_MC_scale_down->SetLineColor(kGreen+1);
  ratio_MC_scale_down->SetMarkerStyle(8);
  dummy_MC_scale_down->SetMarkerStyle(8);
  ratio_MC_scale_down->SetMarkerSize(1.3);
  dummy_MC_scale_down->SetMarkerSize(1.3);
  ratio_MC_scale_down->SetLineWidth(2);
  dummy_MC_scale_down->SetLineWidth(2);
  ratio_MC_scale_down->Draw("e1p same");

  leg_scale_MC = new TLegend(0.639098,0.810105,0.803258,0.874564);
  leg_scale_MC->AddEntry(dummy_MC_scale_up,"Diboson scaled with +50% (MC)","lep");
  leg_scale_MC->AddEntry(dummy_MC_scale_down,"Diboson scaled with -50% (MC)","lep");
  leg_scale_MC->Draw();
  
  c6->cd();

  TGraphAsymmErrors* ratio_DATA_scale_up = new TGraphAsymmErrors(h_DATA_fele_pt_scale_up, h_DATA_jet_pt_scale_up);
  TGraphAsymmErrors* ratio_DATA_scale_down = new TGraphAsymmErrors(h_DATA_fele_pt_scale_down, h_DATA_jet_pt_scale_down);
  h_DATA_fele_pt_scale_up->Divide(h_DATA_jet_pt_scale_up);
  h_DATA_fele_pt_scale_down->Divide(h_DATA_jet_pt_scale_down);
  TH1D* dummy_DATA_scale_up = (TH1D*) h_DATA_fele_pt_scale_up->Clone("dummy_DATA_scale_up");
  TH1D* dummy_DATA_scale_down = (TH1D*) h_DATA_fele_pt_scale_down->Clone("dummy_DATA_scale_down");
  ratio_DATA_scale_up->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_scale_up->GetYaxis()->SetRangeUser(0,0.0014);
  ratio_DATA_scale_up->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_scale_up->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_scale_up->SetMarkerColor(kBlue+1);
  ratio_DATA_scale_up->SetLineColor(kBlue+1);
  dummy_DATA_scale_up->SetMarkerColor(kBlue+1);
  dummy_DATA_scale_up->SetLineColor(kBlue+1);
  ratio_DATA_scale_up->SetMarkerStyle(8);
  dummy_DATA_scale_up->SetMarkerStyle(8);
  ratio_DATA_scale_up->SetMarkerSize(1.3);
  dummy_DATA_scale_up->SetMarkerSize(1.3);
  ratio_DATA_scale_up->SetLineWidth(2);
  dummy_DATA_scale_up->SetLineWidth(2);
  ratio_DATA_scale_up->Draw("ap");

  ratio_DATA_scale_down->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_scale_down->GetYaxis()->SetRangeUser(0,0.0014);
  ratio_DATA_scale_down->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_scale_down->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_scale_down->SetMarkerColor(kGreen+1);
  ratio_DATA_scale_down->SetLineColor(kGreen+1);
  dummy_DATA_scale_down->SetMarkerColor(kGreen+1);
  dummy_DATA_scale_down->SetLineColor(kGreen+1);
  ratio_DATA_scale_down->SetMarkerStyle(8);
  dummy_DATA_scale_down->SetMarkerStyle(8);
  ratio_DATA_scale_down->SetMarkerSize(1.3);
  dummy_DATA_scale_down->SetMarkerSize(1.3);
  ratio_DATA_scale_down->SetLineWidth(2);
  dummy_DATA_scale_down->SetLineWidth(2);
  ratio_DATA_scale_down->Draw("e1p same");

  leg_scale_DATA = new TLegend(0.639098,0.810105,0.803258,0.874564);
  leg_scale_DATA->AddEntry(dummy_DATA_scale_up,"Diboson scaled with +50% (DATA)","lep");
  leg_scale_DATA->AddEntry(dummy_DATA_scale_down,"Diboson scaled with -50% (DATA)","lep");
  leg_scale_DATA->Draw();

  c7->cd();

  TGraphAsymmErrors* ratio_DATA_MC_up_down = new TGraphAsymmErrors(dummy_DATA_scale_up, dummy_MC_scale_up, "pois");
  TGraphAsymmErrors* ratio_DATA_MC_down_up = new TGraphAsymmErrors(dummy_DATA_scale_down, dummy_MC_scale_down, "pois");
  TH1D* dummy_DATA_MC_up_down = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_MC_up_down");
  TH1D* dummy_DATA_MC_down_up = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_MC_down_up");
  cout << endl;
  cout << "Used option 'pois' for the DATA/MC plot, can you use it?" << endl;
  cout << endl;

  ratio_DATA_MC_down_up->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC_down_up->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC_down_up->GetYaxis()->SetRangeUser(0,3.5);
  ratio_DATA_MC_down_up->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_MC_down_up->SetMarkerColor(kGreen+1);
  ratio_DATA_MC_down_up->SetLineColor(kGreen+1);
  dummy_DATA_MC_down_up->SetMarkerColor(kGreen+1);
  dummy_DATA_MC_down_up->SetLineColor(kGreen+1);
  ratio_DATA_MC_down_up->SetMarkerStyle(8);
  dummy_DATA_MC_down_up->SetMarkerStyle(8);
  ratio_DATA_MC_down_up->SetMarkerSize(1.3);
  dummy_DATA_MC_down_up->SetMarkerSize(1.3);
  ratio_DATA_MC_down_up->SetLineWidth(2);
  dummy_DATA_MC_down_up->SetLineWidth(2);
  ratio_DATA_MC_down_up->Draw("ap");

  ratio_DATA_MC_up_down->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC_up_down->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC_up_down->GetYaxis()->SetRangeUser(0,3.5);
  ratio_DATA_MC_up_down->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_MC_up_down->SetMarkerColor(kBlue+1);
  ratio_DATA_MC_up_down->SetLineColor(kBlue+1);
  dummy_DATA_MC_up_down->SetMarkerColor(kBlue+1);
  dummy_DATA_MC_up_down->SetLineColor(kBlue+1);
  ratio_DATA_MC_up_down->SetMarkerStyle(8);
  dummy_DATA_MC_up_down->SetMarkerStyle(8);
  ratio_DATA_MC_up_down->SetMarkerSize(1.3);
  dummy_DATA_MC_up_down->SetMarkerSize(1.3);
  ratio_DATA_MC_up_down->SetLineWidth(2);
  dummy_DATA_MC_up_down->SetLineWidth(2);
  ratio_DATA_MC_up_down->Draw("e1psame");

  leg_up_down = new TLegend(0.714286,0.824042,0.884712,0.881533);
  leg_up_down->AddEntry(dummy_DATA_MC_up_down,"DATA/MC (Data up, MC up)","lep");
  leg_up_down->AddEntry(dummy_DATA_MC_down_up,"DATA/MC (Data down, MC down)","lep");
  leg_up_down->Draw();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

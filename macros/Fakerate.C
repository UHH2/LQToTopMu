#include "TFile.h"
#include "TTree.h"
#include <TH1F.h>
#include <cmath>
using namespace std;

void Fakerate()
{
  /*
    one canvas for the nice plots
    one canvas for real sf with sys+stat error
    one canvas for wqs up+down + real sf
    one canvas for btag up+down + real sf
*/

  //Fakerate of pt
  //setup canvas
  TCanvas *sf = new TCanvas("sf","Scalefactor",200,10,800,600);
  TCanvas *wqs = new TCanvas("wqs","wqs variation",200,10,800,600);
  TCanvas *btag = new TCanvas("btag","btag variation",200,10,800,600);

  //load DATA
  TFile *DY = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.MC.DYJets.root", "READ");
  TFile *Diboson = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.MC.Diboson.root", "READ");
  TFile *QCD = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.MC.QCD.root", "READ");
  TFile *SingleTop = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.MC.SingleTop.root", "READ");
  TFile *TTbar = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.MC.TTbar.root", "READ");
  TFile *WJets = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.MC.WJets.root", "READ");
  TFile *DATA = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Fake_dR04_newPre_MET60_NJet2_NoHT_with_SF_2bJet_btagSF_HT350_scaled_Diboson/uhh2.AnalysisModuleRunner.DATA.DATA.root", "READ");

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
  h_DATA_fele_pt_rebin3->Add(Diboson_real_ele_rebin3, -1);
  h_DATA_fele_pt_rebin3->Add(DY_real_ele_rebin3, -1);
  h_DATA_fele_pt_rebin3->Add(QCD_real_ele_rebin3, -1);
  h_DATA_fele_pt_rebin3->Add(SingleTop_real_ele_rebin3, -1);
  h_DATA_fele_pt_rebin3->Add(TTbar_real_ele_rebin3, -1);
  h_DATA_fele_pt_rebin3->Add(WJets_real_ele_rebin3, -1);

  //divide fake ele pt and jet pt
  TGraphAsymmErrors* ratio_MC = new TGraphAsymmErrors(h_sum_MC_fele_pt_rebin3, h_sum_MC_jet_pt_rebin3);
  h_sum_MC_fele_pt_rebin3->Divide(h_sum_MC_jet_pt_rebin3);
  TH1D* dummy_MC = (TH1D*) h_sum_MC_fele_pt_rebin3->Clone("dummy_MC");
  ratio_MC->GetXaxis()->SetRangeUser(20,800);
  ratio_MC->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_MC->GetYaxis()->SetRangeUser(0.000001,0.0014);
  ratio_MC->GetYaxis()->SetTitle("#epsilon");
  ratio_MC->SetMarkerColor(kGreen+1);
  ratio_MC->SetLineColor(kGreen+1);
  dummy_MC->SetMarkerColor(kGreen+1);
  dummy_MC->SetLineColor(kGreen+1);
  ratio_MC->SetMarkerStyle(25);
  dummy_MC->SetMarkerStyle(25);
  ratio_MC->SetMarkerSize(0.6);
  dummy_MC->SetMarkerSize(0.6);
  ratio_MC->SetLineWidth(2);
  dummy_MC->SetLineWidth(2);
  ratio_MC->SetTitle("p_{T} fakerate [MC & DATA]");

  //h_DATA_fele_pt_rebin3->Divide(h_DATA_jet_pt_rebin3);
  TGraphAsymmErrors* ratio_DATA = new TGraphAsymmErrors(h_DATA_fele_pt_rebin3, h_DATA_jet_pt_rebin3);
  h_DATA_fele_pt_rebin3->Divide(h_DATA_jet_pt_rebin3);
  TH1D* dummy_DATA = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA");
  ratio_DATA->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA->GetYaxis()->SetRangeUser(0.000001,0.0014);
  ratio_DATA->SetMarkerColor(kBlack);
  ratio_DATA->SetLineColor(kBlack);
  dummy_DATA->SetMarkerColor(kBlack);
  dummy_DATA->SetLineColor(kBlack);
  ratio_DATA->SetMarkerStyle(8);
  dummy_DATA->SetMarkerStyle(8);
  ratio_DATA->SetMarkerSize(0.6);
  dummy_DATA->SetMarkerSize(0.6);
  ratio_DATA->SetLineWidth(2);
  dummy_DATA->SetLineWidth(2);
  
  //h_DATA_fele_pt_rebin3->Divide(h_sum_MC_fele_pt_rebin3);
  

  TGraphAsymmErrors* ratio_DATA_MC = new TGraphAsymmErrors(dummy_DATA, dummy_MC, "pois");
  TH1D* dummy_DATA_MC = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_MC");
  ratio_DATA_MC->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC->GetYaxis()->SetRangeUser(0,3.5);
  ratio_DATA_MC->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_MC->SetMarkerColor(kBlue);
  ratio_DATA_MC->SetLineColor(kBlue);
  dummy_DATA_MC->SetMarkerColor(kBlue);
  dummy_DATA_MC->SetLineColor(kBlue);
  ratio_DATA_MC->SetMarkerStyle(8);
  dummy_DATA_MC->SetMarkerStyle(8);
  ratio_DATA_MC->SetMarkerSize(0.6);
  dummy_DATA_MC->SetMarkerSize(0.6);
  ratio_DATA_MC->SetLineWidth(2);
  dummy_DATA_MC->SetLineWidth(2);
  ratio_DATA_MC->SetTitle("p_{T} fakerate [DATA/MC]");
  //ratio_DATA_MC->Draw("ap");

  double x_SF,y_SF;
  ratio_DATA_MC->GetPoint(0,x_SF,y_SF);
  ratio_DATA_MC->GetPoint(1,x_SF,y_SF);
  ratio_DATA_MC->GetPoint(2,x_SF,y_SF);

  
  
  ratio_DATA_MC_nice = (TGraphErrors*)ratio_DATA_MC->Clone("ratio_DATA_MC_nice");

  ratio_DATA_nice = (TGraphErrors*)ratio_DATA->Clone("ratio_DATA_nice");
  ratio_MC_nice = (TGraphErrors*)ratio_MC->Clone("ratio_MC_nice");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* ***********  Nice Ratio plots from fits ********** */

  // TFile *down_fakerate = new TFile("fakerate_down.root", "READ");
  // TFile *down_sf = new TFile("sf_down.root", "READ");

  // TH1D* up_error = (TH1D*)up->Get("cfitsigold");
  // TH1D* down_error = (TH1D*)down->Get("cfitsigold");

  // TGraphAsymmErrors* ratio_DATA = new TGraphAsymmErrors(h_DATA_fele_pt_rebin3, h_DATA_jet_pt_rebin3);

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

  /*TFile* out = new TFile("/nfs/dust/cms/user/skottkej/LQToTopMu/Run2_80X_v2/Optimization/SF_FakeEle.root","RECREATE");
  ratio_DATA_MC_nice->Write();
  out->Close();
  delete out;*/

  TLine *line = new TLine(20,1,800,1);
  line->SetLineColor(14);
  line->SetLineStyle(2);
  line->Draw();

  blub = (TGraphErrors*)ratio_DATA_MC->Clone("blub");

  sf->cd();
  
  double sf_wqs_up_1, sf_wqs_down_1, sf_btag_up_1, sf_btag_down_1, sf_1;
  double sf_wqs_up_2, sf_wqs_down_2, sf_btag_up_2, sf_btag_down_2, sf_2;
  double sf_wqs_up_3, sf_wqs_down_3, sf_btag_up_3, sf_btag_down_3, sf_3;
  double d_wqs_up_1, d_wqs_down_1, d_btag_up_1, d_btag_down_1;
  double d_wqs_up_2, d_wqs_down_2, d_btag_up_2, d_btag_down_2;
  double d_wqs_up_3, d_wqs_down_3, d_btag_up_3, d_btag_down_3;
  double sigma_syst_up_1, sigma_syst_down_1;
  double sigma_syst_up_2, sigma_syst_down_2;
  double sigma_syst_up_3, sigma_syst_down_3;
  double sigma_stat_up_1, sigma_stat_down_1;
  double sigma_stat_up_2, sigma_stat_down_2;
  double sigma_stat_up_3, sigma_stat_down_3;
  double sigma_total_up_1, sigma_total_down_1;
  double sigma_total_up_2, sigma_total_down_2;
  double sigma_total_up_3, sigma_total_down_3;

  sf_1 = 1.68716;
  sf_2 = 1.21191;
  sf_3 = 1.04207;

  sf_wqs_up_1 = 1.48704;
  sf_wqs_up_2 = 0.985038;
  sf_wqs_up_3 = 0.908532;

  sf_wqs_down_1 = 1.89142;
  sf_wqs_down_2 = 1.45017;
  sf_wqs_down_3 = 1.17408;

  sf_btag_up_1 = 1.48648;
  sf_btag_up_2 = 0.971895;
  sf_btag_up_3 = 0.8622;

  sf_btag_down_1 = 1.89058;
  sf_btag_down_2 = 1.46863;
  sf_btag_down_3 = 1.21813;

  sigma_stat_up_1 = ratio_DATA_MC_nice->GetErrorYhigh(0);
  sigma_stat_up_2 = ratio_DATA_MC_nice->GetErrorYhigh(1);
  sigma_stat_up_3 = ratio_DATA_MC_nice->GetErrorYhigh(2);

  sigma_stat_down_1 = ratio_DATA_MC_nice->GetErrorYlow(0);
  sigma_stat_down_2 = ratio_DATA_MC_nice->GetErrorYlow(1);
  sigma_stat_down_3 = ratio_DATA_MC_nice->GetErrorYlow(2);


  d_wqs_up_1 = sf_wqs_up_1 - sf_1;
  d_wqs_down_1 = sf_wqs_down_1 - sf_1;
  d_btag_up_1 = sf_btag_up_1 - sf_1;
  d_btag_down_1 = sf_btag_down_1 - sf_1;

  d_wqs_up_2 = sf_wqs_up_2 - sf_2;
  d_wqs_down_2 = sf_wqs_down_2 - sf_2;
  d_btag_up_2 = sf_btag_up_2 - sf_2;
  d_btag_down_2 = sf_btag_down_2 - sf_2;

  d_wqs_up_3 = sf_wqs_up_3 - sf_3;
  d_wqs_down_3 = sf_wqs_down_3 - sf_3;
  d_btag_up_3 = sf_btag_up_3 - sf_3;
  d_btag_down_3 = sf_btag_down_3 - sf_3;


  sigma_syst_up_1 = sqrt(pow(d_wqs_down_1, 2) + pow(d_btag_down_1, 2));
  sigma_syst_down_1 = sqrt(pow(d_wqs_up_1, 2) + pow(d_btag_up_1, 2));

  sigma_syst_up_2 = sqrt(pow(d_wqs_down_2, 2) + pow(d_btag_down_2, 2));
  sigma_syst_down_2 = sqrt(pow(d_wqs_up_2, 2) + pow(d_btag_up_2, 2));

  sigma_syst_up_3 = sqrt(pow(d_wqs_down_3, 2) + pow(d_btag_down_3, 2));
  sigma_syst_down_3 = sqrt(pow(d_wqs_up_3, 2) + pow(d_btag_up_3, 2));

  cout << endl;
  cout << endl;
  cout << "syst error up bin 1:" << sigma_syst_up_1 << endl;
  cout << "syst error up bin 2:" << sigma_syst_up_2 << endl;
  cout << "syst error up bin 3:" << sigma_syst_up_3 << endl;
  cout << endl;
  cout << endl;

  cout << "syst error down bin 1:" << sigma_syst_down_1 << endl;
  cout << "syst error down bin 2:" << sigma_syst_down_2 << endl;
  cout << "syst error down bin 3:" << sigma_syst_down_3 << endl;
  cout << endl;
  cout << endl;



  sigma_total_up_1 = sqrt(pow(sigma_stat_up_1, 2) + pow(sigma_syst_up_1, 2));
  sigma_total_up_2 = sqrt(pow(sigma_stat_up_2, 2) + pow(sigma_syst_up_2, 2));
  sigma_total_up_3 = sqrt(pow(sigma_stat_up_3, 2) + pow(sigma_syst_up_3, 2));

  sigma_total_down_1 = sqrt(pow(sigma_stat_down_1, 2) + pow(sigma_syst_down_1, 2));
  sigma_total_down_2 = sqrt(pow(sigma_stat_down_2, 2) + pow(sigma_syst_down_2, 2));
  sigma_total_down_3 = sqrt(pow(sigma_stat_down_3, 2) + pow(sigma_syst_down_3, 2));

  const Int_t n = 3;
  Double_t x[n]   = {60, 150, 500};
  Double_t y[n]   = {sf_1, sf_2, sf_3};
  Double_t exl[n] = {0, 0, 0};
  Double_t eyl[n] = {sigma_total_down_1, sigma_total_down_2, sigma_total_down_3};
  Double_t exh[n] = {0, 0, 0};
  Double_t eyh[n] = {sigma_total_up_1, sigma_total_up_2, sigma_total_up_3};
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


  cout << "stat Upper Error Bin 1: " << blub->GetErrorYhigh(0) << endl;
  cout << "stat Upper Error Bin 2: " << blub->GetErrorYhigh(1) << endl;
  cout << "stat Upper Error Bin 3: " << blub->GetErrorYhigh(2) << endl;
  cout << "stat Lower Error Bin 1: " << blub->GetErrorYlow(0) << endl;
  cout << "stat Lower Error Bin 2: " << blub->GetErrorYlow(1) << endl;
  cout << "stat Lower Error Bin 3: " << blub->GetErrorYlow(2) << endl;


  cout << "sys Upper Error Bin 1: " << gr->GetErrorYhigh(0) - blub->GetErrorYhigh(0) << endl;
  cout << "sys Upper Error Bin 2: " << gr->GetErrorYhigh(1) - blub->GetErrorYhigh(1) << endl;
  cout << "sys Upper Error Bin 3: " << gr->GetErrorYhigh(2) - blub->GetErrorYhigh(2) << endl;
  cout << "sys Lower Error Bin 1: " << gr->GetErrorYhigh(0) - blub->GetErrorYlow(0) << endl;
  cout << "sys Lower Error Bin 2: " << gr->GetErrorYhigh(1) - blub->GetErrorYlow(1) << endl;
  cout << "sys Lower Error Bin 3: " << gr->GetErrorYhigh(2) - blub->GetErrorYlow(2) << endl;



  m_rp1_top->cd();

  ratio_DATA_nice->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_nice->GetYaxis()->SetLabelSize(0.04);
  ratio_MC_nice->GetYaxis()->SetLabelSize(0.04);
  ratio_DATA_nice->GetYaxis()->SetTitleSize(0.07);
  ratio_DATA_nice->GetYaxis()->SetTitleOffset(0.9);
  ratio_DATA_nice->Draw("ap");
  ratio_MC_nice->Draw("e1p same");

  double x_eff_DATA,y_eff_DATA;
  ratio_DATA_nice->GetPoint(0,x_eff_DATA,y_eff_DATA);
  ratio_DATA_nice->GetPoint(1,x_eff_DATA,y_eff_DATA);
  ratio_DATA_nice->GetPoint(2,x_eff_DATA,y_eff_DATA);

  double x_eff_MC,y_eff_MC;
  ratio_MC_nice->GetPoint(0,x_eff_MC,y_eff_MC);
  ratio_MC_nice->GetPoint(1,x_eff_MC,y_eff_MC);
  ratio_MC_nice->GetPoint(2,x_eff_MC,y_eff_MC);


  double sf1, sf2, sf3;
  double error_x1, error_x2, error_x3;
  ratio_DATA_MC_nice->GetPoint(0, error_x1, sf1);
  ratio_DATA_MC_nice->GetPoint(1, error_x2, sf2);
  ratio_DATA_MC_nice->GetPoint(2, error_x3, sf3);

  leg_nice = new TLegend(0.746429,0.727759,0.90998,0.888294);
  leg_nice->AddEntry(dummy_MC_nice,"MC","lep");
  leg_nice->AddEntry(dummy_DATA_nice,"DATA","lep");
  leg_nice->SetLineWidth(0);
  leg_nice->Draw();




  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // scale Diboson up and down
  gStyle->SetOptTitle(0);

  gStyle->SetOptTitle(0);

  //scale the Diboson plots
  Diboson_real_ele_scale_up->Scale(1.1861);
  h_Diboson_fele_pt_scale_up->Scale(1.1861);
  h_Diboson_jet_pt_scale_up->Scale(1.1861);

  Diboson_real_ele_scale_down->Scale(0.8139);
  h_Diboson_fele_pt_scale_down->Scale(0.8139);
  h_Diboson_jet_pt_scale_down->Scale(0.8139);


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
  ratio_MC_scale_up->GetYaxis()->SetTickLength(0.015);
  ratio_MC_scale_up->GetYaxis()->SetTitle("DATA/MC");
  ratio_MC_scale_up->GetYaxis()->SetTitleSize(0.07);
  ratio_MC_scale_up->GetYaxis()->SetTitleOffset(0.65);
  ratio_MC_scale_up->SetMarkerColor(kBlue);
  ratio_MC_scale_up->SetLineColor(kBlue);
  dummy_MC_scale_up->SetMarkerColor(kBlue);
  dummy_MC_scale_up->SetLineColor(kBlue);
  ratio_MC_scale_up->SetMarkerStyle(8);
  dummy_MC_scale_up->SetMarkerStyle(8);
  ratio_MC_scale_up->SetMarkerSize(1.3);
  dummy_MC_scale_up->SetMarkerSize(1.3);
  ratio_MC_scale_up->SetLineWidth(2);
  dummy_MC_scale_up->SetLineWidth(2);

  ratio_MC_scale_down->GetXaxis()->SetRangeUser(20,800);
  ratio_MC_scale_down->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_MC_scale_down->GetYaxis()->SetTitle("DATA/MC");
  ratio_MC_scale_down->GetYaxis()->SetTickLength(0.015);
  ratio_MC_scale_down->GetYaxis()->SetTitleSize(0.07);
  ratio_MC_scale_down->GetYaxis()->SetTitleOffset(0.65);
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
  

  TGraphAsymmErrors* ratio_DATA_scale_up = new TGraphAsymmErrors(h_DATA_fele_pt_scale_up, h_DATA_jet_pt_scale_up);
  TGraphAsymmErrors* ratio_DATA_scale_down = new TGraphAsymmErrors(h_DATA_fele_pt_scale_down, h_DATA_jet_pt_scale_down);
  h_DATA_fele_pt_scale_up->Divide(h_DATA_jet_pt_scale_up);
  h_DATA_fele_pt_scale_down->Divide(h_DATA_jet_pt_scale_down);
  TH1D* dummy_DATA_scale_up = (TH1D*) h_DATA_fele_pt_scale_up->Clone("dummy_DATA_scale_up");
  TH1D* dummy_DATA_scale_down = (TH1D*) h_DATA_fele_pt_scale_down->Clone("dummy_DATA_scale_down");
  ratio_DATA_scale_up->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_scale_up->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_scale_up->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_scale_up->GetYaxis()->SetTickLength(0.015);
  ratio_DATA_scale_up->GetYaxis()->SetTitleSize(0.07);
  ratio_DATA_scale_up->GetYaxis()->SetTitleOffset(0.65);
  ratio_DATA_scale_up->SetMarkerColor(kBlue);
  ratio_DATA_scale_up->SetLineColor(kBlue);
  dummy_DATA_scale_up->SetMarkerColor(kBlue);
  dummy_DATA_scale_up->SetLineColor(kBlue);
  ratio_DATA_scale_up->SetMarkerStyle(8);
  dummy_DATA_scale_up->SetMarkerStyle(8);
  ratio_DATA_scale_up->SetMarkerSize(1.3);
  dummy_DATA_scale_up->SetMarkerSize(1.3);
  ratio_DATA_scale_up->SetLineWidth(2);
  dummy_DATA_scale_up->SetLineWidth(2);

  ratio_DATA_scale_down->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_scale_down->GetYaxis()->SetTitle("#epsilon");
  ratio_DATA_scale_down->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_scale_down->GetYaxis()->SetTickLength(0.015);
  ratio_DATA_scale_down->GetYaxis()->SetTitleSize(0.07);
  ratio_DATA_scale_down->GetYaxis()->SetTitleOffset(0.65);
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


  wqs->cd();

  TGraphAsymmErrors* ratio_DATA_MC_up_down = new TGraphAsymmErrors(dummy_DATA_scale_up, dummy_MC_scale_up, "pois");
  TGraphAsymmErrors* ratio_DATA_MC_down_up = new TGraphAsymmErrors(dummy_DATA_scale_down, dummy_MC_scale_down, "pois");
  TH1D* dummy_DATA_MC_up_down = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_MC_up_down");
  TH1D* dummy_DATA_MC_down_up = (TH1D*) h_DATA_fele_pt_rebin3->Clone("dummy_DATA_MC_down_up");

  ratio_DATA_MC_down_up->GetXaxis()->SetTickLength(0.055);
  ratio_DATA_MC_down_up->GetYaxis()->SetTickLength(0.015);
  ratio_DATA_MC_down_up->GetYaxis()->SetRangeUser(0,3.4);  
  ratio_DATA_MC_down_up->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC_down_up->GetYaxis()->SetTitleOffset(1);
  ratio_DATA_MC_down_up->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC_down_up->GetYaxis()->SetTitle("DATA/MC");
  ratio_DATA_MC_down_up->SetMarkerColor(kBlue);
  ratio_DATA_MC_down_up->SetLineColor(kBlue);
  dummy_DATA_MC_down_up->SetMarkerColor(kBlue);
  dummy_DATA_MC_down_up->SetLineColor(kBlue);
  gStyle->SetEndErrorSize(5);
  gPad->SetTicks(1,1);
  ratio_DATA_MC_down_up->SetMarkerStyle(8);
  dummy_DATA_MC_down_up->SetMarkerStyle(8);
  ratio_DATA_MC_down_up->SetLineStyle(kDashed);
  dummy_DATA_MC_down_up->SetLineStyle(kDashed);
  ratio_DATA_MC_down_up->SetMarkerSize(1.2);
  dummy_DATA_MC_down_up->SetMarkerSize(1.2);
  ratio_DATA_MC_down_up->SetLineWidth(2);
  dummy_DATA_MC_down_up->SetLineWidth(2);
  ratio_DATA_MC_down_up->Draw("ap");
  

  ratio_DATA_MC_up_down->SetMarkerColor(kGreen+1);
  ratio_DATA_MC_up_down->SetLineColor(kGreen+1);
  dummy_DATA_MC_up_down->SetMarkerColor(kGreen+1);
  dummy_DATA_MC_up_down->SetLineColor(kGreen+1);
  ratio_DATA_MC_up_down->SetMarkerStyle(8);
  dummy_DATA_MC_up_down->SetMarkerStyle(8);
  ratio_DATA_MC_up_down->SetLineStyle(kDashed);
  dummy_DATA_MC_up_down->SetLineStyle(kDashed);
  ratio_DATA_MC_up_down->SetMarkerSize(1.2);
  dummy_DATA_MC_up_down->SetMarkerSize(1.2);
  ratio_DATA_MC_up_down->SetLineWidth(2);
  dummy_DATA_MC_up_down->SetLineWidth(2);
  ratio_DATA_MC_up_down->Draw("e1psame");

  dummy =  (TGraphErrors*)ratio_DATA_MC_nice->Clone();

  dummy->Draw("e1psame");
  dummy->SetMarkerSize(1.2);
  line->Draw();

  leg_wqs = new TLegend(0.70912,0.682609,0.872639,0.843478);
  leg_wqs->AddEntry(dummy,"Skalierungsfaktor","lep");
  leg_wqs->AddEntry(ratio_DATA_MC_up_down,"+#sigma","lep");
  leg_wqs->AddEntry(ratio_DATA_MC_down_up,"-#sigma","lep");
  leg_wqs->SetLineWidth(0);
  leg_wqs->Draw();

  btag->cd();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //load DATA
  TFile *DY_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.MC.DYJets.root", "READ");
  TFile *Diboson_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.MC.Diboson.root", "READ");
  TFile *QCD_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.MC.QCD.root", "READ");
  TFile *SingleTop_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.MC.SingleTop.root", "READ");
  TFile *TTbar_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.MC.TTbar.root", "READ");
  TFile *WJets_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.MC.WJets.root", "READ");
  TFile *DATA_up = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_up/uhh2.AnalysisModuleRunner.DATA.DATA.root", "READ");

  //Rebin 3  
  //create hists
  // get pt hists of the jets which are the closest to the real electrons which are not part of the best pair, this needs to be substracted from DATA and MC
  TH1D* DY_real_ele_rebin3_up = (TH1D*)DY_up->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* Diboson_real_ele_rebin3_up = (TH1D*)Diboson_up->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* QCD_real_ele_rebin3_up = (TH1D*)QCD_up->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* SingleTop_real_ele_rebin3_up = (TH1D*)SingleTop_up->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* TTbar_real_ele_rebin3_up = (TH1D*)TTbar_up->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* WJets_real_ele_rebin3_up = (TH1D*)WJets_up->Get("FinalSelection/pt_jet_notclose_rebin_3");

  // get pt hists of all jets for Nele >= 3
  TH1D* h_DY_fele_pt_rebin3_up = (TH1D*)DY_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_Diboson_fele_pt_rebin3_up = (TH1D*)Diboson_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_QCD_fele_pt_rebin3_up = (TH1D*)QCD_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_SingleTop_fele_pt_rebin3_up = (TH1D*)SingleTop_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_TTbar_fele_pt_rebin3_up = (TH1D*)TTbar_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_WJets_fele_pt_rebin3_up = (TH1D*)WJets_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_DATA_fele_pt_rebin3_up = (TH1D*)DATA_up->Get("FinalSelection/pt_jet_data_rebin_3");

  // get pt hists of all jets for Nele >= 2
  TH1D* h_DY_jet_pt_rebin3_up = (TH1D*)DY_up->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_Diboson_jet_pt_rebin3_up = (TH1D*)Diboson_up->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_QCD_jet_pt_rebin3_up = (TH1D*)QCD_up->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_SingleTop_jet_pt_rebin3_up = (TH1D*)SingleTop_up->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_TTbar_jet_pt_rebin3_up = (TH1D*)TTbar_up->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_WJets_jet_pt_rebin3_up = (TH1D*)WJets_up->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_DATA_jet_pt_rebin3_up = (TH1D*)DATA_up->Get("0bJetLoose/pt_jets_rebin_3");

  //for adding hists
  TH1D* h_sum_MC_fele_pt_rebin3_up = (TH1D*)DY_up->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_sum_MC_jet_pt_rebin3_up = (TH1D*)DY_up->Get("0bJetLoose/pt_jets_rebin_3");


  gStyle->SetOptStat(0);

  //Add pt hists for the jets of the electrons w/o a gen particle (MC) hists for Nele >= 3
  h_sum_MC_fele_pt_rebin3_up->Add(h_QCD_fele_pt_rebin3_up);
  h_sum_MC_fele_pt_rebin3_up->Add(h_Diboson_fele_pt_rebin3_up);
  h_sum_MC_fele_pt_rebin3_up->Add(h_SingleTop_fele_pt_rebin3_up);
  h_sum_MC_fele_pt_rebin3_up->Add(h_TTbar_fele_pt_rebin3_up);
  h_sum_MC_fele_pt_rebin3_up->Add(h_WJets_fele_pt_rebin3_up);

  //Add all MC hists for Nele >= 2
  h_sum_MC_jet_pt_rebin3_up->Add(h_QCD_jet_pt_rebin3_up);
  h_sum_MC_jet_pt_rebin3_up->Add(h_Diboson_jet_pt_rebin3_up);
  h_sum_MC_jet_pt_rebin3_up->Add(h_SingleTop_jet_pt_rebin3_up);
  h_sum_MC_jet_pt_rebin3_up->Add(h_TTbar_jet_pt_rebin3_up);
  h_sum_MC_jet_pt_rebin3_up->Add(h_WJets_jet_pt_rebin3_up);

  //Substract the pt of jets which belong to real electrons, Nele >= 3
  h_DATA_fele_pt_rebin3_up->Add(Diboson_real_ele_rebin3_up, -1);
  h_DATA_fele_pt_rebin3_up->Add(DY_real_ele_rebin3_up, -1);
  h_DATA_fele_pt_rebin3_up->Add(QCD_real_ele_rebin3_up, -1);
  h_DATA_fele_pt_rebin3_up->Add(SingleTop_real_ele_rebin3_up, -1);
  h_DATA_fele_pt_rebin3_up->Add(TTbar_real_ele_rebin3_up, -1);
  h_DATA_fele_pt_rebin3_up->Add(WJets_real_ele_rebin3_up, -1);

  //divide fake ele pt and jet pt
  TGraphAsymmErrors* ratio_MC_up = new TGraphAsymmErrors(h_sum_MC_fele_pt_rebin3_up, h_sum_MC_jet_pt_rebin3_up);
  h_sum_MC_fele_pt_rebin3_up->Divide(h_sum_MC_jet_pt_rebin3_up);
  TH1D* dummy_MC_up = (TH1D*) h_sum_MC_fele_pt_rebin3_up->Clone("dummy_MC_up");

  TGraphAsymmErrors* ratio_DATA_up = new TGraphAsymmErrors(h_DATA_fele_pt_rebin3_up, h_DATA_jet_pt_rebin3_up);
  h_DATA_fele_pt_rebin3_up->Divide(h_DATA_jet_pt_rebin3_up);
  TH1D* dummy_DATA_up = (TH1D*) h_DATA_fele_pt_rebin3_up->Clone("dummy_DATA_up");

  TGraphAsymmErrors* ratio_DATA_MC_up = new TGraphAsymmErrors(dummy_DATA_up, dummy_MC_up, "pois");
  TH1D* dummy_DATA_MC_up = (TH1D*) h_DATA_fele_pt_rebin3_up->Clone("dummy_DATA_MC_up");

  ratio_DATA_MC_up->GetXaxis()->SetTickLength(0.055);
  ratio_DATA_MC_up->GetYaxis()->SetTickLength(0.015);
  ratio_DATA_MC_up->GetYaxis()->SetTitleOffset(1);
  gStyle->SetEndErrorSize(5);
  ratio_DATA_MC_up->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC_up->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_DATA_MC_up->GetYaxis()->SetRangeUser(0,3.4);
  ratio_DATA_MC_up->GetYaxis()->SetTitle("DATA/MC");
  ratio_DATA_MC_up->SetMarkerColor(kGreen+1);
  ratio_DATA_MC_up->SetLineColor(kGreen+1);
  dummy_DATA_MC_up->SetMarkerColor(kGreen+1);
  dummy_DATA_MC_up->SetLineColor(kGreen+1);
  ratio_DATA_MC_up->SetMarkerStyle(8);
  dummy_DATA_MC_up->SetMarkerStyle(8);
  ratio_DATA_MC_up->SetLineStyle(kDashed);
  dummy_DATA_MC_up->SetLineStyle(kDashed);
  ratio_DATA_MC_up->SetMarkerSize(1.2);
  dummy_DATA_MC_up->SetMarkerSize(1.2);
  ratio_DATA_MC_up->SetLineWidth(2);
  dummy_DATA_MC_up->SetLineWidth(2);
  ratio_DATA_MC_up->Draw("ap");
  gPad->SetTicks(1,1);
  line->Draw();






  //load DATA
  TFile *DY_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.MC.DYJets.root", "READ");
  TFile *Diboson_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.MC.Diboson.root", "READ");
  TFile *QCD_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.MC.QCD.root", "READ");
  TFile *SingleTop_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.MC.SingleTop.root", "READ");
  TFile *TTbar_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.MC.TTbar.root", "READ");
  TFile *WJets_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.MC.WJets.root", "READ");
  TFile *DATA_down = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_sf_down/uhh2.AnalysisModuleRunner.DATA.DATA.root", "READ");

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Rebin 3  
  //create hists
  // get pt hists of the jets which are the closest to the real electrons which are not part of the best pair, this needs to be substracted from DATA and MC
  TH1D* DY_real_ele_rebin3_down = (TH1D*)DY_down->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* Diboson_real_ele_rebin3_down = (TH1D*)Diboson_down->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* QCD_real_ele_rebin3_down = (TH1D*)QCD_down->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* SingleTop_real_ele_rebin3_down = (TH1D*)SingleTop_down->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* TTbar_real_ele_rebin3_down = (TH1D*)TTbar_down->Get("FinalSelection/pt_jet_notclose_rebin_3");
  TH1D* WJets_real_ele_rebin3_down = (TH1D*)WJets_down->Get("FinalSelection/pt_jet_notclose_rebin_3");

  // get pt hists of all jets for Nele >= 3
  TH1D* h_DY_fele_pt_rebin3_down = (TH1D*)DY_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_Diboson_fele_pt_rebin3_down = (TH1D*)Diboson_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_QCD_fele_pt_rebin3_down = (TH1D*)QCD_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_SingleTop_fele_pt_rebin3_down = (TH1D*)SingleTop_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_TTbar_fele_pt_rebin3_down = (TH1D*)TTbar_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_WJets_fele_pt_rebin3_down = (TH1D*)WJets_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_DATA_fele_pt_rebin3_down = (TH1D*)DATA_down->Get("FinalSelection/pt_jet_data_rebin_3");

  // get pt hists of all jets for Nele >= 2
  TH1D* h_DY_jet_pt_rebin3_down = (TH1D*)DY_down->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_Diboson_jet_pt_rebin3_down = (TH1D*)Diboson_down->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_QCD_jet_pt_rebin3_down = (TH1D*)QCD_down->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_SingleTop_jet_pt_rebin3_down = (TH1D*)SingleTop_down->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_TTbar_jet_pt_rebin3_down = (TH1D*)TTbar_down->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_WJets_jet_pt_rebin3_down = (TH1D*)WJets_down->Get("0bJetLoose/pt_jets_rebin_3");
  TH1D* h_DATA_jet_pt_rebin3_down = (TH1D*)DATA_down->Get("0bJetLoose/pt_jets_rebin_3");

  //for adding hists
  TH1D* h_sum_MC_fele_pt_rebin3_down = (TH1D*)DY_down->Get("FinalSelection/pt_fjet_nogen_notclose_rebin_3");
  TH1D* h_sum_MC_jet_pt_rebin3_down = (TH1D*)DY_down->Get("0bJetLoose/pt_jets_rebin_3");


  gStyle->SetOptStat(0);

  //Add pt hists for the jets of the electrons w/o a gen particle (MC) hists for Nele >= 3
  h_sum_MC_fele_pt_rebin3_down->Add(h_QCD_fele_pt_rebin3_down);
  h_sum_MC_fele_pt_rebin3_down->Add(h_Diboson_fele_pt_rebin3_down);
  h_sum_MC_fele_pt_rebin3_down->Add(h_SingleTop_fele_pt_rebin3_down);
  h_sum_MC_fele_pt_rebin3_down->Add(h_TTbar_fele_pt_rebin3_down);
  h_sum_MC_fele_pt_rebin3_down->Add(h_WJets_fele_pt_rebin3_down);

  //Add all MC hists for Nele >= 2
  h_sum_MC_jet_pt_rebin3_down->Add(h_QCD_jet_pt_rebin3_down);
  h_sum_MC_jet_pt_rebin3_down->Add(h_Diboson_jet_pt_rebin3_down);
  h_sum_MC_jet_pt_rebin3_down->Add(h_SingleTop_jet_pt_rebin3_down);
  h_sum_MC_jet_pt_rebin3_down->Add(h_TTbar_jet_pt_rebin3_down);
  h_sum_MC_jet_pt_rebin3_down->Add(h_WJets_jet_pt_rebin3_down);

  //Substract the pt of jets which belong to real electrons, Nele >= 3
  h_DATA_fele_pt_rebin3_down->Add(Diboson_real_ele_rebin3_down, -1);
  h_DATA_fele_pt_rebin3_down->Add(DY_real_ele_rebin3_down, -1);
  h_DATA_fele_pt_rebin3_down->Add(QCD_real_ele_rebin3_down, -1);
  h_DATA_fele_pt_rebin3_down->Add(SingleTop_real_ele_rebin3_down, -1);
  h_DATA_fele_pt_rebin3_down->Add(TTbar_real_ele_rebin3_down, -1);
  h_DATA_fele_pt_rebin3_down->Add(WJets_real_ele_rebin3_down, -1);
  //divide fake ele pt and jet pt
  TGraphAsymmErrors* ratio_MC_down = new TGraphAsymmErrors(h_sum_MC_fele_pt_rebin3_down, h_sum_MC_jet_pt_rebin3_down);
  h_sum_MC_fele_pt_rebin3_down->Divide(h_sum_MC_jet_pt_rebin3_down);
  TH1D* dummy_MC_down = (TH1D*) h_sum_MC_fele_pt_rebin3_down->Clone("dummy_MC_down");

  //h_DATA_fele_pt_rebin3->Divide(h_DATA_jet_pt_rebin3);
  TGraphAsymmErrors* ratio_DATA_down = new TGraphAsymmErrors(h_DATA_fele_pt_rebin3_down, h_DATA_jet_pt_rebin3_down);
  h_DATA_fele_pt_rebin3_down->Divide(h_DATA_jet_pt_rebin3_down);
  TH1D* dummy_DATA_down = (TH1D*) h_DATA_fele_pt_rebin3_down->Clone("dummy_DATA_down");
  

  TGraphAsymmErrors* ratio_DATA_MC_down = new TGraphAsymmErrors(dummy_DATA_down, dummy_MC_down, "pois");
  TH1D* dummy_DATA_MC_down = (TH1D*) h_DATA_fele_pt_rebin3_down->Clone("dummy_DATA_MC_down");

  ratio_DATA_MC_down->GetXaxis()->SetTickLength(0.055);
  ratio_DATA_MC_down->GetYaxis()->SetTickLength(0.015);
  ratio_DATA_MC_down->GetYaxis()->SetTitleOffset(1);
  gStyle->SetEndErrorSize(5);
  ratio_DATA_MC_down->GetXaxis()->SetRangeUser(20,800);
  ratio_DATA_MC_down->GetYaxis()->SetRangeUser(0,3.4);
  ratio_DATA_MC_down->SetMarkerColor(kBlue);
  ratio_DATA_MC_down->SetLineColor(kBlue);
  dummy_DATA_MC_down->SetMarkerColor(kBlue);
  dummy_DATA_MC_down->SetLineColor(kBlue);
  ratio_DATA_MC_down->SetMarkerStyle(8);
  dummy_DATA_MC_down->SetMarkerStyle(8);
  ratio_DATA_MC_down->SetLineStyle(kDashed);
  dummy_DATA_MC_down->SetLineStyle(kDashed);
  ratio_DATA_MC_down->SetMarkerSize(1.2);
  dummy_DATA_MC_down->SetMarkerSize(1.2);
  ratio_DATA_MC_down->SetLineWidth(2);
  dummy_DATA_MC_down->SetLineWidth(2);
  ratio_DATA_MC_down->Draw("e1psame");
  dummy->Draw("e1psame");

  leg_btag = new TLegend(0.70912,0.682609,0.872639,0.843478);
  leg_btag->SetLineWidth(0);
  leg_btag->AddEntry(dummy,"Skalierungsfaktor","lep");
  leg_btag->AddEntry(dummy_DATA_MC_up,"+#sigma","lep");
  leg_btag->AddEntry(dummy_DATA_MC_down,"-#sigma","lep");
  leg_btag->Draw();




  cfitsigold->SaveAs("fakerate.eps");
  sf->SaveAs("Scalefactor.eps");
  wqs->SaveAs("wqs.eps");
  btag->SaveAs("btag.eps");
}

#include "TFile.h"
#include "TTree.h"
#include <TH1F.h>
#include <cmath>

using namespace std;

void Diboson_WQS(){
  //Load data
  TFile *DY = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.DYJets.root", "READ");
  TFile *Diboson = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.Diboson.root", "READ");
  TFile *Diboson_WW = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.Diboson_WW.root", "READ");
  TFile *Diboson_WZ = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.Diboson_WZ.root", "READ");
  TFile *Diboson_ZZ = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.Diboson_ZZ.root", "READ");
  TFile *QCD = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.QCD.root", "READ");
  TFile *SingleTop = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.SingleTop.root", "READ");
  TFile *TTbar = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.TTbar.root", "READ");
  TFile *WJets = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.WJets.root", "READ");

  TFile *DATA = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.DATA.DATA.root", "READ");

  //Get sum event weights
  TH1D* DY_N = (TH1D*)DY->Get("FinalSelection/sum_event_weights");
  TH1D* Diboson_N = (TH1D*)Diboson->Get("FinalSelection/sum_event_weights");
  TH1D* Diboson_WW_N = (TH1D*)Diboson_WW->Get("FinalSelection/sum_event_weights");
  TH1D* Diboson_WZ_N = (TH1D*)Diboson_WZ->Get("FinalSelection/sum_event_weights");
  TH1D* Diboson_ZZ_N = (TH1D*)Diboson_ZZ->Get("FinalSelection/sum_event_weights");
  TH1D* QCD_N = (TH1D*)QCD->Get("FinalSelection/sum_event_weights");
  TH1D* SingleTop_N = (TH1D*)SingleTop->Get("FinalSelection/sum_event_weights");
  TH1D* TTbar_N = (TH1D*)TTbar->Get("FinalSelection/sum_event_weights");
  TH1D* WJets_N = (TH1D*)WJets->Get("FinalSelection/sum_event_weights");

  TH1D* DATA_N = (TH1D*)DATA->Get("FinalSelection/sum_event_weights");

  //clone the histograms
  TH1D* DY_sf = (TH1D*) DY_N->Clone("DY_sf");
  TH1D* Diboson_sf = (TH1D*) Diboson_N->Clone("Diboson_sf");
  TH1D* QCD_sf = (TH1D*) QCD_N->Clone("QCD_sf");
  TH1D* SingleTop_sf = (TH1D*) SingleTop_N->Clone("SingelTop_sf");
  TH1D* TTbar_sf = (TH1D*) TTbar_N->Clone("TTbar_sf");
  TH1D* WJets_sf = (TH1D*) WJets_N->Clone("WJets_sf");
  TH1D* DATA_sf = (TH1D*) DATA_N->Clone("DATA_sf");
  TH1D* blub = (TH1D*) DY_N->Clone("blub");


  blub->Add(QCD_sf,1);
  blub->Add(SingleTop_sf,1);
  blub->Add(TTbar_sf,1);
  blub->Add(WJets_sf,1);
  blub->Add(Diboson_sf,1);

  cout << "reinheit: " << Diboson_sf->GetBinContent(1)/blub->GetBinContent(1) << endl;
  cout << "Anteil WW: " << Diboson_WW_N->GetBinContent(1)/Diboson_sf->GetBinContent(1) << endl;
  cout << "Anteil WZ: " << Diboson_WZ_N->GetBinContent(1)/Diboson_sf->GetBinContent(1) << endl;
  cout << "Anteil ZZ: " << Diboson_ZZ_N->GetBinContent(1)/Diboson_sf->GetBinContent(1) << endl;


  DATA_sf->Add(DY_sf,-1);
  DATA_sf->Add(QCD_sf,-1);
  DATA_sf->Add(SingleTop_sf,-1);
  DATA_sf->Add(TTbar_sf,-1);
  DATA_sf->Add(WJets_sf,-1);

  TH1D* sf = (TH1D*) DATA_sf->Clone("sf");

  sf->Divide(Diboson_sf);
  
  double sigma_D = DATA_sf->GetBinError(1);
  cout << "Error of Data: " << sigma_D << endl;

  double sigma_MC = Diboson_sf->GetBinError(1);
  cout << "Error of Diboson: " << sigma_MC << endl;

  double content_DATA = DATA_sf->GetBinContent(1);
  cout << "content of Data: " << content_DATA << endl;

  double content_Diboson = Diboson_sf->GetBinContent(1);
  cout << "content of Diboson: " << content_Diboson << endl;
  cout << endl;

  double content_sf = sf->GetBinContent(1);
  cout << "Aplly SF on Diboson: " << content_sf << endl;

  double sigma_sf = sqrt(pow((1/(content_Diboson)*sigma_D),2) + (pow((content_DATA*sigma_MC)/pow(content_Diboson,2),2)));

  cout << "Vary Diboson by: " << sigma_sf/content_sf * 100 << "%" << endl;

}

#include "TFile.h"
#include "TTree.h"
#include <TH1F.h>
#include <cmath>

using namespace std;

void SF_Diboson(){
  //Load data
  TFile *DY = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.DYJets.root", "READ");
  TFile *Diboson = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.Diboson.root", "READ");
  TFile *QCD = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.QCD.root", "READ");
  TFile *SingleTop = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.SingleTop.root", "READ");
  TFile *TTbar = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.TTbar.root", "READ");
  TFile *WJets = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.MC.WJets.root", "READ");

  TFile *DATA = new TFile("../../../../LQToTopMu/Run2_80X_v2/Optimization/27200fb_NoMuonSF_Diboson_newPre_0Jets_NoHT_deleteJets_with_SF_2bJet_btagSF/uhh2.AnalysisModuleRunner.DATA.DATA.root", "READ");

  //Get sum event weights
  TH1D* DY_bJet = (TH1D*)DY->Get("FinalSelection/N_bJets_loose");
  TH1D* Diboson_bJet = (TH1D*)Diboson->Get("FinalSelection/N_bJets_loose");
  TH1D* QCD_bJet = (TH1D*)QCD->Get("FinalSelection/N_bJets_loose");
  TH1D* SingleTop_bJet = (TH1D*)SingleTop->Get("FinalSelection/N_bJets_loose");
  TH1D* TTbar_bJet = (TH1D*)TTbar->Get("FinalSelection/N_bJets_loose");
  TH1D* WJets_bJet = (TH1D*)WJets->Get("FinalSelection/N_bJets_loose");
  TH1D* DATA_bJet = (TH1D*)DATA->Get("FinalSelection/N_bJets_loose");

  TH1D* DY_clone = (TH1D*) DY_bJet->Clone("DY_clone");
  TH1D* Diboson_clone = (TH1D*) Diboson_bJet->Clone("Diboson_clone");
  TH1D* QCD_clone = (TH1D*) QCD_bJet->Clone("QCD_clone");
  TH1D* SingleTop_clone = (TH1D*) SingleTop_bJet->Clone("SingleTop_clone");
  TH1D* TTbar_clone = (TH1D*) TTbar_bJet->Clone("TTbar_clone");
  TH1D* WJets_clone = (TH1D*) WJets_bJet->Clone("WJets_clone");
  TH1D* DATA_clone = (TH1D*) DATA_bJet->Clone("DATA_clone");

  double DY_Nbjet_loose_0, Diboson_Nbjet_loose_0, QCD_Nbjet_loose_0, SingleTop_Nbjet_loose_0, TTbar_Nbjet_loose_0, WJets_Nbjet_loose_0, DATA_Nbjet_loose_0;
  double DY_Nbjet_loose_1, Diboson_Nbjet_loose_1, QCD_Nbjet_loose_1, SingleTop_Nbjet_loose_1, TTbar_Nbjet_loose_1, WJets_Nbjet_loose_1, DATA_Nbjet_loose_1;
  double DY_Nbjet_loose_2, Diboson_Nbjet_loose_2, QCD_Nbjet_loose_2, SingleTop_Nbjet_loose_2, TTbar_Nbjet_loose_2, WJets_Nbjet_loose_2, DATA_Nbjet_loose_2;
  double SF_Nbjet_loose_0, SF_Nbjet_loose_1, SF_Nbjet_loose_2;

  DY_Nbjet_loose_0 = DY_clone->GetBinContent(1);
  Diboson_Nbjet_loose_0 = Diboson_clone->GetBinContent(1);
  QCD_Nbjet_loose_0 = QCD_clone->GetBinContent(1);
  SingleTop_Nbjet_loose_0 = SingleTop_clone->GetBinContent(1);
  TTbar_Nbjet_loose_0 = TTbar_clone->GetBinContent(1);
  WJets_Nbjet_loose_0 = WJets_clone->GetBinContent(1);
  DATA_Nbjet_loose_0 = DATA_clone->GetBinContent(1);

  DY_Nbjet_loose_1 = DY_clone->GetBinContent(2);
  Diboson_Nbjet_loose_1 = Diboson_clone->GetBinContent(2);
  QCD_Nbjet_loose_1 = QCD_clone->GetBinContent(2);
  SingleTop_Nbjet_loose_1 = SingleTop_clone->GetBinContent(2);
  TTbar_Nbjet_loose_1 = TTbar_clone->GetBinContent(2);
  WJets_Nbjet_loose_1 = WJets_clone->GetBinContent(2);
  DATA_Nbjet_loose_1 = DATA_clone->GetBinContent(2);

  DY_Nbjet_loose_2 = DY_clone->GetBinContent(3);
  Diboson_Nbjet_loose_2 = Diboson_clone->GetBinContent(3);
  QCD_Nbjet_loose_2 = QCD_clone->GetBinContent(3);
  SingleTop_Nbjet_loose_2 = SingleTop_clone->GetBinContent(3);
  TTbar_Nbjet_loose_2 = TTbar_clone->GetBinContent(3);
  WJets_Nbjet_loose_2 = WJets_clone->GetBinContent(3);
  DATA_Nbjet_loose_2 = DATA_clone->GetBinContent(3);

  SF_Nbjet_loose_0 = (DATA_Nbjet_loose_0 - WJets_Nbjet_loose_0 - TTbar_Nbjet_loose_0 - SingleTop_Nbjet_loose_0 - QCD_Nbjet_loose_0 - DY_Nbjet_loose_0)/Diboson_Nbjet_loose_0;

  cout << endl;
  cout << "Scale factor for Diboson for 0 bJets_loose: " << SF_Nbjet_loose_0 << endl;

  SF_Nbjet_loose_1 = (DATA_Nbjet_loose_1 - WJets_Nbjet_loose_1 - TTbar_Nbjet_loose_1 - SingleTop_Nbjet_loose_1 - QCD_Nbjet_loose_1 - DY_Nbjet_loose_1)/Diboson_Nbjet_loose_1;

  cout << "Scale factor for Diboson for 1 bJets_loose: " << SF_Nbjet_loose_1 << endl;

  SF_Nbjet_loose_2 = (DATA_Nbjet_loose_2 - WJets_Nbjet_loose_2 - TTbar_Nbjet_loose_2 - SingleTop_Nbjet_loose_2 - QCD_Nbjet_loose_2 - DY_Nbjet_loose_2)/Diboson_Nbjet_loose_2;

  cout << "Scale factor for Diboson for 2 bJets_loose: " << SF_Nbjet_loose_2 << endl;
  cout << endl;





  DATA_clone->Add(DY_clone,-1);
  DATA_clone->Add(QCD_clone,-1);
  DATA_clone->Add(SingleTop_clone,-1);
  DATA_clone->Add(TTbar_clone,-1);
  DATA_clone->Add(WJets_clone,-1);

  TH1D* sf_1 = (TH1D*) DATA_clone->Clone("sf_1");
  TH1D* sf_2 = (TH1D*) DATA_clone->Clone("sf_2");

  sf_1->Divide(Diboson_clone);
  
  double sigma_D_1 = DATA_clone->GetBinError(2);
  double sigma_MC_1 = Diboson_clone->GetBinError(2);
  double content_DATA_1 = DATA_clone->GetBinContent(2);
  double content_Diboson_1 = Diboson_clone->GetBinContent(2);

  double content_clone_1 = sf_1->GetBinContent(2);
  cout << "Scale factor for Diboson for 1 bJets_loose: " << content_clone_1 << endl;

  double sigma_clone_1 = sqrt(pow((1/(content_Diboson_1)*sigma_D_1),2) + (pow((content_DATA_1*sigma_MC_1)/pow(content_Diboson_1,2),2)));
  cout << "Vary Diboson by: " << sigma_clone_1/content_clone_1 * 100 << "%" << endl;
  cout << "Scale Factor up: " << (1+sigma_clone_1/content_clone_1) * content_clone_1 << endl;
  cout << "Scale Factor down: " << (1-sigma_clone_1/content_clone_1) * content_clone_1 << endl;
  cout << endl;



  sf_2->Divide(Diboson_clone);
  
  double sigma_D = DATA_clone->GetBinError(3);
  double sigma_MC = Diboson_clone->GetBinError(3);
  double content_DATA = DATA_clone->GetBinContent(3);
  double content_Diboson = Diboson_clone->GetBinContent(3);


  double content_clone = sf_2->GetBinContent(3);
  cout << "Scale factor for Diboson for 2 bJets_loose: " << content_clone << endl;
  cout << "Upper error bin 3: " << sf_2->GetBinErrorUp(3) << endl;
  cout << "Lower error bin 3: " << sf_2->GetBinErrorLow(3) << endl;
  cout << "Upper error bin 2: " << sf_2->GetBinErrorUp(2) << endl;
  cout << "Lower error bin 2: " << sf_2->GetBinErrorLow(2) << endl;
  cout << "Upper error bin 1: " << sf_2->GetBinErrorUp(1) << endl;
  cout << "Lower error bin 1: " << sf_2->GetBinErrorLow(1) << endl;

  double sigma_clone = sqrt(pow((1/(content_Diboson)*sigma_D),2) + (pow((content_DATA*sigma_MC)/pow(content_Diboson,2),2)));
  cout << "Vary Diboson by: " << sigma_clone/content_clone * 100 << "%" << endl;
  cout << "Scale Factor up: " << (1+sigma_clone/content_clone) * content_clone << endl;
  cout << "Scale Factor down: " << (1-sigma_clone/content_clone) * content_clone << endl;



}

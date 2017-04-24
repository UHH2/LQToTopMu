#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstruction.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuAnalysisModule: public AnalysisModule {
  public:

    explicit LQToTopMuAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    unique_ptr<CommonModules> common;

    unique_ptr<AnalysisModule> Ele_SF_Reco, Ele_SF_TightID, Ele_SF_LooseID, SF_btag;
    
    unique_ptr<JetCleaner> jetcleaner;
    
    // declare the Selections to use.
    unique_ptr<Selection>  ele_sel, nele_sel, njet_sel, m_mumu_veto, htlept_sel, met_sel, m_eleele_veto, nbtag_loose_sel, ht_sel;
    double dR_max, dR_gen_reco, Z_mass, sf_0bjet, sf_1bjet, sf_2bjet;
    int Nele;
    /*unique_ptr<AnalysisModule> Electron_printer, Jet_printer, GenParticles_printer;*/
    
    // store the Hists collection as member variables.
    unique_ptr<Hists> h_2ele, h_jets_2ele, h_ele_2ele, h_mu_2ele, h_event_2ele, h_topjets_2ele;
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts;
    unique_ptr<Hists> h_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto;
    unique_ptr<Hists> h_htlept0, h_jets_htlept0, h_ele_htlept0, h_mu_htlept0, h_event_htlept0, h_topjets_htlept0;
    unique_ptr<Hists> h_MET60, h_jets_MET60, h_ele_MET60, h_mu_MET60, h_event_MET60, h_topjets_MET60;
    unique_ptr<Hists> h_NJet_2, h_jets_NJet_2, h_ele_NJet_2, h_mu_NJet_2, h_event_NJet_2, h_topjets_NJet_2;
    unique_ptr<Hists> h_3ele, h_jets_3ele, h_ele_3ele, h_mu_3ele, h_event_3ele, h_topjets_3ele;
    unique_ptr<Hists> h_dR, h_jets_dR, h_ele_dR, h_mu_dR, h_event_dR, h_topjets_dR;
    unique_ptr<Hists> h_ht350, h_jets_ht350, h_ele_ht350, h_mu_ht350, h_event_ht350, h_topjets_ht350;
    unique_ptr<Hists> h_0bJetLoose, h_jets_0bJetLoose, h_ele_0bJetLoose, h_mu_0bJetLoose, h_event_0bJetLoose, h_topjets_0bJetLoose, h_btageff_0bJetLoose;
    unique_ptr<Hists> h_real_ele, h_jets_real_ele, h_ele_real_ele, h_mu_real_ele, h_event_real_ele, h_topjets_real_ele;
    unique_ptr<Hists> h_Npvs15, h_jets_Npvs15, h_ele_Npvs15, h_mu_Npvs15, h_event_Npvs15, h_topjets_Npvs15;
    unique_ptr<Hists> h_Npvs30, h_jets_Npvs30, h_ele_Npvs30, h_mu_Npvs30, h_event_Npvs30, h_topjets_Npvs30;
    unique_ptr<Hists> h_ht700, h_jets_ht700, h_ele_ht700, h_mu_ht700, h_event_ht700, h_topjets_ht700;
    unique_ptr<Hists> h_ht1050, h_jets_ht1050, h_ele_ht1050, h_mu_ht1050, h_event_ht1050, h_topjets_ht1050;
    unique_ptr<Hists> h_htmax, h_jets_htmax, h_ele_htmax, h_mu_htmax, h_event_htmax, h_topjets_htmax;
    unique_ptr<Hists> h_htjets500, h_jets_htjets500, h_ele_htjets500, h_mu_htjets500, h_event_htjets500, h_topjets_htjets500;
    unique_ptr<Hists> h_htjetsmax, h_jets_htjetsmax, h_ele_htjetsmax, h_mu_htjetsmax, h_event_htjetsmax, h_topjets_htjetsmax;
    unique_ptr<Hists> h_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection;

    MuonId MuId;
    ElectronId EleId;

    bool is_mc, isSR, isDY, isDiboson, is_Diboson_SF, is_apply_SF;
    string Sys_PU, Sys_BTag;

    JetId Btag_loose;

    // create handles to pass  information to LQToTopMuHists.cxx
    vector<unique_ptr<AnalysisModule>> recomodules;
    uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps;

    uhh2::Event::Handle<vector<LorentzVector>> h_jet_ele_ele;
    uhh2::Event::Handle<vector<LorentzVector>> h_jet_ele_fele;
    uhh2::Event::Handle<vector<LorentzVector>> h_jet_fele_fele;

    uhh2::Event::Handle<vector<LorentzVector>> h_ele;
    uhh2::Event::Handle<vector<LorentzVector>> h_jet;

    uhh2::Event::Handle<vector<LorentzVector>> h_fele_nogen;
    uhh2::Event::Handle<vector<LorentzVector>> h_fjet_nogen;
    uhh2::Event::Handle<vector<LorentzVector>> h_fele_nogen_notZ;
    uhh2::Event::Handle<vector<LorentzVector>> h_fjet_nogen_notZ;
    uhh2::Event::Handle<vector<LorentzVector>> h_fele_nogen_inZ;
    uhh2::Event::Handle<vector<LorentzVector>> h_fjet_nogen_inZ;

    uhh2::Event::Handle<vector<LorentzVector>> h_ele_close;
    uhh2::Event::Handle<vector<LorentzVector>> h_fele_nogen_close;
    uhh2::Event::Handle<vector<LorentzVector>> h_ele_notclose;
    uhh2::Event::Handle<vector<LorentzVector>> h_fele_nogen_notclose;

    uhh2::Event::Handle<vector<LorentzVector>> h_jet_close;
    uhh2::Event::Handle<vector<LorentzVector>> h_fjet_nogen_close;
    uhh2::Event::Handle<vector<LorentzVector>> h_jet_notclose;
    uhh2::Event::Handle<vector<LorentzVector>> h_fjet_nogen_notclose;

    uhh2::Event::Handle<vector<double>> h_M_fele_ele;
    uhh2::Event::Handle<vector<double>> h_M_ele_ele_wofele;

    uhh2::Event::Handle<vector<double>> h_weights;
    uhh2::Event::Handle<bool> h_is_pt_scale;


    //Get histo of scalefactors for fakerate
    TFile *scaling = new TFile("/nfs/dust/cms/user/skottkej/LQToTopMu/Run2_80X_v2/Optimization/SF_FakeEle.root", "READ");
    TGraphAsymmErrors* scalefactor= (TGraphAsymmErrors*)scaling->Get("ratio_DATA_MC_nice");
    bool berror;
  };


  LQToTopMuAnalysisModule::LQToTopMuAnalysisModule(Context & ctx){
    berror=false;
    if(berror)  cout << "Hello World from LQToTopMuAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()){
      if(berror) cout << " " << kv.first << " = " << kv.second << endl;
    }

    is_mc = ctx.get("dataset_type") == "MC";

    if(ctx.get("dataset_version") == "Diboson_WW" || ctx.get("dataset_version") == "Diboson_WZ" || ctx.get("dataset_version") == "Diboson_ZZ") {
      is_Diboson_SF = true;
    }
    else {
      is_Diboson_SF = false; 
    }

    Sys_PU = ctx.get("Systematic_PU");
    Sys_BTag = ctx.get("Systematic_BTag");
    isSR = ctx.get("use_SR") == "true";
    isDY = (ctx.get("channel") == "DY" || ctx.get("channel") == "dy");
    isDiboson = (ctx.get("channel") == "Diboson" || ctx.get("channel") == "diboson");
    if(!isDY && !isDiboson) throw runtime_error("In AnalysisModulce.cxx: Invalid channel speficied, must be: 'DY', 'dy', 'Diboson', or 'diboson'.");
    is_apply_SF = ctx.get("apply_SF") == "true";
    

    // 1. setup other modules. CommonModules and the JetCleaner:

    /* Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
       Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
       GenParticles_printer.reset(new GenParticlesPrinter(ctx));*/

    MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4),MuonIso(0.15));
    if(isDY) EleId = AndId<Electron>(ElectronID_Spring15_25ns_loose,PtEtaCut(30.0, 2.4));
    if(isDiboson) EleId = AndId<Electron>(ElectronID_Spring15_25ns_tight,PtEtaCut(30.0, 2.4));

    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);

    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.5));

    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->init(ctx,Sys_PU);


    h_hyps                = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");

    h_jet_ele_ele         = ctx.declare_event_output<vector<LorentzVector>>("h_jet_ele_ele");
    h_jet_ele_fele        = ctx.declare_event_output<vector<LorentzVector>>("h_jet_ele_fele");
    h_jet_fele_fele       = ctx.declare_event_output<vector<LorentzVector>>("h_jet_fele_fele");

    h_ele                 = ctx.declare_event_output<vector<LorentzVector>>("h_ele");
    h_jet                 = ctx.declare_event_output<vector<LorentzVector>>("h_jet");

    h_fele_nogen          = ctx.declare_event_output<vector<LorentzVector>>("h_fele_nogen");
    h_fjet_nogen          = ctx.declare_event_output<vector<LorentzVector>>("h_fjet_nogen");
    h_fele_nogen_notZ     = ctx.declare_event_output<vector<LorentzVector>>("h_fele_nogen_notZ");
    h_fjet_nogen_notZ     = ctx.declare_event_output<vector<LorentzVector>>("h_fjet_nogen_notZ");
    h_fele_nogen_inZ      = ctx.declare_event_output<vector<LorentzVector>>("h_fele_nogen_inZ");
    h_fjet_nogen_inZ      = ctx.declare_event_output<vector<LorentzVector>>("h_fjet_nogen_inZ");

    h_ele_close           = ctx.declare_event_output<vector<LorentzVector>>("h_ele_close");
    h_fele_nogen_close    = ctx.declare_event_output<vector<LorentzVector>>("h_fele_nogen_close");
    h_ele_notclose        = ctx.declare_event_output<vector<LorentzVector>>("h_ele_notclose");
    h_fele_nogen_notclose = ctx.declare_event_output<vector<LorentzVector>>("h_fele_nogen_notclose");

    h_jet_close           = ctx.declare_event_output<vector<LorentzVector>>("h_jet_close");
    h_fjet_nogen_close    = ctx.declare_event_output<vector<LorentzVector>>("h_fjet_nogen_close");
    h_jet_notclose        = ctx.declare_event_output<vector<LorentzVector>>("h_jet_notclose");
    h_fjet_nogen_notclose = ctx.declare_event_output<vector<LorentzVector>>("h_fjet_nogen_notclose");

    h_M_fele_ele          = ctx.declare_event_output<vector<double>>("h_M_fele_ele");
    h_M_ele_ele_wofele    = ctx.declare_event_output<vector<double>>("h_M_ele_ele_wofele");

    h_weights             = ctx.declare_event_output<vector<double>>("h_weights");
    h_is_pt_scale         = ctx.declare_event_output<bool>("h_is_pt_scale");



    // 2. set up selections
    //Selection
    m_eleele_veto.reset(new InvMass2EleVeto_Inverted(71, 111));
    ele_sel.reset(new NElectronSelection(2, -1));
    nele_sel.reset(new NElectronSelection(3, 3));
    if(isDY) met_sel.reset(new METSelection(0, 60));
    if(isDiboson) met_sel.reset(new METSelection(0, 30));
    htlept_sel.reset(new HTLeptSelection(0, -1));
    ht_sel.reset(new HtSelection(350));
    nbtag_loose_sel.reset(new NJetSelection(0, 2, Btag_loose)); 


    dR_max = 0.4;
    dR_gen_reco = 0.1;
    Z_mass = 91.2;

    //scalefactors for different numbers of bJets
    sf_0bjet=1.12482;
    sf_1bjet=3.13195;
    sf_2bjet=2.10099;

    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));

    //SF modules
    Ele_SF_Reco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/skottkej/CMSSW_8_0_20/src/UHH2/common/data/egammaEffi.txt_EGM2D_Reconstruction_scale_factor.root", 1));

    if(isDY) Ele_SF_LooseID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/skottkej/CMSSW_8_0_20/src/UHH2/common/data/egammaEffi.txt_EGM2D_LooseCutBased_ID_WP_scale_factor.root", 0));
    if(isDiboson) Ele_SF_TightID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/skottkej/CMSSW_8_0_20/src/UHH2/common/data/egammaEffi.txt_EGM2D_TightCutBased_ID_WP_scale_factor.root", 0));

    if(isDY) SF_btag.reset(new MCBTagScaleFactor(ctx,CSVBTag::WP_LOOSE,"jets",Sys_BTag,"mujets", "incl", "MCBtagEfficiencies_DY"));
    if(isDiboson) SF_btag.reset(new MCBTagScaleFactor(ctx,CSVBTag::WP_LOOSE,"jets",Sys_BTag,"mujets", "incl", "MCBtagEfficiencies_Diboson"));


    // 3. Set up Hists classes:
    h_2ele.reset(new LQToTopMuHists(ctx, "2ele", isSR));
    h_jets_2ele.reset(new JetHists(ctx, "Jets_2ele"));
    h_ele_2ele.reset(new ElectronHists(ctx, "Ele_2ele"));
    h_mu_2ele.reset(new MuonHists(ctx, "Mu_2ele"));
    h_event_2ele.reset(new EventHists(ctx, "Event_2ele"));
    h_topjets_2ele.reset(new TopJetHists(ctx, "TopJets_2ele"));

    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts", isSR));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));

    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto", isSR));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
    h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));

    h_htlept0.reset(new LQToTopMuHists(ctx, "HTLept0", isSR));
    h_jets_htlept0.reset(new JetHists(ctx, "Jets_HTLept0"));
    h_ele_htlept0.reset(new ElectronHists(ctx, "Ele_HTLept0"));
    h_mu_htlept0.reset(new MuonHists(ctx, "Mu_HTLept0"));
    h_topjets_htlept0.reset(new TopJetHists(ctx, "TopJets_HTLept0"));
    h_event_htlept0.reset(new EventHists(ctx, "Event_HTLept0"));

    h_MET60.reset(new LQToTopMuHists(ctx, "MET60", isSR));
    h_jets_MET60.reset(new JetHists(ctx, "Jets_MET60"));
    h_ele_MET60.reset(new ElectronHists(ctx, "Ele_MET60"));
    h_mu_MET60.reset(new MuonHists(ctx, "Mu_MET60"));
    h_event_MET60.reset(new EventHists(ctx, "Event_MET60"));
    h_topjets_MET60.reset(new TopJetHists(ctx, "TopJets_MET60"));

    h_NJet_2.reset(new LQToTopMuHists(ctx, "NJet_2", isSR));
    h_jets_NJet_2.reset(new JetHists(ctx, "Jets_NJet_2"));
    h_ele_NJet_2.reset(new ElectronHists(ctx, "Ele_NJet_2"));
    h_mu_NJet_2.reset(new MuonHists(ctx, "Mu_NJet_2"));
    h_event_NJet_2.reset(new EventHists(ctx, "Event_NJet_2"));
    h_topjets_NJet_2.reset(new TopJetHists(ctx, "TopJets_NJet_2"));

    h_ht350.reset(new LQToTopMuHists(ctx, "ht350", isSR));
    h_jets_ht350.reset(new JetHists(ctx, "Jets_ht350"));
    h_ele_ht350.reset(new ElectronHists(ctx, "Ele_ht350"));
    h_mu_ht350.reset(new MuonHists(ctx, "Mu_ht350"));
    h_event_ht350.reset(new EventHists(ctx, "Event_ht350"));
    h_topjets_ht350.reset(new TopJetHists(ctx, "TopJets_ht350"));

    h_0bJetLoose.reset(new LQToTopMuHists(ctx, "0bJetLoose", isSR));
    h_jets_0bJetLoose.reset(new JetHists(ctx, "Jets_0bJetLoose"));
    h_ele_0bJetLoose.reset(new ElectronHists(ctx, "Ele_0bJetLoose"));
    h_mu_0bJetLoose.reset(new MuonHists(ctx, "Mu_0bJetLoose"));
    h_event_0bJetLoose.reset(new EventHists(ctx, "Event_0bJetLoose"));
    h_topjets_0bJetLoose.reset(new TopJetHists(ctx, "TopJets_0bJetLoose"));    

    h_btageff_0bJetLoose.reset(new BTagMCEfficiencyHists(ctx, "BTagEff_0bJetLoose",CSVBTag::WP_LOOSE));

    h_3ele.reset(new LQToTopMuHists(ctx, "3Ele", isSR));
    h_jets_3ele.reset(new JetHists(ctx, "Jets_3Ele"));
    h_ele_3ele.reset(new ElectronHists(ctx, "Ele_3Ele"));
    h_mu_3ele.reset(new MuonHists(ctx, "Mu_3Ele"));
    h_event_3ele.reset(new EventHists(ctx, "Event_3Ele"));
    h_topjets_3ele.reset(new TopJetHists(ctx, "TopJets_3Ele"));

    h_dR.reset(new LQToTopMuHists(ctx, "dR", isSR));
    h_jets_dR.reset(new JetHists(ctx, "Jets_dR"));
    h_ele_dR.reset(new ElectronHists(ctx, "Ele_dR"));
    h_mu_dR.reset(new MuonHists(ctx, "Mu_dR"));
    h_event_dR.reset(new EventHists(ctx, "Event_dR"));
    h_topjets_dR.reset(new TopJetHists(ctx, "TopJets_dR"));

    h_real_ele.reset(new LQToTopMuHists(ctx, "real_ele", isSR));
    h_jets_real_ele.reset(new JetHists(ctx, "Jets_real_ele"));
    h_ele_real_ele.reset(new ElectronHists(ctx, "Ele_real_ele"));
    h_mu_real_ele.reset(new MuonHists(ctx, "Mu_real_ele"));
    h_event_real_ele.reset(new EventHists(ctx, "Event_real_ele"));
    h_topjets_real_ele.reset(new TopJetHists(ctx, "TopJets_real_ele"));

    h_Npvs15.reset(new LQToTopMuHists(ctx, "Npvs15", isSR));
    h_jets_Npvs15.reset(new JetHists(ctx, "Jets_Npvs15"));
    h_ele_Npvs15.reset(new ElectronHists(ctx, "Ele_Npvs15"));
    h_mu_Npvs15.reset(new MuonHists(ctx, "Mu_Npvs15"));
    h_event_Npvs15.reset(new EventHists(ctx, "Event_Npvs15"));
    h_topjets_Npvs15.reset(new TopJetHists(ctx, "TopJets_Npvs15"));

    h_Npvs30.reset(new LQToTopMuHists(ctx, "Npvs30", isSR));
    h_jets_Npvs30.reset(new JetHists(ctx, "Jets_Npvs30"));
    h_ele_Npvs30.reset(new ElectronHists(ctx, "Ele_Npvs30"));
    h_mu_Npvs30.reset(new MuonHists(ctx, "Mu_Npvs30"));
    h_event_Npvs30.reset(new EventHists(ctx, "Event_Npvs30"));
    h_topjets_Npvs30.reset(new TopJetHists(ctx, "TopJets_Npvs30"));

    h_ht700.reset(new LQToTopMuHists(ctx, "ht700", isSR));
    h_jets_ht700.reset(new JetHists(ctx, "Jets_ht700"));
    h_ele_ht700.reset(new ElectronHists(ctx, "Ele_ht700"));
    h_mu_ht700.reset(new MuonHists(ctx, "Mu_ht700"));
    h_event_ht700.reset(new EventHists(ctx, "Event_ht700"));
    h_topjets_ht700.reset(new TopJetHists(ctx, "TopJets_ht700"));

    h_ht1050.reset(new LQToTopMuHists(ctx, "ht1050", isSR));
    h_jets_ht1050.reset(new JetHists(ctx, "Jets_ht1050"));
    h_ele_ht1050.reset(new ElectronHists(ctx, "Ele_ht1050"));
    h_mu_ht1050.reset(new MuonHists(ctx, "Mu_ht1050"));
    h_event_ht1050.reset(new EventHists(ctx, "Event_ht1050"));
    h_topjets_ht1050.reset(new TopJetHists(ctx, "TopJets_ht1050"));

    h_htmax.reset(new LQToTopMuHists(ctx, "htmax", isSR));
    h_jets_htmax.reset(new JetHists(ctx, "Jets_htmax"));
    h_ele_htmax.reset(new ElectronHists(ctx, "Ele_htmax"));
    h_mu_htmax.reset(new MuonHists(ctx, "Mu_htmax"));
    h_event_htmax.reset(new EventHists(ctx, "Event_htmax"));
    h_topjets_htmax.reset(new TopJetHists(ctx, "TopJets_htmax"));

    h_htjets500.reset(new LQToTopMuHists(ctx, "htjets500", isSR));
    h_jets_htjets500.reset(new JetHists(ctx, "Jets_htjets500"));
    h_ele_htjets500.reset(new ElectronHists(ctx, "Ele_htjets500"));
    h_mu_htjets500.reset(new MuonHists(ctx, "Mu_htjets500"));
    h_event_htjets500.reset(new EventHists(ctx, "Event_htjets500"));
    h_topjets_htjets500.reset(new TopJetHists(ctx, "TopJets_htjets500"));

    h_htjetsmax.reset(new LQToTopMuHists(ctx, "htjetsmax", isSR));
    h_jets_htjetsmax.reset(new JetHists(ctx, "Jets_htjetsmax"));
    h_ele_htjetsmax.reset(new ElectronHists(ctx, "Ele_htjetsmax"));
    h_mu_htjetsmax.reset(new MuonHists(ctx, "Mu_htjetsmax"));
    h_event_htjetsmax.reset(new EventHists(ctx, "Event_htjetsmax"));
    h_topjets_htjetsmax.reset(new TopJetHists(ctx, "TopJets_htjetsmax"));



    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection", isSR));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection"));
  }


  bool LQToTopMuAnalysisModule::process(Event & event) {
    if(is_Diboson_SF) event.weight *= sf_0bjet;

    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);

    Nele = event.electrons->size();

    // MLQ reco
    if(nele_sel->passes(event)){
      for(auto & m : recomodules){
	m->process(event);
      }
    }

    Ele_SF_Reco->process(event);
    if(isDY) Ele_SF_LooseID->process(event);
    if(isDiboson) Ele_SF_TightID->process(event);


    double real_weight = event.weight;
    vector<double> scalefactors;
    bool is_pt_scale = false;

    if(is_mc && is_apply_SF){
      double x, y;
      scalefactor->GetPoint(0,x,y);
      scalefactors.push_back(y);
      scalefactor->GetPoint(1,x,y);
      scalefactors.push_back(y);
      scalefactor->GetPoint(2,x,y);
      scalefactors.push_back(y);
      scalefactors.push_back(real_weight);
      event.set(h_weights, scalefactors);

      event.set(h_is_pt_scale, is_pt_scale);
    }


    if(!ele_sel->passes(event)) return false;
    h_2ele->fill(event);
    h_jets_2ele->fill(event);
    h_ele_2ele->fill(event);
    h_mu_2ele->fill(event);
    h_event_2ele->fill(event);
    h_topjets_2ele->fill(event);

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);

    //InvMassVeto
    if(!m_eleele_veto->passes(event)) return false;
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    h_event_InvMassVeto->fill(event);
    h_topjets_InvMassVeto->fill(event);

    //HT Lept
    if(!htlept_sel->passes(event)) return false;
    h_jets_htlept0->fill(event);
    h_htlept0->fill(event);
    h_ele_htlept0->fill(event);
    h_mu_htlept0->fill(event);
    h_event_htlept0->fill(event);
    h_topjets_htlept0->fill(event);

    //MET 60
    if(isDY){
      if(!met_sel->passes(event)) return false;
    }
    if(isDiboson){
      if(met_sel->passes(event)) return false;
    }
    h_jets_MET60->fill(event);
    h_MET60->fill(event);
    h_ele_MET60->fill(event);
    h_mu_MET60->fill(event);
    h_event_MET60->fill(event);
    h_topjets_MET60->fill(event);


    vector<bool> is_real;
    vector<Electron> fake_electrons;
    if(is_mc){
      /*	Jet_printer->process(event);
		Electron_printer->process(event);
		GenParticles_printer->process(event); */

      //set up array containing information wether a specific genparticle has been matched to a reco-ele
      const int Ngenp = event.genparticles->size();
      bool genp_used[Ngenp];
      for(int i=0; i<Ngenp; i++){
	genp_used[i] = false;
      }
      //set up array containing information wether a specific electron is real (has a gen particle) or not (has NO gen particle)
      for(int i=0; i<Nele; i++){
	is_real.push_back(true);
      }
      if(is_real.size() != event.electrons->size()) throw runtime_error("In AnalysisModule.cxx: Vector containing information on (real,fake) has not the same size as event.electrons.");

      //match reco-ele with a gen-ele within R=0.1  (R is variable, you can change it at "2. set up selections")
      for(int i=0; i<Nele; i++) {
	Electron reco_ele = event.electrons->at(i);
	double dr_min = 999999;
	int idx_closest_genele = -1;
	for(unsigned int j=0; j<event.genparticles->size(); j++) {
	  if(genp_used[j]) continue;
	  GenParticle gen_ele = event.genparticles->at(j);
	  if(abs(gen_ele.pdgId()) == 11){
	    double dr_tmp = deltaR(reco_ele,gen_ele);
	    // cout << "dr_tmp for reco-ele no " << i << " and gen-part no " << j << ": " << dr_tmp << endl;
	    if(dr_tmp < dr_min){
	      dr_min = dr_tmp;
	      if(dr_min < dR_gen_reco){
		idx_closest_genele = j;
	      }
	    }
	  }
	}
	genp_used[idx_closest_genele] = true;
	//save fake-eles
	// cout << "dR min when matching reco with gen eles for ele no " << i  << ": " << dr_min << endl;
	if(dr_min > dR_gen_reco){
	  fake_electrons.push_back(reco_ele);
	  is_real.at(i) = false;
	}
      }
    }

    //get the pair of electrons which are the closest to the Z-mass
    double diff_min = 999999;
    int Z_ele1=-1;
    int Z_ele2=-1;

    for(int i=0; i<Nele; i++){
      for(int j=0; j<Nele; j++){
	if(j>i){
	  if(event.electrons->at(i).charge() != event.electrons->at(j).charge()){
	    double M_ee = (event.electrons->at(i).v4() + event.electrons->at(j).v4()).M();
	    double diff = abs(M_ee-Z_mass);
	    if(diff < diff_min){
	      diff_min = diff;
	      Z_ele1=i;
	      Z_ele2=j;
	    }
	  }
	}
      }
    }

    //set up array containing information whether a specific electron is part of the pair closest to the Z-mass
    bool ele_Z[Nele];
    for(int i=0; i<Nele; i++){
      ele_Z[i] = false;
    }
    
    ele_Z[Z_ele1] = true;
    ele_Z[Z_ele2] = true;
    //set up array containing information whether a specific jet got matched with an electron which is part of the pair closest to the Z-mass
    bool jet_hypo[event.jets->size()];
    for(unsigned int i=0; i<event.jets->size(); i++){
      jet_hypo[i] = false;
    }

    //get the jets which are the closest to the respective electrons of the closest pair to the Z-mass
    vector<double> dr_mins;
    vector<int> idx_closest_jet;
    vector<LorentzVector> closest_jet;
    if(event.jets->size() > 0){
      for(int i=0; i<Nele; i++){
	double dr_min = 99999999;
	int idx_jet = -1;
	for(unsigned int j=0; j<event.jets->size(); j++){
	  double dr_tmp = deltaR(event.electrons->at(i),event.jets->at(j));
	  if(dr_tmp < dr_min){
	    dr_min = dr_tmp;
	    idx_jet = j;
	  }
	}
	dr_mins.push_back(dr_min);
	idx_closest_jet.push_back(idx_jet);
	if(idx_jet != -1) closest_jet.push_back(event.jets->at(idx_jet).v4());
	else throw runtime_error("AnalysisModule.cxx: Not working as intended!");
	if((isDY && ele_Z[i] && dr_min <= 0.1) || (!isDY && dr_min <= 0.1)){
	  jet_hypo[idx_jet] = true;
	}
      }  
    }

    if(is_mc && event.jets->size()>0){
      vector<LorentzVector> jet_ele_ele, jet_ele_fele, jet_fele_fele;
      if(Nele == 2){
	for(int i=0; i<Nele; i++){
	  if(is_real.at(0) && is_real.at(1)){
	    if(jet_hypo[idx_closest_jet.at(i)]){
	      jet_ele_ele.push_back(closest_jet.at(i));
	      jet_ele_fele.push_back({0,0,0,0});
	      jet_fele_fele.push_back({0,0,0,0});
	    }
	  }
	  if((!is_real.at(0) && is_real.at(1)) || (is_real.at(0) && !is_real.at(1))){
	    if(jet_hypo[idx_closest_jet.at(i)]){
	      jet_ele_ele.push_back({0,0,0,0});
	      jet_ele_fele.push_back(closest_jet.at(i));
	      jet_fele_fele.push_back({0,0,0,0});
	    }
	  }
	  if(!is_real.at(0) && !is_real.at(1)){
	    if(jet_hypo[idx_closest_jet.at(i)]){
	      jet_ele_ele.push_back({0,0,0,0});
	      jet_ele_fele.push_back({0,0,0,0});
	      jet_fele_fele.push_back(closest_jet.at(i));
	    }
	  }
	}
      }
      event.set(h_jet_ele_ele, jet_ele_ele);
      event.set(h_jet_ele_fele, jet_ele_fele);
      event.set(h_jet_fele_fele, jet_fele_fele);
    }
    //Delete the jets which are closest to the electrons in the best hypothesis and set new event.jets
    vector<Jet> event_jets_new;
    if(event.jets->size()>0){
      for(int j=0; j<Nele; j++){
	// cout << "Distance between Electron no. " << j << " and closest jet: " << dr_mins.at(j) << ", which is jet no. " << idx_closest_jet.at(j) << ", is in Z-pair?: " <<ele_Z[j] <<  endl;
      }
      for(unsigned int i=0; i<event.jets->size(); i++){
	if(!jet_hypo[i]){
	  event_jets_new.push_back(event.jets->at(i));
	}
      }
      swap(event_jets_new, *event.jets);
    }



    if(is_Diboson_SF){
      std::vector<Jet> bjets_loose;
      for (unsigned int i =0; i<event.jets->size(); ++i) {
      	if(Btag_loose(event.jets->at(i),event)) {
    	  bjets_loose.push_back(event.jets->at(i));
    	}
      }

      int NbJets_loose = bjets_loose.size();

      if(NbJets_loose == 1) event.weight *= sf_1bjet;
      if(NbJets_loose == 2) event.weight *= sf_2bjet;
    }



    //NJet >= 2
    if(event.jets->size() < 2 && isDY) return false;
    h_jets_NJet_2->fill(event);
    h_NJet_2->fill(event);
    h_ele_NJet_2->fill(event);
    h_mu_NJet_2->fill(event);
    h_event_NJet_2->fill(event);
    h_topjets_NJet_2->fill(event);

    //HT >= 350
    if(!ht_sel->passes(event) && isDY) return false;
    h_jets_ht350->fill(event);
    h_ht350->fill(event);
    h_ele_ht350->fill(event);
    h_mu_ht350->fill(event);
    h_event_ht350->fill(event);
    h_topjets_ht350->fill(event);

    //NbJet
    if(!nbtag_loose_sel->passes(event)) return false;
    SF_btag->process(event);

    if(is_mc && is_apply_SF){      
      scalefactors.at(3) = event.weight;
      event.set(h_weights, scalefactors);

      is_pt_scale = true;
      event.set(h_is_pt_scale, is_pt_scale);
    }

    h_jets_0bJetLoose->fill(event);
    h_0bJetLoose->fill(event);
    h_ele_0bJetLoose->fill(event);
    h_mu_0bJetLoose->fill(event);
    h_event_0bJetLoose->fill(event);
    h_topjets_0bJetLoose->fill(event);
    h_btageff_0bJetLoose->fill(event);

    //Min 3 Ele
    if(!nele_sel->passes(event)) return false;
    h_jets_3ele->fill(event);
    h_3ele->fill(event);
    h_ele_3ele->fill(event);
    h_mu_3ele->fill(event);
    h_event_3ele->fill(event);
    h_topjets_3ele->fill(event);

    int fele_size=0;
    if(is_mc) {
      //match each fake-ele (determined on gen-lvl) with a jet within R=0.4 (this jet is responsible for faking that electron)   (R is variable, you can change it at "2. set up selections")
      vector<Jet> faking_jets;
      vector<int> idx_faking_jet;
      for(unsigned int i=0; i<fake_electrons.size(); i++){
	double dr_min = 999999;
	int idx_matching_jet = -1;
	for(unsigned int j=0; j<event.jets->size(); j++){
	  double dr_tmp = deltaR(event.jets->at(j), fake_electrons.at(i));
	  if(dr_tmp < dr_min){
	    dr_min = dr_tmp;
	    idx_matching_jet = j;
	  }
	}
	//save jets corresponding to fake electrons
	if(dr_min < dR_max){
	  faking_jets.push_back(event.jets->at(idx_matching_jet));
	  idx_faking_jet.push_back(idx_matching_jet);
	}
	else{
	  idx_faking_jet.push_back(-1);
	}
      }
      if(idx_faking_jet.size() != fake_electrons.size()) throw runtime_error("In AnalysisModule.cxx: faking_jets does not have the same size as fake_electrons.");

      //look for invariant mass 71 < M < 111
      vector<LorentzVector> fele_nogen_inZ, fjet_nogen_inZ, fele_nogen_notZ, fjet_nogen_notZ;
      vector<double> M_fele_ele, M_ele_ele_wofele;
      bool is_Z[fake_electrons.size()];
      for(unsigned int i=0; i<fake_electrons.size(); i++){
	is_Z[i] = false;
      }
      // consider pairs of 1 fake and 1 real electron
      for(unsigned int i=0; i<fake_electrons.size(); i++){
	double M_feleele;
	for(int j=0; j<Nele; j++){
	  if(fake_electrons[i].v4() != event.electrons->at(j).v4()){
	    M_feleele = (fake_electrons[i].v4() + event.electrons->at(j).v4()).M();
	    M_fele_ele.push_back(M_feleele);
	    if(71 < M_feleele && M_feleele < 111){
	      is_Z[i]=true;
	    }
	  }
	}
      }
      //consider only pairs of 2 real electrons
      for(int j=0; j<Nele; j++){
	for(int k=0; k<Nele; k++){
	  if(is_real[j] && is_real[k]){
	    if(k>j){
	      double M_eleele = (event.electrons->at(j).v4() + event.electrons->at(k).v4()).M();
	      M_ele_ele_wofele.push_back(M_eleele);
	    }
	  }
	}
      }

      //set handle for electrons and jets
      for(unsigned int i=0; i<fake_electrons.size(); i++) {
	if(is_Z[i]){
	  fele_nogen_inZ.push_back(fake_electrons[i].v4());
	  if(idx_faking_jet[i] != -1) fjet_nogen_inZ.push_back(event.jets->at(idx_faking_jet.at(i)).v4());
	  else fjet_nogen_inZ.push_back({-1,-1,-1,-1});
	  fele_nogen_notZ.push_back({0,0,0,0});
	  fjet_nogen_notZ.push_back({0,0,0,0});
	}
	else{
	  fele_nogen_notZ.push_back(fake_electrons[i].v4());
	  if(idx_faking_jet[i] != -1) fjet_nogen_notZ.push_back(event.jets->at(idx_faking_jet.at(i)).v4());
	  else fjet_nogen_notZ.push_back({-1,-1,-1,-1});
	  fele_nogen_inZ.push_back({0,0,0,0});
	  fjet_nogen_inZ.push_back({0,0,0,0});
	}
      }

      event.set(h_fele_nogen_inZ, fele_nogen_inZ);
      event.set(h_fjet_nogen_inZ, fjet_nogen_inZ);
      event.set(h_fele_nogen_notZ, fele_nogen_notZ);
      event.set(h_fjet_nogen_notZ, fjet_nogen_notZ);


      //set handle for eles w/o gen-ele
      vector<LorentzVector> fele_nogen, fjet_nogen;
      // cout << "size fake-ele: " << fake_electrons.size() << ", faking jets size: " << faking_jets.size() << endl;
      for(unsigned int i=0; i<fake_electrons.size(); i++){
	fele_nogen.push_back(fake_electrons[i].v4());
	if(idx_faking_jet[i] != -1) fjet_nogen.push_back(event.jets->at(idx_faking_jet.at(i)).v4());
	else fjet_nogen.push_back({-1,-1,-1,-1});
	// cout << "dR fake-ele -- faking jet for fake-ele no. " << i << ": " << deltaR(fake_electrons[i], faking_jets[i]) << endl;
      }
      event.set(h_fele_nogen, fele_nogen);
      event.set(h_fjet_nogen, fjet_nogen);

      event.set(h_M_fele_ele, M_fele_ele);
      event.set(h_M_ele_ele_wofele, M_ele_ele_wofele);
      fele_size = fake_electrons.size();
    }





    //look for invariant mass closest to Z mass (91.2)
    vector<LorentzVector> fele_nogen_notclose, fjet_nogen_notclose, fele_nogen_close, fjet_nogen_close;
    vector<LorentzVector> ele_notclose, jet_notclose, ele_close, jet_close;
    vector<LorentzVector> ele, jet;

    if(is_mc){
      for(int i=0; i<Nele; i++){
	if(!ele_Z[i]){
	  if(!is_real[i]){
	    fele_nogen_notclose.push_back(event.electrons->at(i).v4());
	    if(event.jets->size() > 0) fjet_nogen_notclose.push_back(closest_jet.at(i));
	    else fjet_nogen_notclose.push_back({0,0,0,0});
	    ele_notclose.push_back({0,0,0,0});
	    jet_notclose.push_back({0,0,0,0});
	    fele_nogen_close.push_back({0,0,0,0});
	    fjet_nogen_close.push_back({0,0,0,0});
	    ele_close.push_back({0,0,0,0});
	    jet_close.push_back({0,0,0,0});

	    //cout << "Ele no. " << i << " is a fake ele and not part of the 'best pair'. The closest jet to this ele is jet no. " << idx_closest_jet.at(i) << " with distance " << dr_mins.at(i) << endl;
	  }
	  else{
	    ele_notclose.push_back(event.electrons->at(i).v4());
	    if(event.jets->size() > 0) jet_notclose.push_back(closest_jet.at(i));
	    else jet_notclose.push_back({0,0,0,0});
	    fele_nogen_notclose.push_back({0,0,0,0});
	    fjet_nogen_notclose.push_back({0,0,0,0});
	    fele_nogen_close.push_back({0,0,0,0});
	    fjet_nogen_close.push_back({0,0,0,0});
	    ele_close.push_back({0,0,0,0});
	    jet_close.push_back({0,0,0,0});

	    // cout << "Ele no. " << i << " is a real ele and not part of the 'best pair'. The mass difference to 91.2 is: " << diff_min << ". The closest jet to this ele is jet no. " << idx_closest_jet.at(i) << " with distance " << dr_mins.at(i) << endl;
	  }
	}
	else{
	  if(!is_real[i]){
	    fele_nogen_close.push_back(event.electrons->at(i).v4());
	    fjet_nogen_close.push_back({0,0,0,0});
	    fele_nogen_notclose.push_back({0,0,0,0});
	    fjet_nogen_notclose.push_back({0,0,0,0});
	    ele_notclose.push_back({0,0,0,0});
	    jet_notclose.push_back({0,0,0,0});
	    ele_close.push_back({0,0,0,0});
	    jet_close.push_back({0,0,0,0});

	    //cout << "Ele no. " << i << " is a fake ele and IS part of the 'best pair'. The closest jet to this ele is jet no. " << idx_closest_jet.at(i) << " with distance " << dr_mins.at(i) << endl;
	  }
	  else{
	    ele_close.push_back(event.electrons->at(i).v4());
	    jet_close.push_back({0,0,0,0});
	    fele_nogen_notclose.push_back({0,0,0,0});
	    fjet_nogen_notclose.push_back({0,0,0,0});
	    ele_notclose.push_back({0,0,0,0});
	    jet_notclose.push_back({0,0,0,0});
	    fele_nogen_close.push_back({0,0,0,0});
	    fjet_nogen_close.push_back({0,0,0,0});

	    //cout << "Ele no. " << i << " is a real ele and IS part of the 'best pair'. The closest jet to this ele is jet no. " << idx_closest_jet.at(i) << " with distance " << dr_mins.at(i) << endl;
	  }
	}
      }

      event.set(h_fele_nogen_notclose, fele_nogen_notclose);
      event.set(h_ele_notclose, ele_notclose);
      event.set(h_fele_nogen_close, fele_nogen_close);
      event.set(h_ele_close, ele_close);

      event.set(h_fjet_nogen_notclose, fjet_nogen_notclose);
      event.set(h_jet_notclose, jet_notclose);
      event.set(h_fjet_nogen_close, fjet_nogen_close);
      event.set(h_jet_close, jet_close);
    }

    for(int i=0; i<Nele; i++){
      if(!ele_Z[i]){
	ele.push_back(event.electrons->at(i).v4());
	if(event.jets->size() > 0) jet.push_back(closest_jet.at(i));
	else jet.push_back({0,0,0,0});
      }
      else {
	ele.push_back({0,0,0,0});
	jet.push_back({0,0,0,0});
      }
    }
    event.set(h_ele, ele);
    event.set(h_jet, jet);
    


    //dR <= 0.4
    h_jets_dR->fill(event);
    h_dR->fill(event);
    h_ele_dR->fill(event);
    h_mu_dR->fill(event);
    h_event_dR->fill(event);
    h_topjets_dR->fill(event);

    if(fele_size == 0){
      //only real ele
      h_jets_real_ele->fill(event);
      h_real_ele->fill(event);
      h_ele_real_ele->fill(event);
      h_mu_real_ele->fill(event);
      h_event_real_ele->fill(event);
      h_topjets_real_ele->fill(event);
    }

    int Npvs = event.pvs->size();
    if(Npvs <= 15){
      h_jets_Npvs15->fill(event);
      h_Npvs15->fill(event);
      h_ele_Npvs15->fill(event);
      h_mu_Npvs15->fill(event);
      h_event_Npvs15->fill(event);
      h_topjets_Npvs15->fill(event); 
    }
    if(Npvs > 15){
      h_jets_Npvs30->fill(event);
      h_Npvs30->fill(event);
      h_ele_Npvs30->fill(event);
      h_mu_Npvs30->fill(event);
      h_event_Npvs30->fill(event);
      h_topjets_Npvs30->fill(event);   
    }

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

    if(isSR){
      for(const auto & muon : *event.muons){
	ht_lep += muon.pt();
      }
    }
    ht = ht_lep + ht_jets + met;

    if(ht >= 350 && ht <=700){
      h_jets_ht700->fill(event);
      h_ht700->fill(event);
      h_ele_ht700->fill(event);
      h_mu_ht700->fill(event);
      h_event_ht700->fill(event);
      h_topjets_ht700->fill(event);
    }
    if(ht > 700 && ht <=1050){
      h_jets_ht1050->fill(event);
      h_ht1050->fill(event);
      h_ele_ht1050->fill(event);
      h_mu_ht1050->fill(event);
      h_event_ht1050->fill(event);
      h_topjets_ht1050->fill(event);
    }
    if(ht > 1050){
      h_jets_htmax->fill(event);
      h_htmax->fill(event);
      h_ele_htmax->fill(event);
      h_mu_htmax->fill(event);
      h_event_htmax->fill(event);
      h_topjets_htmax->fill(event);
    }

    if(ht_jets <= 500){
      h_jets_htjets500->fill(event);
      h_htjets500->fill(event);
      h_ele_htjets500->fill(event);
      h_mu_htjets500->fill(event);
      h_event_htjets500->fill(event);
      h_topjets_htjets500->fill(event);
    }
    if(ht_jets > 500){
      h_jets_htjetsmax->fill(event);
      h_htjetsmax->fill(event);
      h_ele_htjetsmax->fill(event);
      h_mu_htjetsmax->fill(event);
      h_event_htjetsmax->fill(event);
      h_topjets_htjetsmax->fill(event);
    }

    //Final Selection
 
    h_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);

    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuAnalysisModule)

}
  

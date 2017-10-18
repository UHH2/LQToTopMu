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
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
#include "UHH2/LQToTopMu/include/MET2dHists.h"
#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstruction.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"


using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    unique_ptr<CommonModules> common;
    unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer, GenParticles_printer, syst_module, SF_muonID, SF_muonIso, SF_btag, SF_eleReco, SF_eleID;
    unique_ptr<MCMuonScaleFactor> SF_muonTrigger;
    unique_ptr<ElectronFakeRateWeights> SF_eleFakeRate;
    unique_ptr<MuonTrkWeights> SF_muonTrk;

    unique_ptr<JetCleaner> jetcleaner;
    
    // declare the Selections to use.
    unique_ptr<Selection>  nele_sel, njet_sel, nbtag_loose_sel, m_mumu_veto, htlept_sel, mttbar_gen_sel, m_muele_veto, m_ee_veto, ht_sel, dr_lepjet_sel, genlvl_TopDilepton_sel, ele_trigger_sel1, ele_trigger_sel2, genlvl_ZMuMu_sel, genlvl_ZEE_sel;
    
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts, h_eff_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_event_1ele, h_topjets_1ele, h_tau_1ele, h_eff_1ele, h_lumi_1ele;
    unique_ptr<Hists> h_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_event_1bJetLoose, h_topjets_1bJetLoose, h_tau_1bJetLoose, h_eff_1bJetLoose, h_lumi_1bJetLoose;
    unique_ptr<Hists> h_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto, h_tau_InvMassVeto, h_eff_InvMassVeto, h_lumi_InvMassVeto;
    unique_ptr<Hists> h_htlept200, h_jets_htlept200, h_ele_htlept200, h_mu_htlept200, h_event_htlept200, h_topjets_htlept200, h_tau_htlept200, h_btageff_htlept200, h_eff_htlept200, h_lumi_htlept200;
    unique_ptr<Hists> h_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection, h_tau_finalSelection, h_eff_finalSelection, h_lumi_finalSelection;
    unique_ptr<Hists> h_ht_InvMassVeto, h_ht_finalSelection;
    unique_ptr<Hists> h_Sideband;
    unique_ptr<Hists> h_PDF_variations;

    
    MuonId MuId;
    ElectronId EleId;
    JetId Btag_loose, Btag_medium, Btag_tight;
    CSVBTag::wp wp_btag_loose;
    int n_misid;


    bool do_scale_variation, is_mc, do_pdf_variations, is_foreff;
    string Sys_MuonID, Sys_BTag, Sys_MuonTrigger, Sys_MuonIso, Sys_PU, Sys_EleID, Sys_EleReco, Sys_WJ, Sys_TTbar, Sys_DY, Sys_ST, Sys_DB, Sys_QCD, Sys_TTV, Sys_EleFake, Sys_MuFake, Sys_MuonTrk;
    TString dataset_version;
 
    
    vector<unique_ptr<AnalysisModule>> recomodules, muonic_recomodules;
    uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps, h_muonic_hyps;

    uhh2::Event::Handle<double> h_FakeRateWeightEle, h_FakeRateWeightMu;
    
  };
  
  
  LQToTopMuAnalysisModule::LQToTopMuAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    n_misid = 0;
    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    dataset_version = ctx.get("dataset_version");
    is_mc = ctx.get("dataset_type") == "MC";
    is_foreff = ctx.get("IsForEff") == "true";
    Sys_MuonID = ctx.get("Systematic_MuonID");
    Sys_MuonTrigger = ctx.get("Systematic_MuonTrigger");
    Sys_MuonIso = ctx.get("Systematic_MuonIso");
    Sys_MuonTrk = ctx.get("Systematic_MuonTrk");
    Sys_MuFake = ctx.get("Systematic_MuFake");
    Sys_EleID = ctx.get("Systematic_EleID");
    Sys_EleReco = ctx.get("Systematic_EleReco");
    Sys_EleFake = ctx.get("Systematic_EleFake");
    Sys_BTag = ctx.get("Systematic_BTag");
    Sys_PU = ctx.get("Systematic_PU");
    Sys_TTbar = ctx.get("Systematic_TTbar");
    Sys_DY = ctx.get("Systematic_DY");
    Sys_ST = ctx.get("Systematic_ST");
    Sys_DB = ctx.get("Systematic_DB");
    Sys_QCD = ctx.get("Systematic_QCD");
    Sys_TTV = ctx.get("Systematic_TTV");
    Sys_WJ = ctx.get("Systematic_WJ");


    // 1. setup other modules. CommonModules and the JetCleaner:
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
    Muon_printer.reset(new MuonPrinter("Muon-Printer"));
    GenParticles_printer.reset(new GenParticlesPrinter(ctx));
    MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4),MuonIso(0.15));
    EleId = AndId<Electron>(ElectronID_Spring16_loose,PtEtaCut(30.0, 2.4));

    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
    wp_btag_loose = CSVBTag::WP_LOOSE;
    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));  // 2.5!!!


    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->init(ctx,Sys_PU);

    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, Sys_MuonID));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, Sys_MuonTrigger));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, Sys_MuonIso));
    SF_muonTrk.reset(new MuonTrkWeights(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", Sys_MuonTrk));

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", Sys_EleReco));
    SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", Sys_EleID));

    SF_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_loose,"jets",Sys_BTag));

    
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
    h_muonic_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassMuonicLQReconstruction");
    if(Sys_EleFake == "nominal")   h_FakeRateWeightEle =ctx.get_handle<double>("FakeRateWeightEle");
    else if(Sys_EleFake == "up")   h_FakeRateWeightEle =ctx.get_handle<double>("FakeRateWeightEleUp");
    else if(Sys_EleFake == "down") h_FakeRateWeightEle =ctx.get_handle<double>("FakeRateWeightEleDown");
    else throw runtime_error("Sys_EleFake is not one of the following: ['up', 'down', 'nominal']");
    if(Sys_MuFake == "nominal")   h_FakeRateWeightMu =ctx.get_handle<double>("FakeRateWeightMu");
    else if(Sys_MuFake == "up")   h_FakeRateWeightMu =ctx.get_handle<double>("FakeRateWeightMuUp");
    else if(Sys_MuFake == "down") h_FakeRateWeightMu =ctx.get_handle<double>("FakeRateWeightMuDown");
    else throw runtime_error("Sys_MuFake is not one of the following: ['up', 'down', 'nominal']");

    
    // 2. set up selections
    //Selection
    njet_sel.reset(new NJetSelection(2, -1));  
    nbtag_loose_sel.reset(new NJetSelection(1, -1, Btag_loose));  //1, -1 
    m_mumu_veto.reset(new InvMass2MuVeto(0, 111));
    m_ee_veto.reset(new InvMassEleEleVeto(0.,111.));
    nele_sel.reset(new NElectronSelection(1, -1));
    htlept_sel.reset(new HTLeptSelection(200., -1));
    ht_sel.reset(new HtSelection(350., -1));
    genlvl_ZMuMu_sel.reset(new GenLvlZMuMuSelection());
    genlvl_ZEE_sel.reset(new GenLvlZEESelection());


    //ele_trigger_sel1.reset(new TriggerSelection("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"));
    ele_trigger_sel1.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    ele_trigger_sel2.reset(new TriggerSelection("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"));
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));
    muonic_recomodules.emplace_back(new HighMassMuonicLQReconstruction(ctx,LQNeutrinoReconstruction));
    muonic_recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassMuonicLQReconstruction"));

    //systematics modules
    syst_module.reset(new MCScaleVariation(ctx));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
    h_eff_nocuts.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));
 
    h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
    h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
    h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
    h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));  
    h_event_1ele.reset(new EventHists(ctx, "Event_1Ele"));
    h_topjets_1ele.reset(new TopJetHists(ctx, "TopJets_1Ele"));
    h_eff_1ele.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_1Ele"));
    h_lumi_1ele.reset(new LuminosityHists(ctx, "Lumi_1Ele"));

    h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
    h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
    h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
    h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));
    h_event_1bJetLoose.reset(new EventHists(ctx, "Event_1bJetLoose"));
    h_topjets_1bJetLoose.reset(new TopJetHists(ctx, "TopJets_1bJetLoose"));
    h_eff_1bJetLoose.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_1bJetLoose"));
    h_lumi_1bJetLoose.reset(new LuminosityHists(ctx, "Lumi_1bJetLoose"));


    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
    h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));
    h_eff_InvMassVeto.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_InvMassVeto"));
    h_lumi_InvMassVeto.reset(new LuminosityHists(ctx, "Lumi_InvMassVeto"));


    h_htlept200.reset(new LQToTopMuHists(ctx, "HTLept200"));
    h_jets_htlept200.reset(new JetHists(ctx, "Jets_HTLept200"));
    h_ele_htlept200.reset(new ElectronHists(ctx, "Ele_HTLept200"));
    h_mu_htlept200.reset(new MuonHists(ctx, "Mu_HTLept200"));
    h_topjets_htlept200.reset(new TopJetHists(ctx, "TopJets_HTLept200"));
    h_event_htlept200.reset(new EventHists(ctx, "Event_HTLept200"));
    h_eff_htlept200.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_HTLept200"));
    h_lumi_htlept200.reset(new LuminosityHists(ctx, "Lumi_HTLept200"));
    h_btageff_htlept200.reset(new BTagMCEfficiencyHists(ctx, "BTagEff_HTLept200",wp_btag_loose));


    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection")); 
    h_eff_finalSelection.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_FinalSelection"));
    h_lumi_finalSelection.reset(new LuminosityHists(ctx, "Lumi_FinalSelection"));
    h_ht_finalSelection.reset(new HT2dHists(ctx, "HT2d_FinalSelection"));
    h_PDF_variations.reset(new LQToTopMuPDFHists(ctx, "PDF_variations", true, do_pdf_variations));
    h_Sideband.reset(new LQToTopMuHists(ctx, "Sideband_weights_applied"));
    
    
  }
  
  
  bool LQToTopMuAnalysisModule::process(Event & event) {
    //cout << endl << endl << "+++NEW EVENT+++" << endl;

    if(is_mc){
      double factor_xsec = -1;
      int control = (dataset_version.Contains("TTbar") && Sys_TTbar != "nominal") + (dataset_version.Contains("DYJets") && Sys_DY != "nominal") + (dataset_version.Contains("SingleTop") && Sys_ST != "nominal") + (dataset_version.Contains("WJets") && Sys_WJ != "nominal") + (dataset_version.Contains("Diboson") && Sys_DB != "nominal") + (dataset_version.Contains("QCD") && Sys_QCD != "nominal")+ (dataset_version.Contains("TTV") && Sys_TTV != "nominal");
      if(!(control == 0 || control == 1)) throw runtime_error("In LQToTopMuSidebandAnalysisModule.cxx: More than one rate systematic is set to something different than 'nominal'");

      if(control == 0) factor_xsec = 0;
      else if(control == 1){
	if(dataset_version.Contains("TTbar") && Sys_TTbar != "nominal")       factor_xsec = 0.056;
	else if(dataset_version.Contains("DYJets") && Sys_DY != "nominal")    factor_xsec = 0.1;
	else if(dataset_version.Contains("SingleTop") && Sys_ST != "nominal") factor_xsec = 0.1;
	else if(dataset_version.Contains("WJets") && Sys_WJ != "nominal")     factor_xsec = 0.1;
	else if(dataset_version.Contains("Diboson") && Sys_DB != "nominal")   factor_xsec = 0.2;
	else if(dataset_version.Contains("QCD") && Sys_QCD != "nominal")      factor_xsec = 1;
	else if(dataset_version.Contains("TTV") && Sys_TTV != "nominal")      factor_xsec = 0.25;
	else if(dataset_version.Contains("LQ"))                               factor_xsec = 0;
      }
      double sf_xsec = 1;
      if(Sys_TTbar == "up" || Sys_DY == "up"|| Sys_ST == "up"|| Sys_DB == "up"|| Sys_WJ == "up"|| Sys_QCD == "up"|| Sys_TTV == "up") sf_xsec += factor_xsec;
      else if(Sys_TTbar == "down" || Sys_DY == "down"|| Sys_ST == "down"|| Sys_DB == "down"|| Sys_WJ == "down"|| Sys_QCD == "down"|| Sys_TTV == "down") sf_xsec -= factor_xsec;
      else if(control != 0) throw runtime_error("In LQToTopMuAnalysisModule.cxx: Invalid direction for 'Sys_Rate_YYY' specified.");

      event.weight *= sf_xsec;
    }


    //Jet_printer->process(event);
    //Muon_printer->process(event);
    //Electron_printer->process(event);
    //GenParticles_printer->process(event);

    //Apply fake-rate scale factors
    double SF_ele = 1, SF_mu = 1;
    double n_matched_to_taus = 0, n_matched_to_muons = 0;
    if(is_mc){
      //cout <<endl << "NEW EVENT" << endl << "Before applying SF: weight = " << event.weight << endl;
      SF_ele = event.get(h_FakeRateWeightEle);
      if(event.electrons->size() == 0 && SF_ele != 1) throw runtime_error("There are no electrons in the event, still the fake-rate SF is != 1...");
      event.weight *= SF_ele;    
      SF_mu = event.get(h_FakeRateWeightMu);
      if(event.muons->size() == 0 && SF_mu != 1) throw runtime_error("There are no muons in the event, still the fake-rate SF is != 1...");     
      event.weight *= SF_mu;
    }

    if(!is_foreff){
      //apply muon SFs as in preselection
      //SF_muonTrigger->process(event);
      SF_muonTrigger->process_onemuon(event,0);
      SF_muonID->process(event);
      SF_muonIso->process(event);
      SF_muonTrk->process(event);

      if(event.electrons->size() >= 1){
	SF_eleReco->process(event);
	SF_eleID->process(event);
      }
    }

    //Only for NoLep selection
    if(is_foreff){
      if(!(ele_trigger_sel1->passes(event) || ele_trigger_sel2->passes(event)) ) return false;
      //if(!ht_sel->passes(event)) return false;
    }
     
    if(do_scale_variation){
      syst_module->process(event);    
    }



    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);
    if(!njet_sel->passes(event)) return false;
    if(!ht_sel->passes(event)) return false;




    /*ttgenprod->process(event);
      LQgenprod->process(event);*/
    // MLQ reco
    if(nele_sel->passes(event) && event.muons->size() >= 2){
      for(auto & m : recomodules){
	m->process(event);
      }
    }


    

    double sum_mu_charge = 0;
    for(const auto & mu : *event.muons) sum_mu_charge += mu.charge();

    
    if(event.muons->size() == 3 && event.electrons->size() == 0 && fabs(sum_mu_charge) == 1){
      for(auto & m : muonic_recomodules){
	m->process(event);
      }
    }

  
    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_eff_nocuts->fill(event);
    h_lumi_nocuts->fill(event);



    //Nele
    // if(nele_sel->passes(event)) {
    //   h_1ele->fill(event);
    //   h_jets_1ele->fill(event);
    //   h_ele_1ele->fill(event);
    //   h_mu_1ele->fill(event);
    //   h_event_1ele->fill(event);
    //   h_topjets_1ele->fill(event);
    //   h_eff_1ele->fill(event);
    //   h_lumi_1ele->fill(event);
    // }

    //1 bTag1 loose
    if(!nbtag_loose_sel->passes(event)) return false; 
    SF_btag->process(event);
    //LOESCHEMICH
    h_jets_1bJetLoose->fill(event);
    h_1bJetLoose->fill(event);
    h_ele_1bJetLoose->fill(event);
    h_mu_1bJetLoose->fill(event);
    h_event_1bJetLoose->fill(event);
    h_topjets_1bJetLoose->fill(event);
    h_eff_1bJetLoose->fill(event);
    h_lumi_1bJetLoose->fill(event);

    //InvMassVeto
    if(!(m_mumu_veto->passes(event) || is_foreff)) return false; 
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    h_event_InvMassVeto->fill(event);
    h_topjets_InvMassVeto->fill(event);
    h_eff_InvMassVeto->fill(event);
    h_lumi_InvMassVeto->fill(event);

    //HT Lept
    if(!htlept_sel->passes(event)) return false;
    h_jets_htlept200->fill(event);
    h_htlept200->fill(event);
    h_ele_htlept200->fill(event);
    h_mu_htlept200->fill(event);
    h_event_htlept200->fill(event);
    h_topjets_htlept200->fill(event);    
    h_eff_htlept200->fill(event);
    h_lumi_htlept200->fill(event);
    h_btageff_htlept200->fill(event);

    //Final Selection
    h_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);
    h_eff_finalSelection->fill(event);
    h_lumi_finalSelection->fill(event);
    h_ht_finalSelection->fill(event);

    h_PDF_variations->fill(event);
    h_Sideband->fill(event);


    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuAnalysisModule)
} 


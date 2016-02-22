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

  class LQToTopMuSidebandAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuSidebandAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    std::unique_ptr<CommonModules> common;
    std::unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer, GenParticles_printer, syst_module, SF_muonID, SF_muonTrigger, SF_muonIso, SF_btag;
    
    std::unique_ptr<JetCleaner> jetcleaner;
    
    // declare the Selections to use.
    std::unique_ptr<Selection>  nele_sel, njet_sel, nbtag_loose_sel, nbtag_med_sel, nbtag_tight_sel, htlept_sel, mttbar_gen_sel, m_muele_veto;
    
    // store the Hists collection as member variables. 
    std::unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts;
    std::unique_ptr<Hists> h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_event_1ele, h_topjets_1ele, h_tau_1ele;
    std::unique_ptr<Hists> h_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_event_1bJetLoose, h_topjets_1bJetLoose, h_tau_1bJetLoose;
    std::unique_ptr<Hists> h_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto, h_tau_InvMassVeto;
    std::unique_ptr<Hists> h_htlept200, h_jets_htlept200, h_ele_htlept200, h_mu_htlept200, h_event_htlept200, h_topjets_htlept200, h_tau_htlept200, h_btageff_htlept200;
    std::unique_ptr<Hists> h_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection, h_tau_finalSelection;
    std::unique_ptr<Hists> h_ht_InvMassVeto, h_ht_finalSelection;
    std::unique_ptr<Hists> h_Sideband, h_Sideband_inclusive;
    std::unique_ptr<Hists> h_PDF_variations;

    
    MuonId MuId;
    ElectronId EleId;
    JetId Btag_loose, Btag_medium, Btag_tight;
    CSVBTag::wp wp_btag_loose;

    bool do_scale_variation, is_mc, do_pdf_variations;


    TFile* file_alpha_inclusive;
    TGraphAsymmErrors* alpha_inclusive;
    TH1D* norm_inclusive;
    string filepath_alpha, Sys_MuonID, Sys_BTag, Sys_MuonTrigger, Sys_MuonIso, Sys_PU;

    std::vector<std::unique_ptr<AnalysisModule>> recomodules;
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_hyps;
    
  };
  
  
  LQToTopMuSidebandAnalysisModule::LQToTopMuSidebandAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuSidebandAnalysisModule!" << endl;

    do_scale_variation = false;
    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    is_mc = ctx.get("dataset_type") == "MC";
    filepath_alpha = ctx.get("filepath_alpha");
    Sys_MuonID = ctx.get("Systematic_MuonID");
    Sys_MuonTrigger = ctx.get("Systematic_MuonTrigger");
    Sys_MuonIso = ctx.get("Systematic_MuonIso");
    Sys_BTag = ctx.get("Systematic_BTag");
    Sys_PU = ctx.get("Systematic_PU");
    const char* c_filepath_alpha = filepath_alpha.c_str();

    file_alpha_inclusive = new TFile(c_filepath_alpha,"READ");
    alpha_inclusive = (TGraphAsymmErrors*)file_alpha_inclusive->Get("Graph");
    norm_inclusive = (TH1D*)file_alpha_inclusive->Get("h_normalization");

    // 1. setup other modules. CommonModules and the JetCleaner:
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
    Muon_printer.reset(new MuonPrinter("Muon-Printer"));
    GenParticles_printer.reset(new GenParticlesPrinter(ctx));
    MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4),MuonIso(0.12));
    EleId = AndId<Electron>(ElectronID_Spring15_25ns_medium,PtEtaCut(30.0, 2.4));

    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
    wp_btag_loose = CSVBTag::WP_LOOSE;
    Btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
    Btag_tight = CSVBTag(CSVBTag::WP_TIGHT); 
    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.5)); 

    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->init(ctx,Sys_PU);
    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_7_4_15_patch1/src/UHH2/common/data/MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1", 1, "tightID", Sys_MuonID));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_7_4_15_patch1/src/UHH2/common/data/SingleMuonTrigger_Z_RunD_Reco74X_Nov20.root", "IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins", 0.5, "trigger", Sys_MuonTrigger));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_7_4_15_patch1/src/UHH2/common/data/MuonIso_Z_RunCD_Reco74X_Dec1.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1", 1, "iso", Sys_MuonIso));
    SF_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_loose,"jets",Sys_BTag));

    
    h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");

    
    // 2. set up selections
    //Selection
    nbtag_loose_sel.reset(new NJetSelection(1, -1, Btag_loose));  
    nbtag_med_sel.reset(new NJetSelection(1, -1, Btag_medium));
    nbtag_tight_sel.reset(new NJetSelection(2, -1, Btag_tight));
    nele_sel.reset(new NElectronSelection(1, 1));
    htlept_sel.reset(new HTLeptSelection(200., -1));
    m_muele_veto.reset(new InvMassMuEleVeto(71.,111.));
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));

    //systematics modules
    syst_module.reset(new MCScaleVariation(ctx));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
 
    h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
    h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
    h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
    h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));  
    h_event_1ele.reset(new EventHists(ctx, "Event_1Ele"));
    h_topjets_1ele.reset(new TopJetHists(ctx, "TopJets_1Ele"));

    h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
    h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
    h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
    h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));
    h_event_1bJetLoose.reset(new EventHists(ctx, "Event_1bJetLoose"));
    h_topjets_1bJetLoose.reset(new TopJetHists(ctx, "TopJets_1bJetLoose"));

    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
    h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));


    h_htlept200.reset(new LQToTopMuHists(ctx, "HTLept200"));
    h_jets_htlept200.reset(new JetHists(ctx, "Jets_HTLept200"));
    h_ele_htlept200.reset(new ElectronHists(ctx, "Ele_HTLept200"));
    h_mu_htlept200.reset(new MuonHists(ctx, "Mu_HTLept200"));
    h_topjets_htlept200.reset(new TopJetHists(ctx, "TopJets_HTLept200"));
    h_event_htlept200.reset(new EventHists(ctx, "Event_HTLept200"));
    h_btageff_htlept200.reset(new BTagMCEfficiencyHists(ctx, "BTagEff_HTLept200",wp_btag_loose));


    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection")); 
    h_ht_finalSelection.reset(new HT2dHists(ctx, "HT2d_FinalSelection"));
    h_Sideband.reset(new LQToTopMuHists(ctx, "Sideband_weights_applied"));
    h_Sideband_inclusive.reset(new LQToTopMuHists(ctx, "Sideband_inclusive_weights_applied"));
    h_PDF_variations.reset(new LQToTopMuPDFHists(ctx, "PDF_variations", do_pdf_variations));
    
  }
  
  
  bool LQToTopMuSidebandAnalysisModule::process(Event & event) {
    
    
    //apply muon SFs as in preselection
    if(is_mc){
      SF_muonTrigger->process(event);
      SF_muonID->process(event);
      //SF_muonIso->process(event);
    }

    if(do_scale_variation){
      syst_module->process(event);    
    }
    
    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);


    //HT
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
    for(const auto & muon : *event.muons){
      ht_lep += muon.pt();
    }
    ht = ht_lep + ht_jets + met;


    //if(is_mc)GenParticles_printer->process(event);

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
 
    //Nele
    if(!nele_sel->passes(event)) return false;
    h_1ele->fill(event);
    h_jets_1ele->fill(event);
    h_ele_1ele->fill(event);
    h_mu_1ele->fill(event);
    h_event_1ele->fill(event);
    h_topjets_1ele->fill(event);


    //1 bTag1 loose
    if(!nbtag_loose_sel->passes(event)) return false;
    SF_btag->process(event);

    h_jets_1bJetLoose->fill(event);
    h_1bJetLoose->fill(event);
    h_ele_1bJetLoose->fill(event);
    h_mu_1bJetLoose->fill(event);
    h_event_1bJetLoose->fill(event);
    h_topjets_1bJetLoose->fill(event);
 

    //InvMassVeto
    if(!m_muele_veto->passes(event)) return false;
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    h_event_InvMassVeto->fill(event);
    h_topjets_InvMassVeto->fill(event);


    if(!htlept_sel->passes(event)) return false;
    h_jets_htlept200->fill(event);
    h_htlept200->fill(event);
    h_ele_htlept200->fill(event);
    h_mu_htlept200->fill(event);
    h_event_htlept200->fill(event);
    h_topjets_htlept200->fill(event);
    h_btageff_htlept200->fill(event);


    h_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);
    h_ht_finalSelection->fill(event);

    h_PDF_variations->fill(event);

    
    double original_weight = event.weight;
    double d_alpha_inclusive = alpha_inclusive->Eval(ht);
    double d_norm_inclusive = norm_inclusive->GetBinContent(1);
    double sideband_weight_inclusive = d_alpha_inclusive * d_norm_inclusive;
    
    //change weights
    event.weight *= sideband_weight_inclusive;
    h_Sideband_inclusive->fill(event);
    //restore weights
    event.weight = original_weight;

    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuSidebandAnalysisModule)
} 


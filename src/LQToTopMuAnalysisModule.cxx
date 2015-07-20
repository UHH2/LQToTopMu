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
#include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
#include "UHH2/LQToTopMu/include/MET2dHists.h"
#include "UHH2/LQToTopMu/include/HT2dHists.h"
//#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
//#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstruction.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQGen.h"


using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    std::unique_ptr<CommonModules> common;
    std::unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer, GenParticles_printer;
    
    std::unique_ptr<JetCleaner> jetcleaner;
    std::unique_ptr<MuonCleaner> muoncleaner;
    std::unique_ptr<MuonCleaner> muoncleaner_iso;
    std::unique_ptr<ElectronCleaner> electroncleaner;
    
    // declare the Selections to use.
    std::unique_ptr<Selection>  nele_sel, njet_sel, nbtag_loose_sel, nbtag_med_sel, nbtag_tight_sel, m_mumu_veto, m_mumu_veto_inverted, ht_sel, pt_lead_mu_sel, pt_2nd_mu_sel, pt_lead_jet_sel, pt_2nd_jet_sel, cmsTopTag_sel, hepTopTag_sel, ptrel_mu1jet_sel, eta_leadjet_sel, htjets_sel, htlept_sel, met_sel, nmuon_sel, sideband_pt_lead_mu_sel, sideband_pt_2nd_mu_sel, sideband_njet_sel, sideband_ptrel_mujet_sel, sideband_pt_lead_jet_sel, sideband_pt_lead_jet_sel2, sideband_pt_lead_jet_sel3, sideband_pt_lead_jet_sel4, sideband_met_sel, ZMuMu_sel;
    
    // store the Hists collection as member variables. 
    std::unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts;
    std::unique_ptr<Hists> h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_event_1ele, h_topjets_1ele, h_tau_1ele;
    std::unique_ptr<Hists> h_3jets, h_jets_3jets, h_ele_3jets, h_mu_3jets, h_event_3jets, h_topjets_3jets, h_tau_3jets;
    std::unique_ptr<Hists> h_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_event_1bJetLoose, h_topjets_1bJetLoose, h_tau_1bJetLoose;
    std::unique_ptr<Hists> h_2bJetsLoose, h_jets_2bJetsLoose, h_mu_2bJetsLoose, h_ele_2bJetsLoose, h_event_2bJetsLoose, h_topjets_2bJetsLoose, h_tau_2bJetsLoose;
    std::unique_ptr<Hists> h_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto, h_tau_InvMassVeto;
    std::unique_ptr<Hists> h_InvMassVetoInverted, h_jets_InvMassVetoInverted, h_mu_InvMassVetoInverted, h_ele_InvMassVetoInverted, h_event_InvMassVetoInverted, h_topjets_InvMassVetoInverted, h_tau_InvMassVetoInverted;
    std::unique_ptr<Hists> h_pt_leadmu300, h_jets_pt_leadmu300, h_ele_pt_leadmu300, h_mu_pt_leadmu300, h_event_pt_leadmu300, h_topjets_pt_leadmu300, h_tau_pt_leadmu300;
    std::unique_ptr<Hists> h_pt_2ndmu90, h_jets_pt_2ndmu90, h_ele_pt_2ndmu90, h_mu_pt_2ndmu90, h_event_pt_2ndmu90, h_topjets_pt_2ndmu90, h_taupt_2ndmu90;
    std::unique_ptr<Hists> h_pt_leadjet200, h_jets_pt_leadjet200, h_ele_pt_leadjet200, h_mu_pt_leadjet200, h_event_pt_leadjet200, h_topjets_pt_leadjet200, h_tau_pt_leadjet200;
    std::unique_ptr<Hists> h_pt_2ndjet120, h_jets_pt_2ndjet120, h_ele_pt_2ndjet120, h_mu_pt_2ndjet120, h_event_pt_2ndjet120, h_topjets_pt_2ndjet120, h_tau_pt_2ndjet120;
    std::unique_ptr<Hists> h_1cmstoptag, h_jets_1cmstoptag, h_ele_1cmstoptag, h_mu_1cmstoptag, h_event_1cmstoptag, h_topjets_1cmstoptag, h_tau_1cmstoptag;
    std::unique_ptr<Hists> h_1heptoptag, h_jets_1heptoptag, h_ele_1heptoptag, h_mu_1heptoptag, h_event_1heptoptag, h_topjets_1heptoptag, h_tau_1heptoptag;
    std::unique_ptr<Hists> h_ptrel_mu1jet200, h_jets_ptrel_mu1jet200, h_ele_ptrel_mu1jet200, h_mu_ptrel_mu1jet200, h_event_ptrel_mu1jet200, h_topjets_ptrel_mu1jet200, h_tau_ptrel_mu1jet200;
    std::unique_ptr<Hists> h_eta_leadjet18, h_jets_eta_leadjet18, h_ele_eta_leadjet18, h_mu_eta_leadjet18, h_event_eta_leadjet18, h_topjets_eta_leadjet18, h_tau_eta_leadjet18;
    std::unique_ptr<Hists> h_htjets250, h_jets_htjets250, h_ele_htjets250, h_mu_htjets250, h_event_htjets250, h_topjets_htjets250, h_tau_htjets250;
    std::unique_ptr<Hists> h_htlept200, h_jets_htlept200, h_ele_htlept200, h_mu_htlept200, h_event_htlept200, h_topjets_htlept200, h_tau_htlept200;
    std::unique_ptr<Hists> h_met40, h_jets_met40, h_ele_met40, h_mu_met40, h_event_met40, h_topjets_met40, h_tau_met40;
    std::unique_ptr<Hists> h_3muons, h_jets_3muons, h_ele_3muons, h_mu_3muons, h_event_3muons, h_topjets_3muons, h_tau_3muons;
    std::unique_ptr<Hists> h_ht450, h_jets_ht450, h_ele_ht450, h_mu_ht450, h_event_ht450, h_topjets_ht450;

    std::unique_ptr<Hists> h_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection, h_tau_finalSelection;

    std::unique_ptr<Hists> h_sideband_zmumu, h_jets_sideband_zmumu, h_ele_sideband_zmumu, h_mu_sideband_zmumu, h_event_sideband_zmumu, h_topjets_sideband_zmumu;
    std::unique_ptr<Hists> h_sideband_ptleadmu, h_jets_sideband_ptleadmu, h_ele_sideband_ptleadmu, h_mu_sideband_ptleadmu, h_event_sideband_ptleadmu, h_topjets_sideband_ptleadmu;
    std::unique_ptr<Hists> h_sideband_pt2ndmu, h_jets_sideband_pt2ndmu, h_ele_sideband_pt2ndmu, h_mu_sideband_pt2ndmu, h_event_sideband_pt2ndmu, h_topjets_sideband_pt2ndmu;
    std::unique_ptr<Hists> h_ht_InvMassVeto, h_ht_finalSelection;

    Event::Handle<TTbarGen> h_ttbargen;
    Event::Handle<LQGen> h_LQLQbargen;
    std::unique_ptr<AnalysisModule> ttgenprod;
    std::unique_ptr<AnalysisModule> LQgenprod;

    
    MuonId MuIso;
    TopJetId cmsTopTag;
    TopJetId hepTopTag;
    JetId Btag_loose, Btag_medium, Btag_tight;
    
    std::vector<std::unique_ptr<AnalysisModule>> recomodules;
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_hyps;
    
  };
  
  
  LQToTopMuAnalysisModule::LQToTopMuAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuAnalysisModule!" << endl;
    
    // 1. setup other modules. CommonModules and the JetCleaner:
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
    Muon_printer.reset(new MuonPrinter("Muon-Printer"));
    GenParticles_printer.reset(new GenParticlesPrinter(ctx));
    common.reset(new CommonModules());
    common->disable_jersmear();
    common->disable_jec();
    common->init(ctx);
    
    MuIso = MuonIso(0.12);
    cmsTopTag = CMSTopTag(50,140,250);
    hepTopTag = HEPTopTag(200, 0.85, 1.15, 0.35, 0.35);
    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
    Btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
    Btag_tight = CSVBTag(CSVBTag::WP_TIGHT); 

    jetcleaner.reset(new JetCleaner(30.0, 2.5)); 
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.1))));
    muoncleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuIso, PtEtaCut(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(30.0, 2.5))));
    
    h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
    /*ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    LQgenprod.reset(new LQGenProducer(ctx, "LQLQbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    h_LQLQbargen = ctx.get_handle<LQGen>("LQLQbargen");*/
    
    // 2. set up selections
    //Selection
    njet_sel.reset(new NJetSelection(3, -1));
    nbtag_loose_sel.reset(new NJetSelection(1, -1, Btag_loose));  //default: (1, -1)
    nbtag_med_sel.reset(new NJetSelection(1, -1, Btag_medium));
    nbtag_tight_sel.reset(new NJetSelection(1, -1, Btag_tight));
    m_mumu_veto.reset(new InvMass2MuVeto(71, 111)); //81, 101
    m_mumu_veto_inverted.reset(new InvMass2MuVetoInverted(71, 111)); //81, 101
    nele_sel.reset(new NElectronSelection(1, -1));
    ht_sel.reset(new HtSelection(450,-1));
    pt_lead_mu_sel.reset(new PtLeadingMuonSelection(300., -1));
    pt_2nd_mu_sel.reset(new Pt2ndMuonSelection(90,-1));
    pt_lead_jet_sel.reset(new PtLeadingJetSelection(200,-1));
    pt_2nd_jet_sel.reset(new Pt2ndJetSelection(120,-1));
    cmsTopTag_sel.reset(new NTopJetSelection(1, -1, cmsTopTag));
    hepTopTag_sel.reset(new NTopJetSelection(1, -1, hepTopTag)); // use right topjet collection!
    ptrel_mu1jet_sel.reset(new PtRelMu1JetSelection(200, -1));
    eta_leadjet_sel.reset(new EtaLeadingJetSelection(1.8));
    htjets_sel.reset(new HTJetsSelection(250., -1));
    htlept_sel.reset(new HTLeptSelection(200., -1));
    met_sel.reset(new METSelection(40, -1));
    nmuon_sel.reset(new NMuonSelection(3, -1));



    //DY Sideband Selection
    sideband_pt_lead_mu_sel.reset(new PtLeadingMuonSelection(200, -1));
    sideband_pt_2nd_mu_sel.reset(new Pt2ndMuonSelection(200, -1)); 
    sideband_njet_sel.reset(new NJetSelection(0, -1));
    sideband_ptrel_mujet_sel.reset(new PtRelMuJetSelection(60, -1));
    sideband_pt_lead_jet_sel.reset(new PtLeadingJetSelection(100, -1));
    sideband_pt_lead_jet_sel2.reset(new PtLeadingJetSelection(250, 300));
    sideband_pt_lead_jet_sel3.reset(new PtLeadingJetSelection(350, 400));
    sideband_pt_lead_jet_sel4.reset(new PtLeadingJetSelection(450, 500));
    sideband_met_sel.reset(new METSelection(100, -1));
    ZMuMu_sel.reset(new GenLvlZMuMuSelection());
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));
    //recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassLQReconstruction"));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
    h_tau_nocuts.reset(new TauHists(ctx, "Tau_NoCuts"));

    /*h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
    h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
    h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
    h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));  
    h_event_1ele.reset(new EventHists(ctx, "Event_1Ele"));
    h_topjets_1ele.reset(new TopJetHists(ctx, "TopJets_1Ele"));
    h_tau_1ele.reset(new TauHists(ctx, "Tau_1Ele"));
    h_hyphists.reset(new HypothesisHistsOwn(ctx, "Chi2_Hists", "HighMassLQReconstruction", "Chi2"));*/
    //h_hyphists.reset(new HypothesisHistsOwn(ctx, "CorrectMatch_Hists", "HighMassLQReconstruction", "CorrectMatch"));

    h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
    h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
    h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
    h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));
    h_event_1bJetLoose.reset(new EventHists(ctx, "Event_1bJetLoose"));
    h_topjets_1bJetLoose.reset(new TopJetHists(ctx, "TopJets_1bJetLoose"));
    //h_tau_1bJetLoose.reset(new TauHists(ctx, "Tau_1bJetLoose"));

    /*h_2bJetsLoose.reset(new LQToTopMuHists(ctx, "2bJetsLoose"));
    h_jets_2bJetsLoose.reset(new JetHists(ctx, "Jets_2bJetsLoose"));
    h_ele_2bJetsLoose.reset(new ElectronHists(ctx, "Ele_2bJetsLoose"));
    h_mu_2bJetsLoose.reset(new MuonHists(ctx, "Mu_2bJetsLoose"));
    h_event_2bJetsLoose.reset(new EventHists(ctx, "Event_2bJetsLoose"));
    h_topjets_2bJetsLoose.reset(new TopJetHists(ctx, "TopJets_2bJetsLoose"));
    h_tau_nocuts.reset(new TauHists(ctx, "Tau_NoCuts"));*/

    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
    h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));
    //h_tau_InvMassVeto.reset(new TauHists(ctx, "Tau_InvMassVeto"));

    /*h_InvMassVetoInverted.reset(new LQToTopMuHists(ctx, "InvMassVetoInverted"));
    h_jets_InvMassVetoInverted.reset(new JetHists(ctx, "Jets_InvMassVetoInverted"));
    h_ele_InvMassVetoInverted.reset(new ElectronHists(ctx, "Ele_InvMassVetoInverted"));
    h_mu_InvMassVetoInverted.reset(new MuonHists(ctx, "Mu_InvMassVetoInverted"));
    h_event_InvMassVetoInverted.reset(new EventHists(ctx, "Event_InvMassVetoInverted"));
    h_topjets_InvMassVetoInverted.reset(new TopJetHists(ctx, "TopJets_InvMassVetoInverted"));
    h_tau_InvMassVetoInverted.reset(new TauHists(ctx, "Tau_InvMassVetoInverted"));*/

    //h_ht_InvMassVeto.reset(new HT2dHists(ctx, "HT2d_InvMassVeto"));


    /*h_3jets.reset(new LQToTopMuHists(ctx, "3Jets"));
    h_jets_3jets.reset(new JetHists(ctx, "Jets_3Jets"));
    h_ele_3jets.reset(new ElectronHists(ctx, "Ele_3Jets"));
    h_mu_3jets.reset(new MuonHists(ctx, "Mu_3Jets"));
    h_event_3jets.reset(new EventHists(ctx, "Event_3Jets"));
    h_topjets_3jets.reset(new TopJetHists(ctx, "TopJets_3Jets"));*/
    //h_tau_3jets.reset(new TauHists(ctx, "Tau_3Jets"));

    /*h_pt_leadmu300.reset(new LQToTopMuHists(ctx, "Pt_LeadMu300"));
    h_jets_pt_leadmu300.reset(new JetHists(ctx, "Jets_Pt_LeadMu300"));
    h_ele_pt_leadmu300.reset(new ElectronHists(ctx, "Ele_Pt_LeadMu300"));
    h_mu_pt_leadmu300.reset(new MuonHists(ctx, "Mu_Pt_LeadMu300"));
    h_topjets_pt_leadmu300.reset(new TopJetHists(ctx, "TopJets_Pt_LeadMu300"));
    h_event_pt_leadmu300.reset(new EventHists(ctx, "Event_Pt_LeadMu300"));*/
    //h_tau_pt_leadmu300.reset(new TauHists(ctx, "Tau_Pt_LeadMu300"));

    /*h_pt_2ndmu90.reset(new LQToTopMuHists(ctx, "Pt_2ndMu90"));
    h_jets_pt_2ndmu90.reset(new JetHists(ctx, "Jets_Pt_2ndMu90"));
    h_ele_pt_2ndmu90.reset(new ElectronHists(ctx, "Ele_Pt_2ndMu90"));
    h_mu_pt_2ndmu90.reset(new MuonHists(ctx, "Mu_Pt_2ndMu90"));
    h_topjets_pt_2ndmu90.reset(new TopJetHists(ctx, "TopJets_Pt_2ndMu90"));
    h_event_pt_2ndmu90.reset(new EventHists(ctx, "Event_Pt_2ndMu90"));
    h_tau_pt_2ndmu90.reset(new TauHists(ctx, "Tau_Pt_2ndMu90"));*/

    /*h_ptrel_mu1jet200.reset(new LQToTopMuHists(ctx, "Ptrel_Mu1Jet200"));
    h_jets_ptrel_mu1jet200.reset(new JetHists(ctx, "Jets_Ptrel_Mu1Jet200"));
    h_ele_ptrel_mu1jet200.reset(new ElectronHists(ctx, "Ele_Ptrel_Mu1Jet200"));
    h_mu_ptrel_mu1jet200.reset(new MuonHists(ctx, "Mu_Ptrel_Mu1Jet200"));
    h_topjets_ptrel_mu1jet200.reset(new TopJetHists(ctx, "TopJets_Ptrel_Mu1Jet200"));
    h_event_ptrel_mu1jet200.reset(new EventHists(ctx, "Event_Ptrel_Mu1Jet200"));
    h_tau_ptrel_mu1jet200.reset(new TauHists(ctx, "Tau_Ptrel_Mu1Jet200"));*/

    /*h_1cmstoptag.reset(new LQToTopMuHists(ctx, "1CMSTopTag"));
    h_jets_1cmstoptag.reset(new JetHists(ctx, "Jets_1CMSTopTag"));
    h_ele_1cmstoptag.reset(new ElectronHists(ctx, "Ele_1CMSTopTag"));
    h_mu_1cmstoptag.reset(new MuonHists(ctx, "Mu_1CMSTopTag"));
    h_topjets_1cmstoptag.reset(new TopJetHists(ctx, "TopJets_1CMSTopTag"));
    h_event_1cmstoptag.reset(new EventHists(ctx, "Event_1CMSTopTag"));
    h_tau_1cmstoptag.reset(new TauHists(ctx, "Tau_1CMSTopTag"));*/

    /*h_1heptoptag.reset(new LQToTopMuHists(ctx, "1HEPTopTag"));
    h_jets_1heptoptag.reset(new JetHists(ctx, "Jets_1HEPTopTag"));
    h_ele_1heptoptag.reset(new ElectronHists(ctx, "Ele_1HEPTopTag"));
    h_mu_1heptoptag.reset(new MuonHists(ctx, "Mu_1HEPTopTag"));
    h_topjets_1heptoptag.reset(new TopJetHists(ctx, "TopJets_1HEPTopTag"));
    h_event_1heptoptag.reset(new EventHists(ctx, "Event_1HEPTopTag"));
    h_tau_1heptoptag.reset(new TauHists(ctx, "Tau_1HEPTopTag"));*/

    /*h_eta_leadjet18.reset(new LQToTopMuHists(ctx, "Eta_LeadJet18"));
    h_jets_eta_leadjet18.reset(new JetHists(ctx, "Jets_Eta_LeadJet18"));
    h_ele_eta_leadjet18.reset(new ElectronHists(ctx, "Ele_Eta_LeadJet18"));
    h_mu_eta_leadjet18.reset(new MuonHists(ctx, "Mu_Eta_LeadJet18"));
    h_topjets_eta_leadjet18.reset(new TopJetHists(ctx, "TopJets_Eta_LeadJet18"));
    h_event_eta_leadjet18.reset(new EventHists(ctx, "Event_Eta_LeadJet18"));
    //h_tau_eta_leadjet18.reset(new TauHists(ctx, "Tau_Eta_LeadJet18"));*/

    /*h_htjets250.reset(new LQToTopMuHists(ctx, "HTJets250"));
    h_jets_htjets250.reset(new JetHists(ctx, "Jets_HTJets250"));
    h_ele_htjets250.reset(new ElectronHists(ctx, "Ele_HTJets250"));
    h_mu_htjets250.reset(new MuonHists(ctx, "Mu_HTJets250"));
    h_topjets_htjets250.reset(new TopJetHists(ctx, "TopJets_HTJets250"));
    h_event_htjets250.reset(new EventHists(ctx, "Event_HTJets250"));*/
    //h_tau_htjets250.reset(new TauHists(ctx, "Tau_HTJets250"));

    h_htlept200.reset(new LQToTopMuHists(ctx, "HTLept200"));
    h_jets_htlept200.reset(new JetHists(ctx, "Jets_HTLept200"));
    h_ele_htlept200.reset(new ElectronHists(ctx, "Ele_HTLept200"));
    h_mu_htlept200.reset(new MuonHists(ctx, "Mu_HTLept200"));
    h_topjets_htlept200.reset(new TopJetHists(ctx, "TopJets_HTLept200"));
    h_event_htlept200.reset(new EventHists(ctx, "Event_HTLept200"));
    //h_tau_htlept200.reset(new TauHists(ctx, "Tau_HTLept200"));

    /*h_met40.reset(new LQToTopMuHists(ctx, "MET40"));
    h_jets_met40.reset(new JetHists(ctx, "Jets_MET40"));
    h_ele_met40.reset(new ElectronHists(ctx, "Ele_MET40"));
    h_mu_met40.reset(new MuonHists(ctx, "Mu_MET40"));
    h_topjets_met40.reset(new TopJetHists(ctx, "TopJets_MET40"));
    h_event_met40.reset(new EventHists(ctx, "Event_MET40"));*/
    //h_tau_met40.reset(new TauHists(ctx, "Tau_MET40"));

    /*h_3muons.reset(new LQToTopMuHists(ctx, "3Muons"));
    h_jets_3muons.reset(new JetHists(ctx, "Jets_3Muons"));
    h_ele_3muons.reset(new ElectronHists(ctx, "Ele_3Muons"));
    h_mu_3muons.reset(new MuonHists(ctx, "Mu_3Muons"));
    h_event_3muons.reset(new EventHists(ctx, "Event_3Muons"));
    h_topjets_3muons.reset(new TopJetHists(ctx, "TopMuons_3Muons"));*/
    //h_tau_3muons.reset(new TauHists(ctx, "Tau_3Muons"));

    /*h_ht450.reset(new LQToTopMuHists(ctx, "HT450"));
    h_jets_ht450.reset(new JetHists(ctx, "Jets_HT450"));
    h_ele_ht450.reset(new ElectronHists(ctx, "Ele_HT450"));
    h_mu_ht450.reset(new MuonHists(ctx, "Mu_HT450"));
    h_topjets_ht450.reset(new TopJetHists(ctx, "TopJets_HT450"));
    h_event_ht450.reset(new EventHists(ctx, "Event_HT450"));*/




    /*h_sideband_zmumu.reset(new LQToTopMuHists(ctx, "Sideband_ZMUMU"));
    h_jets_sideband_zmumu.reset(new JetHists(ctx, "Jets_Sideband_ZMUMU"));
    h_ele_sideband_zmumu.reset(new ElectronHists(ctx, "Ele_Sideband_ZMUMU"));
    h_mu_sideband_zmumu.reset(new MuonHists(ctx, "Mu_Sideband_ZMUMU"));
    h_topjets_sideband_zmumu.reset(new TopJetHists(ctx, "TopJets_Sideband_ZMUMU"));
    h_event_sideband_zmumu.reset(new EventHists(ctx, "Event_Sideband_ZMUMU"));*/

    /*h_sideband_ptleadmu.reset(new LQToTopMuHists(ctx, "Sideband_PtLeadMu"));
    h_jets_sideband_ptleadmu.reset(new JetHists(ctx, "Jets_Sideband_PtLeadMu"));
    h_ele_sideband_ptleadmu.reset(new ElectronHists(ctx, "Ele_Sideband_PtLeadMu"));
    h_mu_sideband_ptleadmu.reset(new MuonHists(ctx, "Mu_Sideband_PtLeadMu"));
    h_topjets_sideband_ptleadmu.reset(new TopJetHists(ctx, "TopJets_Sideband_PtLeadMu"));
    h_event_sideband_ptleadmu.reset(new EventHists(ctx, "Event_Sideband_PtLeadMu"));*/

    /*h_sideband_pt2ndmu.reset(new LQToTopMuHists(ctx, "Sideband_Pt2ndMu"));
    h_jets_sideband_pt2ndmu.reset(new JetHists(ctx, "Jets_Sideband_Pt2ndMu"));
    h_ele_sideband_pt2ndmu.reset(new ElectronHists(ctx, "Ele_Sideband_Pt2ndMu"));
    h_mu_sideband_pt2ndmu.reset(new MuonHists(ctx, "Mu_Sideband_Pt2ndMu"));
    h_topjets_sideband_pt2ndmu.reset(new TopJetHists(ctx, "TopJets_Sideband_Pt2ndMu"));
    h_event_sideband_pt2ndmu.reset(new EventHists(ctx, "Event_Sideband_Pt2ndMu"));*/



    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection")); 
    h_ht_finalSelection.reset(new HT2dHists(ctx, "HT2d_FinalSelection"));
    //h_tau_finalSelection.reset(new TauHists(ctx, "Tau_FinalSelection"));

    
  }
  
  
  bool LQToTopMuAnalysisModule::process(Event & event) {

    //for sideband selection, reweight certain events
    auto ptmu1 = event.muons->at(0).pt();
    // good after InvMassVeto
    double corrfactor_ptmu1reweight = (0.030397+0.00059467*ptmu1-0.0000007976*ptmu1*ptmu1+0.0000000003271*ptmu1*ptmu1*ptmu1)*(1.77298+0.00061177*ptmu1-0.000001238*ptmu1*ptmu1);

    //good after HTLEP
    //double corrfactor_ptmu1reweight = 1.1337*(0.10969+0.00006686*ptmu1+0.000000175*ptmu1*ptmu1-0.0000000001767*ptmu1*ptmu1*ptmu1);
    //event.weight *= corrfactor_ptmu1reweight;

    //cout << "eventweight: " << event.weight << endl;

    common->process(event);
    
    muoncleaner->process(event);
    muoncleaner_iso->process(event);
    electroncleaner->process(event);
    jetcleaner->process(event);


    /*Electron_printer->process(event);
    Muon_printer->process(event);
    Jet_printer->process(event);*/
    //GenParticles_printer->process(event);

    /*ttgenprod->process(event);
      LQgenprod->process(event);*/

    // MLQ reco
    for(auto & m : recomodules){
      m->process(event);
      }
    
    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    //h_tau_nocuts->fill(event);

    //Nele
    /*if (!nele_sel->passes(event)) return false;
    h_1ele->fill(event);
    h_jets_1ele->fill(event);
    h_ele_1ele->fill(event);
    h_mu_1ele->fill(event);
    h_hyphists->fill(event);
    h_event_1ele->fill(event);
    h_topjets_1ele->fill(event);*/
    //h_tau_1ele->fill(event);

    //1 bTag1 loose
    if(!nbtag_loose_sel->passes(event)) return false;
    h_jets_1bJetLoose->fill(event);
    h_1bJetLoose->fill(event);
    h_ele_1bJetLoose->fill(event);
    h_mu_1bJetLoose->fill(event);
    h_event_1bJetLoose->fill(event);
    h_topjets_1bJetLoose->fill(event);
    //h_tau_1bJetLoose->fill(event);

    //InvMassVeto
    if(m_mumu_veto->passes(event)) return false; // mit !: Signal, ohne !: Sideband
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    h_event_InvMassVeto->fill(event);
    h_topjets_InvMassVeto->fill(event);
    //h_tau_InvMassVeto->fill(event);

    //InvMassVetoInverted
    /*if(!m_mumu_veto_inverted->passes(event)) return false;
    h_jets_InvMassVetoInverted->fill(event);
    h_InvMassVetoInverted->fill(event);
    h_ele_InvMassVetoInverted->fill(event);
    h_mu_InvMassVetoInverted->fill(event);
    h_event_InvMassVetoInverted->fill(event);
    h_topjets_InvMassVetoInverted->fill(event);
    //h_tau_InvMassVetoInverted->fill(event);*/
    //h_ht_InvMassVeto->fill(event);
  
    //Njets >= 3
    /*if(!njet_sel->passes(event)) return false;
    h_3jets->fill(event);
    h_jets_3jets->fill(event);
    h_ele_3jets->fill(event);
    h_mu_3jets->fill(event);
    h_event_3jets->fill(event);
    h_topjets_3jets->fill(event);*/
    //h_tau_3jets->fill(event);

    //pt leading mu < 300
    /*if(!pt_lead_mu_sel->passes(event)) return false; 
    h_jets_pt_leadmu300->fill(event);
    h_pt_leadmu300->fill(event);
    h_ele_pt_leadmu300->fill(event);
    h_mu_pt_leadmu300->fill(event);
    h_event_pt_leadmu300->fill(event);
    h_topjets_pt_leadmu300->fill(event);*/
    //h_tau_pt_leadmu300->fill(event);

    //pt 2nd mu > 90
    /*if(!pt_2nd_mu_sel->passes(event)) return false;
    h_jets_pt_2ndmu90->fill(event);
    h_pt_2ndmu90->fill(event);
    h_ele_pt_2ndmu90->fill(event);
    h_mu_pt_2ndmu90->fill(event);
    h_event_pt_2ndmu90->fill(event);
    h_topjets_pt_2ndmu90->fill(event);
    //h_tau_pt_2ndmu90->fill(event);*/

    /*if(!ptrel_mu1jet_sel->passes(event)) return false;
    h_jets_ptrel_mu1jet200->fill(event);
    h_ptrel_mu1jet200->fill(event);
    h_ele_ptrel_mu1jet200->fill(event);
    h_mu_ptrel_mu1jet200->fill(event);
    h_event_ptrel_mu1jet200->fill(event);
    h_topjets_ptrel_mu1jet200->fill(event);
    //h_tau_ptrel_mu1jet200->fill(event);*/

    //CMSTopTag
    /*if(!cmsTopTag_sel->passes(event)) return false;
    h_jets_1cmstoptag->fill(event);
    h_1cmstoptag->fill(event);
    h_ele_1cmstoptag->fill(event);
    h_mu_1cmstoptag->fill(event);
    h_event_1cmstoptag->fill(event);
    h_topjets_1cmstoptag->fill(event);
    //h_tau_1cmstoptag->fill(event);*/

    //HEPTopTag
    /*if(!hepTopTag_sel->passes(event)) return false;
    h_jets_1heptoptag->fill(event);
    h_1heptoptag->fill(event);
    h_ele_1heptoptag->fill(event);
    h_mu_1heptoptag->fill(event);
    h_event_1heptoptag->fill(event);
    h_topjets_1heptoptag->fill(event);
    //h_tau_1heptoptag->fill(event);*/

    /*if(!eta_leadjet_sel->passes(event)) return false;
    h_jets_eta_leadjet18->fill(event);
    h_eta_leadjet18->fill(event);
    h_ele_eta_leadjet18->fill(event);
    h_mu_eta_leadjet18->fill(event);
    h_event_eta_leadjet18->fill(event);
    h_topjets_eta_leadjet18->fill(event);
    //h_tau_eta_leadjet18->fill(event);*/

    /*if(!htjets_sel->passes(event)) return false;
    h_jets_htjets250->fill(event);
    h_htjets250->fill(event);
    h_ele_htjets250->fill(event);
    h_mu_htjets250->fill(event);
    h_event_htjets250->fill(event);
    h_topjets_htjets250->fill(event);*/
    //h_tau_htjets250->fill(event);

    if(!htlept_sel->passes(event)) return false;
    h_jets_htlept200->fill(event);
    h_htlept200->fill(event);
    h_ele_htlept200->fill(event);
    h_mu_htlept200->fill(event);
    h_event_htlept200->fill(event);
    h_topjets_htlept200->fill(event);
    //h_tau_htlept200->fill(event);

    /*if(!met_sel->passes(event)) return false;
    h_jets_met40->fill(event);
    h_met40->fill(event);
    h_ele_met40->fill(event);
    h_mu_met40->fill(event);
    h_event_met40->fill(event);
    h_topjets_met40->fill(event);*/
    //h_tau_met40->fill(event);

    //Nmuons >= 3
    /*if(!nmuon_sel->passes(event)) return false;
    h_3muons->fill(event);
    h_jets_3muons->fill(event);
    h_ele_3muons->fill(event);
    h_mu_3muons->fill(event);
    h_event_3muons->fill(event);
    h_topjets_3muons->fill(event);*/
    //h_tau_3muons->fill(event);

    // HT > 450
    /*if(!ht_sel->passes(event)) return false;
    h_ht450->fill(event);
    h_jets_ht450->fill(event);
    h_ele_ht450->fill(event);
    h_mu_ht450->fill(event);
    h_event_ht450->fill(event);
    h_topjets_ht450->fill(event);*/


    //SIDEBAND: Pt Lead Mu 200
    /*if(!sideband_pt_lead_mu_sel->passes(event)) return false;
    h_jets_sideband_ptleadmu->fill(event);
    h_sideband_ptleadmu->fill(event);
    h_ele_sideband_ptleadmu->fill(event);
    h_mu_sideband_ptleadmu->fill(event);
    h_event_sideband_ptleadmu->fill(event);
    h_topjets_sideband_ptleadmu->fill(event);*/

    //SIDEBAND: Pt 2nd Mu 200
    /*if(!sideband_pt_2nd_mu_sel->passes(event)) return false;
    h_jets_sideband_pt2ndmu->fill(event);
    h_sideband_pt2ndmu->fill(event);
    h_ele_sideband_pt2ndmu->fill(event);
    h_mu_sideband_pt2ndmu->fill(event);
    h_event_sideband_pt2ndmu->fill(event);
    h_topjets_sideband_pt2ndmu->fill(event);*/


    // SIDEBAND: Z -> MuMu on Gen Lvl
    /*if(!ZMuMu_sel->passes(event)) return false;
    h_jets_sideband_zmumu->fill(event);
    h_sideband_zmumu->fill(event);
    h_ele_sideband_zmumu->fill(event);
    h_mu_sideband_zmumu->fill(event);
    h_event_sideband_zmumu->fill(event);
    h_topjets_sideband_zmumu->fill(event);*/





    //Final Selection
    h_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);
    //h_tau_finalSelection->fill(event);
    h_ht_finalSelection->fill(event);


    return true;
  }
  
  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the LQToTopMuAnalysisModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuAnalysisModule)
} 


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
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/PrintingModules.h"
//#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    std::unique_ptr<CommonModules> common;
    std::unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer;
    
    std::unique_ptr<JetCleaner> jetcleaner;
    std::unique_ptr<MuonCleaner> muoncleaner;
    std::unique_ptr<MuonCleaner> muoncleaner_iso;
    std::unique_ptr<ElectronCleaner> electroncleaner;
    
    // declare the Selections to use.
    std::unique_ptr<Selection>  nele_sel, njet_sel, nbtag_loose_sel, nbtag_med_sel, nbtag_tight_sel, m_mumu_veto, ht_sel, pt_lead_mu_sel, pt_2nd_mu_sel, pt_lead_jet_sel, pt_2nd_jet_sel, cmsTopTag_sel, hepTopTag_sel;
    
    // store the Hists collection as member variables. 
    std::unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts;
    std::unique_ptr<Hists> h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_event_1ele, h_topjets_1ele;
    std::unique_ptr<Hists> h_4jets, h_jets_4jets, h_ele_4jets, h_mu_4jets, h_event_4jets, h_topjets_4jets;
    std::unique_ptr<Hists> h_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_event_1bJetLoose, h_topjets_1bJetLoose;
    std::unique_ptr<Hists> h_2bJetsLoose, h_jets_2bJetsLoose, h_mu_2bJetsLoose, h_ele_2bJetsLoose, h_event_2bJetsLoose, h_topjets_2bJetsLoose;
    std::unique_ptr<Hists> h_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto;
    std::unique_ptr<Hists> h_pt_leadmu90, h_jets_pt_leadmu90, h_ele_pt_leadmu90, h_mu_pt_leadmu90, h_event_pt_leadmu90, h_topjets_pt_leadmu90;
    std::unique_ptr<Hists> h_pt_2ndmu90, h_jets_pt_2ndmu90, h_ele_pt_2ndmu90, h_mu_pt_2ndmu90, h_event_pt_2ndmu90, h_topjets_pt_2ndmu90;
    std::unique_ptr<Hists> h_pt_leadjet200, h_jets_pt_leadjet200, h_ele_pt_leadjet200, h_mu_pt_leadjet200, h_event_pt_leadjet200, h_topjets_pt_leadjet200;
    std::unique_ptr<Hists> h_pt_2ndjet120, h_jets_pt_2ndjet120, h_ele_pt_2ndjet120, h_mu_pt_2ndjet120, h_event_pt_2ndjet120, h_topjets_pt_2ndjet120;
    std::unique_ptr<Hists> h_1cmstoptag, h_jets_1cmstoptag, h_ele_1cmstoptag, h_mu_1cmstoptag, h_event_1cmstoptag, h_topjets_1cmstoptag;
    std::unique_ptr<Hists> h_1heptoptag, h_jets_1heptoptag, h_ele_1heptoptag, h_mu_1heptoptag, h_event_1heptoptag, h_topjets_1heptoptag;
    std::unique_ptr<Hists> h_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection;
    
    MuonId MuIso;
    TopJetId cmsTopTag;
    TopJetId hepTopTag;
    JetId Btag_loose, Btag_medium, Btag_tight;
    
    std::vector<std::unique_ptr<AnalysisModule>> recomodules;
    uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps;
    
  };
  
  
  LQToTopMuAnalysisModule::LQToTopMuAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuAnalysisModule!" << endl;
    
    // 1. setup other modules. CommonModules and the JetCleaner:
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
    Muon_printer.reset(new MuonPrinter("Muon-Printer"));
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
    
    h_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>("HighMassReconstruction");
    
    // 2. set up selections
    //Selection
    njet_sel.reset(new NJetSelection(4, -1));
    nbtag_loose_sel.reset(new NJetSelection(2, -1, Btag_loose));  //default: (1, -1)
    nbtag_med_sel.reset(new NJetSelection(1, -1, Btag_medium));
    nbtag_tight_sel.reset(new NJetSelection(1, -1, Btag_tight));
    m_mumu_veto.reset(new InvMass2MuVeto(81, 101));
    nele_sel.reset(new NElectronSelection(1, -1));
    ht_sel.reset(new HtSelection(900, -1));
    pt_lead_mu_sel.reset(new PtLeadingMuonSelection(90,-1));
    pt_2nd_mu_sel.reset(new Pt2ndMuonSelection(90,-1));
    pt_lead_jet_sel.reset(new PtLeadingJetSelection(200,-1));
    pt_2nd_jet_sel.reset(new Pt2ndJetSelection(120,-1));
    cmsTopTag_sel.reset(new NTopJetSelection(1, -1, cmsTopTag));
    hepTopTag_sel.reset(new NTopJetSelection(1, -1, hepTopTag)); // use right topjet collection!
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new PrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassTTbarReconstruction(ctx,NeutrinoReconstruction));
    recomodules.emplace_back(new Chi2Discriminator(ctx,"HighMassReconstruction"));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));

    /*h_4jets.reset(new LQToTopMuHists(ctx, "4Jets"));
    h_jets_4jets.reset(new JetHists(ctx, "Jets_4Jets"));
    h_ele_4jets.reset(new ElectronHists(ctx, "Ele_4Jets"));
    h_mu_4jets.reset(new MuonHists(ctx, "Mu_4Jets"));
    h_event_4jets.reset(new EventHists(ctx, "Event_4Jets"));
    h_topjets_4jets.reset(new TopJetHists(ctx, "TopJets_4Jets"));*/

    /*h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
    h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
    h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
    h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));
    h_event_1bJetLoose.reset(new EventHists(ctx, "Event_1bJetLoose"));
    h_topjets_1bJetLoose.reset(new TopJetHists(ctx, "TopJets_1bJetLoose"));*/

    /*h_2bJetsLoose.reset(new LQToTopMuHists(ctx, "2bJetsLoose"));
    h_jets_2bJetsLoose.reset(new JetHists(ctx, "Jets_2bJetsLoose"));
    h_ele_2bJetsLoose.reset(new ElectronHists(ctx, "Ele_2bJetsLoose"));
    h_mu_2bJetsLoose.reset(new MuonHists(ctx, "Mu_2bJetsLoose"));
    h_event_2bJetsLoose.reset(new EventHists(ctx, "Event_2bJetsLoose"));
    h_topjets_2bJetsLoose.reset(new TopJetHists(ctx, "TopJets_2bJetsLoose"));

    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
    h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));*/
    
    /*h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
    h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
    h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
    h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));  
    h_event_1ele.reset(new EventHists(ctx, "Event_1Ele"));
    h_topjets_1ele.reset(new TopJetHists(ctx, "TopJets_1Ele"));
    h_hyphists.reset(new HypothesisHistsOwn(ctx, "Chi2_Hists", "HighMassReconstruction", "Chi2"));*/

    /*h_pt_leadmu90.reset(new LQToTopMuHists(ctx, "Pt_LeadMu90"));
    h_jets_pt_leadmu90.reset(new JetHists(ctx, "Jets_Pt_LeadMu90"));
    h_ele_pt_leadmu90.reset(new ElectronHists(ctx, "Ele_Pt_LeadMu90"));
    h_mu_pt_leadmu90.reset(new MuonHists(ctx, "Mu_Pt_LeadMu90"));
    h_topjets_pt_leadmu90.reset(new TopJetHists(ctx, "TopJets_Pt_LeadMu90"));
    h_event_pt_leadmu90.reset(new EventHists(ctx, "Event_Pt_LeadMu90"));

    h_pt_2ndmu90.reset(new LQToTopMuHists(ctx, "Pt_2ndMu90"));
    h_jets_pt_2ndmu90.reset(new JetHists(ctx, "Jets_Pt_2ndMu90"));
    h_ele_pt_2ndmu90.reset(new ElectronHists(ctx, "Ele_Pt_2ndMu90"));
    h_mu_pt_2ndmu90.reset(new MuonHists(ctx, "Mu_Pt_2ndMu90"));
    h_topjets_pt_2ndmu90.reset(new TopJetHists(ctx, "TopJets_Pt_2ndMu90"));
    h_event_pt_2ndmu90.reset(new EventHists(ctx, "Event_Pt_2ndMu90"));*/

    /*h_pt_leadjet200.reset(new LQToTopMuHists(ctx, "Pt_LeadJet200"));
    h_jets_pt_leadjet200.reset(new JetHists(ctx, "Jets_Pt_LeadJet200"));
    h_ele_pt_leadjet200.reset(new ElectronHists(ctx, "Ele_Pt_LeadJet200"));
    h_mu_pt_leadjet200.reset(new MuonHists(ctx, "Mu_Pt_LeadJet200"));
    h_topjets_pt_leadjet200.reset(new TopJetHists(ctx, "TopJets_Pt_LeadJet200"));
    h_event_pt_leadjet200.reset(new EventHists(ctx, "Event_Pt_LeadJet200"));*/

    /*h_pt_2ndjet120.reset(new LQToTopMuHists(ctx, "Pt_2ndJet120"));
    h_jets_pt_2ndjet120.reset(new JetHists(ctx, "Jets_Pt_2ndJet120"));
    h_ele_pt_2ndjet120.reset(new ElectronHists(ctx, "Ele_Pt_2ndJet120"));
    h_mu_pt_2ndjet120.reset(new MuonHists(ctx, "Mu_Pt_2ndJet120"));
    h_topjets_pt_2ndjet120.reset(new TopJetHists(ctx, "TopJets_Pt_2ndJet120"));
    h_event_pt_2ndjet120.reset(new EventHists(ctx, "Event_Pt_2ndJet120"));*/

    /*h_1cmstoptag.reset(new LQToTopMuHists(ctx, "1CMSTopTag"));
    h_jets_1cmstoptag.reset(new JetHists(ctx, "Jets_1CMSTopTag"));
    h_ele_1cmstoptag.reset(new ElectronHists(ctx, "Ele_1CMSTopTag"));
    h_mu_1cmstoptag.reset(new MuonHists(ctx, "Mu_1CMSTopTag"));
    h_topjets_1cmstoptag.reset(new TopJetHists(ctx, "TopJets_1CMSTopTag"));
    h_event_1cmstoptag.reset(new EventHists(ctx, "Event_1CMSTopTag"));*/

    /*h_1heptoptag.reset(new LQToTopMuHists(ctx, "1HEPTopTag"));
    h_jets_1heptoptag.reset(new JetHists(ctx, "Jets_1HEPTopTag"));
    h_ele_1heptoptag.reset(new ElectronHists(ctx, "Ele_1HEPTopTag"));
    h_mu_1heptoptag.reset(new MuonHists(ctx, "Mu_1HEPTopTag"));
    h_topjets_1heptoptag.reset(new TopJetHists(ctx, "TopJets_1HEPTopTag"));
    h_event_1heptoptag.reset(new EventHists(ctx, "Event_1HEPTopTag"));*/

    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection"));

    
  }
  
  
  bool LQToTopMuAnalysisModule::process(Event & event) {
    
    // 1. run all modules other modules.
    common->process(event);
    
    muoncleaner->process(event);
    muoncleaner_iso->process(event);
    electroncleaner->process(event);
    jetcleaner->process(event);

    //um MLQ_HT_Mix zu berechnen
    for(auto & m : recomodules){
      m->process(event);
      }
    
    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    
    // Selection   
    //Njets
    /*if(!njet_sel->passes(event)) return false;
    h_4jets->fill(event);
    h_jets_4jets->fill(event);
    h_ele_4jets->fill(event);
    h_mu_4jets->fill(event);
    h_event_4jets->fill(event);
    h_topjets_4jets->fill(event);*/

    //3x bTag: loose, medium, tight
    /*if(!nbtag_loose_sel->passes(event)) return false;
    h_jets_2bJetsLoose->fill(event);
    h_2bJetsLoose->fill(event);
    h_ele_2bJetsLoose->fill(event);
    h_mu_2bJetsLoose->fill(event);
    h_event_2bJetsLoose->fill(event);
    h_topjets_2bJetsLoose->fill(event);*/
    /*h_jets_1bJetLoose->fill(event);
    h_1bJetLoose->fill(event);
    h_ele_1bJetLoose->fill(event);
    h_mu_1bJetLoose->fill(event);
    h_event_1bJetLoose->fill(event);
    h_topjets_1bJetLoose->fill(event);*/


    //InvMassVeto
    /*if(!m_mumu_veto->passes(event)) return false;
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    h_event_InvMassVeto->fill(event);
    h_topjets_InvMassVeto->fill(event);*/
    
    //Nele
    /*if (!nele_sel->passes(event)) return false;
    h_1ele->fill(event);
    h_jets_1ele->fill(event);
    h_ele_1ele->fill(event);
    h_mu_1ele->fill(event);
    h_hyphists->fill(event);
    h_event_1ele->fill(event);
    h_topjets_1ele->fill(event);*/

    //pt leading mu
    /*if(!pt_lead_mu_sel->passes(event)) return false;
    h_jets_pt_leadmu90->fill(event);
    h_pt_leadmu90->fill(event);
    h_ele_pt_leadmu90->fill(event);
    h_mu_pt_leadmu90->fill(event);
    h_event_pt_leadmu90->fill(event);
    h_topjets_pt_leadmu90->fill(event);*/

    //pt 2nd mu
    /*if(!pt_2nd_mu_sel->passes(event)) return false;
    h_jets_pt_2ndmu90->fill(event);
    h_pt_2ndmu90->fill(event);
    h_ele_pt_2ndmu90->fill(event);
    h_mu_pt_2ndmu90->fill(event);
    h_event_pt_2ndmu90->fill(event);
    h_topjets_pt_2ndmu90->fill(event);*/

    //pt leading jet
    /*if(!pt_lead_jet_sel->passes(event)) return false;
    h_jets_pt_leadjet200->fill(event);
    h_pt_leadjet200->fill(event);
    h_ele_pt_leadjet200->fill(event);
    h_mu_pt_leadjet200->fill(event);
    h_event_pt_leadjet200->fill(event);
    h_topjets_pt_leadjet200->fill(event);*/

    //pt 2nd jet
    /*if(!pt_2nd_jet_sel->passes(event)) return false;
    h_jets_pt_2ndjet120->fill(event);
    h_pt_2ndjet120->fill(event);
    h_ele_pt_2ndjet120->fill(event);
    h_mu_pt_2ndjet120->fill(event);
    h_event_pt_2ndjet120->fill(event);
    h_topjets_pt_2ndjet120->fill(event);*/

    //CMSTopTag
    /*if(!cmsTopTag_sel->passes(event)) return false;
    h_jets_1cmstoptag->fill(event);
    h_1cmstoptag->fill(event);
    h_ele_1cmstoptag->fill(event);
    h_mu_1cmstoptag->fill(event);
    h_event_1cmstoptag->fill(event);
    h_topjets_1cmstoptag->fill(event);*/

    //HEPTopTag
    /*if(!hepTopTag_sel->passes(event)) return false;
    h_jets_1heptoptag->fill(event);
    h_1heptoptag->fill(event);
    h_ele_1heptoptag->fill(event);
    h_mu_1heptoptag->fill(event);
    h_event_1heptoptag->fill(event);
    h_topjets_1heptoptag->fill(event);*/

    //Final Selection
    h_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);

    return true;
  }
  
  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the LQToTopMuAnalysisModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuAnalysisModule)
} 


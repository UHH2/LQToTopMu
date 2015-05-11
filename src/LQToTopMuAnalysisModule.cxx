#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/PrintingModules.h"

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
    std::unique_ptr<Selection>  nele_sel, njet_sel, nbtag_loose_sel, nbtag_med_sel, nbtag_tight_sel, m_mumu_veto, ht_sel, pt_mu_sel;
    
    // store the Hists collection as member variables. 
    std::unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_4jets, h_jets_4jets, h_ele_4jets, h_mu_4jets, h_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_2bJetsLoose, h_jets_2bJetsLoose, h_mu_2bJetsLoose, h_ele_2bJetsLoose, h_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_ht900, h_jets_ht900, h_ele_ht900, h_mu_ht900, h_ptmu240, h_jets_ptmu240, h_ele_ptmu240, h_mu_ptmu240;
    std::unique_ptr<Hists> h_1bJetMedium, h_jets_1bJetMedium, h_mu_1bJetMedium, h_ele_1bJetMedium, h_1bJetTight, h_jets_1bJetTight, h_mu_1bJetTight, h_ele_1bJetTight;
    
    MuonId MuIso;
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
    nbtag_med_sel.reset(new NJetSelection(1, -1, Btag_medium));
    nbtag_tight_sel.reset(new NJetSelection(1, -1, Btag_tight));
    nbtag_loose_sel.reset(new NJetSelection(2, -1, Btag_loose));  //default: (1, -1)
    m_mumu_veto.reset(new InvMass2MuVeto(81, 101));
    nele_sel.reset(new NElectronSelection(1, -1));
    ht_sel.reset(new HtSelection(900, -1));
    pt_mu_sel.reset(new PtLeadingMuonSelection(240,-1));

    
    //make reconstruction hypotheses
    recomodules.emplace_back(new PrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassTTbarReconstruction(ctx,NeutrinoReconstruction));
    recomodules.emplace_back(new Chi2Discriminator(ctx,"HighMassReconstruction"));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));

    /*h_4jets.reset(new LQToTopMuHists(ctx, "4Jets"));
    h_jets_4jets.reset(new JetHists(ctx, "Jets_4Jets"));
    h_ele_4jets.reset(new ElectronHists(ctx, "Ele_4Jets"));
    h_mu_4jets.reset(new MuonHists(ctx, "Mu_4Jets"));*/

    /*h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
    h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
    h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
    h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));*/

    h_2bJetsLoose.reset(new LQToTopMuHists(ctx, "2bJetsLoose"));
    h_jets_2bJetsLoose.reset(new JetHists(ctx, "Jets_2bJetsLoose"));
    h_ele_2bJetsLoose.reset(new ElectronHists(ctx, "Ele_2bJetsLoose"));
    h_mu_2bJetsLoose.reset(new MuonHists(ctx, "Mu_2bJetsLoose"));
    
    /*h_1bJetMedium.reset(new LQToTopMuHists(ctx, "1bJetMedium"));
    h_jets_1bJetMedium.reset(new JetHists(ctx, "Jets_1bJetMedium"));
    h_ele_1bJetMedium.reset(new ElectronHists(ctx, "Ele_1bJetMedium"));
    h_mu_1bJetMedium.reset(new MuonHists(ctx, "Mu_1bJetMedium"));*/
    
    /*h_1bJetTight.reset(new LQToTopMuHists(ctx, "1bJetTight"));
    h_jets_1bJetTight.reset(new JetHists(ctx, "Jets_1bJetTight"));
    h_ele_1bJetTight.reset(new ElectronHists(ctx, "Ele_1bJetTight"));
    h_mu_1bJetTight.reset(new MuonHists(ctx, "Mu_1bJetTight"));*/

    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    
    h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
    h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
    h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
    h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));
    
    h_hyphists.reset(new HypothesisHistsOwn(ctx, "Chi2_Hists", "HighMassReconstruction", "Chi2"));

    /* h_ht900.reset(new LQToTopMuHists(ctx, "HT900"));
    h_jets_ht900.reset(new JetHists(ctx, "Jets_HT900"));
    h_ele_ht900.reset(new ElectronHists(ctx, "Ele_HT900"));
    h_mu_ht900.reset(new MuonHists(ctx, "Mu_HT900"));*/

    /*h_ptmu240.reset(new LQToTopMuHists(ctx, "PtMu240"));
    h_jets_ptmu240.reset(new JetHists(ctx, "Jets_PtMu240"));
    h_ele_ptmu240.reset(new ElectronHists(ctx, "Ele_PtMu240"));
    h_mu_ptmu240.reset(new MuonHists(ctx, "Mu_PtMu240"));*/

    
  }
  
  
  bool LQToTopMuAnalysisModule::process(Event & event) {
    
    // 1. run all modules other modules.
    common->process(event);
    
    muoncleaner->process(event);
    muoncleaner_iso->process(event);
    electroncleaner->process(event);
    jetcleaner->process(event);

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    
    // Selection   
    //Njets
    /*if(!njet_sel->passes(event)) return false;
    h_4jets->fill(event);
    h_jets_4jets->fill(event);
    h_ele_4jets->fill(event);
    h_mu_4jets->fill(event);*/

    //3x bTag: loose, medium, tight
    if(!nbtag_loose_sel->passes(event)) return false;
    h_jets_2bJetsLoose->fill(event);
    h_2bJetsLoose->fill(event);
    h_ele_2bJetsLoose->fill(event);
    h_mu_2bJetsLoose->fill(event); 
    /*h_jets_1bJetLoose->fill(event);
    h_1bJetLoose->fill(event);
    h_ele_1bJetLoose->fill(event);
    h_mu_1bJetLoose->fill(event);*/

    /*if(!nbtag_med_sel->passes(event)) return false;
    h_jets_1bJetMedium->fill(event);
    h_1bJetMedium->fill(event);
    h_ele_1bJetMedium->fill(event);
    h_mu_1bJetMedium->fill(event);*/
    
    /*if(!nbtag_tight_sel->passes(event)) return false;
    h_jets_1bJetTight->fill(event);
    h_1bJetTight->fill(event);
    h_ele_1bJetTight->fill(event);
    h_mu_1bJetTight->fill(event);*/

    //InvMassVeto
    if(!m_mumu_veto->passes(event)) return false;
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    
    //Nele
    if (!nele_sel->passes(event)) return false;
    h_1ele->fill(event);
    h_jets_1ele->fill(event);
    h_ele_1ele->fill(event);
    h_mu_1ele->fill(event);
    
    for(auto & m : recomodules){
      m->process(event);
    }
    h_hyphists->fill(event);

    //HT > 500 GeV
    /*  if(!ht_sel->passes(event)) return false;
    h_ht900->fill(event);
    h_jets_ht900->fill(event);
    h_ele_ht900->fill(event);
    h_mu_ht900->fill(event);*/

    //PT leading muon > 150 GeV/c
    /*if(!pt_mu_sel->passes(event)) return false;
    h_ptmu240->fill(event);
    h_jets_ptmu240->fill(event);
    h_ele_ptmu240->fill(event);
    h_mu_ptmu240->fill(event);*/
    
    return true;
  }
  
  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the LQToTopMuAnalysisModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuAnalysisModule)
} 


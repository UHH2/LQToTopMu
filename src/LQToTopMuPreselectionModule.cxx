#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class LQToTopMuPreselectionModule: public AnalysisModule {
public:
    
    explicit LQToTopMuPreselectionModule(Context & ctx);
    virtual bool process(Event & event) override;
  
private:
  
  std::unique_ptr<CommonModules> common;
  
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<JetLeptonCleaner> jetleptoncleaner;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<MuonCleaner> muoncleaner_iso;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  
  // declare the Selections to use.
  std::unique_ptr<Selection> njet_sel, nmuon_sel, ht_sel;
  
  // store the Hists collection as member variables. 
  std::unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets,
    h_ht350, h_jets_ht350, h_ele_ht350, h_mu_ht350, h_2mu, h_jets_2mu, h_ele_2mu, h_mu_2mu;

  MuonId MuIso;
};


LQToTopMuPreselectionModule::LQToTopMuPreselectionModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuPreselectionModule!" << endl;
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;

    



    // 1. setup other modules. CommonModules and the JetCleaner:
    common.reset(new CommonModules());
    common->init(ctx);
    MuIso = MuonIso(0.12);
    jetleptoncleaner.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));
    jetcleaner.reset(new JetCleaner(30.0, 2.5)); 
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.1))));
    muoncleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuIso, PtEtaCut(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(30.0, 2.5))));

    
    // 2. set up selections

    //Preselection
    njet_sel.reset(new NJetSelection(2, -1)); // see common/include/NSelections.h
    nmuon_sel.reset(new NMuonSelection(2, -1));
    ht_sel.reset(new HtSelection(350));


    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));

    h_cleaner.reset(new LQToTopMuHists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));

    h_2jets.reset(new LQToTopMuHists(ctx, "2Jets"));
    h_jets_2jets.reset(new JetHists(ctx, "Jets_2Jets"));
    h_ele_2jets.reset(new ElectronHists(ctx, "Ele_2Jets"));
    h_mu_2jets.reset(new MuonHists(ctx, "Mu_2Jets"));

    h_ht350.reset(new LQToTopMuHists(ctx, "HT350"));
    h_jets_ht350.reset(new JetHists(ctx, "Jets_HT350"));
    h_ele_ht350.reset(new ElectronHists(ctx, "Ele_HT350"));
    h_mu_ht350.reset(new MuonHists(ctx, "Mu_HT350"));

    h_2mu.reset(new LQToTopMuHists(ctx, "2Mu"));
    h_jets_2mu.reset(new JetHists(ctx, "Jets_2Mu"));
    h_ele_2mu.reset(new ElectronHists(ctx, "Ele_2Mu"));
    h_mu_2mu.reset(new MuonHists(ctx, "Mu_2Mu"));
}


bool LQToTopMuPreselectionModule::process(Event & event) {

    
  // 1. run all modules other modules.
  common->process(event);

  h_nocuts->fill(event);
  h_jets_nocuts->fill(event);
  h_ele_nocuts->fill(event);
  h_mu_nocuts->fill(event);

  
  muoncleaner->process(event);
  muoncleaner_iso->process(event);
  electroncleaner->process(event);
  jetleptoncleaner->process(event);
  jetcleaner->process(event);
  
  h_cleaner->fill(event);
  h_jets_cleaner->fill(event);
  h_ele_cleaner->fill(event);
  h_mu_cleaner->fill(event);
  
  // 2. test selections and fill histograms
  
  //--  Preselection  
  //auto jets = event.jets;
  if (!njet_sel->passes(event)) return false;
  h_2jets->fill(event);
  h_jets_2jets->fill(event);
  h_ele_2jets->fill(event);
  h_mu_2jets->fill(event);

  if (!ht_sel->passes(event)) return false;
  h_ht350->fill(event);
  h_jets_ht350->fill(event);
  h_ele_ht350->fill(event);
  h_mu_ht350->fill(event);

  if (!nmuon_sel->passes(event)) return false;
  h_2mu->fill(event);
  h_jets_2mu->fill(event);
  h_ele_2mu->fill(event);
  h_mu_2mu->fill(event);
  
  
  // 3. decide whether or not to keep the current event in the output:
  return true;
}
  
  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the LQToTopMuPreselectionModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuPreselectionModule)
  
}

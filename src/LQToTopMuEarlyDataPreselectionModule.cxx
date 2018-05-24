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
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuPreselectionHists.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuEarlyDataPreselectionModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuEarlyDataPreselectionModule(Context & ctx);
    virtual bool process(Event & event) override;
  
  private:
  
    unique_ptr<CommonModules> common;  
    unique_ptr<JetCleaner> jetcleaner;
  
    // declare the Selections to use.
    unique_ptr<Selection> njet_sel, nmuon_sel, nele_sel, trigger_sel1, trigger_sel2;
  
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_lumi_nocuts, 
      h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_lumi_cleaner, 
      h_1mu, h_jets_1mu, h_ele_1mu, h_mu_1mu, h_event_1mu, h_lumi_1mu, 
      h_muontrigger, h_jets_muontrigger, h_ele_muontrigger, h_mu_muontrigger, h_event_muontrigger, h_lumi_muontrigger, 
      h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets, h_event_2jets, h_lumi_2jets;

    MuonId MuId;
    ElectronId EleId;
    bool is_mc;
    TString triggermode;

  };


  LQToTopMuEarlyDataPreselectionModule::LQToTopMuEarlyDataPreselectionModule(Context & ctx){
    
    cout << "Hello from LQToTopMuPreselectionModule!" << endl;
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    is_mc = ctx.get("dataset_type") == "MC";
    triggermode = ctx.get("Trigger");
    if(triggermode != "Iso" && triggermode != "NonIso") throw runtime_error("Invalid triggermode specified. Allowed are 'Iso' and 'NonIso' in the PreSelection.");
    if(triggermode == "Iso")          MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4),MuonIso(0.15));
    else if (triggermode == "NonIso") MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(55.0, 2.4));

    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    //common->set_electron_id(EleId);
    common->disable_lumisel();
    common->disable_mcpileupreweight();
    common->set_muon_id(MuId);
    common->init(ctx);

    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));


    // 2. set up selections

    //Preselection

    trigger_sel1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    //trigger_sel2.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    njet_sel.reset(new NJetSelection(2, -1));
    nmuon_sel.reset(new NMuonSelection(1, 1)); 

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuPreselectionHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_cleaner.reset(new LQToTopMuPreselectionHists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_Cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_Cleaner"));

    h_1mu.reset(new LQToTopMuPreselectionHists(ctx, "1Mu"));
    h_jets_1mu.reset(new JetHists(ctx, "Jets_1Mu"));
    h_ele_1mu.reset(new ElectronHists(ctx, "Ele_1Mu"));
    h_mu_1mu.reset(new MuonHists(ctx, "Mu_1Mu"));
    h_event_1mu.reset(new EventHists(ctx, "Event_1Mu"));
    h_lumi_1mu.reset(new LuminosityHists(ctx, "Lumi_1Mu"));

    h_muontrigger.reset(new LQToTopMuPreselectionHists(ctx, "MuonTrigger"));
    h_jets_muontrigger.reset(new JetHists(ctx, "Jets_MuonTrigger"));
    h_ele_muontrigger.reset(new ElectronHists(ctx, "Ele_MuonTrigger"));
    h_mu_muontrigger.reset(new MuonHists(ctx, "Mu_MuonTrigger"));
    h_event_muontrigger.reset(new EventHists(ctx, "Event_MuonTrigger"));
    h_lumi_muontrigger.reset(new LuminosityHists(ctx, "Lumi_MuonTrigger"));

    h_2jets.reset(new LQToTopMuPreselectionHists(ctx, "2Jets"));
    h_jets_2jets.reset(new JetHists(ctx, "Jets_2Jets"));
    h_ele_2jets.reset(new ElectronHists(ctx, "Ele_2Jets"));
    h_mu_2jets.reset(new MuonHists(ctx, "Mu_2Jets"));
    h_event_2jets.reset(new EventHists(ctx, "Event_2Jets"));
    h_lumi_2jets.reset(new LuminosityHists(ctx, "Lumi_2Jets"));


  }


  bool LQToTopMuEarlyDataPreselectionModule::process(Event & event) {



    //non-iso trigger is active starting from RunC
    if(triggermode == "NonIso") if(event.run < 299368) return false;

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);

    h_cleaner->fill(event);
    h_jets_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_event_cleaner->fill(event);
    h_lumi_cleaner->fill(event);

    if(!nmuon_sel->passes(event)) return false;
    h_1mu->fill(event);
    h_jets_1mu->fill(event);
    h_ele_1mu->fill(event);
    h_mu_1mu->fill(event);
    h_event_1mu->fill(event);
    h_lumi_1mu->fill(event);

    if(trigger_sel1->passes(event)){
      h_muontrigger->fill(event);
      h_jets_muontrigger->fill(event);
      h_ele_muontrigger->fill(event);
      h_mu_muontrigger->fill(event);
      h_event_muontrigger->fill(event);
      h_lumi_muontrigger->fill(event);
    }
  
    if(!njet_sel->passes(event)) return false;
    h_2jets->fill(event);
    h_jets_2jets->fill(event);
    h_ele_2jets->fill(event);
    h_mu_2jets->fill(event);
    h_event_2jets->fill(event);
    h_lumi_2jets->fill(event);

    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuEarlyDataPreselectionModule)
  
}

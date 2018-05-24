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
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/LQToTopMuTagProbeHists.h"
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

  class LQToTopMuEarlyDataAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuEarlyDataAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    unique_ptr<CommonModules> common;
    unique_ptr<JetCleaner> jetcleaner;
    unique_ptr<ElectronCleaner> elecleaner;
    
    // declare the Selections to use.
    unique_ptr<Selection>  trigger_iso_sel, trigger_noniso_sel, trigger_mu_sel, nele_sel_iso_lowpt, nele_sel_iso_highpt, nele_sel_noniso_highpt, nele_sel;
    
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_controlhists_ele_input, h_controlhists_mu_input, h_controlhists_jets_input, h_controlhists_event_input, h_controlhists_lumi_input, h_tagprobe_base, h_controlhists_ele_base, h_controlhists_mu_base, h_controlhists_jets_base, h_controlhists_event_base, h_controlhists_lumi_base, h_tagprobe_plateau_iso_inclusive_before, h_tagprobe_plateau_noniso_inclusive_before, h_tagprobe_plateau_iso_lowpt_before, h_tagprobe_plateau_iso_highpt_before,
                      h_tagprobe_plateau_iso_50_before, h_tagprobe_trigger, h_tagprobe_plateau_iso_inclusive_after, h_tagprobe_plateau_noniso_inclusive_after, h_tagprobe_plateau_iso_lowpt_after, h_tagprobe_plateau_iso_highpt_after, h_tagprobe_plateau_iso_50_after;
    
    MuonId MuId_iso, MuId_noniso;
    ElectronId EleId_base, EleId_noniso_base, EleId_iso_lowpt, EleId_iso_highpt, EleId_noniso_highpt;


    bool is_mc;
    string Sys_PU;
    TString triggermode;
    
  };
  
  
  LQToTopMuEarlyDataAnalysisModule::LQToTopMuEarlyDataAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuEarlyDataAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
    is_mc = ctx.get("dataset_type") == "MC";
    Sys_PU = ctx.get("Systematic_PU");
    triggermode = ctx.get("Trigger");
    if(triggermode != "Iso" && triggermode != "NonIso" && triggermode != "Comb") throw runtime_error("Invalid triggermode 'Trigger' specified in xml file.");
    string triggerpath_e1 = "HLT_Ele27_WPTight_Gsf_v*";
    //string triggerpath_e1 = "HLT_Ele35_WPTight_Gsf_v*";
    //string triggerpath_e1 = "HLT_Ele38_WPTight_Gsf_v*";
    //string triggerpath_e1 = "HLT_Ele40_WPTight_Gsf_v*";
    string triggerpath_e2 = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*TEST";

    string triggerpath_mu;
    if(triggermode == "Iso" || triggermode == "Comb") triggerpath_mu  = "HLT_IsoMu24_v*";
    else if(triggermode == "NonIso") triggerpath_mu  = "HLT_Mu50_v*";


    // 1. setup other modules. CommonModules and the JetCleaner:
    MuId_iso            = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4),MuonIso(0.15));
    MuId_noniso         = AndId<Muon>(MuonIDTight(),PtEtaCut(55.0, 2.4));
    EleId_base          = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(10.0, 2.4));          //10!!!
    EleId_noniso_base   = AndId<Electron>(ElectronID_Spring16_tight_noIso,PtEtaCut(10.0, 2.4));
    EleId_iso_lowpt     = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(30.0, 2.4)); 
    EleId_iso_highpt    = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(120.0, 2.4)); 
    EleId_noniso_highpt = AndId<Electron>(ElectronID_Spring16_tight_noIso,PtEtaCut(120.0, 2.4)); 
    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));


    common.reset(new CommonModules());
    //common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->disable_mcpileupreweight();
    if(triggermode == "Iso" || triggermode == "Comb"){
      common->set_electron_id(EleId_base);
      common->set_muon_id(MuId_iso);
    }
    else if(triggermode == "NonIso"){
      common->set_electron_id(EleId_noniso_base);
      common->set_muon_id(MuId_noniso);  
    }
    common->init(ctx,Sys_PU);

    
    // 2. set up selections
    //Selection
    trigger_mu_sel.reset(new TriggerSelection(triggerpath_mu));
    trigger_iso_sel.reset(new TriggerSelection(triggerpath_e1));
    trigger_noniso_sel.reset(new TriggerSelection(triggerpath_e2)); 

    nele_sel.reset(new NElectronSelection(1,1));
    nele_sel_iso_lowpt    .reset(new NElectronSelection(1, 1, EleId_iso_lowpt));
    nele_sel_iso_highpt   .reset(new NElectronSelection(1, 1, EleId_iso_highpt));
    nele_sel_noniso_highpt.reset(new NElectronSelection(1, 1, EleId_noniso_highpt));
    
    // 3. Set up Hists classes:    
    h_tagprobe_base.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Base"));  
    h_controlhists_ele_input.reset(new ElectronHists(ctx, "Control_Ele_Input"));     
    h_controlhists_mu_input.reset(new MuonHists(ctx, "Control_Mu_Input"));     
    h_controlhists_jets_input.reset(new JetHists(ctx, "Control_Jets_Input"));     
    h_controlhists_event_input.reset(new EventHists(ctx, "Control_Event_Input"));     
    h_controlhists_lumi_input.reset(new LuminosityHists(ctx, "Control_Lumi_Input"));     
    h_controlhists_ele_base.reset(new ElectronHists(ctx, "Control_Ele_Base"));     
    h_controlhists_mu_base.reset(new MuonHists(ctx, "Control_Mu_Base"));     
    h_controlhists_jets_base.reset(new JetHists(ctx, "Control_Jets_Base"));     
    h_controlhists_event_base.reset(new EventHists(ctx, "Control_Event_Base"));     
    h_controlhists_lumi_base.reset(new LuminosityHists(ctx, "Control_Lumi_Base"));     
    h_tagprobe_plateau_iso_inclusive_before.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_Inclusive_Before"));  
    h_tagprobe_plateau_noniso_inclusive_before.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_NonIso_Inclusive_Before"));   
    h_tagprobe_plateau_iso_lowpt_before.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_Lowpt_Before")); 
    h_tagprobe_plateau_iso_highpt_before.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_Highpt_Before")); 
    h_tagprobe_plateau_iso_50_before.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_50_Before")); 

    h_tagprobe_trigger.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Trigger"));    
    h_tagprobe_plateau_iso_inclusive_after.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_Inclusive_After"));   
    h_tagprobe_plateau_noniso_inclusive_after.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_NonIso_Inclusive_After"));   
    h_tagprobe_plateau_iso_lowpt_after.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_Lowpt_After")); 
    h_tagprobe_plateau_iso_highpt_after.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_Highpt_After")); 
    h_tagprobe_plateau_iso_50_after.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau_Iso_50_After")); 

    
    
  }
  
  
  bool LQToTopMuEarlyDataAnalysisModule::process(Event & event) {    
    //cout << endl << endl << "---------- New Event ----------" << endl << endl;;

    //Non-iso trigger present only starting from RunC, which starts at 299368
    if(triggermode == "NonIso" || triggermode == "Comb"){
      if(event.run < 299368) return false;
    }

    //if(event.run > 299329) return false;
    //if(event.run < 299368 || event.run > 302019) return false;
    //if(event.run < 302026) return false;

    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);
    
    h_controlhists_ele_input->fill(event);
    h_controlhists_mu_input->fill(event);
    h_controlhists_jets_input->fill(event);
    h_controlhists_event_input->fill(event);
    h_controlhists_lumi_input->fill(event);
    
    if(!trigger_mu_sel->passes(event)) return false;
    if(!nele_sel->passes(event)) return false;
    h_controlhists_ele_base->fill(event);
    h_controlhists_mu_base->fill(event);
    h_controlhists_jets_base->fill(event);
    h_controlhists_event_base->fill(event);
    h_controlhists_lumi_base->fill(event);

    h_tagprobe_base->fill(event);

    
    if(nele_sel_iso_lowpt->passes(event)) 
      h_tagprobe_plateau_iso_inclusive_before->fill(event);

    if(nele_sel_noniso_highpt->passes(event)) 
      h_tagprobe_plateau_noniso_inclusive_before->fill(event);

    if(triggermode == "Comb"){
      if(nele_sel_iso_lowpt->passes(event) && !nele_sel_iso_highpt->passes(event)) 
	h_tagprobe_plateau_iso_lowpt_before->fill(event);
    

      if(nele_sel_iso_highpt->passes(event)) 
	h_tagprobe_plateau_iso_highpt_before->fill(event);
    }

    if(event.electrons->at(0).pt() >= 50) 
      h_tagprobe_plateau_iso_50_before->fill(event);

    
    /* ---------------------- Electron Trigger --------------------*/
    if(triggermode == "Iso"){
      if(!trigger_iso_sel->passes(event)) return false;
    }
    else if(triggermode == "NonIso"){
      if(!trigger_noniso_sel->passes(event)) return false;
    }
    else if(triggermode == "Comb"){
      if(!(trigger_iso_sel->passes(event) || trigger_noniso_sel->passes(event))) return false;
    }
    /* ------------------------------------------------------------*/
    
    h_tagprobe_trigger->fill(event);

    if(nele_sel_iso_lowpt->passes(event))
      h_tagprobe_plateau_iso_inclusive_after->fill(event);


    if(nele_sel_noniso_highpt->passes(event))
      h_tagprobe_plateau_noniso_inclusive_after->fill(event);
    

    if(triggermode == "Comb"){
      if(nele_sel_iso_lowpt->passes(event) && !nele_sel_iso_highpt->passes(event))
	h_tagprobe_plateau_iso_lowpt_after->fill(event);
    

      if(nele_sel_iso_highpt->passes(event))
	h_tagprobe_plateau_iso_highpt_after->fill(event);
    }

    if(event.electrons->at(0).pt() >= 50)
      h_tagprobe_plateau_iso_50_after->fill(event);
    


    return false;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuEarlyDataAnalysisModule)
} 


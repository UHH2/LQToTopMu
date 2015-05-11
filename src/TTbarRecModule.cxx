#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/HypothesisHists.h"
//#include "UHH2/LQToTopMu/include/HypothesisHists.h"
#include "UHH2/common/include/ElectronHists.h"

using namespace std;
using namespace uhh2;

/** 
 * Example Module that performs a selection of two jets and one muon and fills the TTbar reconstruction hypotheses for selected events into the ouptut tree
 */
class TTbarRecModule: public AnalysisModule {
public:
    
  explicit TTbarRecModule(Context & ctx);
  virtual bool process(Event & event);
  
private:
  JetId jet_kinematic;
  TopJetId topjet_kinematic;
  MuonId muid;
  
  std::vector<std::unique_ptr<AnalysisModule>> cleanermodules;
  std::vector<std::unique_ptr<AnalysisModule>> recomodules;
  std::unique_ptr<Hists> h_hyphists, h_ele;
  AndSelection selection;
  
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps;
  
};


TTbarRecModule::TTbarRecModule(Context & ctx): selection(ctx, "selection") {
  jet_kinematic = PtEtaCut(50.0, 2.4);
  //jet_kinematic = PtEtaCut(30.0, 2.5);
  //topjet_kinematic = PtEtaCut(150.0,2.4);
  muid = AndId<Muon>(MuonIDTight(), PtEtaCut(20.0, 2.4));
  
  // clean the objects:
  cleanermodules.emplace_back(new JetCleaner(jet_kinematic));
  cleanermodules.emplace_back(new MuonCleaner(muid));
  cout << "test3" << endl;
  
  // make the selection, step-by-step. Note that in most cases, no explicit
  // object id is passed, as the cleaners have removed the objects not passing the id already.
  selection.add<NMuonSelection>("nm >= 1", 1);
  //selection.add<NElectronSelection>("ne >= 1",1); //
  selection.add<NJetSelection>("nj >= 2", 2);
  selection.add<NTopJetSelection>("ntj >= 1", 1);
  cout << "test4" << endl;
  
  //make reconstruction hypotheses
  recomodules.emplace_back(new PrimaryLepton(ctx));
  recomodules.emplace_back(new HighMassTTbarReconstruction(ctx,NeutrinoReconstruction));
  recomodules.emplace_back(new Chi2Discriminator(ctx,"HighMassReconstruction"));
  cout << "test5" << endl;
  
  //set up the histograms
  h_hyphists.reset(new HypothesisHists(ctx,"Chi2Hists","HighMassReconstruction","Chi2"));
  cout << "test6" << endl;
  h_ele.reset(new ElectronHists(ctx, "ElectronHists"));
  cout << "test7" << endl;
  
  h_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>("HighMassReconstruction");
  cout << "test8" << endl;
  
}


bool TTbarRecModule::process(Event & event) {
  for(auto & m : cleanermodules){
    m->process(event);
  }
  if(selection.passes(event)){
    for(auto & m : recomodules){
      m->process(event);
    }
    
    h_hyphists->fill(event);
    //h_ele->fill(event);
    
    return true;
  }
  return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the ExampleModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TTbarRecModule)

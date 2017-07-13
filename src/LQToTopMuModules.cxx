#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include <TH1D.h>


using namespace uhh2;
using namespace std;

JetLeptonOverlapCleaner::JetLeptonOverlapCleaner(double RJet_): RJet(RJet_){}
  
bool JetLeptonOverlapCleaner::process(Event & event){

   vector<Jet> result_jets;
   vector<Muon> result_muons;
   vector<Electron> result_electrons;

   cout << "############# Start of a new event ##################" << endl << endl;

   cout << "------------- first deal with the jets --------------" << endl << endl;
   //Check, if jets are identical to leptons
   //if so, kick them out
   for(const Jet & thisjet : *event.jets){
     auto thisjet_v4_raw = thisjet.v4() * thisjet.JEC_factor_raw();
     bool keep_thisjet = true;
     cout << "new jet: " << endl;
     for(const Muon & thismu : *event.muons){
       //if true, jet is considered identical to the muon
       cout << "deltaR between this muon and the jet: " << deltaR(thismu,thisjet) << ", mupt/jetptraw: " << thismu.pt()/thisjet_v4_raw.pt() << endl;
       if (deltaR(thismu,thisjet) <= 0.05 && (thismu.pt() <= thisjet_v4_raw.pt()*1.1 && thismu.pt() >= thisjet_v4_raw.pt()*0.8)){
	 keep_thisjet = false;
	 cout << "at least bc/ of this, this jet will not be kept" << endl;
       }
     }//end muons
     for(const Electron & thisele : *event.electrons){
       //if true, jet is considered identical to the electron
      cout << "deltaR between this electron and the jet: " << deltaR(thisele,thisjet) << ", elept/jetptraw: " << thisele.pt()/thisjet_v4_raw.pt() << endl;
       if (deltaR(thisele,thisjet) <= 0.05 && (thisele.pt() <= thisjet_v4_raw.pt()*1.1 && thisele.pt() >= thisjet_v4_raw.pt()*0.8)){
	 keep_thisjet = false;
	 cout << "at least bc/ of this, this jet will not be kept" << endl;
       }
     }//end electrons
     if(keep_thisjet) result_jets.push_back(thisjet);
   }
   //final collection of jets that are not identical to leptons
   swap(*event.jets,result_jets);



   //now check, if leptons are lying inside a jet
   //if so, kick the lepton out

   //first for muons
   cout << "----------------- coming to the muons ------------" << endl << endl;
   for(const Muon & thismu : *event.muons){
     bool keep_thismu = true;
     cout << "new muon: " << endl;
     //already kicked out jets that actually were leptons
     for(const Jet & thisjet : *event.jets){
       cout << "deltaR between this muon and a jet: " <<  deltaR(thismu,thisjet) << endl;
       if(deltaR(thismu,thisjet) <= RJet) {
	 keep_thismu = false;
	 cout << "at least bc/ of this, this muon will not be kept" << endl;
       }
     }
     if(keep_thismu) result_muons.push_back(thismu);
   }
   
   //the same for electrons
   cout << "----------------- coming to the electrons ------------" << endl << endl;
   for(const Electron & thisele : *event.electrons){
     bool keep_thisele = true;
     cout << "new electron: " << endl;
     //already kicked out jets that actually were leptons
     for(const Jet & thisjet : *event.jets){
       cout << "deltaR between this electron and a jet: " <<  deltaR(thisele,thisjet) << endl;
       if(deltaR(thisele,thisjet) <= RJet){
	 keep_thisele = false;
	 cout << "at least bc/ of this, this muon will not be kept" << endl;
       }
     }
     if(keep_thisele) result_electrons.push_back(thisele);
   }

   swap(*event.muons,result_muons);
   swap(*event.electrons,result_electrons);

   cout << "############### End of the event ###################" << endl << endl;
   return true;
}


ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  
  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  
  Eff_lowpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_TTbar_eff"));					 
  Eff_highpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_TTbar_eff"));					 
  Eff_lowpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_DATA_eff"));					 
  Eff_highpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_DATA_eff"));					 
  
  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
  
}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData) return true;

  //cout << "Weight before SF: " << event.weight << endl;
  double prob_notrig_mc = 1, prob_notrig_data = 1;
  for(const auto & ele : *event.electrons){
    double eta = ele.eta();
    if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Ele-|eta| > 2.4 is not supported at the moment.");


    //find right bin in eta
    int idx = 0;
    bool lowpt = false;
    if(30 <= ele.pt() && ele.pt() < 120){
      lowpt = true;
      //lowpt trigger
      bool keep_going = true;
      while(keep_going){
	double x,y;
	Eff_lowpt_MC->GetPoint(idx,x,y);
	keep_going = eta > x + Eff_lowpt_MC->GetErrorXhigh(idx);
	if(keep_going) idx++;
      }
    }
    else if(ele.pt() >= 120){
     //highpt trigger
      bool keep_going = true;
      while(keep_going){
	double x,y;
	Eff_highpt_MC->GetPoint(idx,x,y);
	keep_going = eta > x + Eff_highpt_MC->GetErrorXhigh(idx);
	if(keep_going) idx++;
      }
    }
    else throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Electron has pt<30. Clean electron collection before applying weights.");

    //access efficiencies for MC and DATA, possibly accout for systematics = statistical + add. 2% up/down
    double eff_data = -1, eff_mc = -1, dummy_x;
    double stat_data = -1, stat_mc = -1, tp = 0.02, total_syst_data = -1, total_syst_mc = -1;
    if(lowpt){
      Eff_lowpt_MC->GetPoint(idx,dummy_x,eff_mc);
      Eff_lowpt_DATA->GetPoint(idx,dummy_x,eff_data);

      if(SysDirection == "up"){		
	stat_mc = Eff_lowpt_MC->GetErrorYlow(idx);	
	stat_data = Eff_lowpt_DATA->GetErrorYhigh(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      	eff_mc -= total_syst_mc;    	
      	eff_data += total_syst_data;	
      }							
      else if(SysDirection == "down"){
	stat_mc = Eff_lowpt_MC->GetErrorYhigh(idx);	
	stat_data = Eff_lowpt_DATA->GetErrorYlow(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
			
      	eff_mc += Eff_lowpt_MC->GetErrorYhigh(idx);    	
      	eff_data -= Eff_lowpt_DATA->GetErrorYlow(idx);	
      }                                                 
    }
    else{
      Eff_highpt_MC->GetPoint(idx,dummy_x,eff_mc);
      Eff_highpt_DATA->GetPoint(idx,dummy_x,eff_data);

      if(SysDirection == "up"){	
	stat_mc = Eff_highpt_MC->GetErrorYlow(idx);	
	stat_data = Eff_highpt_DATA->GetErrorYhigh(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
			
	eff_mc -= Eff_highpt_MC->GetErrorYlow(idx);    	
	eff_data += Eff_highpt_DATA->GetErrorYhigh(idx);	
      }							
      else if(SysDirection == "down"){	
	stat_mc = Eff_highpt_MC->GetErrorYhigh(idx);	
	stat_data = Eff_highpt_DATA->GetErrorYlow(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
					
	eff_mc += Eff_highpt_MC->GetErrorYhigh(idx);    	
	eff_data -= Eff_highpt_DATA->GetErrorYlow(idx);	
      }                                                 
    }

    //multiply to the efficiency for not triggering
    prob_notrig_mc *= 1-eff_mc;
    prob_notrig_data *= 1-eff_data;

    //cout << "Efficiency for this ele -- MC: " << eff_mc << ", DATA" << eff_data << endl;
    //cout << "prob for not triggering -- MC: " << prob_notrig_mc << ", DATA: " << prob_notrig_data << endl;

  }

  //Scale weight by (1-prob_notrig_data) / (1-prob_notrig_mc)
  double SF = (1-prob_notrig_data)/(1-prob_notrig_mc);
  event.weight *= SF;

  //cout << "Weight with SF: " << event.weight << endl; 
  return true;
}

JetCorrectorVariable::JetCorrectorVariable(uhh2::Context & ctx, const std::vector<std::string> & JEC_files): JetCorrector(ctx, JEC_files){}

bool JetCorrectorVariable::correct_collection(uhh2::Event & event, std::vector<Jet> & jets){  
  
    //apply jet corrections
    for(auto & jet : jets){
      correct_jet(*corrector, jet, event, jec_uncertainty, direction);
    }
    return true;
};

/*
JetSmearerVariable::JetSmearerVariable(uhh2::Context & ctx, const string label_jets, const string label_genjets, const JERSmearing::SFtype1& JER_sf) : GenericJetResolutionSmearer(ctx,label_jets,label_genjets,true,JER_sf){}

bool JetResolutionSmearer::process(uhh2::Event & event) {

  m_gjrs->process(event);
  return true;
}

JetResolutionSmearer::~JetResolutionSmearer(){}
*/


/*
JetSmearerVariable::JetSmearerVariable(uhh2::Context & ctx, const JERSmearing::SFtype1& JER_sf): JetResolutionSmearer(ctx, JER_sf){}

bool JetSmearerVariable::smear_collection(uhh2::Event & event, std::vector<Jet> & jets, std::vector<Particle> &genjets){
    
    LorentzVector met;
    if(event.met) {
      met = event.met->v4();
    }
    for(unsigned int i=0; i<jets.size(); ++i) {
      auto & jet = jets.at(i);
      // find next genjet:
      auto closest_genjet = closestParticle(jet, genjets);
      // ignore unmatched jets (=no genjets at all or large DeltaR), or jets with very low genjet pt:
      if(closest_genjet == nullptr || deltaR(*closest_genjet, jet) > 0.3) continue;
      auto genpt = closest_genjet->pt();
      if(genpt < 15.0f) {
	continue;
      }
      LorentzVector jet_v4 = jet.v4();
      float recopt = jet_v4.pt();
      float abseta = fabs(jet_v4.eta());

      int ieta(-1);

      for(unsigned int idx=0; idx<JER_SFs_.size(); ++idx){

	const float min_eta = idx ? JER_SFs_.at(idx-1).at(0) : 0.;
	const float max_eta =       JER_SFs_.at(idx)  .at(0);

	if(min_eta <= abseta && abseta < max_eta){ ieta = idx; break; }
      }
      if(ieta < 0) {
	cout << "WARNING: JetSmearerVariable: index for JER-smearing SF not found for jet with |eta| = " << abseta << endl;
	cout << "         no JER smearing is applied." << endl;
	continue;
      }

      float c;
      if(direction == 0){
	c = JER_SFs_.at(ieta).at(1);
      }
      else if(direction == 1){
	c = JER_SFs_.at(ieta).at(2);
      }
      else{
	c = JER_SFs_.at(ieta).at(3);
      }
      float new_pt = std::max(0.0f, genpt + c * (recopt - genpt));
      jet_v4 *= new_pt / recopt;

      //update JEC_factor_raw needed for smearing MET
      float factor_raw = jet.JEC_factor_raw();
      factor_raw *= recopt/new_pt;

      jet.set_JEC_factor_raw(factor_raw);
      jet.set_v4(jet_v4);
    }

    return true;
}
*/
ElectronFakeRateWeights::ElectronFakeRateWeights(Context & ctx, const std::vector<std::string> & JEC_files, TString path_, TString SysDirection_, const string label_jets, const string label_genjets):  path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronFakeRateWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  
  SF.reset((TGraphAsymmErrors*)file->Get("ScaleFactors"));						 
  
  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronFakeRateWeights.process(): Invalid SysDirection specified."); 

  jet_corrector.reset(new JetCorrectorVariable(ctx, JEC_files));

  bool jer_was_applied = ctx.get("meta_jer_applied", "") == "true";
  if(jer_was_applied) ctx.set_metadata("jer_applied", "false", true);
  //jet_smearer.reset(new JetSmearerVariable(ctx, label_jets, label_genjets, JERSmearing::SF_13TeV_2016));
  jet_smearer.reset(new GenericJetResolutionSmearer(ctx, label_jets, label_genjets, true, JERSmearing::SF_13TeV_2016));
  if(!jer_was_applied) ctx.set_metadata("jer_applied", "false", true);

  jet_id = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.5));

  FakeRateWeightEle = ctx.get_handle<double>("FakeRateWeightEle");
  FakeRateWeightEleUp = ctx.get_handle<double>("FakeRateWeightEleUp");
  FakeRateWeightEleDown = ctx.get_handle<double>("FakeRateWeightEleDown");
  h_jets = ctx.get_handle<std::vector<Jet>>(label_jets);

}

bool ElectronFakeRateWeights::process(Event & event){


  vector<Jet> *jets = &event.get(h_jets);
  jet_corrector->correct_collection(event, *jets);
  jet_smearer->process(event);

  if(event.isRealData || event.electrons->size() < 1){
    event.set(FakeRateWeightEle,1.);
    event.set(FakeRateWeightEleUp,1.);
    event.set(FakeRateWeightEleDown,1.);
    return false;
  }


    
  
  //find fake-electrons
  vector<bool> is_fake;




  //if the number of gen and reco-muons are the same, assume muons are real
  unsigned int n_genele = 0;
  for(const auto & gp : *event.genparticles){
    if(fabs(gp.pdgId()) != 11) continue;
    n_genele++;
  }
  //no fakes if number of gen and reco electrons is the same
  if(n_genele == event.electrons->size()){
    event.set(FakeRateWeightEle,1.);
    event.set(FakeRateWeightEleUp,1.);
    event.set(FakeRateWeightEleDown,1.);
    return false;
  } 
  else{
    //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
    unsigned int n_matched_to_ele = 0, n_matched_to_tau = 0;
    for(const auto & ele : *event.electrons){
      bool is_matched = false;
      for(const auto & gp : *event.genparticles){
	if(fabs(gp.pdgId()) == 11){
	  if(deltaR(gp,ele) < 0.1 && !is_matched){
	    is_matched = true;
	    n_matched_to_ele++;
	  }
	}
	else if(fabs(gp.pdgId()) == 15){ 
	  if(deltaR(gp,ele) < 0.2 && !is_matched){
	    is_matched = true;
	    n_matched_to_tau++;
	  }
	}
      }
      is_fake.push_back(!is_matched);
    }
    if(n_matched_to_tau + n_matched_to_ele == event.electrons->size()){
      event.set(FakeRateWeightEle,1.);
      event.set(FakeRateWeightEleUp,1.);
      event.set(FakeRateWeightEleDown,1.);
      return false;
    }
  }

  //find jets matching fake-electrons
  vector<bool> faking_jet;
  for(auto & jet : *jets){

    bool faking = false;
    double dr_min = 999;
    for(unsigned int i=0; i<event.electrons->size(); i++){

      if(!is_fake[i]) continue;
      //consider only jets to be eligible for faking the electron that fulfill the jetId used when deriving the scale factors
      if(!jet_id(jet,event)) continue;

      if(deltaR(event.electrons->at(i), jet) < 0.4) faking = true;
      if(deltaR(event.electrons->at(i), jet) < dr_min) dr_min = deltaR(event.electrons->at(i), jet);
    }
    //cout << "jet passes jetid: " << jet_id(jet,event) <<  endl;
    //cout << "Jet is faking an ele? : " << faking << ", pT: " << jet.pt() << ", dRmin to fake_ele: " << dr_min << endl;
    faking_jet.push_back(faking);
  }

  //apply SF to the eventweight for each faking jet
  double SF_final = 1, SF_final_up = 1, SF_final_down = 1;
  for(unsigned int i=0; i<jets->size(); i++){
    if(!faking_jet[i]) continue;
    int bin = -1;
    if(jets->at(i).pt() < 100) bin = 0;
    else if(jets->at(i).pt() < 200) bin = 1;
    else bin = 2;

    //possibly accout for systematics = statistical, inflate stat unc. by 50% if >800 GeV 
    double x,scale_factor, scale_factor_up, scale_factor_down, mult = 1;
    if(jets->at(i).pt() > 800) mult = 1.5;
    SF->GetPoint(bin,x,scale_factor);
    scale_factor_up   = scale_factor + mult * SF->GetErrorYhigh(bin);
    scale_factor_down = scale_factor - mult * SF->GetErrorYlow(bin);

    if(SysDirection == "up") event.weight *= scale_factor_up;
    else if(SysDirection == "down") event.weight *= scale_factor_down;
    else if(SysDirection == "nominal") event.weight *= scale_factor; 
    else throw runtime_error("In ElectronFakeRateWeights::process(): SysDirection is not one of the following: ['up', 'down', 'nominal']");

    SF_final      *= scale_factor;
    SF_final_up   *= scale_factor_up;
    SF_final_down *= scale_factor_down;

  }

  event.set(FakeRateWeightEle,SF_final);
  event.set(FakeRateWeightEleUp,SF_final_up);
  event.set(FakeRateWeightEleDown,SF_final_down);
    

  return true;
}


MuonFakeRateWeights::MuonFakeRateWeights(Context & ctx, TString path_, TString SysDirection_):  path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: MuonFakeRateWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  SF.reset((TGraphAsymmErrors*)file->Get("ScaleFactors"));						 
  
  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, MuonFakeRateWeights.process(): Invalid SysDirection specified."); 


  FakeRateWeightMu = ctx.get_handle<double>("FakeRateWeightMu");
  FakeRateWeightMuUp = ctx.get_handle<double>("FakeRateWeightMuUp");
  FakeRateWeightMuDown = ctx.get_handle<double>("FakeRateWeightMuDown");

}

bool MuonFakeRateWeights::process(Event & event){

  if(event.isRealData || event.muons->size() < 1){
    event.set(FakeRateWeightMu,1.);
    event.set(FakeRateWeightMuUp,1.);
    event.set(FakeRateWeightMuDown,1.);
    return false;
  }
  
  //find fake-muons
  vector<bool> is_fake;
  unsigned int n_genmu = 0;
  for(const auto & gp : *event.genparticles){
    if(fabs(gp.pdgId()) != 13) continue;
    n_genmu++;
  }

  if(n_genmu == event.muons->size()){
    event.set(FakeRateWeightMu,1.);
    event.set(FakeRateWeightMuUp,1.);
    event.set(FakeRateWeightMuDown,1.);
    return false;
  }
  else{
    //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
    unsigned int n_matched_to_muons = 0, n_matched_to_taus = 0;
    for(const auto & mu : *event.muons){
      bool is_matched = false;
      for(const auto & gp : *event.genparticles){
	if(fabs(gp.pdgId()) == 13){
	  if(deltaR(gp,mu) < 0.1 && !is_matched){
	    is_matched = true;
	    n_matched_to_muons++;
	  }
	}
	else if(fabs(gp.pdgId()) == 15){ 
	  if(deltaR(gp,mu) < 0.2 && !is_matched){
	    is_matched = true;
	    n_matched_to_taus++;
	  }
	}
      }
      is_fake.push_back(!is_matched);
    }
    if(n_matched_to_taus + n_matched_to_muons == event.muons->size()){
    event.set(FakeRateWeightMu,1.);
    event.set(FakeRateWeightMuUp,1.);
    event.set(FakeRateWeightMuDown,1.);
    return false;
    }
  }


  //apply SF to the eventweight for each fake muon
  double SF_final = 1, SF_final_up = 1, SF_final_down = 1;
  for(unsigned int i=0; i<event.muons->size(); i++){
    if(!is_fake[i]) continue;
    int bin = 0;

    double x,scale_factor, scale_factor_up, scale_factor_down;
    SF->GetPoint(bin,x,scale_factor);
    scale_factor_up   = scale_factor + SF->GetErrorYhigh(bin);
    scale_factor_down = scale_factor - SF->GetErrorYlow(bin);

    if(SysDirection == "up") event.weight *= scale_factor_up;
    else if(SysDirection == "down") event.weight *= scale_factor_down;
    else if(SysDirection == "nominal") event.weight *= scale_factor; 
    else throw runtime_error("In MuonFakeRateWeights::process(): SysDirection is not one of the following: ['up', 'down', 'nominal']");

    SF_final      *= scale_factor;
    SF_final_up   *= scale_factor_up;
    SF_final_down *= scale_factor_down;

  }

  event.set(FakeRateWeightMu,SF_final);
  event.set(FakeRateWeightMuUp,SF_final_up);
  event.set(FakeRateWeightMuDown,SF_final_down);


  return true;
}


ZEEFinder::ZEEFinder(){}

pair<int,int> ZEEFinder::search(Event & event){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >= 2 electrons, no pair can be found. Going to return (-1,-1)." << endl; 
      pair<int,int> dummy(-1,-1);
      return dummy;
    }
  }

  const int Nele = event.electrons->size();
  int idx_best_ele1 = -1;
  int idx_best_ele2 = -1;
  double Mee_best = 99999999;
  for(int i=0; i<Nele; i++){
    for(int j=0; j<Nele; j++){
      if(j>i){
	double Mee = (event.electrons->at(i).v4()+event.electrons->at(j).v4()).M();
	if(fabs(Mee - 91.2) < fabs(Mee_best - 91.2)){
	  Mee_best = Mee;
	  idx_best_ele1 = i;
	  idx_best_ele2 = j;
	}
      }
    }
  }
  pair<int,int> dummy(idx_best_ele1,idx_best_ele2);
  return dummy;
}

pair<int,int> ZEEFinder::search(const Event & event){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >= 2 electrons, no pair can be found. Going to return (-1,-1)." << endl; 
      pair<int,int> dummy(-1,-1);
      return dummy;
    }
  }

  const int Nele = event.electrons->size();
  int idx_best_ele1 = -1;
  int idx_best_ele2 = -1;
  double Mee_best = 99999999;
  for(int i=0; i<Nele; i++){
    for(int j=0; j<Nele; j++){
      if(j>i){
	double Mee = (event.electrons->at(i).v4()+event.electrons->at(j).v4()).M();
	if(fabs(Mee - 91.2) < fabs(Mee_best - 91.2)){
	  Mee_best = Mee;
	  idx_best_ele1 = i;
	  idx_best_ele2 = j;
	}
      }
    }
  }
  pair<int,int> dummy(idx_best_ele1,idx_best_ele2);
  return dummy;
}


ElectronJetOverlapCleaner::ElectronJetOverlapCleaner(){}

bool ElectronJetOverlapCleaner::process(Event & event, int idx1, int idx2){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >= 2 electrons, no ElectronJetOverlapCleaning is applied." << endl; 
      return false;
    }
  }

  if(idx1 < 0 || idx2 < 0){
    throw runtime_error("ElectronJetOverlapCleaner was not given two valid electron indices.");
    return false;
  }

  vector<Jet> event_jets_new;
  int Njets = event.jets->size();
  for(int i=0; i<Njets; i++){
    double dR1 = deltaR(event.jets->at(i),event.electrons->at(idx1));
    double dR2 = deltaR(event.jets->at(i),event.electrons->at(idx2));
    if(dR1 > 0.1 && dR2 > 0.1) event_jets_new.push_back(event.jets->at(i));
    //cout << "In ElectronJet-Cleaner: This jet has (dR1, dR2): (" << dR1 << ", " << dR2 << "), it has electron-multiplicity: " << event.jets->at(i).electronMultiplicity() << endl;
  }
  if(fabs(Njets - event_jets_new.size()) > 2) throw runtime_error("In LQToTopMuModules.cxx::ElectronJetOverlapCleaner::process(): When deleting jets overlapping with best electrons, more than 2 were deleted. This must not be.");
  swap(event_jets_new, *event.jets);
  return true;
}

DibosonScaleFactors::DibosonScaleFactors(Context & ctx, TString path_, TString SysDirectionXSec_, TString SysDirectionBTag_) : path(path_), SysDirectionXSec(SysDirectionXSec_), SysDirectionBTag(SysDirectionBTag_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  TString dataset_version = ctx.get("dataset_version");
  if(!dataset_version.Contains("Diboson")){
    is_diboson = false;
    cout << "DibosonScaleFactors will not have an effect on this non-Diboson sample (dataset_version = '" + dataset_version + "')" << endl;
    return;
  }
  else is_diboson = true;
}

bool DibosonScaleFactors::process(Event & event){
  if(event.isRealData) return true;
  if(!is_diboson) return true;
  
  unique_ptr<TFile> file;
  file.reset(new TFile(path,"READ"));

  unique_ptr<TH1D> XSecSF, BTagSF;
  XSecSF.reset((TH1D*)file->Get("Diboson_XSec_SF"));
  BTagSF.reset((TH1D*)file->Get("Diboson_BTag_SF"));
  
  double xsec_SF, btag1_SF, btag2_SF;
  if(SysDirectionXSec == "nominal")   xsec_SF = XSecSF->GetBinContent(1);
  else if(SysDirectionXSec == "up")   xsec_SF = XSecSF->GetBinContent(1) + XSecSF->GetBinError(1);
  else if(SysDirectionXSec == "down") xsec_SF = XSecSF->GetBinContent(1) - XSecSF->GetBinError(1);
  else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionXSec specified.");

  if(SysDirectionBTag == "nominal"){
    btag1_SF = BTagSF->GetBinContent(1);
    btag2_SF = BTagSF->GetBinContent(2);
  }
  else if(SysDirectionBTag == "up"){
    btag1_SF = BTagSF->GetBinContent(1) + BTagSF->GetBinError(1);
    btag2_SF = BTagSF->GetBinContent(2) + BTagSF->GetBinError(2);
  }
  else if(SysDirectionBTag == "down"){
    btag1_SF = BTagSF->GetBinContent(1) - BTagSF->GetBinError(1);
    btag2_SF = BTagSF->GetBinContent(2) - BTagSF->GetBinError(2);
  }
  else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionBTag specified.");

  //First apply XSec SF to all Diboson events before applying 1&2-btag SF
  //cout << "Weight before xsec SF: " << event.weight << ", SF: " << xsec_SF << endl;
  event.weight *= xsec_SF;
  //cout << "Weight after xsec SF: " << event.weight << endl;

  //Count number of loose b-jets in the event
  int n_bjets = 0;
  CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  for (unsigned int i =0; i<event.jets->size(); ++i) {
    if(Btag_loose(event.jets->at(i),event)){
      n_bjets++;
    } 
  }

  //cout << "Number of btags in the event: " << n_bjets << endl;
  //cout << "SF1: " << btag1_SF << ", SF2: " << btag2_SF << endl;
  //cout << "Weight before applying BTag SF: " << event.weight << endl;
  //apply 1&2-btag SF
  if(n_bjets > 2) throw runtime_error("In DibosonScaleFactors::process(): More than 2 b-jets present in the event. Scale factors have only been derived for Nbjets <= 2.");
  else if(n_bjets == 1) event.weight *= btag1_SF;
  else if(n_bjets == 2) event.weight *= btag2_SF; 
  //cout << "Weight after applying BTag SF: " << event.weight << endl;

  return true;
}


MuonTrkWeights::MuonTrkWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: MuonTrkWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  
  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  
  Trk_SF.reset((TGraphAsymmErrors*)file->Get("ratio_eff_eta3_dr030e030_corr"));					 
  
  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
  
}


bool MuonTrkWeights::process(Event & event){

  if(event.isRealData) return true;

  double SF = 1.0;
  for(const auto & mu : *event.muons){
    double eta = mu.eta();
    if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, MuonTrkWeights::process(): Mu-|eta| > 2.4 is not supported at the moment.");
    
    //find right bin in eta
    int idx = 0;
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Trk_SF->GetPoint(idx,x,y);
      keep_going = eta > x + Trk_SF->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
    
    double eff_fnl = 1., dummy_x;
    Trk_SF->GetPoint(idx,dummy_x,eff_fnl);
    
    SF *= eff_fnl;
    if(SysDirection == "up"){
      SF *= 1.005;
    }
    if(SysDirection == "down"){
      SF *= 0.995;
    }
  
  }

  event.weight *= SF;

  return true;

}




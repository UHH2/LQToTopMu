#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"


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

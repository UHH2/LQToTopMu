#include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include <math.h>

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

//!only fill histograms of this class after all nessecary recomodules have processed the event in the analysis cycle and if all requirements for reconstrucions (e.g. >= 1 electron etc.) are met, if not required explicitly before filling!

LQToTopMuRecoHists::LQToTopMuRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  

book<TH1F>("MLQ_HT_Mix", "M_{LQ,mean} & H_{T}", 100, 0, 5000);


//For MLQ reconstruction
 h_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>("HighMassReconstruction");
 m_discriminator_name ="Chi2";

}

void LQToTopMuRecoHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  //HT
  auto met = event.met->pt();
  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    ht += jet.pt();
  }
  //Bedeutung der for-Schleife
  /*const auto jets = event.jets;
    for(unsigned int i=0; i<jets.size();i++){
    Jet jet=jets[i];
    ht +=jet.pt();
    }*/
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
  }
  ht = ht_lep + ht_jets + met;


  //HT / MLQ Mix

  //Fill HT, if Nele = 0, else
  //reconstruct MLQ and fill MLQmean
  int Nele = event.electrons->size();
  if(Nele == 0){hist("MLQ_HT_Mix")->Fill(ht, weight);}
  if(Nele >= 1){   
    std::vector<ReconstructionHypothesis> hyps = event.get(h_hyps);
    const ReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
    
    /*double mttbar_rec = 0;
    if( (hyp->top_v4()+hyp->antitop_v4()).isTimelike() )
      mttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).M();
    else{
      mttbar_rec = sqrt( -(hyp->top_v4()+hyp->antitop_v4()).mass2());
      }*/
    
    /*double mtoplep=0;
    double mtophad=0;
    if(hyp->toplep_v4().isTimelike()) mtoplep = hyp->toplep_v4().M();
    if(hyp->tophad_v4().isTimelike()) mtophad = hyp->tophad_v4().M();*/
    
    //Get Muons and Electrons
    std::vector<Muon>*my_muons = event.muons;
    std::vector<Electron>*my_electrons = event.electrons;
    
    LorentzVector Muon1 = (my_muons->at(0).v4());
    LorentzVector Muon2 = (my_muons->at(1).v4());
    
    double charge_M1 = (my_muons->at(0).charge());
    //double charge_M2 = (my_muons->at(1).charge());
    double charge_E1 = (my_electrons->at(0).charge());
    
    //Combine Top and Muon (Electron and Muon Charge have to be opposite)
    double mLQlep_rec = 0;
    double mLQhad_rec = 0;
    //double mLQmax_rec = 0;
    double mLQmed_rec = 0;

    if(charge_M1 != charge_E1){
      if( (hyp->toplep_v4()+Muon1).isTimelike() )
	mLQlep_rec = (hyp->toplep_v4()+Muon1).M();
      else{
	mLQlep_rec = sqrt( -(hyp->toplep_v4()+Muon1).mass2());
      }
      if( (hyp->tophad_v4()+Muon2).isTimelike() )
	mLQhad_rec = (hyp->tophad_v4()+Muon2).M();
      else{
	mLQhad_rec = sqrt( -(hyp->tophad_v4()+Muon2).mass2());
      }
    }
    else{
      if( (hyp->toplep_v4()+Muon2).isTimelike() )
	mLQlep_rec = (hyp->toplep_v4()+Muon2).M();
      else{
	mLQlep_rec = sqrt( -(hyp->toplep_v4()+Muon2).mass2());
      } 
      if( (hyp->tophad_v4()+Muon1).isTimelike() )
	mLQhad_rec = (hyp->tophad_v4()+Muon1).M();
      else{
	mLQhad_rec = sqrt( -(hyp->tophad_v4()+Muon1).mass2());
      }
    }

    /*if(mLQhad_rec>mLQlep_rec){
      mLQmax_rec = mLQhad_rec;
    }
    else{
      mLQmax_rec = mLQlep_rec;
      }*/

    mLQmed_rec = (mLQhad_rec + mLQlep_rec)/2;
    hist("MLQ_HT_Mix")->Fill(mLQmed_rec, weight);

  }
}




LQToTopMuRecoHists::~LQToTopMuRecoHists(){}




#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;

HypothesisHistsOwn::HypothesisHistsOwn(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  TString name = discriminator_name;
    if(discriminator_name=="Chi2"){
      name = "#Chi^{2}";
    }
    else{
      name += " discriminator";
    }


 
    M_LQlep_rec  = book<TH1F>("M_LQlep_rec", "M_{LQ,lep}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQhad_rec  = book<TH1F>("M_LQhad_rec", "M_{LQ,had}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmax_rec  = book<TH1F>("M_LQmax_rec", "M_{LQmax}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmean_rec = book<TH1F>("M_LQmean_rec", "M_{LQmean}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    //rebinned for theta-analysis
    //double binsxLQmean[3] = {0,350,2000};
    //double binsxLQmean[8] = {0,200,250,300,350,400,600,1200};
    double binsxLQmean[19] = {0,50,100, 150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
    M_LQmean_rec_rebin = book<TH1F>("M_LQmean_rec_rebin", "M_{LQmean}^{rec} [GeV/c^{2}]", 18, binsxLQmean );

    M_ttbar_rec = book<TH1F>( "M_ttbar_rec", "M_{t#bar{t}}^{rec} [GeV/c^{2}]", 100, 0, 5000 ) ;
 
    M_toplep_rec = book<TH1F>( "M_toplep_rec", "M^{top,lep} [GeV/c^{2}]", 70, 0, 700 ) ;
    M_tophad_rec = book<TH1F>( "M_tophad_rec", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;
  
    /*  M_tophad_rec_1jet = book<TH1F>( "M_tophad_rec_1jet", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;
    M_tophad_rec_2jet = book<TH1F>( "M_tophad_rec_2jet", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;
    M_tophad_rec_3jet = book<TH1F>( "M_tophad_rec_3jet", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;
    
    Pt_toplep_rec = book<TH1F>( "Pt_toplep_rec", "P_{T}^{top,lep} [GeV/c]", 60, 0, 1200 ) ;
    Pt_tophad_rec = book<TH1F>( "Pt_tophad_rec", "P_{T}^{top,had} [GeV/c]", 60, 0, 1200 ) ;*/
    
    Pt_ttbar_rec = book<TH1F>( "Pt_ttbar_rec", "P_{T,t#bar{t}}^{rec} [GeV/c]", 60, 0, 600 ) ;
    
    
    h_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>(hyps_name);
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    m_discriminator_name = discriminator_name;
}


void HypothesisHistsOwn::fill(const uhh2::Event & e){



  std::vector<ReconstructionHypothesis> hyps = e.get(h_hyps);
  const ReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
  double weight = e.weight;

   double mttbar_rec = 0;


    if( (hyp->top_v4()+hyp->antitop_v4()).isTimelike() )
    mttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).M();
  else{
    mttbar_rec = sqrt( -(hyp->top_v4()+hyp->antitop_v4()).mass2());
  }
  double ptttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).Pt();

  M_ttbar_rec->Fill(mttbar_rec, weight);
  Pt_ttbar_rec->Fill ( ptttbar_rec, weight);
  
  double mtoplep=0;
  double mtophad=0;

  if(hyp->toplep_v4().isTimelike()) mtoplep = hyp->toplep_v4().M();
  if(hyp->tophad_v4().isTimelike()) mtophad = hyp->tophad_v4().M();
  M_toplep_rec->Fill(mtoplep,weight);
  M_tophad_rec->Fill(mtophad,weight);


  //Get Muons and Electrons
  std::vector<Muon>*my_muons = e.muons;
  std::vector<Electron>*my_electrons = e.electrons;

  LorentzVector Muon1 = (my_muons->at(0).v4());
  LorentzVector Muon2 = (my_muons->at(1).v4());

  double charge_M1 = (my_muons->at(0).charge());
  //double charge_M2 = (my_muons->at(1).charge());
  double charge_E1 = (my_electrons->at(0).charge());
  
  //Combine Top and Muon (Electron and Muon Charge have to be opposite)
  double mLQlep_rec = 0;
  double mLQhad_rec = 0;
  double mLQmax_rec = 0;
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
  
  if(mLQhad_rec>mLQlep_rec){
    mLQmax_rec = mLQhad_rec;
  }
  else{
    mLQmax_rec = mLQlep_rec;
    }
  
  mLQmed_rec = (mLQhad_rec + mLQlep_rec)/2;
  
  M_LQlep_rec->Fill(mLQlep_rec, weight);
  M_LQhad_rec->Fill(mLQhad_rec, weight);
  M_LQmax_rec->Fill(mLQmax_rec, weight);
  M_LQmean_rec->Fill(mLQmed_rec, weight);
  M_LQmean_rec_rebin->Fill(mLQmed_rec, weight);

  
}

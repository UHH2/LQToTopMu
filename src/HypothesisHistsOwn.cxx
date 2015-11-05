#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;
using namespace std;

HypothesisHistsOwn::HypothesisHistsOwn(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  TString name = discriminator_name;
  double min=0;
  double max=500;
  
  if(discriminator_name=="Chi2"){
    name = "#Chi^{2}";
  }
  else{
    name += " discriminator";
  }

  if( discriminator_name=="CorrectMatch"){
    min=0;
    max=2;
  }


    Discriminator = book<TH1F>("Discriminator",name,100,min,max);
    Discriminator_2 = book<TH1F>("Discriminator_2",name,50,0,10);
    Discriminator_3 = book<TH1F>("Discriminator_3",name,300,0,30); 
 
    M_LQlep_rec  = book<TH1F>("M_LQlep_rec", "M_{LQ,lep}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQhad_rec  = book<TH1F>("M_LQhad_rec", "M_{LQ,had}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmax_rec  = book<TH1F>("M_LQmax_rec", "M_{LQmax}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmean_rec = book<TH1F>("M_LQmean_rec", "M_{LQmean}^{rec} [GeV/c^{2}]", 40, 0, 2000 );

    M_LQ_rec_diff = book<TH1F>("M_LQ_rec_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV/c^2]", 50, -500, 500);
    M_LQ_rec_diff_rel = book<TH1F>("M_LQ_rec_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep}) / M_{LQ}^{mean} [GeV/c^]", 50, -0.5, 0.5);

    //rebinned for theta-analysis
    double binsxLQmean[19] = {0,50,100, 150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
    M_LQmean_rec_rebin = book<TH1F>("M_LQmean_rec_rebin", "M_{LQmean}^{rec} [GeV/c^{2}]", 18, binsxLQmean );

    M_ttbar_rec = book<TH1F>( "M_ttbar_rec", "M_{t#bar{t}}^{rec} [GeV/c^{2}]", 100, 0, 5000 ) ;
 
    M_toplep_rec = book<TH1F>( "M_toplep_rec", "M^{top,lep} [GeV/c^{2}]", 70, 0, 700 ) ;
    M_tophad_rec = book<TH1F>( "M_tophad_rec", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;

    
    Pt_ttbar_rec = book<TH1F>( "Pt_ttbar_rec", "P_{T,t#bar{t}}^{rec} [GeV/c]", 60, 0, 600 ) ;
    
    
    h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>(hyps_name);
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    m_discriminator_name = discriminator_name;
}


void HypothesisHistsOwn::fill(const uhh2::Event & e){



  std::vector<LQReconstructionHypothesis> hyps = e.get(h_hyps);
  const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
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
  else mtoplep = sqrt( -(hyp->toplep_v4()).mass2());
  if(hyp->tophad_v4().isTimelike()) mtophad = hyp->tophad_v4().M();
  else mtophad = sqrt( -(hyp->tophad_v4()).mass2());
  M_toplep_rec->Fill(mtoplep,weight);
  M_tophad_rec->Fill(mtophad,weight);

  /*cout << "in hyphists: MTopLep: " << mtoplep << ", MTopHad: " << mtophad << endl;
  cout << "Timelike: " << hyp->toplep_v4().isTimelike() << endl;
  cout << "2 Masses: " << hyp->toplep_v4().M() << " and " << sqrt( -(hyp->toplep_v4()).mass2()) << endl;
  cout << "4 vec:    " << hyp->toplep_v4() << endl;*/

  Discriminator->Fill(hyp->discriminator(m_discriminator_name) ,weight);
  Discriminator_2->Fill(hyp->discriminator(m_discriminator_name) ,weight); 
  Discriminator_3->Fill(hyp->discriminator(m_discriminator_name) ,weight);

  
  //Combine Top and Muon (Electron and Muon Charge have to be opposite)
  double mLQlep_rec = 0;
  double mLQhad_rec = 0;
  double mLQmax_rec = 0;
  double mLQmed_rec = 0;
  double mLQ_rec_diff = 0;
  double mLQ_rec_diff_rel = 0;

  if( (hyp->LQlep_v4()).isTimelike() ) {mLQlep_rec = (hyp->LQlep_v4()).M();}
  else {mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());}
  if( (hyp->LQhad_v4()).isTimelike() ) {mLQhad_rec = (hyp->LQhad_v4()).M();}
  else {mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());}
  
  if(mLQhad_rec>mLQlep_rec){
    mLQmax_rec = mLQhad_rec;
  }
  else{
    mLQmax_rec = mLQlep_rec;
    }
  
  mLQmed_rec = (mLQhad_rec + mLQlep_rec)/2;
  mLQ_rec_diff = mLQhad_rec - mLQlep_rec;
  mLQ_rec_diff_rel = mLQ_rec_diff/mLQmed_rec;
  //cout << "in hyphists: MLQmean: " << mLQmed_rec << ", mLQdiff: " << mLQ_rec_diff << endl << endl;
  
  M_LQlep_rec->Fill(mLQlep_rec, weight);
  M_LQhad_rec->Fill(mLQhad_rec, weight);
  M_LQmax_rec->Fill(mLQmax_rec, weight);
  M_LQmean_rec->Fill(mLQmed_rec, weight);
  M_LQmean_rec_rebin->Fill(mLQmed_rec, weight);
  if(hyp->discriminator(m_discriminator_name) < 20){ // 999999 is set as discr. value on CorrectMatchDiscr, if one of the required conditions is not matched
    M_LQ_rec_diff->Fill(mLQ_rec_diff, weight);
    M_LQ_rec_diff_rel->Fill(mLQ_rec_diff_rel,weight);
  }
  
}

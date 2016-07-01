#include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <sstream>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;


LQToTopMuPDFHists::LQToTopMuPDFHists(Context & ctx, const string & dirname, bool use_pdf_weights_): Hists(ctx, dirname), use_pdf_weights(use_pdf_weights_){  

is_mc = ctx.get("dataset_type") == "MC";
 //For MLQ reconstruction
  h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
  m_discriminator_name ="Chi2";

  double bins_low_NoEle2[11] = {0,350,500,650,800,950,1100,1250,1450,1750,2050};
  double bins_mlq_low2[6] = {0,200,400,600,800,1000};

  for(int i=0; i<100; i++){
    stringstream ss_name;
    ss_name << "H_T_rebin3_PDF_"  << i+1 ;
    stringstream ss_name2;
    ss_name2 << "H_T_from350_rebin_PDF_"  << i+1 ;
    stringstream ss_name3;
    ss_name3 << "M_LQ_comb_rebin2_PDF_"  << i+1 ;
    stringstream ss_name4;
    ss_name4 << "H_T_comb_NoEle_rebin2_PDF_"  << i+1 ;
    stringstream ss_title;
    ss_title << "H_{T} [GeV] for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title2;
    ss_title2 << "H_{T} [GeV] (from 350) for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title3;
    ss_title3 << "M_{LQ} [GeV] for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title4;
    ss_title4 << "H_{T} [GeV] No Ele for PDF No. "  << i+1 << " out of 100" ;

    string s_name = ss_name.str();
    string s_title = ss_title.str();
    string s_name2 = ss_name2.str();
    string s_title2 = ss_title2.str();
    string s_name3 = ss_name3.str();
    string s_title3 = ss_title3.str();
    string s_name4 = ss_name4.str();
    string s_title4 = ss_title4.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();
    const char* char_name2 = s_name2.c_str();
    const char* char_title2 = s_title2.c_str();
    const char* char_name3 = s_name3.c_str();
    const char* char_title3 = s_title3.c_str();
    const char* char_name4 = s_name4.c_str();
    const char* char_title4 = s_title4.c_str();
    histo_names[i] = s_name;
    histo_names2[i] = s_name2;
    histo_names3[i] = s_name3;
    histo_names4[i] = s_name4;

    book<TH1F>(char_name, char_title, 10,bins_low_NoEle2);
    book<TH1F>(char_name2, char_title2, 48, 0,4200);
    book<TH1F>(char_name3, char_title3,5, bins_mlq_low2);
    book<TH1F>(char_name4, char_title4, 10,bins_low_NoEle2);
  }

}

void LQToTopMuPDFHists::fill(const Event & event){
  double weight = event.weight;

  if(is_mc){
    //HT
    auto met = event.met->pt();
    double ht = 0.0;
    double ht_jets = 0.0;
    double ht_lep = 0.0;
    for(const auto & jet : *event.jets){
      ht_jets += jet.pt();
    }
    for(const auto & electron : *event.electrons){
      ht_lep += electron.pt();
    }
    for(const auto & muon : *event.muons){
      ht_lep += muon.pt();
    }
    
    ht = ht_lep + ht_jets + met;
    
    if(event.genInfo->systweights().size()){
      for(int i=0; i<100; i++){
	if(use_pdf_weights){
	  double pdf_weight = event.genInfo->systweights().at(i+9);
	  double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	  const char* name = histo_names[i].c_str();
	  const char* name2 = histo_names2[i].c_str();
	  const char* name4 = histo_names4[i].c_str();
	  if(ht <= 2000) hist(name)->Fill(ht,fillweight);
	  else hist(name)->Fill(2000,fillweight);
	  hist(name2)->Fill(ht,fillweight);
	  if(event.electrons->size() == 0){
	    if(ht <= 2000) hist(name4)->Fill(ht,fillweight);
	    else hist(name4)->Fill(2000,fillweight);
	  }
	}
      }
    }

    bool charge_opposite = false;
    for(unsigned int i=0; i<event.muons->size(); i++){
      for(unsigned int j=0; j<event.muons->size(); j++){
	if(j>i){
	  if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
	    charge_opposite = true;
	  }
	}
      }
    }
    if(event.electrons->size() >= 1 && event.muons->size() >= 2 && charge_opposite){   
      std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps); 
      const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );

      double mLQlep_rec = 0;
      double mLQhad_rec = 0;
      double mLQmed_rec = 0;


      if( (hyp->LQlep_v4()).isTimelike() ) mLQlep_rec = (hyp->LQlep_v4()).M();
      else mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());
      if( (hyp->LQhad_v4()).isTimelike() ) mLQhad_rec = (hyp->LQhad_v4()).M();
      else mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());

      mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;

      if(event.genInfo->systweights().size()){
	for(int i=0; i<100; i++){
	  if(use_pdf_weights){
	    double pdf_weight = event.genInfo->systweights().at(i+9);
	    double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	    const char* name3 = histo_names3[i].c_str();
	    hist(name3)->Fill(mLQmed_rec,fillweight);
	  }
	}
      } //systweights filled
    } //LQ reconstructable

  } //is_mc
}

LQToTopMuPDFHists::~LQToTopMuPDFHists(){}















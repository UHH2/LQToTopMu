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

  double bins_low_NoEle2[11] = {0,350,500,650,800,950,1100,1250,1450,1750,2050};

  for(int i=0; i<102; i++){
    stringstream ss_name;
    ss_name << "H_T_rebin3_PDF_"  << i+1 ;
    stringstream ss_name2;
    ss_name2 << "H_T_from350_rebin_PDF_"  << i+1 ;
    stringstream ss_title;
    ss_title << "H_{T} [GeV] for PDF No. "  << i+1 << " out of 102" ;
    stringstream ss_title2;
    ss_title2 << "H_{T} [GeV] (from 350) for PDF No. "  << i+1 << " out of 102" ;

    string s_name = ss_name.str();
    string s_title = ss_title.str();
    string s_name2 = ss_name2.str();
    string s_title2 = ss_title2.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();
    const char* char_name2 = s_name2.c_str();
    const char* char_title2 = s_title2.c_str();
    histo_names[i] = s_name;
    histo_names2[i] = s_name2;

    book<TH1F>(char_name, char_title, 10,bins_low_NoEle2);
    book<TH1F>(char_name2, char_title2, 80, 0,7000);
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
      for(int i=0; i<102; i++){
	if(use_pdf_weights){
	  double pdf_weight = event.genInfo->systweights().at(i+9);
	  double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	  const char* name = histo_names[i].c_str();
	  const char* name2 = histo_names2[i].c_str();
	  if(ht <= 2000) hist(name)->Fill(ht,fillweight);
	  else hist(name)->Fill(2000,fillweight);
	  hist(name2)->Fill(ht,fillweight);
	}
      }
    }
    //else cout << "systweights not filled for this sample." << endl;
  }
}

LQToTopMuPDFHists::~LQToTopMuPDFHists(){}















#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/core/include/Event.h"
#include <math.h>

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuHists::LQToTopMuHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);  
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);
  book<TH1F>("pt_jet1", "#p_{T}^{jet 1}", 100, 0, 3000);
  book<TH1F>("pt_jet2", "#p_{T}^{jet 2}", 100, 0, 3000);
  book<TH1F>("pt_jet3", "#p_{T}^{jet 3}", 100, 0, 3000);
  book<TH1F>("N_bJets_loose", "#N_{Bjets}^{loose}", 10, 0, 10);
  book<TH1F>("N_bJets_med", "#N_{Bjets}^{medium}", 10, 0, 10);
  book<TH1F>("N_bJets_tight", "#N_{Bjets}^{tight}", 10, 0, 10);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1F>("M_mumu", "M_(#mu#mu)",50 , 0, 1000);
  double bins_pt_low[13] = {0,30,60,90,120,150,180,210,240,270,300,350,800};
  book<TH1F>("Pt_mu1", "P_{T}^{leading #mu}",12 ,bins_pt_low);

  // general
  book<TH1F>("N_pv", "N_{PV}", 50, 0, 50);
  book<TH1F>("H_T", "H_{T}", 100, 0, 3000);
  double bins_low[9] = {0,350,500,700,900,1100,1300,1500,3000};
  book<TH1F>("H_T_rebin", "H_{T}", 8, bins_low);

  //book <TH1F>("M_t_had", "M_{t,had}", 50, 0, 500);
  //book <TH1F>("M_t_lep", "M_{t,lep}", 70, 0, 700);
  //book <TH1F>("M_ttbar", "M_{t#bar{t}}", 100, 0, 5000);
}


void LQToTopMuHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);
  
  if(Njets>=1){
    hist("eta_jet1")->Fill(jets->at(0).eta(), weight);
  }
  if(Njets>=2){
    hist("eta_jet2")->Fill(jets->at(1).eta(), weight);
  }
  if(Njets>=3){
    hist("eta_jet3")->Fill(jets->at(2).eta(), weight);
  }
  if(Njets>=4){
    hist("eta_jet4")->Fill(jets->at(3).eta(), weight);
  }
  if(Njets>=1){
    hist("pt_jet1")->Fill(jets->at(0).pt(), weight);
  }
  if(Njets>=2){
    hist("pt_jet2")->Fill(jets->at(1).pt(), weight);
  }
  if(Njets>=3){
    hist("pt_jet3")->Fill(jets->at(2).pt(), weight);
  }

  //# b-jets
  std::vector<Jet> bjets_loose, bjets_med, bjets_tight;
  for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.423) { //loose: >0.423, medium: >0.814, tight: >0.914
      bjets_loose.push_back(jets->at(i));
    }
    if(jets->at(i).btag_combinedSecondaryVertex()>0.814) { //loose: >0.423, medium: >0.814, tight: >0.914
      bjets_med.push_back(jets->at(i));
    }
    if(jets->at(i).btag_combinedSecondaryVertex()>0.914) { //loose: >0.423, medium: >0.814, tight: >0.914
      bjets_tight.push_back(jets->at(i));
    }
  }
  int NbJets_loose = bjets_loose.size();
  hist("N_bJets_loose")->Fill(NbJets_loose,weight);
  int NbJets_med = bjets_med.size();
  hist("N_bJets_med")->Fill(NbJets_med,weight);
  int NbJets_tight = bjets_tight.size();
  hist("N_bJets_tight")->Fill(NbJets_tight,weight);


  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      hist("pt_mu")->Fill(thismu.pt(), weight);
      hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }
  hist("Pt_mu1")->Fill(event.muons->at(0).pt(), weight);
  
  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

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
  hist("H_T")->Fill(ht, weight);
  hist("H_T_rebin")->Fill(ht, weight);

  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j > i){
	M_mumu = (muons[i] + muons[j]).M();
	hist("M_mumu")->Fill(M_mumu, weight);	
      }
    }
  }
  

  //M_LQLQbar
  /*  
  int Nele = event.electrons->size();
  if(Nele >= 1){


 //M_top,lep reconstruction
 
  //step 1: calculate pz_neutrino following B2G-12-006 AN
    std::vector<Electron>* electrons = event.electrons;
    LorentzVector Electron = electrons->at(0).v4();
    LorentzVectorXYZE ElectronXYZE = toXYZ(Electron);
    
    double mu = ((80.385*80.385/2) + Electron.pt() * met * cos(Electron.phi() - event.met->phi()));
    double A = - (Electron.pt() * Electron.pt() );
    double B = mu * ElectronXYZE.pz();
    double C = mu * mu - Electron.energy() * Electron.energy() * met * met;
    double discriminant = B * B - A * C;

    double p_z1_neutrino, p_z2_neutrino;
    if( discriminant < 0){
      p_z1_neutrino = -B / A;}
    else{
      p_z1_neutrino = (-B - sqrt(discriminant)) / A;
      p_z2_neutrino = (-B + sqrt(discriminant)) / A;
    }
    
    //step 2: calculate eta = y (m_neutrino = 0)
    double energy1_neutrino = sqrt(pow(met,2)+pow(p_z1_neutrino,2));
    double energy2_neutrino = sqrt(pow(met,2)+pow(p_z2_neutrino,2));
    double eta1_neutrino = 0;
    double eta2_neutrino = 0;
    eta1_neutrino = 0.5 * log((energy1_neutrino + p_z1_neutrino)/(energy1_neutrino - p_z1_neutrino));   
    if( discriminant >=0){ 
      eta2_neutrino = 0.5 * log((energy2_neutrino + p_z2_neutrino)/(energy2_neutrino - p_z2_neutrino));
    }
    
    //step 3: build LorentzVectors
    LorentzVector Neutrino1 = {met, eta1_neutrino, event.met->phi(), energy1_neutrino};
    LorentzVector Neutrino2 = {0,0,0,0};
    if(discriminant >= 0 ){
      Neutrino2 = {met, eta2_neutrino, event.met->phi(), energy2_neutrino};
    }
    
    //--now bring the hadronic part into the game
    //Anzahl moeglicher hadr.  Kombinationen herausfinden - fuer Arraygroesse
    int Nkomb = 0;

    //Top aus 1 Jet rekonstruiert (N ueber 1 = N)
    if(Njets >= 2){
    for(int i=0; i<Njets; i++){Nkomb++;}}
    //Top aus 2 Jets rekonstruiert - +(N ueber 2)
    if(Njets >= 3){
      for(int i=0; i<Njets-1; i++){
	for(int j=0; j<Njets; j++){
	    if(j > i){
	      Nkomb++;}}}}
    //Top aus 3 Jets rekonstruiert - +(N ueber 3)
    if(Njets >=4){
      for(int i=0; i<Njets-2; i++){
	for(int j=0; j<Njets-1; j++){
	  for(int k=0; k<Njets; k++){
	    if(k>j && j>i){
	      Nkomb++;}}}}}
    //Top aus 4 Jets rekonstruiert - +(N ueber 4)
    if(Njets >=5){
      for(int i=0; i<Njets-3; i++){
	for(int j=0; j<Njets-2; j++){
	  for(int k=0; k<Njets-1; k++){
	    for(int l=0;l<Njets; l++){
	      if(l>k && k>j && j>i){
		Nkomb++;}}}}}}

    Nkomb = 2*Nkomb; // Neutrino kann 2 Loesungen haben-->doppelt so viele  Kombinationen

    //step 4: declare variables and set up arrays
    double chi2Hadronic[Nkomb*Njets]; //Nkomb * Njets: # of had Hyps per lept Hyp * # of lept Hyps
    for(int i=0; i<Nkomb*Njets; i++){chi2Hadronic[i] = 1000000;}
    double residualHadronic[Nkomb*Njets];
    for(int i=0; i<Nkomb*Njets; i++){residualHadronic[i] = 1000000;}
    LorentzVector jet[Njets];
    for(int i=0; i<Njets; i++){jet[i] = jets->at(i).v4();}

    double residualLeptonic[2*Njets];
    for(int i=0; i<2*Njets; i++){residualLeptonic[i] = 1000000;}
    double chi2Leptonic[2*Njets];
    for(int i=0; i<2*Njets; i++){chi2Leptonic[i] = 1000000;}
 

    double chi2ges[Nkomb*Njets][2*Njets];
    for(int i=0; i<Nkomb*Njets; i++){
      for(int j=0; j<2*Njets; j++){chi2ges[i][j] = 1000000;}}


    LorentzVector HypTopLep[2*Njets];
    for(int i=0; i<2*Njets; i++){
      HypTopLep[i] = {0,0,0,0};
    }
    LorentzVector HypTopHad[Nkomb*Njets];
    for(int i=0; i<Nkomb*Njets; i++){
      HypTopHad[i] = {0,0,0,0};
      }
    
        Nkomb = 0;

    //1'st real solution - calculate all possible combinations of jets assigned to top_lep and top_had and the corresponding chi2ges
	for (int jetlep=0; jetlep<Njets; jetlep++){ //loop over Jets for leptonic Hyp.
      residualLeptonic[jetlep] = (Neutrino1 + Electron + jet[jetlep]).M() - 174;
      chi2Leptonic[jetlep] = (residualLeptonic[jetlep]/18)*(residualLeptonic[jetlep]/18);

      //loop over combinations of Jets for hadronic Hyps
      //Top aus 1 Jet
      if(Njets >= 2){
	for(int i=0; i<Njets; i++){
	  if(i != jetlep){
	    residualHadronic[Nkomb] = jet[i].M() - 181;
	    chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
	    chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
	     HypTopLep[jetlep] = Neutrino1 + Electron + jet[jetlep];
	     HypTopHad[Nkomb] = jet[i];}
	  Nkomb++;}}

      //Top aus 2 Jets
      if(Njets >= 3){
	for(int i=0; i<Njets-1; i++){
	  for(int j=0; j<Njets; j++){
	    if(j > i){
	      if(i != jetlep && j != jetlep){
		residualHadronic[Nkomb] = (jet[i] + jet[j]).M() - 181;
		chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
		chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
		HypTopLep[jetlep] = Neutrino1 + Electron + jet[jetlep];
		HypTopHad[Nkomb] = jet[i] + jet[j];}
	      Nkomb++;}}}}

      //Top aus 3 Jets
      if(Njets >=4){
	for(int i=0; i<Njets-2; i++){
	  for(int j=0; j<Njets-1; j++){
	    for(int k=0; k<Njets; k++){
	      if(k>j && j>i){
		if(i!=jetlep && j!=jetlep && k!=jetlep){
		  residualHadronic[Nkomb] = (jet[i] + jet[j] + jet[k]).M() - 181;
		  chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
		  chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
		  HypTopLep[jetlep] = Neutrino1 + Electron + jet[jetlep];
		  HypTopHad[Nkomb] = jet[i] + jet[j] + jet[k];}
		Nkomb++;}}}}}

      //Top aus 4 Jets
      if(Njets >=5){
	for(int i=0; i<Njets-3; i++){
	  for(int j=0; j<Njets-2; j++){
	    for(int k=0; k<Njets-1; k++){
	      for(int l=0;l<Njets; l++){
		if(l>k && k>j && j>i){
		  if(i!=jetlep && j!=jetlep && k!=jetlep && l!=jetlep){
		    residualHadronic[Nkomb] = (jet[i] + jet[j] + jet[k] + jet[l]).M() - 181;
		    chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
		    chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
		    HypTopLep[jetlep] = Neutrino1 + Electron + jet[jetlep];
		    HypTopHad[Nkomb] = jet[i] + jet[j] + jet[k] + jet[l];}
		  Nkomb++;}}}}}}

    }
    
    //2'nd real solution - calculate all possible combinations of jets assigned to top_lep and top_had
    if(0 != p_z2_neutrino){	
      for (int jetlep=Njets; jetlep<(2*Njets); jetlep++){ //same as above, but only if there are 2 solutions
	residualLeptonic[jetlep] = (Neutrino2 + Electron + jet[jetlep-Njets]).M() - 174;
	chi2Leptonic[jetlep] = (residualLeptonic[jetlep]/18)*(residualLeptonic[jetlep]/18);

	//Top aus 1 Jet
	if(Njets >= 2){
	  for(int i=0; i<Njets; i++){
	    if(i != (jetlep-Njets)){
	      residualHadronic[Nkomb] = jet[i].M() - 181;
	      chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
	      chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
	      HypTopLep[jetlep] = Neutrino2 + Electron + jet[jetlep-Njets];
	      HypTopHad[Nkomb] = jet[i];}
	    Nkomb++;}}
	
	//Top aus 2 Jets
	if(Njets >= 3){
	  for(int i=0; i<Njets-1; i++){
	    for(int j=0; j<Njets; j++){
	      if(j > i){
		if(i != (jetlep-Njets) && j != (jetlep-Njets)){
		  residualHadronic[Nkomb] = (jet[i] + jet[j]).M() - 181;
		  chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
		  chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
		  HypTopLep[jetlep] = Neutrino2 + Electron + jet[jetlep-Njets];
		  HypTopHad[Nkomb] = jet[i] + jet[j];}
		Nkomb++;}}}}
	
	//Top aus 3 Jets
	if(Njets >=4){
	  for(int i=0; i<Njets-2; i++){
	    for(int j=0; j<Njets-1; j++){
	      for(int k=0; k<Njets; k++){
		if(k>j && j>i){
		  if(i!=(jetlep-Njets) && j!=(jetlep-Njets) && k!=(jetlep-Njets)){
		    residualHadronic[Nkomb] = (jet[i] + jet[j] + jet[k]).M() - 181;
		    chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
		    chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
		    HypTopLep[jetlep] = Neutrino2 + Electron + jet[jetlep-Njets];
		    HypTopHad[Nkomb] = jet[i] + jet[j] + jet[k];}
		  Nkomb++;}}}}}
	
	//Top aus 4 Jets
	if(Njets >=5){
	  for(int i=0; i<Njets-3; i++){
	    for(int j=0; j<Njets-2; j++){
	      for(int k=0; k<Njets-1; k++){
		for(int l=0;l<Njets; l++){
		  if(l>k && k>j && j>i){
		    if(i!=(jetlep-Njets) && j!=(jetlep-Njets) && k!=(jetlep-Njets) && l!=(jetlep-Njets)){
		      residualHadronic[Nkomb] = (jet[i] + jet[j] + jet[k] + jet[l]).M() - 181;
		      chi2Hadronic[Nkomb] = (residualHadronic[Nkomb]/15)*(residualHadronic[Nkomb]/15);
		      chi2ges[Nkomb][jetlep] = chi2Hadronic[Nkomb] + chi2Leptonic[jetlep];
		      HypTopLep[jetlep] = Neutrino2 + Electron + jet[jetlep-Njets];
		      HypTopHad[Nkomb] = jet[i] + jet[j] + jet[k] + jet[l];}
		    Nkomb++;}}}}}}
      }
    }
 
    //find minimal chi2_ges

    double chi2gesmin = 1000000;
    int FinalKombHad = 0;
    int FinalKombLep = 0;
  
   
    for(int i=0;i<Nkomb ;i++){// # Jetkombinationen fuer hadr. Top (= "frueheres Nkomb * Njets")
      for(int j=0; j<2*Njets; j++){//Jet und Neutrinokombinationen fuer lept. Top
	if (chi2ges[i][j] < chi2gesmin){
	  chi2gesmin = chi2ges[i][j];
	  FinalKombHad = i;
	  FinalKombLep = j;
	}
      }
    }

    double TopmassRecoLeptonic = residualLeptonic[FinalKombLep] + 174;  
    //cout << "chi2Leptonic minimal value: " <<  chi2Leptonic[FinalKombLep] << endl;
    //cout << "# of combination: " << FinalKombLep << endl;
    //cout << "reconstructed Topmass leptonic: " << TopmassRecoLeptonic << endl;   
    //cout << "gleiches Ergebnis?: " << HypTopLep[FinalKombLep].M() << endl;
    hist("M_t_lep")->Fill(TopmassRecoLeptonic, weight); //only for testing

    double TopmassRecoHadronic = residualHadronic[FinalKombHad] + 181;
    //cout << "chi2Hadronic minimal value: " <<  chi2Hadronic[FinalKombHad] << endl;
    //cout << "# of combination: " << FinalKombHad << endl;
    //cout << "reconstructed Topmass hadronic: " << TopmassRecoHadronic << endl;
    //cout << "gleiches Ergebnis?: " << HypTopHad[FinalKombHad].M() << endl;
    hist("M_t_had")->Fill(TopmassRecoHadronic, weight);

    //cout << "chi2_ges minimal value: " << chi2gesmin << endl << endl;

    double Mttbar_reco = (HypTopLep[FinalKombLep] + HypTopHad[FinalKombHad]).M();
    //cout << "reconstructed Mttbar: " << Mttbar_reco << endl;
    hist("M_ttbar")->Fill(Mttbar_reco, weight);

  } //Nele = 1
    */
} //Methode

LQToTopMuHists::~LQToTopMuHists(){}

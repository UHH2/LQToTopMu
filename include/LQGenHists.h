#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/LQToTopMu/include/LQGen.h"

/** \brief Histograms for LQLQbar quantities on generator (parton) level
 * 
 * LQGen container has to be filled before calling this histogram class
 */
class LQGenHists: public uhh2::Hists {
public:
    LQGenHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F* MLQLQbar_gen, *Pt_LQLQbar_gen, *shat, *M_LQ, *M_antiLQ, *Pt_LQ, *Pt_antiLQ, *Pt_LQ_over_shat, *Pt_antiLQ_over_shat, *Pt_LQ_over_M_LQLQbar, *Pt_antiLQ_over_M_LQLQbar, *eta_LQ, *eta_antiLQ, *y_LQ, *y_antiLQ, *phi_LQ, *phi_antiLQ, *diffabseta, *diffabsy, *deltaR_LQ_decays, *deltaR_antiLQ_decays, *deltaR_top_decays, *deltaR_antitop_decays, *Pt_mu, *Pt_antimu, *eta_mu, *eta_antimu, *y_mu, *y_antimu, *phi_mu, *phi_antimu, *M_mu, *M_antimu, *Pt_top, *Pt_antitop, *eta_top, *eta_antitop, *y_top, *y_antitop, *M_top, *M_antitop, *phi_top, *phi_antitop
      //, *cosThetastar_LQ_LQLQframe, *cosThetastar_antiLQ_LQLQframe, *Pt_LQ_LQLQframe, *Pt_antiLQ_LQLQframe, *DecayChannel
;

    TH2F* M_LQLQbar_vs_shat, *M_LQLQbar_vs_deltaR_LQ, *M_LQLQbar_vs_deltaR_antiLQ, *shat_vs_deltaR_LQ, *shat_vs_deltaR_antiLQ, *Pt_LQ_vs_deltaR_LQ, *Pt_antiLQ_vs_deltaR_antiLQ, *M_LQLQbar_vs_deltaR_top, *M_LQLQbar_vs_deltaR_antitop, *M_LQLQbar_vs_Pt_LQ, *M_LQLQbar_vs_Pt_antiLQ, *shat_vs_Pt_LQ, *shat_vs_Pt_antiLQ, *Pt_LQ_vs_Pt_antiLQ,/* *M_LQLQbar_vs_Pt_LQ_LQLQframe, *M_LQLQbar_vs_Pt_antiLQ_LQLQframe,*/ *M_LQLQbar_vs_eta_LQ, *M_LQLQbar_vs_eta_antiLQ;
    uhh2::Event::Handle<LQGen> h_LQLQbargen;
};

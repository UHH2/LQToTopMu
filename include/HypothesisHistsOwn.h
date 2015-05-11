#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/TTbarGen.h"

/** \brief Common histograms for reconstruction hypotheses
 *
 * hyps_name is the name of the reconstruction hypothesis collection, for instance "HighMassReconstruction"
 * discriminator_name is the name of the discriminator used to choose the best reconstruction hypothesis, for instance "Chi2"
 */
class HypothesisHistsOwn: public uhh2::Hists {
public:
    HypothesisHistsOwn(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F *M_LQlep_rec, *M_LQhad_rec, *M_LQmax_rec, *M_LQmean_rec, *M_LQmean_rec_rebin;
    TH1F *M_ttbar_rec, *M_toplep_rec, *M_tophad_rec, *M_tophad_rec_1jet, *M_tophad_rec_2jet, *M_tophad_rec_3jet;
    TH1F *Pt_toplep_rec, *Pt_tophad_rec, *Pt_ttbar_rec;

    uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    std::string m_discriminator_name;
    
};

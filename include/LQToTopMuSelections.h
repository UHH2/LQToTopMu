#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/GenParticle.h"

#include <vector>

namespace uhh2examples {
    
/* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
 * below 20% of the average of the leading two jets, where the minimum deltaphi and
 * maximum third jet pt fraction can be changed in the constructor.
 * The jets are assumed to be sorted in pt.
 */


  
  class HtSelection: public uhh2::Selection {
  public:
    explicit HtSelection(double ht_min=0., double ht_max=-1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };


  class InvMass2MuVeto: public uhh2::Selection {
  public:
    explicit InvMass2MuVeto(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };


  class InvMass2MuVeto_Inverted: public uhh2::Selection {
  public:
    explicit InvMass2MuVeto_Inverted(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };


  class InvMass2EleVeto_Inverted: public uhh2::Selection {
  public:
    explicit InvMass2EleVeto_Inverted(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };


  class HTLeptSelection : public uhh2::Selection{
  public:
    explicit HTLeptSelection(double ht_min = 0., double ht_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  
  class METSelection : public uhh2::Selection{
  public:
    explicit METSelection(double MET_min = 0., double MET_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double MET_min, MET_max;
    };


   class dRSelection : public uhh2::Selection{
  public:
    explicit dRSelection(double dR_min = 0., double dR_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double dR_min, dR_max;
  };

}

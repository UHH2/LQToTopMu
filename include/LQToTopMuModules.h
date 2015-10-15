#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"

class JetLeptonOverlapCleaner: public uhh2::AnalysisModule {

 public:
  explicit JetLeptonOverlapCleaner(double RJet = 0.4);
  virtual bool process(uhh2::Event & event) override;

 private:
  double RJet;

};

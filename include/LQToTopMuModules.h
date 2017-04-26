#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>

class JetLeptonOverlapCleaner: public uhh2::AnalysisModule {

 public:
  explicit JetLeptonOverlapCleaner(double RJet = 0.4);
  virtual bool process(uhh2::Event & event) override;

 private:
  double RJet;

};



class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirection;

};

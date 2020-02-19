//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// Generator includes
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

// Reweight includes
#include "RwCalculators/GReWeightNuXSecCOHuB.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCOHuB::GReWeightNuXSecCOHuB() : GReWeightModel("COHuB")
{
  this->Reset();
}
//_______________________________________________________________________________________
GReWeightNuXSecCOHuB::~GReWeightNuXSecCOHuB()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecCOHuB::IsHandled(GSyst_t syst) const
{
  if ( syst == kXSecTwkDial_NormCCCOH
    || syst == kXSecTwkDial_NormNCCOH )
  {
    return true;
  }

  return false;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCOHuB::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
  if ( type == kScCoherent ) return true;
  else return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOHuB::SetSystematic(GSyst_t syst, double twk_dial)
{
  if ( !this->IsHandled(syst) ) return;

  switch(syst) {
    case ( kXSecTwkDial_NormCCCOH ) :
      fCCNormTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_NormNCCOH ) :
      fNCNormTwkDial = twk_dial;
      break;
    default:
      return;
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOHuB::Reset(void)
{
  fCCNormTwkDial = 0.;
  fNCNormTwkDial = 0.;

  fCurNormCC = 1.;
  fCurNormNC = 1.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOHuB::Reconfigure(void)
{
  GSystUncertainty* gsu = GSystUncertainty::Instance();

  // Note: this assumes that the error is symmetric.
  // TODO: consider changing this to handle asymmetric errors on the normalization
  double frac_err_cc_norm = gsu->OneSigmaErr( kXSecTwkDial_NormCCCOH );
  double frac_err_nc_norm = gsu->OneSigmaErr( kXSecTwkDial_NormNCCOH );

  fCurNormCC = std::max(0., 1. + fCCNormTwkDial * frac_err_cc_norm);
  fCurNormNC = std::max(0., 1. + fNCNormTwkDial * frac_err_nc_norm);
}
//_______________________________________________________________________________________
double GReWeightNuXSecCOHuB::CalcWeight(const genie::EventRecord& event)
{
  Interaction* interaction = event.Summary();

  bool is_coh = interaction->ProcInfo().IsCoherent();
  if ( !is_coh ) return 1.;

  bool is_cc = interaction->ProcInfo().IsWeakCC();
  bool is_nc = interaction->ProcInfo().IsWeakNC();

  double weight = 1.;
  if ( is_cc ) weight = fCurNormCC;
  else if ( is_nc ) weight = fCurNormNC;

  return weight;
}

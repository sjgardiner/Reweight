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
    || syst == kXSecTwkDial_NormNCCOH
    || syst == kProdTwkDial_NormCCNonCOHPi0
    || syst == kProdTwkDial_NormNCNonCOHPi0 )
  {
    return true;
  }

  return false;
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCOHuB::AppliesTo(ScatteringType_t /*type*/, bool /*is_cc*/) const
{
  // Accept all events, since we can get pi0 from FSIs as well as from the
  // vertex
  return true;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOHuB::SetSystematic(GSyst_t syst, double twk_dial)
{
  if ( !this->IsHandled(syst) ) return;

  switch(syst) {
    case ( kXSecTwkDial_NormCCCOH ) :
      fCCCOHNormTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_NormNCCOH ) :
      fNCCOHNormTwkDial = twk_dial;
      break;
    case ( kProdTwkDial_NormCCNonCOHPi0 ) :
      fCCNonCOHPi0NormTwkDial = twk_dial;
      break;
    case ( kProdTwkDial_NormNCNonCOHPi0 ) :
      fNCNonCOHPi0NormTwkDial = twk_dial;
      break;

    default:
      return;
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOHuB::Reset(void)
{
  fCCCOHNormTwkDial = 0.;
  fNCCOHNormTwkDial = 0.;


  fCCNonCOHPi0NormTwkDial = 0.;
  fNCNonCOHPi0NormTwkDial = 0.;

  fCurCOHNormCC = 1.;
  fCurCOHNormNC = 1.;

  fCurNonCOHPi0NormCC = 1.;
  fCurNonCOHPi0NormNC = 1.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOHuB::Reconfigure(void)
{
  GSystUncertainty* gsu = GSystUncertainty::Instance();

  // Note: this assumes that the error is symmetric.
  // TODO: consider changing this to handle asymmetric errors on the normalization
  double frac_err_coh_cc_norm = gsu->OneSigmaErr( kXSecTwkDial_NormCCCOH );
  double frac_err_coh_nc_norm = gsu->OneSigmaErr( kXSecTwkDial_NormNCCOH );

  double frac_err_noncoh_pi0_cc_norm = gsu->OneSigmaErr( kProdTwkDial_NormCCNonCOHPi0 );
  double frac_err_noncoh_pi0_nc_norm = gsu->OneSigmaErr( kProdTwkDial_NormNCNonCOHPi0 );

  fCurCOHNormCC = std::max(0., 1. + fCCCOHNormTwkDial * frac_err_coh_cc_norm);
  fCurCOHNormNC = std::max(0., 1. + fNCCOHNormTwkDial * frac_err_coh_nc_norm);

  fCurNonCOHPi0NormCC = std::max(0., 1. + fCCNonCOHPi0NormTwkDial
    * frac_err_noncoh_pi0_cc_norm);
  fCurNonCOHPi0NormNC = std::max(0., 1. + fNCNonCOHPi0NormTwkDial
    * frac_err_noncoh_pi0_nc_norm);
}
//_______________________________________________________________________________________
double GReWeightNuXSecCOHuB::CalcWeight(const genie::EventRecord& event)
{
  Interaction* interaction = event.Summary();

  bool is_coh = interaction->ProcInfo().IsCoherent();

  // Check whether or not a pi0 was produced in the final state.
  // Note that a null pointer will evaluate to false and any other value will
  // evaluate to true when implicitly cast to bool below.
  GHepParticle* has_fs_pi0 = event.FindParticle( kPdgPi0,
    kIStStableFinalState, 0 );

  bool is_cc = interaction->ProcInfo().IsWeakCC();
  bool is_nc = interaction->ProcInfo().IsWeakNC();

  double weight = 1.;
  if ( is_coh ) {
    if ( is_cc ) weight = fCurCOHNormCC;
    else if ( is_nc ) weight = fCurCOHNormNC;
  }
  else if ( has_fs_pi0 && !is_coh ) {
    if ( is_cc ) weight = fCurNonCOHPi0NormCC;
    else if ( is_nc ) weight = fCurNonCOHPi0NormNC;
  }

  return weight;
}

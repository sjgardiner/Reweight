//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// GENIE/Generator includes
#include "Framework/Interaction/Interaction.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/PDGCodes.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightRESBugFix.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

namespace {
  const int ROOTINO = 0;
}

//_______________________________________________________________________________________
GReWeightRESBugFix::GReWeightRESBugFix() : GReWeightModel("RESbugfix")
{
}
//_______________________________________________________________________________________
GReWeightRESBugFix::~GReWeightRESBugFix()
{
}
//_______________________________________________________________________________________
bool GReWeightRESBugFix::IsHandled(GSyst_t syst) const
{
  // Handles only the dummy knob that implements the bug fix
  if ( syst == kBugFixDial_RESRootino ) return true;
  else return false;
}
//_______________________________________________________________________________________
bool GReWeightRESBugFix::AppliesTo(ScatteringType_t type, bool /*is_cc*/) const
{
  // This calculator is relevant for all RES events
  if ( type == kScResonant ) return true;
  return false;
}
//_______________________________________________________________________________________
void GReWeightRESBugFix::SetSystematic(GSyst_t /*syst*/, double /*twk_dial*/)
{
  // The tweak dial value doesn't actually matter for this calculator,
  // so this function is a no-op.
}
//_______________________________________________________________________________________
void GReWeightRESBugFix::Reset(void)
{
}
//_______________________________________________________________________________________
void GReWeightRESBugFix::Reconfigure(void)
{
}
//_______________________________________________________________________________________
double GReWeightRESBugFix::CalcWeight(const genie::EventRecord& event)
{
  Interaction* interaction = event.Summary();

  // The bug only occurs for resonant events
  bool is_res = interaction->ProcInfo().IsResonant();
  if ( !is_res ) return 1.;

  // Only events that produce a P33(1600) or F17(1970) are affected
  Resonance_t res_code = interaction->ExclTag().Resonance();
  if ( res_code != kP33_1600 && res_code != kF17_1970 ) return 1.;

  // Either of these two resonances are mislabeled by the bug
  // as a "rootino" (PDG code == 0) in the event record. If we find
  // this rootino (confirming the presence of the bug during event
  // generation), then zero out the weight for the current event.
  double weight = 1.;
  genie::GHepParticle* bad_resonance = event.FindParticle( ROOTINO,
    genie::kIStPreDecayResonantState, 0 );
  if ( bad_resonance ) weight = 0.;

  return weight;
}
//_______________________________________________________________________________________

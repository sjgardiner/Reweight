//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCOHuB

\brief    Alternate COH cross section weight calculator for MicroBooNE.
          Temporary workaround for problems uncovered with GReWeightNuXSecCOH.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Feb 18, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_COH_UB_H_
#define _G_REWEIGHT_NU_XSEC_COH_UB_H_

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

namespace genie {
namespace rew   {

 class GReWeightNuXSecCOHuB : public GReWeightModel
 {
 public:
   GReWeightNuXSecCOHuB();
  ~GReWeightNuXSecCOHuB();

   // implement the GReWeightI interface
   bool   AppliesTo      (ScatteringType_t type, bool is_cc) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 private:

   // Separate tweak dials for COH-CC and COH-NC
   double fCCCOHNormTwkDial;
   double fNCCOHNormTwkDial;

   // Separate tweak dials for CC/NC non-COH pi0 production
   double fCCNonCOHPi0NormTwkDial;
   double fNCNonCOHPi0NormTwkDial;

   // Configured normalization factors
   double fCurCOHNormCC;
   double fCurCOHNormNC;

   double fCurNonCOHPi0NormCC;
   double fCurNonCOHPi0NormNC;
 };

} // rew   namespace
} // genie namespace

#endif

//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightRESBugFix

\brief    Zeros out RES events that produced a P33(1600) or F17(1970) resonance
          and contain a "rootino" (PDG code == 0) in the particle list. This
          is one way of fixing a bug in GENIE v3.0.6 that was first reported
          during the October 2019 User Forum (see GENIE docDB #153,
          https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=153).

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Dec 27, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_RES_BUG_FIX_H_
#define _G_REWEIGHT_RES_BUG_FIX_H_

// standard library includes
#include <string>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"

namespace genie {

namespace rew   {

 class GReWeightRESBugFix : public GReWeightModel {

 public:

   GReWeightRESBugFix();
  ~GReWeightRESBugFix();

   // Implement the GReWeightI interface
   bool   AppliesTo      (ScatteringType_t type, bool is_cc) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

 };

} // rew   namespace
} // genie namespace

#endif // _G_REWEIGHT_RES_BUG_FIX_H_

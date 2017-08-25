/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCROTATECOM_H
#define DCROTATECOM_H
#include "DCComponent.h"
#include "CBMC.h"

namespace mol_setup { class MolKind; }

namespace cbmc {
   class DCData;

   class DCRotateCOM : public DCComponent
   {  
   public:
      DCRotateCOM(DCData* data);
      void PrepareNew(TrialMol& newMol, uint molIndex);
      void PrepareOld(TrialMol& oldMol, uint molIndex);
      void BuildOld(TrialMol& oldMol, uint molIndex);
      void BuildNew(TrialMol& newMol, uint molIndex);
      DCComponent* Clone() { return new DCRotateCOM(*this); };

   private:
      DCData* data;
      XYZ COM;
      uint atomNumber;
   };
}

#endif

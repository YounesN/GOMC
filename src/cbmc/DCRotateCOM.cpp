#define _USE_MATH_DEFINES 
#include <math.h> 
#include "DCRotateCOM.h" 
#include "DCData.h" 
#include "TrialMol.h" 
#include "MolSetup.h" 
#include "Forcefield.h" 
#include "PRNG.h" 
#include "NumLib.h" 
 
namespace cbmc 
{ 
 
   DCRotateCOM::DCRotateCOM(DCData* data) 
     : data(data) {} 
 
 
   void DCRotateCOM::PrepareNew(TrialMol& newMol, uint molIndex) 
   { 
     newMol.SetWeight(1.0);
     atomNumber = newMol.GetCoords().Count();
     //old center of mass
     XYZ oldCOM = newMol.GetCOM();
     //new center of mass that need to be transfered
     COM = newMol.GetSeed();
     XYZ diff = COM - oldCOM;
     
     for(uint p = 0; p < atomNumber; p++)
     {
       newMol.SetAtomCoords(p, newMol.AtomPosition(p) + diff);
     }
   } 
 
   void DCRotateCOM::PrepareOld(TrialMol& oldMol, uint molIndex) 
   { 
     oldMol.SetWeight(1.0);
     atomNumber = oldMol.GetCoords().Count();
     //old center of mass
     XYZ oldCOM = oldMol.GetCOM();
     //new center of mass that need to be transfered
     COM = oldMol.GetSeed();
     XYZ diff = COM - oldCOM;
     
     for(uint p = 0; p < atomNumber; p++)
     {
       oldMol.SetAtomCoords(p, oldMol.AtomPosition(p) + diff);
     }
   } 
 
 
   void DCRotateCOM::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const CalculateEnergy& calc = data->calc; 
      const EwaldCached *calcEwald = data->calcEwald; 
      const Forcefield& ff = data->ff; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* ljWeights = data->ljWeights; 
      double* inter = data->inter; 
      double* real = data->real; 
 
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0); 
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      //get info about existing geometry 
      newMol.ShiftBasis(COM); 
      const XYZ center = COM; 
      XYZArray* positions = data->multiPositions; 

      if(atomNumber == 1)
      {
	nLJTrials = 1;
      }
 
      for (uint p = 0; p < atomNumber; ++p) 
      { 
	positions[p].Set(0, newMol.AtomPosition(p)); 
	//data->axes.UnwrapPBC(positions[p], 0, 1, newMol.GetBox(), center); 
	positions[p].Add(0, -center);
      } 
 
      //counting backward to preserve prototype 
      for (uint lj = nLJTrials; lj-- > 0;) 
      { 
         //convert chosen torsion to 3D positions 
         RotationMatrix spin = 
            RotationMatrix::UniformRandom(prng(), prng(), prng()); 
         for (uint p = 0; p < atomNumber; ++p) 
         { 
               //find positions 
               positions[p].Set(lj, spin.Apply(positions[p][0])); 
               positions[p].Add(lj, center); 
         } 
      } 
 
      for (uint p = 0; p < atomNumber; ++p) 
      { 
         data->axes.WrapPBC(positions[p], newMol.GetBox());  
      } 
 
 
      for (uint p = 0; p < atomNumber; ++p) 
      { 
	 calc.ParticleInter(inter, real, positions[p], p,  
			    molIndex, newMol.GetBox(), nLJTrials); 
      }  
 
      double stepWeight = 0; 
      for (uint lj = 0; lj < nLJTrials; ++lj) 
      { 
	ljWeights[lj] = exp(-ff.beta * (inter[lj] + real[lj])); 
         stepWeight += ljWeights[lj]; 
      } 
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight); 

      for(uint p = 0; p < atomNumber; ++p) 
      { 
         newMol.AddAtom(p, positions[p][winner]); 
      } 
 
      newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner], 0.0, 0.0,
			      0.0)); 
      newMol.MultWeight(stepWeight); 
   } 
 
   void DCRotateCOM::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const CalculateEnergy& calc = data->calc; 
      const EwaldCached * calcEwald = data->calcEwald; 
      const Forcefield& ff = data->ff; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* ljWeights = data->ljWeights; 
      double* inter = data->inter; 
      double* real = data->real; 
 
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0);  
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      if(atomNumber == 1)
      {
	nLJTrials = 1;
      }

      //get info about existing geometry 
      oldMol.ShiftBasis(COM); 
      const XYZ center = COM; 
      XYZArray* positions = data->multiPositions; 
 
      for (uint p = 0; p < atomNumber; ++p) 
      { 
         //get position and shift to origin 
         positions[p].Set(0, oldMol.AtomPosition(p)); 
         //data->axes.UnwrapPBC(positions[p], 0, 1, oldMol.GetBox(), center); 
         positions[p].Add(0, -center); 
      }  
 
      //counting backward to preserve prototype 
      for (uint lj = nLJTrials; lj-- > 1;) 
      { 
         //convert chosen torsion to 3D positions 
         RotationMatrix spin = 
            RotationMatrix::UniformRandom(prng(), prng(), prng()); 
         for (uint p = 0; p < atomNumber; ++p) 
         { 
            //find positions 
            positions[p].Set(lj, spin.Apply(positions[p][0])); 
            positions[p].Add(lj, center); 
         } 
      } 
 
      for (uint p = 0; p < atomNumber; ++p) 
      { 
         positions[p].Add(0, center); 
         data->axes.WrapPBC(positions[p], oldMol.GetBox()); 
      } 
 
      for (uint p = 0; p < atomNumber; ++p)
      { 
	calc.ParticleInter(inter, real, positions[p], p, 
                            molIndex, oldMol.GetBox(), nLJTrials); 
      } 

      double stepWeight = 0;  
      for (uint lj = 0; lj < nLJTrials; ++lj) 
      { 
	stepWeight += exp(-ff.beta * (inter[lj] + real[lj])); 
      } 

      for (uint p = 0; p < atomNumber; ++p)
      { 
         oldMol.AddAtom(p, positions[p][0]); 
      } 
 
      oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0], 0.0, 0.0, 0.0)); 
      oldMol.MultWeight(stepWeight); 
   } 
 
 
}              

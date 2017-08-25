#ifndef IDENTITYEXCHANGE_H
#define IDENTITYEXCHANGE_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "cbmc/TrialMol.h"

// Identity Exchange Move:
//
// insert the molecule A inside the cavity that has the size of the molecule B
// and vice versa. Only 1 trial must be perform for the first seed.
// Mohammad Soroush Barhaghi. July 2017

class IdentityExchange : public MoveBase
{
 public:

   IdentityExchange(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
      MoveBase(sys, statV), rmax(statV.forcefield.rmax) {}

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint earlyReject, const uint step);

 private:

   void FindRmax();
   double GetCoeff() const;
   uint GetBoxPairAndMol(const double subDraw, const double movPerc);
   MolPick molPick;
   uint sourceBox, destBox;
   uint pStartA, pStartB, pLenA, pLenB;
   uint molIndexA, molIndexB, kindIndexA, kindIndexB;
   bool hasMol;
#if ENSEMBLE == GEMC
   bool subVSourceBox;
   uint subVBox;
#endif
   double rmax;
   XYZ center;
   double W_tc, W_recip;
   double correct_oldA, correct_newA, self_oldA, self_newA;
   double correct_oldB, correct_newB, self_oldB, self_newB;
   cbmc::TrialMol oldMolA, oldMolB, newMolA, newMolB;
   Intermolecular tcNew[BOX_TOTAL], recipLoseA, recipGainA;
   Intermolecular recipLoseB, recipGainB;
   std::vector<uint> molInCav;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline uint IdentityExchange::GetBoxPairAndMol
(const double subDraw, const double movPerc)
{
   //Set the cavity diameter
   FindRmax();
   //Pick box at random
   uint state;

#if ENSEMBLE == GEMC
   double density;
   double maxDens = 0.0;
   uint densB;
   //choose the sourceBox to be the dense phase
   for(uint b = 0; b < BOX_TOTAL; b++)
   {
     density = 0.0;
     for(uint k = 0; k < molLookRef.GetNumKind(); k++)
     {
       density += molLookRef.NumKindInBox(k, b) * boxDimRef.volInv[b] *
	 molRef.kinds[k].molMass;
     }

     if(density > maxDens)
     {
       maxDens = density;
       densB = b;
     }
   }

   //Pick box at random
   prng.PickBox(sourceBox, subDraw, movPerc); 
   //Pick the destination box
   prng.SetOtherBox(destBox, sourceBox);
   //pick the box for searching the molecule in cavity
   prng.PickBox(subVBox, subDraw, movPerc);

   if(subVBox == sourceBox)
     subVSourceBox = true;
   else
     subVSourceBox = false;

   //pick a random location in dense phase
   XYZ axis = boxDimRef.GetAxis(subVBox);
   XYZ temp(prng.randExc(axis.x), prng.randExc(axis.y), prng.randExc(axis.z));
   center = temp;

   //Find the molecule in the cavity
   molInCav.clear();
   hasMol = calcEnRef.FindMolInCavity(molInCav, center, rmax, subVBox);

   if(hasMol)
   {
     if(subVSourceBox)
     {
       //Find a molecule at random
       uint i = prng.randIntExc(molInCav.size());
       molIndexA = molInCav[i];
       kindIndexA = molRef.GetMolKind(molIndexA);
       //shift the center to COM of MoleculeA
       //center = comCurrRef.Get(molIndexA);
       //pick a molecule from other kind in less dense phase.
       state = prng.PickMol(kindIndexA, kindIndexB, molIndexB, destBox);
     }
     else
     {
       //Find a molecule at random
       uint i = prng.randIntExc(molInCav.size());
       molIndexB = molInCav[i];
       kindIndexB = molRef.GetMolKind(molIndexB);
       //shift the center to COM of MoleculeB
       //center = comCurrRef.Get(molIndexB);
       //pick a molecule from other kind in less dense phase.
       state = prng.PickMol(kindIndexB, kindIndexA, molIndexA, sourceBox);
     }
   }
   else
   {
     //reject the move
     state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
   }

#elif ENSEMBLE == GCMC
   sourceBox = mv::BOX0;
   destBox = mv::BOX1;

   //for GCMC pick a random location only in Box0
   XYZ axis = boxDimRef.GetAxis(sourceBox);

   XYZ temp(prng.randExc(axis.x), prng.randExc(axis.y), prng.randExc(axis.z));
   center = temp;

   //Find the molecule in the cavity
   molInCav.clear();
   hasMol = calcEnRef.FindMolInCavity(molInCav, center, rmax, sourceBox);

   if(hasMol)
   {
     //Find a molecule at random
     uint i = prng.randIntExc(molInCav.size());
     molIndexA = molInCav[i];
     kindIndexA = molRef.GetMolKind(molIndexA);
     //shift the center to COM of MoleculeA
     //center = comCurrRef.Get(molIndexA);
     //pick a molecule from other kind in resv.
     state = prng.PickMol(kindIndexA, kindIndexB, molIndexB, destBox);
   }
   else
   {
     //reject the move
     state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
   }
#endif
   
   if(state == mv::fail_state::NO_FAIL)
   {
     pStartA = pLenA = pStartB = pLenB = 0;
     molRef.GetRangeStartLength(pStartA, pLenA, molIndexA);
     molRef.GetRangeStartLength(pStartB, pLenB, molIndexB);
   }

   return state;
}

inline void IdentityExchange::FindRmax()
{
  //rmax = 2.0;
}

inline uint IdentityExchange::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   if(state == mv::fail_state::NO_FAIL)
   {
     //transfering type A from source to dest
     newMolA = cbmc::TrialMol(molRef.kinds[kindIndexA], boxDimRef, destBox);
     oldMolA = cbmc::TrialMol(molRef.kinds[kindIndexA], boxDimRef, sourceBox);
     //transfering type B from dest to source
     newMolB = cbmc::TrialMol(molRef.kinds[kindIndexB], boxDimRef, sourceBox);
     oldMolB = cbmc::TrialMol(molRef.kinds[kindIndexB], boxDimRef, destBox);
     //set the old coordinate
     XYZArray molA(pLenA);
     XYZArray molB(pLenB);
     coordCurrRef.CopyRange(molA, pStartA, 0, pLenA);
     coordCurrRef.CopyRange(molB, pStartB, 0, pLenB);
     boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndexA));
     boxDimRef.UnwrapPBC(molB, destBox, comCurrRef.Get(molIndexB));
     oldMolA.SetCoords(molA, 0);
     oldMolB.SetCoords(molB, 0);

     XYZ axisD = boxDimRef.GetAxis(destBox);        
     XYZ axisS = boxDimRef.GetAxis(sourceBox);

#if ENSEMBLE == GEMC     

     if(subVSourceBox)
     {
       //Pick moleculeA from the cavity in sourceBox and insert it to destBox
       //Pick moleculeB from destBox and insert it in the cavity in sourceBox
      
       //set coordinate of old B to newMolB
       newMolB.SetCoords(molB, 0);
       //rotate COM of MolB around center
       newMolB.SetSeed(center, rmax);   
       //use the COM of oldB to rotate around
       oldMolB.SetSeed(comCurrRef.Get(molIndexB), rmax);

       //set coordinate of oldA to newMolA
       newMolA.SetCoords(molA, 0);
       XYZ tempD(prng.randExc(axisD.x), prng.randExc(axisD.y),
		 prng.randExc(axisD.z));
       //randomly pick a locatoin in destBox to rotate A around that
       newMolA.SetSeed(tempD, rmax);
       //use the COM of oldA to rotate around
       oldMolA.SetSeed(comCurrRef.Get(molIndexA), rmax);
     }
     else
     {
       //Pick moleculeB from the cavity in destBox and insert it to sourceBox
       //Pick moleculeA from sourceBox and insert it in the cavity in destBox

       //set coordinate of old A to newMolA
       newMolA.SetCoords(molA, 0);
       //rotate COM of MolA around center
       newMolA.SetSeed(center, rmax);   
       //use the COM of oldA to rotate around
       oldMolA.SetSeed(comCurrRef.Get(molIndexA), rmax);

       //set coordinate of oldB to newMolB
       newMolB.SetCoords(molB, 0);
       XYZ tempS(prng.randExc(axisS.x), prng.randExc(axisS.y),
		 prng.randExc(axisS.z));
       //randomly pick a locatoin in source to rotate B around that
       newMolB.SetSeed(tempS, rmax);
       //use the COM of oldB to rotate around
       oldMolB.SetSeed(comCurrRef.Get(molIndexB), rmax);
     }
#elif ENSEMBLE == GCMC 
     //Inserting molB from Resv to sourceBox
     //set coordinate of old B to newMolB
     newMolB.SetCoords(molB, 0);
     //rotate COM of MolB around center
     newMolB.SetSeed(center, rmax);   
     //use the COM of oldB to rotate around
     oldMolB.SetSeed(comCurrRef.Get(molIndexB), rmax);
     
     ////Inserting molA from sourceBox to Resv
     //set coordinate of oldA to newMolA
     newMolA.SetCoords(molA, 0);
     XYZ tempD(prng.randExc(axisD.x), prng.randExc(axisD.y),
		prng.randExc(axisD.z));
     //randomly pick a locatoin in destBox to rotate A around that
     newMolA.SetSeed(tempD, rmax);
     //use the COM of oldA to rotate around
     oldMolA.SetSeed(comCurrRef.Get(molIndexA), rmax);
#endif
   }

   return state;
}


inline uint IdentityExchange::Transform()
{
  //Deleting Molecule from cellList to avoid overlaping
  cellList.RemoveMol(molIndexA, sourceBox, coordCurrRef);
  cellList.RemoveMol(molIndexB, destBox, coordCurrRef);
  //Transfer Type A to resv
  molRef.kinds[kindIndexA].BuildID(oldMolA, newMolA, molIndexA);
  //Transfer Type B to box0
  molRef.kinds[kindIndexB].BuildID(oldMolB, newMolB, molIndexB);
     
  return mv::fail_state::NO_FAIL;
}

inline void IdentityExchange::CalcEn()
{   
   W_tc = 1.0, W_recip = 1.0;
   correct_oldA = 0.0, correct_newA = 0.0;
   self_oldA = 0.0, self_newA = 0.0;
   correct_oldB = 0.0, correct_newB = 0.0;
   self_oldB = 0.0, self_newB = 0.0;

   if (ffRef.useLRC)
   {
      double delTC = 0.0;
      for (uint b = 0; b < BOX_TOTAL; ++b)
      {
	 uint kCount[molRef.kindsCount];
	 for (uint k = 0; k < molRef.kindsCount; ++k)
	 {
	    kCount[k] = molLookRef.NumKindInBox(k, b);
	 }

	 if (b == sourceBox)
	 {
	   --kCount[kindIndexA];
	   ++kCount[kindIndexB];	   
	 }
	 else if (b == destBox)
	 {
	   ++kCount[kindIndexA];
	   --kCount[kindIndexB];
	 }
	 tcNew[b].energy = calcEnRef.EnergyCorrection(b, kCount);
	 delTC += tcNew[b].energy - sysPotRef.boxEnergy[b].tc;
      }
     W_tc = exp(-1.0 * ffRef.beta * delTC); 
   }
   
   if (newMolA.GetWeight() != 0.0 && newMolB.GetWeight() != 0.0)
   {
      correct_newA = calcEwald->SwapCorrection(newMolA);
      correct_oldA = calcEwald->SwapCorrection(oldMolA);
      correct_newB = calcEwald->SwapCorrection(newMolB);
      correct_oldB = calcEwald->SwapCorrection(oldMolB);

      self_newA = calcEwald->SwapSelf(newMolA);
      self_oldA = calcEwald->SwapSelf(oldMolA);
      self_newB = calcEwald->SwapSelf(newMolB);
      self_oldB = calcEwald->SwapSelf(oldMolB);

      recipGainA.energy =
	calcEwald->SwapDestRecip(newMolA, destBox, sourceBox, molIndexA);
      recipLoseA.energy =
	calcEwald->SwapSourceRecip(oldMolA, sourceBox, molIndexA);

      recipGainB.energy =
	calcEwald->SwapDestRecip(newMolB, sourceBox, destBox, molIndexB);
      recipLoseB.energy =
	calcEwald->SwapSourceRecip(oldMolB, destBox, molIndexB);
      //need to contribute the self and correction energy 
      W_recip = exp(-1.0 * ffRef.beta * (recipGainA.energy + recipLoseA.energy +
					 recipGainB.energy + recipLoseB.energy +
					 correct_newA - correct_oldA +
					 correct_newB - correct_oldB +
					 self_newA - self_oldA +
					 self_newB - self_oldB));
   }
}

inline double IdentityExchange::GetCoeff() const
{
  double numTypeASource = molLookRef.NumKindInBox(kindIndexA, sourceBox);
  double numTypeADest = molLookRef.NumKindInBox(kindIndexA, destBox);
  double numTypeBSource = molLookRef.NumKindInBox(kindIndexB, sourceBox);
  double numTypeBDest = molLookRef.NumKindInBox(kindIndexB, destBox);
#if ENSEMBLE == GEMC
  
  if(subVSourceBox)
    return numTypeBDest / (numTypeADest + 1);
  else
    return numTypeASource / (numTypeBSource + 1);
  

#elif ENSEMBLE == GCMC
  if(ffRef.isFugacity)
  {
    double delA = BETA * molRef.kinds[kindIndexA].chemPot;
    double insB = BETA * molRef.kinds[kindIndexB].chemPot;
    return insB / delA;
  }
  else
  {
    double delA = exp(-BETA * molRef.kinds[kindIndexA].chemPot);
    double insB = exp(BETA * molRef.kinds[kindIndexB].chemPot);
    return insB * delA;
  }
#endif
}

inline void IdentityExchange::Accept(const uint rejectState, const uint step)
{
   bool result;
   //If we didn't skip the move calculation
   if(rejectState == mv::fail_state::NO_FAIL)
   {
      double molTransCoeff = GetCoeff();
      double WoA = oldMolA.GetWeight();
      double WoB = oldMolB.GetWeight();
      double WnA = newMolA.GetWeight();
      double WnB = newMolB.GetWeight();
     
      double Wrat = (WnA * WnB) * W_tc * W_recip / (WoA * WoB);
      //std::cout << "Wrat: " << Wrat << std::endl;
      result = prng() < molTransCoeff * Wrat;
     
      if(result)
      {
         //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc = tcNew[sourceBox].energy;
         sysPotRef.boxEnergy[destBox].tc = tcNew[destBox].energy;

         //Add rest of energy.
	 sysPotRef.boxEnergy[sourceBox] += newMolB.GetEnergy();
	 sysPotRef.boxEnergy[sourceBox] -= oldMolA.GetEnergy();
	 sysPotRef.boxEnergy[destBox] -= oldMolB.GetEnergy();
	 sysPotRef.boxEnergy[destBox] += newMolA.GetEnergy();
	 

	 //Add Reciprocal energy
	 sysPotRef.boxEnergy[sourceBox].recip += recipLoseA.energy;
	 sysPotRef.boxEnergy[sourceBox].recip += recipGainB.energy;
	 sysPotRef.boxEnergy[destBox].recip += recipGainA.energy;
	 sysPotRef.boxEnergy[destBox].recip += recipLoseB.energy;	 
	 //Add correction energy
	 sysPotRef.boxEnergy[sourceBox].correction -= correct_oldA;
	 sysPotRef.boxEnergy[sourceBox].correction += correct_newB;
	 sysPotRef.boxEnergy[destBox].correction += correct_newA;
	 sysPotRef.boxEnergy[destBox].correction -= correct_oldB;	 
	 //Add self energy
	 sysPotRef.boxEnergy[sourceBox].self -= self_oldA;
	 sysPotRef.boxEnergy[sourceBox].self += self_newB;
	 sysPotRef.boxEnergy[destBox].self += self_newA;
	 sysPotRef.boxEnergy[destBox].self -= self_oldB;
	 
	 for (uint b = 0; b < BOX_TOTAL; b++)
	 {
	    calcEwald->UpdateRecip(b);
	 }

	 //Add type A to dest box
	 newMolA.GetCoords().CopyRange(coordCurrRef, 0, pStartA, pLenA);
	 comCurrRef.SetNew(molIndexA, destBox);
	 molLookRef.ShiftMolBox(molIndexA, sourceBox, destBox, kindIndexA);
	 //Add type B to source box
	 newMolB.GetCoords().CopyRange(coordCurrRef, 0, pStartB, pLenB);
	 comCurrRef.SetNew(molIndexB, sourceBox);
	 molLookRef.ShiftMolBox(molIndexB, destBox, sourceBox, kindIndexB);
	   
	 cellList.AddMol(molIndexA, destBox, coordCurrRef);
	 cellList.AddMol(molIndexB, sourceBox, coordCurrRef);

	 //Retotal
         sysPotRef.Total();
      }
      else
      {
        cellList.AddMol(molIndexA, sourceBox, coordCurrRef);
	cellList.AddMol(molIndexB, destBox, coordCurrRef);
	  
	//when weight is 0, MolDestSwap() will not be executed, thus cos/sin
	//molRef will not be changed. Also since no memcpy, doing restore
	//results in memory overwrite
	if (newMolA.GetWeight() != 0.0 && newMolB.GetWeight() != 0.0)
	  calcEwald->RestoreMol(molIndexA);
	//Need the way to handle the reciprocal issue for cache version
      }
   }
   else  //else we didn't even try because we knew it would fail
      result = false;

#if ENSEMBLE == GEMC
   subPick = mv::GetMoveSubIndex(mv::ID_EXCHANGE, sourceBox);
#elif ENSEMBLE == GCMC
   subPick = mv::GetMoveSubIndex(mv::ID_EXCHANGE);
#endif
   moveSetRef.Update(result, subPick, step);
}

#endif

#endif

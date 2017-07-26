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
    MoveBase(sys, statV) {}

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

   double rmaxA, rmaxB;
   double W_tc, W_recip;
   double correct_oldA, correct_newA, self_oldA, self_newA;
   double correct_oldB, correct_newB, self_oldB, self_newB;
   cbmc::TrialMol oldMolA, oldMolB, newMolA, newMolB;
   Intermolecular tcNew[BOX_TOTAL], recipLoseA, recipGainA;
   Intermolecular recipLoseB, recipGainB;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline uint IdentityExchange::GetBoxPairAndMol
(const double subDraw, const double movPerc)
{
   //Pick two molecule with different ID from two different box
   uint state = prng.PickMolAndBoxPair(molIndexA, molIndexB, kindIndexA,
				       kindIndexB, sourceBox, destBox, subDraw,
				       movPerc);

 
   if ( state != mv::fail_state::NO_TWO_MOLECULE_KIND)
   {
      pStartA = pLenA = pStartB = pLenB = 0;
      molRef.GetRangeStartLength(pStartA, pLenA, molIndexA);
      molRef.GetRangeStartLength(pStartB, pLenB, molIndexB);
   }
   return state;
}

inline void IdentityExchange::FindRmax()
{
  rmaxA = 1.0;
  rmaxB = 1.0;
}

inline uint IdentityExchange::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   //transfering type A from source to dest
   newMolA = cbmc::TrialMol(molRef.kinds[kindIndexA], boxDimRef, destBox);
   oldMolA = cbmc::TrialMol(molRef.kinds[kindIndexA], boxDimRef, sourceBox);
   //transfering type B from dest to source
   newMolB = cbmc::TrialMol(molRef.kinds[kindIndexB], boxDimRef, sourceBox);
   oldMolB = cbmc::TrialMol(molRef.kinds[kindIndexB], boxDimRef, destBox);
   
   oldMolA.SetCoords(coordCurrRef, pStartA);
   oldMolB.SetCoords(coordCurrRef, pStartB);
   //set center of cavity and radius.
   FindRmax();
   //pick the first position around COM of other molecule.
   newMolA.SetSeed(comCurrRef.Get(molIndexB), rmaxB);
   newMolB.SetSeed(comCurrRef.Get(molIndexA), rmaxA);
   //pick the first position around COM of itself.
   oldMolA.SetSeed(comCurrRef.Get(molIndexA), rmaxB);
   oldMolB.SetSeed(comCurrRef.Get(molIndexB), rmaxA);
   W_tc = 1.0;
   return state;
}


inline uint IdentityExchange::Transform()
{
   //Deleting Molecule from cellList to avoid overlaping
   cellList.RemoveMol(molIndexA, sourceBox, coordCurrRef);
   cellList.RemoveMol(molIndexB, destBox, coordCurrRef);

   //Transfer Type A to dest Box
   molRef.kinds[kindIndexA].Build(oldMolA, newMolA, molIndexA);
   //Transfer Type B to source Box
   molRef.kinds[kindIndexB].Build(oldMolB, newMolB, molIndexB);

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
  return numTypeBDest * numTypeASource /
    ((numTypeBSource + 1.0) * (numTypeADest + 1.0));
#elif ENSEMBLE == GCMC
  if (sourceBox == mv::BOX0) //Removing A from Box 0
  {
    if(ffRef.isFugacity)
    {
    }
    else
    {
      double delA = exp(-BETA *  molRef.kinds[kindIndexA].chemPot) *
	numTypeASource;
      double instB = exp(BETA * molRef.kinds[kindIndexB].chemPot) /
	(numTypeBSource + 1) ;
      return delA * instB;
    }
  }
  else //Insert A to Box 0
  {
    if(ffRef.isFugacity)
    {
    }
    else
    {
      double instA =  exp(BETA * molRef.kinds[kindIndexA].chemPot)/
	(numTypeADest + 1);
      double delB = exp(-BETA * molRef.kinds[kindIndexB].chemPot) *
	numTypeBDest;
      return instA * delB;
    }
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
      
      result = prng() < molTransCoeff * Wrat;
      if (result)
      {
         //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc = tcNew[sourceBox].energy;
         sysPotRef.boxEnergy[destBox].tc = tcNew[destBox].energy;

         //Add rest of energy.
         sysPotRef.boxEnergy[sourceBox] -= oldMolA.GetEnergy();
         sysPotRef.boxEnergy[sourceBox] += newMolB.GetEnergy();
         sysPotRef.boxEnergy[destBox] += newMolA.GetEnergy();
         sysPotRef.boxEnergy[destBox] -= oldMolB.GetEnergy();
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

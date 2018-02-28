/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCOnSphere.h"
#include "TrialMol.h"
#include "DCData.h"
#include "XYZArray.h"
#include "PRNG.h"
#include "Forcefield.h"
#include "MolSetup.h"

namespace cbmc
{
DCOnSphere::DCOnSphere(DCData* data, const mol_setup::MolKind kind,
                       uint atom, uint focus) :
  data(data), atom(atom),
  focus(focus)
{
  using namespace mol_setup;
  std::vector<Bond> bonds = AtomBonds(kind, atom);
  for(uint i = 0; i < bonds.size(); ++i) {
    if(bonds[i].a0 == focus || bonds[i].a1 == focus) {
      bondKind = bonds[i].kind;
      break;
    }
  }
}

void DCOnSphere::SetBondLengthNew(TrialMol& newMol)
{
  bondLength = data->ff.bonds.Length(bondKind);
}

void DCOnSphere::SetBondLengthOld(TrialMol& oldMol)
{
  bondLengthOld = sqrt(oldMol.OldDistSq(focus, atom));
}

void DCOnSphere::BondEnergy(TrialMol& mol)
{
  bondEnergy = data->ff.bonds.Calc(bondKind, bondLength);
}

void DCOnSphere::BuildOld(TrialMol& oldMol, uint molIndex)
{
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsNth;
  double* inter = data->inter;
  double* real = data->real;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  //calculate bond energy for old molecule.
  bondEnergy = 0.0;
  BondEnergy(oldMol);

  data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
                                    oldMol.AtomPosition(focus));
  positions.Set(0, oldMol.AtomPosition(atom));
  data->axes.WrapPBC(positions, oldMol.GetBox());

  data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                           oldMol.GetBox(), nLJTrials);


  for (uint trial = 0; trial < nLJTrials; trial++) {
    stepWeight += exp(-1 * data->ff.beta *
                      (inter[trial] + real[trial]));
  }
  oldMol.MultWeight(stepWeight);
  oldMol.AddEnergy(Energy(bondEnergy, 0.0, inter[0], real[0], 0.0,
                          0.0, 0.0));
  oldMol.ConfirmOldAtom(atom);
}

void DCOnSphere::BuildNew(TrialMol& newMol, uint molIndex)
{
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsNth;
  double* inter = data->inter;
  double* real = data->real;
  double* ljWeights = data->ljWeights;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  //calculate bond energy for old molecule.
  bondEnergy = 0.0;
  BondEnergy(newMol);

  data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
                                    newMol.AtomPosition(focus));
  data->axes.WrapPBC(positions, newMol.GetBox());

  data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                           newMol.GetBox(), nLJTrials);


  for (uint trial = 0; trial < nLJTrials; trial++) {
    ljWeights[trial] = exp(-1 * data->ff.beta *
                           (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = data->prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.MultWeight(stepWeight);
  newMol.AddEnergy(Energy(bondEnergy, 0, inter[winner], real[winner], 0.0,
                          0.0, 0.0));
  newMol.AddAtom(atom, positions[winner]);
}
}

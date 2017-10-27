// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>> Helper for W+c event selection >>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Consult analysis documentation (papers, description-wcharm.pdf) for
// better description of applied cuts etc.

#ifndef WCHARM_SELECTION_H
#define WCHARM_SELECTION_H

// additional files from this analysis
#include "wcharm_tree.h"
// ROOT header files and C++ libraries
#include <TLorentzVector.h>
#include <exception>

// constants: electron and muon masses
// (not the best practice to make them global variables, be aware)
const double massEl = 0.000511;
const double massMu = 0.105658;

// Routine for electron selection
// Arguments:
//   const ZTreeWcharm* preselTree: input tree (see wcharm_tree.h), GetEntry() should be done already
//   vec: (output) vector of selected electron.
// Returns:
//   +1 if selected e-,
//   -1 if selected e+,
//   0 if nothing selected.
// See wcharm_tree.h for ZTree variables description.
int SelectEl(const ZTreeWcharm* preselTree, TLorentzVector& vec)
{
  int best = -1; // index of selected electron, initialised with impossible value
  double bestpt = 0.0; // pT of the selected electron, initialised with 0, because electron with highest pT will be selected)

  // loop over electrons
  for(int el = 0; el < preselTree->Nel; el++)
  {
    // skip current electron, if another electron with higher pT was found
    if(TMath::Abs(preselTree->el_pt[el]) < bestpt)
      continue;
    // require pT(el) > 35 GeV
    if(TMath::Abs(preselTree->el_pt[el]) < 35.0)
      continue;
    // require |eta(el)| < 2.1 but exclude 1.44 < |eta(el)| < 1.57
    if(TMath::Abs(preselTree->el_eta[el]) > 2.1)
      continue;
    if(TMath::Abs(preselTree->el_eta[el]) > 1.44 && TMath::Abs(preselTree->el_eta[el]) < 1.57)
      continue;
    // require isolation < 0.05
    if(preselTree->el_isolat[el] > 0.05)
      continue;
    // require transverse mass > 55 GeV
    if(preselTree->el_mt[el] < 55.0)
      continue;
    // require no missing hits
    if(preselTree->el_mh[el] > 0)
      continue;
    // cuts on conversiopn variables not applied: study them if you want
    //if(preselTree->elConvDist[el] < 0.02 && preselTree->elConvDcot[el] < 0.02 && preselTree->elConvDist[el] >= 0.0 && preselTree->elConvDcot[el] >=0.0)
    //  continue;

    // store index and pT of this electron
    best = el;
    bestpt = TMath::Abs(preselTree->el_pt[el]);
  }

  // if no electron found, return 0
  if(best < 0)
    return 0;

  // loop over electrons again to skip events with e+e- pairs
  for(int el = 0; el < preselTree->Nel; el++)
  {
    // skip selected electron
    if(el == best)
      continue;
    // skip if selected electron and current one have the same sign
    if(preselTree->el_pt[el] * preselTree->el_pt[best] > 0)
      continue;
    // skip if current electron has pT < 20 GeV or |eta| > 2.5
    if(TMath::Abs(preselTree->el_pt[el]) < 20.0)
      continue;
    if(TMath::Abs(preselTree->el_eta[el]) > 2.5)
      continue;
    // if we are here, than e+e- pair was found: return 0
    return 0;
  }

  // all selection cuts passed, no e+e- pair: assign output vector to selected electron and return proper value
  vec.SetPtEtaPhiM(TMath::Abs(preselTree->el_pt[best]), preselTree->el_eta[best], preselTree->el_phi[best], massEl);
  return (preselTree->el_pt[best] > 0) ? +1 : -1;
}

// Routine for muon selection
// Arguments:
//   const ZTreeWcharm* preselTree: input tree (see wcharm_tree.h), GetEntry() should be done already
//   vec: (output) vector of selected muon.
// Returns:
//   +1 if selected mu-,
//   -1 if selected mu+,
//   0 if nothing selected.
// See wcharm_tree.h for ZTree variables description.
int SelectMu(const ZTreeWcharm* preselTree, TLorentzVector& vec)
{
  int best = -1; // index of selected muon, initialised with impossible value
  double bestpt = 0.0; // pT of the selected muon, initialised with 0, because muon with highest pT will be selected)

  // loop over muons
  for(int mu = 0; mu < preselTree->Nmu; mu++)
  {
    // skip current muon, if another muon with higher pT was found
    if(TMath::Abs(preselTree->mu_pt[mu]) < bestpt)
      continue;
    // require pT(el) > 25 GeV
    if(TMath::Abs(preselTree->mu_pt[mu]) < 25.0)
      continue;
    // require |eta(el)| < 2.1
    if(TMath::Abs(preselTree->mu_eta[mu]) > 2.1)
      continue;
    // require isolation < 0.12
    if(preselTree->mu_isolat[mu] > 0.12)
      continue;
    // require transverse mass > 40 GeV
    if(preselTree->mu_mt[mu] < 40.0)
      continue;
    // require impact parameter < 0.2
    if(preselTree->mu_imp[mu] > 0.2)
      continue;
    // require muon track number of degrees of greater than 10
    if(preselTree->mu_nChi2[mu] > 10.0)
      continue;
    // require at least 11 tracker hits and at least 1 pixel hit
    if(preselTree->mu_ValidHits[mu] < 12 || preselTree->mu_PixelHits[mu] < 2)
      continue;

    // store index and pT of this muon
    best = mu;
    bestpt = TMath::Abs(preselTree->mu_pt[mu]);
  }

  // if no muon found, return 0
  if(best < 0)
    return 0;

  // loop over muons again to skip events with mu+mu- pairs
  for(int mu = 0; mu < preselTree->Nmu; mu++)
  {
    // skip selected muon
    if(mu == best)
      continue;
    // skip if selected muon and current one have the same sign
    if(preselTree->mu_pt[mu] * preselTree->mu_pt[best] > 0)
      continue;
    // skip if current muon has pT < 25 GeV or |eta| > 2.4
    if(TMath::Abs(preselTree->mu_pt[mu]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->mu_eta[mu]) > 2.4)
      continue;
    // if we are here, than e+e- pair was found: return 0
    return 0;
  }

  // all selection cuts passed, no mu+mu- pair: assign output vector to selected electron and return proper value
  vec.SetPtEtaPhiM(TMath::Abs(preselTree->mu_pt[best]), preselTree->mu_eta[best], preselTree->mu_phi[best], massMu);
  return (preselTree->mu_pt[best] > 0) ? +1 : -1;
}

// Routine for D* selection
// Arguments:
//   const ZTreeWcharm* preselTree: input tree (see wcharm_tree.h), GetEntry() should be done already,
//   vecJ: (output) vector of the jet matched to selected D*,
//   vecD: (output) vector of selected D*,
//   flagInMassWindow: (output) boolean flag: true if mass difference M(D*) - m(D0) of selected D* is inside 145 MeV mass window, false otherwise,
//   deltaMass: (output) mass difference M(D*) - m(D0) of selected D*.
// Returns:
//   +1 if selected D*+,
//   -1 if selected D*-,
//   0 if nothing selected.
// See wcharm_tree.h for ZTree variables description.
int SelectDstar(const ZTreeWcharm* preselTree, TLorentzVector& vecJ, TLorentzVector& vecD, bool& flagInMassWindow, double& deltaMass)
{
  int bestJ = -1; // index of matched jet, initialised with impossible value
  int bestD = -1; // index of selected meson, initialised with impossible value
  double bestpt = 0.0; // pT of the selected meson, initialised with 0, because meson with highest pT will be selected)

  // loop over D*
  for(int d = 0; d < preselTree->DsCand; d++)
  {
    // skip current meson, if another meson with higher pT was found
    if(TMath::Abs(preselTree->ds_pt[d]) < bestpt)
      continue;
    // require D0 decay length significance > 3.0
    if(preselTree->D0dls[d] < 3.0)
      continue;
    // require D0 decay length < 2.0 cm
    if(preselTree->D0dl[d] > 2.0)
      continue;
    // require z component of distance between D0 and primary vertex < 2.0 cm
    if(TMath::Abs(preselTree->pv_z - preselTree->D0sv_z[d]) > 2.0)
      continue;

    // get macthed jet
    int j = 0;
    int nd = 0;
    // loop over all stored jets looking at meson indices assigned to them
    for(j = 0; j < preselTree->Njet; j++)
    {
      nd += preselTree->j_DsInJ[j];
      // if the current index is above the current meson, there is no matched jet:  break loop over jets
      if(nd >= d)
        break;
    }
    // consistency check
    if(j == preselTree->Njet)
    {
      //throw std::logic_error("Wrong j <---> D finding\n");
      //printf("Wrong j <---> D finding\n");
      //return 0;
    }

    // check if matched jet satisfies pT and eta requirements
    if(TMath::Abs(preselTree->j_pt[j]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->j_eta[j]) > 2.5)
      continue;

    // store this meson and matched jet
    bestpt = TMath::Abs(preselTree->ds_pt[d]);
    bestD = d;
    bestJ = j;
  }

  // if no meson found, return 0
  if(bestD < 0)
    return 0;

  // all selection cuts passed: assign output variables and return proper value
  vecJ.SetPtEtaPhiM(TMath::Abs(preselTree->j_pt[bestJ]), preselTree->j_eta[bestJ], preselTree->j_phi[bestJ], preselTree->j_m[bestJ]);
  vecD.SetPtEtaPhiM(TMath::Abs(preselTree->ds_pt[bestD]), preselTree->ds_eta[bestD], preselTree->ds_phi[bestD], preselTree->ds_m[bestD]);
  deltaMass = preselTree->dmDs[bestD];
  flagInMassWindow = (TMath::Abs(deltaMass - 0.145) < 0.005) ? true : false;
  return (preselTree->ds_pt[bestD] > 0) ? +1 : -1;
}

// Routine for D+ selection (both D+ and D- are referred to as D+, unless explicitely stated otherwise)
// Arguments:
//   const ZTreeWcharm* preselTree: input tree (see wcharm_tree.h), GetEntry() should be done already,
//   vecJ: (output) vector of the jet matched to selected D*,
//   vecD: (output) vector of selected D*,
//   flagInMassWindow: (output) boolean flag: true if mass difference (M(D+) - 1870 MeV) is inside 50 MeV mass window, false otherwise (1870 MeV is world average D+ mass).
// Returns:
//   +1 if selected D+,
//   -1 if selected D-,
//   0 if nothing selected.
// See wcharm_tree.h for ZTree variables description.
int SelectDch(const ZTreeWcharm* preselTree, TLorentzVector& vecJ, TLorentzVector& vecD, bool& flagInMassWindow)
{
  int bestJ = -1; // index of matched jet, initialised with impossible value
  int bestD = -1; // index of selected meson, initialised with impossible value
  double bestpt = 0.0; // pT of the selected meson, initialised with 0, because meson with highest pT will be selected)

  // loop over D+
  for(int d = 0; d < preselTree->DCand; d++)
  {
    // skip current meson, if another meson with higher pT was found
    if(TMath::Abs(preselTree->d_pt[d]) < bestpt)
      continue;
    // require D+ decay length significance > 3.0
    if(preselTree->Ddls[d] < 3.0)
      continue;
    // require D+ decay length < 2.0 cm
    if(preselTree->Ddl[d] > 2.0)
      continue;
    // require z component of distance between D+ and primary vertex < 2.0 cm
    if(TMath::Abs(preselTree->pv_z - preselTree->Dsv_z[d]) > 2.0)
      continue;

    // get macthed jet
    int j = 0;
    int nd = 0;
    // loop over all stored jets looking at meson indices assigned to them
    for(j = 0; j < preselTree->Njet; j++)
    {
      nd += preselTree->j_DInJ[j];
      // if the current index is above the current meson, there is no matched jet:  break loop over jets
      if(nd >= d)
        break;
    }
    // consistency check
    if(j == preselTree->Njet)
    {
      //throw std::logic_error("Wrong j <---> D finding\n");
      //printf("Wrong j <---> D finding\n");
      //return 0;
    }

    // check if matched jet satisfies pT and eta requirements
    if(TMath::Abs(preselTree->j_pt[j]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->j_eta[j]) > 2.5)
      continue;

    // store this meson and matched jet
    bestpt = TMath::Abs(preselTree->d_pt[d]);
    bestD = d;
    bestJ = j;
  }

  // if no meson found, return 0
  if(bestD < 0)
    return 0;

  // all selection cuts passed: assign output variables and return proper value
  vecJ.SetPtEtaPhiM(TMath::Abs(preselTree->j_pt[bestJ]), preselTree->j_eta[bestJ], preselTree->j_phi[bestJ], preselTree->j_m[bestJ]);
  vecD.SetPtEtaPhiM(TMath::Abs(preselTree->d_pt[bestD]), preselTree->d_eta[bestD], preselTree->d_phi[bestD], preselTree->d_m[bestD]);
  flagInMassWindow = (TMath::Abs(preselTree->d_m[bestD] - 1.870) < 0.05) ? true : false;
  return (preselTree->d_pt[bestD] > 0) ? +1 : -1;
}

// Routine to select muons from charm.
// Arguments:
//   const ZTreeWcharm* preselTree: input tree (see wcharm_tree.h), GetEntry() should be done already,
//   vecJet: jet vector,
//   vecLep: lepton vector (electron or muon produced with W),
//   flagMu: boolean flag which is true if lepton is muon, otherwise false,
//   vecMu: (output) vector of selected muon.
// Returns:
//   +1 if selected mu-,
//   -1 if selected mu+,
//   0 if nothing selected.
// See wcharm_tree.h for ZTree variables description.
int SelectCharmMu(const ZTreeWcharm* preselTree, const TLorentzVector& vecJet, const TLorentzVector& vecLep, const int flagMu, TLorentzVector& vecMu)
{
  int best = -1; // index of selected muon, initialised with impossible value
  double bestDeltaR = 1000.0; // eta-phi separation of selected muon from jet, initialised with impossible large value, because meson with min separation will be selected

  // loop over muons
  for(int mu = 0; mu < preselTree->Nmu; mu++)
  {
    // create muon 4-vector
    vecMu.SetPtEtaPhiM(TMath::Abs(preselTree->mu_pt[mu]), preselTree->mu_eta[mu], preselTree->mu_phi[mu], massMu);
    // calculate eta-phi separation between muon and jet
    double deltaR = vecJet.DeltaR(vecMu);
    // skip current muon, if another muon with smaller deltaR was found
    if(deltaR > bestDeltaR)
      continue;      
    // require fraction of muon pT over jer pT below 0.6
    if((vecMu.Pt() / vecJet.Pt()) > 0.6)
      continue;
    // require relative pT below 2.5 GeV
    if((vecMu.P() * TMath::Sin(vecMu.Angle(vecJet.Vect()))) > 2.5)
      continue;
    // require pT below 25 GeV
    if(TMath::Abs(preselTree->mu_pt[mu]) > 25.0)
      continue;
    // require |eta| below 2.1
    if(TMath::Abs(preselTree->mu_eta[mu]) > 2.1)
      continue;
    // require impact parameter < 0.2
    if(preselTree->mu_imp[mu] > 0.2)
      continue;
    // require muon track number of degrees of greater than 10
    if(preselTree->mu_nChi2[mu] > 10.0)
      continue;
    // require at least 11 tracker hits and at least 1 pixel hit
    if(preselTree->mu_ValidHits[mu] < 12 || preselTree->mu_PixelHits[mu] < 2)
      continue;

    // calculate dilepton vector
    TLorentzVector dilep = vecLep + vecMu;
    // check that its mass is not below 12 GeV
    if(dilep.M() < 12.0)
      continue;
    // check that its mass is not above 85 GeV, if provided lepton is muon and selected selected muon has opposite sign
    if(flagMu && (preselTree->mu_pt[best] * flagMu) < 0 && dilep.M() > 85.0)
      continue;

    // store index and pT of this muon
    best = mu;
    bestDeltaR = deltaR;
  }

  // check that the best eta-phi separation is not too large (at least one muon passd selection criteria), otherwise return 0
  if(bestDeltaR > 100.0)
    return 0;

  // all selection cuts passed, no mu+mu- pair: assign output vector to selected electron and return proper value
  vecMu.SetPtEtaPhiM(TMath::Abs(preselTree->mu_pt[best]), preselTree->mu_eta[best], preselTree->mu_phi[best], massMu);
  return (preselTree->mu_pt[best] > 0) ? +1 : -1;
}

#endif

#include "tree.h"
#include <TLorentzVector.h>


const double massEl = 0.000511;
const double massMu = 0.105658;

// electron selection
int SelectEl(const ZTree* preselTree, TLorentzVector& vec)
{
  int best = -1;
  double bestpt = 0.0;
  // loop over electrons
  for(int el = 0; el < preselTree->Nel; el++)
  {
    if(TMath::Abs(preselTree->el_pt[el]) < bestpt)
      continue;
    if(TMath::Abs(preselTree->el_pt[el]) < 35.0)
      continue;
    if(TMath::Abs(preselTree->el_eta[el]) > 2.1)
      continue;
    if(TMath::Abs(preselTree->el_eta[el]) > 1.44 && TMath::Abs(preselTree->el_eta[el]) < 1.57)
      continue;
    if(preselTree->el_isolat[el] > 0.05)
      continue;
    if(preselTree->el_mt[el] < 55.0)
      continue;
    if(preselTree->el_mh[el] > 0)
      continue;
    //if(preselTree->elConvDist[el] < 0.02 && preselTree->elConvDcot[el] < 0.02 && preselTree->elConvDist[el] >= 0.0 && preselTree->elConvDcot[el] >=0.0)
    //  continue;
    best = el;
    bestpt = TMath::Abs(preselTree->el_pt[el]);
  }
  if(best < 0)
    return 0;
  for(int el = 0; el < preselTree->Nel; el++)
  {
    if(el == best)
      continue;
    if(preselTree->el_pt[el] * bestpt > 0)
      continue;
    if(TMath::Abs(preselTree->el_pt[el]) < 20.0)
      continue;
    if(TMath::Abs(preselTree->el_eta[el]) > 2.5)
      continue;
    // pair found
    return 0;
  }
  // check that there is no pair
  // all cuts passed: return true
  vec.SetPtEtaPhiM(TMath::Abs(preselTree->el_pt[best]), preselTree->el_eta[best], preselTree->el_phi[best], massEl);
  return (preselTree->el_pt[best] > 0) ? +1 : -1;
}

// muon selection
int SelectMu(const ZTree* preselTree, TLorentzVector& vec)
{
  int best = -1;
  double bestpt = 0.0;
  // loop over muons
  for(int mu = 0; mu < preselTree->Nmu; mu++)
  {
    if(TMath::Abs(preselTree->mu_pt[mu]) < bestpt)
      continue;
    if(TMath::Abs(preselTree->mu_pt[mu]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->mu_eta[mu]) > 2.1)
      continue;
    if(preselTree->mu_isolat[mu] > 0.12)
      continue;
    if(preselTree->mu_mt[mu] < 40.0)
      continue;
    if(preselTree->mu_imp[mu] > 0.2)
      continue;
    if(preselTree->mu_nChi2[mu] > 10.0)
      continue;
    if(preselTree->mu_ValidHits[mu] < 12 || preselTree->mu_PixelHits[mu] < 2)
      continue;
    best = mu;
    bestpt = TMath::Abs(preselTree->mu_pt[mu]);
  }
  if(best < 0)
    return 0;
  for(int mu = 0; mu < preselTree->Nmu; mu++)
  {
    if(mu == best)
      continue;
    if(preselTree->mu_pt[mu] * bestpt > 0)
      continue;
    if(TMath::Abs(preselTree->mu_pt[mu]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->mu_eta[mu]) > 2.4)
      continue;
    // pair found
    return 0;
  }
  // check that there is no pair
  // all cuts passed: return true
  vec.SetPtEtaPhiM(TMath::Abs(preselTree->mu_pt[best]), preselTree->mu_eta[best], preselTree->mu_phi[best], massMu);
  return (preselTree->mu_pt[best] > 0) ? +1 : -1;
}

int SelectDstar(const ZTree* preselTree, TLorentzVector& vecJ, TLorentzVector& vecD, bool& massW, double& deltaM)
{
  int bestJ = -1;
  int bestD = -1;
  double bestpt = 0.0;
  //printf("dupa?\n");
  // loop over dstars
  for(int d = 0; d < preselTree->DsCand; d++)
  {
    if(TMath::Abs(preselTree->ds_pt[d]) < bestpt)
      continue;
    if(preselTree->D0dls[d] < 3.0)
      continue;
    if(preselTree->D0dl[d] > 2.0)
      continue;
    if(TMath::Abs(preselTree->pv_z - preselTree->D0sv_z[d]) > 2.0)
      continue;
    // get jet
    int j = 0;
    int nd = 0;
    for(j = 0; j < preselTree->Njet; j++)
    {
      nd += preselTree->j_DsInJ[j];
      if(nd >= d)
        break;
    }
    if(j == preselTree->Njet)
    {
      printf("Wrong j <---> D finding\n");
      return 0;
    }
    // check jet
    if(TMath::Abs(preselTree->j_pt[j]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->j_eta[j]) > 2.5)
      continue;
    // all passed
    bestpt = TMath::Abs(preselTree->ds_pt[d]);
    bestD = d;
    bestJ = j;
  }
  if(bestD < 0)
    return 0;
  //printf("dupa!\n");
  vecJ.SetPtEtaPhiM(TMath::Abs(preselTree->j_pt[bestJ]), preselTree->j_eta[bestJ], preselTree->j_phi[bestJ], preselTree->j_m[bestJ]);
  vecD.SetPtEtaPhiM(TMath::Abs(preselTree->ds_pt[bestD]), preselTree->ds_eta[bestD], preselTree->ds_phi[bestD], preselTree->ds_m[bestD]);
  deltaM = preselTree->dmDs[bestD];
  massW = (TMath::Abs(deltaM - 0.145) < 0.005) ? true : false;
  return (preselTree->ds_pt[bestD] > 0) ? +1 : -1;
}

int SelectDch(const ZTree* preselTree, TLorentzVector& vecJ, TLorentzVector& vecD, bool& massW)
{
  int bestJ = -1;
  int bestD = -1;
  double bestpt = 0.0;
  //printf("dupa?\n");
  // loop over dstars
  for(int d = 0; d < preselTree->DCand; d++)
  {
    if(TMath::Abs(preselTree->d_pt[d]) < bestpt)
      continue;
    if(preselTree->Ddls[d] < 3.0)
      continue;
    if(preselTree->Ddl[d] > 2.0)
      continue;
    if(TMath::Abs(preselTree->pv_z - preselTree->Dsv_z[d]) > 2.0)
      continue;
    // get jet
    int j = 0;
    int nd = 0;
    for(j = 0; j < preselTree->Njet; j++)
    {
      nd += preselTree->j_DInJ[j];
      if(nd >= d)
        break;
    }
    if(j == preselTree->Njet)
    {
      printf("Wrong j <---> D finding\n");
      return 0;
    }
    // check jet
    if(TMath::Abs(preselTree->j_pt[j]) < 25.0)
      continue;
    if(TMath::Abs(preselTree->j_eta[j]) > 2.5)
      continue;
    // all passed
    bestpt = TMath::Abs(preselTree->d_pt[d]);
    bestD = d;
    bestJ = j;
  }
  if(bestD < 0)
    return 0;
  //printf("dupa!\n");
  vecJ.SetPtEtaPhiM(TMath::Abs(preselTree->j_pt[bestJ]), preselTree->j_eta[bestJ], preselTree->j_phi[bestJ], preselTree->j_m[bestJ]);
  vecD.SetPtEtaPhiM(TMath::Abs(preselTree->d_pt[bestD]), preselTree->d_eta[bestD], preselTree->d_phi[bestD], preselTree->d_m[bestD]);
  massW = (TMath::Abs(preselTree->d_m[bestD] - 1.870) < 0.05) ? true : false;
  return (preselTree->d_pt[bestD] > 0) ? +1 : -1;
}

// muon selection
int SelectCharmMu(const ZTree* preselTree, const TLorentzVector& vecJet, const TLorentzVector& vecLep, TLorentzVector& vecMu, const int signMu = 0)
{
  int best = -1;
  double bestdeltar = 1000.0;
  // loop over muons
  for(int mu = 0; mu < preselTree->Nmu; mu++)
  {
    vecMu.SetPtEtaPhiM(TMath::Abs(preselTree->mu_pt[mu]), preselTree->mu_eta[mu], preselTree->mu_phi[mu], massMu);
    double deltar = vecJet.DeltaR(vecMu);
    if(deltar > bestdeltar)
      continue;      
    if((vecMu.Pt() / vecJet.Pt()) > 0.6)
      continue;
    if((vecMu.P() * TMath::Sin(vecMu.Angle(vecJet.Vect()))) > 2.5)
      continue;
    if(TMath::Abs(preselTree->mu_pt[mu]) > 25.0)
      continue;
    if(TMath::Abs(preselTree->mu_eta[mu]) > 2.1)
      continue;
    //if(preselTree->mu_mt[mu] < 40.0)
    //  continue;
    if(preselTree->mu_imp[mu] > 0.2)
      continue;
    if(preselTree->mu_nChi2[mu] > 10.0)
      continue;
    if(preselTree->mu_ValidHits[mu] < 12 || preselTree->mu_PixelHits[mu] < 2)
      continue;
    TLorentzVector dilep = vecLep + vecMu;
    if(dilep.M() < 12.0)
      continue;
    if(signMu && (preselTree->mu_pt[best] * signMu) < 0 && dilep.M() > 85.0)
      continue;
    // all passed
    best = mu;
    bestdeltar = deltar;
  }
  if(bestdeltar > 100.0)
    return 0;
  // all cuts passed: return true
  vecMu.SetPtEtaPhiM(TMath::Abs(preselTree->mu_pt[best]), preselTree->mu_eta[best], preselTree->mu_phi[best], massMu);
  return (preselTree->mu_pt[best] > 0) ? +1 : -1;
}

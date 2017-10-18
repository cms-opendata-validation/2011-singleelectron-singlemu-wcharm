// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>> Helper for W+c event reconstruction >>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifndef WCHARM_EVENTRECO_H
#define WCHARM_EVENTRECO_H

// additional files from this analysis
#include "wcharm_tree.h"
#include "wcharm_settings.h"
#include "wcharm_selection.h"
// C++ libraries or ROOT header files
#include <map>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>> ZVarHisto class >>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Class for control plot and cross section histograms:
// it stores variable name and histogram.
// A bunch of histograms can be filled using proper ttbar kinematics
// input with just one line (see void FillHistos() below).
// Also see void StoreHistos() for histogram storage.
//
class ZVarHisto
{
  private:
    TH1D* zHisto;
    TString zVar;
    int zCharmFinalState; // charm final state
    bool zMassWindow; //bool flag to determine if mass window is setup when filling histograms, see FillHistos()
    int zSign; // sign of Wc (opposite or same)
  
  public:
    // constructor
    ZVarHisto(const TString& str, TH1D* h, int fs, int sign, bool massWindow)
    {
      zHisto = h;
      zVar = str;
      zCharmFinalState = fs;
      zMassWindow = massWindow;
      zSign = sign;
    }

    // copy constructor
    ZVarHisto(const ZVarHisto& old)
    {
      zHisto = new TH1D(*(old.zHisto));
      zVar = old.zVar;
      zCharmFinalState = old.zCharmFinalState;
      zMassWindow = old.zMassWindow;
      zSign = old.zSign;
    }

    // access histogram
    TH1* H() {return zHisto;}

    // access variable name
    TString V() {return zVar;}

    // access charm final state
    int CharmFinalState() {return zCharmFinalState;}

    // access mass window flag
    bool MassWindow() {return zMassWindow;}

    // access sign of Wc
    int Sign() {return zSign;}
};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>> Fill histogram from ZVarHisto class >>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Arguments:
//   std::vector<ZVarHisto>& VecVarHisto: vector of objects to be filled
//   const double w: weight
//   const int fs: charm final state
//   TLorentzVector* vecLep: lepton momentum
//   TLorentzVector* vecJ: jet momentum
//   TLorentzVector* vecD: meson (or muon from charm) momentum
//   const int sign: sign of Wc
//   const bool massW: flag to determine if mass window should be used
//   const double deltaM: mass window width
//
void FillHistosWCharm(std::vector<ZVarHisto>& VecVarHisto, const double w, const int fs, TLorentzVector* vecLep, TLorentzVector* vecJ, TLorentzVector* vecD, const int sign, const bool massW, const double deltaM)
{
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    // check final state
    if(VecVarHisto[h].CharmFinalState() != fs)
      continue;
    // check mass window
    if(VecVarHisto[h].MassWindow() && !massW)
      continue;
    // check sign
    if(VecVarHisto[h].Sign() && VecVarHisto[h].Sign() != sign)
      continue;

    // retrieve variable name
    TString var = VecVarHisto[h].V();
    // retrieve histogram
    TH1* histo = VecVarHisto[h].H();
    // now fill histograms depending on the variable:
    // jet eta
    if(var == "etaj")
      histo->Fill(TMath::Abs(vecJ->Eta()), w);
    // momentum fraction
    else if(var == "z")
      histo->Fill(TMath::Abs(vecD->P() / vecJ->P()), w);
    // D+ mass
    else if(var == "md")
      histo->Fill(vecD->M(), w);
    // mass difference M(D*) - M(D0)
    else if(var == "dmd")
      histo->Fill(deltaM, w);
    // histogram with one bin used for total number of events at generator level
    else if(var == "gen") 
      histo->Fill(0.0, w);
    // uknown (not implemented) variable
    // (you can implement more variables here if needed)
    else
    {
      //printf("Error: unknown variable %s\n", var.Data());
      //exit(1);
      continue;
    } // end variable for histo
  } // end loop over histos
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>> Store hostogram from ZVarHisto class >>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Store bunch of histogram (argument std::vector<ZVarHisto>& VecVarHisto)
// for specified final state (stores a copy of histogram)
//
// Arguments:
//   const int fs: charm final state
//
void StoreHistos(std::vector<ZVarHisto>& VecVarHisto, const int fs)
{
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    // check final state
    if(VecVarHisto[h].CharmFinalState() != fs)
      continue;
    TH1* histo = VecVarHisto[h].H();
    TString name = histo->GetName();
    TString title = histo->GetTitle();
    histo->SetNameTitle(name, title);
    histo->Write();
  }
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>> Input parameters for eventreco routine (see below)  >>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class ZEventWcharmRecoInput
{
  public:
    TString Name; // name pattern (to be used in output histograms)
    std::vector<ZVarHisto> VecVarHisto; // container with needed histograms
    int WChannel; // 1 mu, 2 e
    int CharmFinalState; // 1 D*, 2 D+, 3 mu
    int Type; // 1 data, 2 MC signal, 3 MC other, 4 MC background (not used)
    bool Gen; // if true, the histogram is filled at true level
    std::vector<TString> VecInFile; // container with input files
    double Weight; // weight for histogram filling
    long MaxNEvents; // maximum number of processed events

    // contstructor
    ZEventWcharmRecoInput()
    {
      // set default values
      Weight = 1.0;
      MaxNEvents = 100e10;
      Gen = false;
    }
    
    // add one more input file (str) to the chain
    void AddToChain(const TString& str)
    {
      VecInFile.push_back(str);
    }
};

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>> Basic routine for ttbar event reconstruction >>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void wcharm_eventreco(ZEventWcharmRecoInput in)
{ 
  printf("****** EVENTRECO ******\n");
  printf("input sample: %s (type %d)\n", in.Name.Data(), in.Type);
  printf("W decaychannel: %d,  charm final state: %d\n", in.WChannel, in.CharmFinalState);
  
  // directory for output ROOT files with histograms
  TString outDir = gHistDir;

  // this flag determines whether generator level information is available
  // (should be available for signal MC)
  bool flagMC = (in.Type == 2 || in.Type == 3);
  
  // output file
  TFile* fout = TFile::Open(TString::Format("%s/%s-c%d-f%d.root", outDir.Data(), in.Name.Data(), in.WChannel, in.CharmFinalState), "recreate");
  
  // input tree
  TChain* chain = new TChain("tree");
  for(int f = 0; f < in.VecInFile.size(); f++)
    chain->Add(in.VecInFile[f]);
  ZTreeWcharm* preselTree = new ZTreeWcharm(flagMC);
  preselTree->Init(chain);

  // disable unneeded branches for faster processing
  chain->SetBranchStatus("metPx", 0);
  chain->SetBranchStatus("metPy", 0);
  chain->SetBranchStatus("evEventNumber", 0);
  chain->SetBranchStatus("evRunNumber", 0);
  chain->SetBranchStatus("el_convDist", 0);
  chain->SetBranchStatus("el_convDcot", 0);
  chain->SetBranchStatus("pv_x", 0);
  chain->SetBranchStatus("pv_y", 0);
  chain->SetBranchStatus("vcx", 0);
  chain->SetBranchStatus("vcy", 0);
  chain->SetBranchStatus("vcz", 0);
  chain->SetBranchStatus("D0dlEr", 0);
  chain->SetBranchStatus("D0sv_x", 0);
  chain->SetBranchStatus("D0sv_y", 0);
  chain->SetBranchStatus("DdlEr", 0);
  chain->SetBranchStatus("Dsv_x", 0);
  chain->SetBranchStatus("Dsv_y", 0);
  // muons can be disabled for W->e channel only if charm final state is not meson
  if(in.WChannel == 2 && in.CharmFinalState != 3)
    chain->SetBranchStatus("mu_*", 0);
  // electrons can be disabled for W->mu channel and any charm final state
  if(in.WChannel == 1)
    chain->SetBranchStatus("el_*", 0);
  // D+ branches can be disabled for c->D+ final state
  if(in.CharmFinalState == 1)
  {
    chain->SetBranchStatus("d_*", 0);
  }
  // D* branches can be disabled for c->D* final state
  if(in.CharmFinalState == 2)
  {
    chain->SetBranchStatus("ds_*", 0);
    chain->SetBranchStatus("D0*", 0);
  }
  // D+ and D* branches can be disabled for c->mu final state
  if(in.CharmFinalState == 3)
  {
    chain->SetBranchStatus("d_*", 0);
    chain->SetBranchStatus("D*", 0);
    chain->SetBranchStatus("ds_*", 0);
  }
  // enable generator level branches only, if generator level should be processed
  if(in.Gen)
  {
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("mcEventType", 1);
    chain->SetBranchStatus("cGen", 1);
  }
    
  // event counters
  long nSel = 0;
  long nGen = 0;

  // determine number of events
  long nEvents = chain->GetEntries();
  if(nEvents > in.MaxNEvents)
    nEvents = in.MaxNEvents;
  printf("nEvents: %ld\n", nEvents);
  // event loop
  for(int e = 0; e < nEvents; e++)
  {
    chain->GetEntry(e);
    if(flagMC)
    {
      // type of event at generator level (described in ZTreeWcharm)
      int ev = TMath::Abs(preselTree->mcEventType);
      // positive mcEventType values for W decaying to muons, negative values for W decaying to electrons
      int sign = (preselTree->mcEventType > 0) ? +1 : -1;

      // process MC signal
      if(in.Type == 2)
      {
        // check W decay
        if(in.WChannel == 1 && sign == -1) continue;
        if(in.WChannel == 2 && sign == +1) continue;
        // check charm final state
        if (!ev || (ev % 3) != (in.CharmFinalState % 3) || (preselTree->cGen % 2) == 0) continue;
      }

      // process MC other than signal
      if(in.Type == 3)
      {
        // check W decay and charm final state
        if(((in.WChannel == 1 && sign == +1) || (in.WChannel == 2 && sign == -1)) && (ev && (ev % 3) == (in.CharmFinalState % 3)) && (preselTree->cGen % 2) != 0) continue;
      }

      // if we are here, than this event is accepted at generator level
      nGen++;
    }

    // analyse generator level
    if(in.Gen)
    {
      // fill histograms
      FillHistosWCharm(in.VecVarHisto, in.Weight, in.CharmFinalState, NULL, NULL, NULL, 0, 0, 0);
      continue;
    }
    
    // analyse reco level
    
    // vectors for lepton, jet and charm meson (or muon from charm)
    TLorentzVector vecLep, vecJ, vecD;
    // lepton sign (here "lepton" refers to the lepton from W decay, unless explicitely stated otherwise)
    int signLep = 0;
    // *****************************************
    // ***************** mu ********************
    // *****************************************
    if(in.WChannel == 1)
    {
      // trigger: 0th bit
      bool trig = ((preselTree->Triggers >> 0) & 1);
      if(!trig) continue;
      signLep = SelectMu(preselTree, vecLep);
      if(!signLep) continue;
    }
    // *****************************************
    // ****************** e ********************
    // *****************************************
    else if(in.WChannel == 2)
    {
      // trigger: 1st or 2nd bits
      bool trig = ((preselTree->Triggers >> 1) & 1) || ((preselTree->Triggers >> 2) & 1);
      if(!trig) continue;
      // call electron selection routine (see wcharm_selection.h for description)
      signLep = SelectEl(preselTree, vecLep);
      if(!signLep) continue;
    }
    // charm selection
    // sign of D meson (or muon from charm)
    int signD = 0;
    // mass window flag will be used for D*
    bool massWindow = false;
    // mass window width will be used for D*
    double deltaM = 0.0;
    // *****************************************
    // ***************** D* ********************
    // *****************************************
    if(in.CharmFinalState == 1)
    {
      // call D* selection routine (see wcharm_selection.h for description)
      signD = SelectDstar(preselTree, vecJ, vecD, massWindow, deltaM);
      if(!signD) continue;
    }
    // *****************************************
    // ***************** D+ ********************
    // *****************************************
    else if(in.CharmFinalState == 2)
    {
      // call D+ selection routine (see wcharm_selection.h for description)
      signD = SelectDch(preselTree, vecJ, vecD, massWindow);
      if(!signD) continue;
    }
    // *****************************************
    // *************** c->mu *******************
    // *****************************************
    else if(in.CharmFinalState == 3)
    {
      // find highest pT jet
      double bestpt = -1.0;
      for(int j = 0; j < preselTree->Njet; j++)
      {
        if(TMath::Abs(preselTree->j_pt[j]) < bestpt)
          continue;
        // check jet pT and eta
        if(TMath::Abs(preselTree->j_pt[j]) < 25.0)
          continue;
        if(TMath::Abs(preselTree->j_eta[j]) > 2.5)
          continue;
        // create jet vector
        vecJ.SetPtEtaPhiM(TMath::Abs(preselTree->j_pt[j]), preselTree->j_eta[j], preselTree->j_phi[j], preselTree->j_m[j]);
        // create vector used for selected muon from charm
        TLorentzVector vecMu;
        // call muon form charm selection routine (see wcharm_selection.h for description)
        int sign = SelectCharmMu(preselTree, vecJ, vecLep, (in.WChannel == 1), vecMu);
        if(!sign) 
          continue;
        // store pT, vector and sign of this jet
        bestpt = TMath::Abs(preselTree->j_pt[j]);
        vecD = vecMu;
        signD = sign;
      }

      // if no muon from charm which satisfy all criteria, skip this event
      if(bestpt < 0.0) continue;
    }

    // if we are here, than this decay is accepted at reco level
    nSel++;
    
    // fill histograms
    int sign = signLep * signD; // sign of Wc
    FillHistosWCharm(in.VecVarHisto, in.Weight, in.CharmFinalState, &vecLep, &vecJ, &vecD, sign, massWindow, deltaM);
  }

  // print numbers of generated (for MC signal) and selected events
  printf("nSel  : %ld\n", nSel);
  if(in.Type == 2) 
  {
    printf("nGen  : %ld\n", nGen);
    // print acceptance (also referred to as efficiency correction factor)
    printf("C = %.2f%%\n", 100. * nSel / nGen);
  }

  // store histograms, close output file
  fout->cd();
  StoreHistos(in.VecVarHisto, in.CharmFinalState);
  fout->Close();
}

#endif

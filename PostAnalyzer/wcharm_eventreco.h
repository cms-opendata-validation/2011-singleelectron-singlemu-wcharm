#ifndef WCHARM_EVENTRECO_H
#define WCHARM_EVENTRECO_H

#include "wcharm_tree.h"
#include "wcharm_selection.h"
#include <map>
#include <TChain.h>
#include <mystaff.cxx>
#include <TCanvas.h>
#include <TFile.h>

class ZVarHisto
{
  private:
    TH1D* zHisto;
    TString zVar;
    int zFS;
    bool zMassWindow;
    int zSign;
  
  public:
    ZVarHisto(const TString& str, TH1D* h, int fs, int sign, bool massWindow)
    {
      zHisto = h;
      zVar = str;
      zFS = fs;
      zMassWindow = massWindow;
      zSign = sign;
    }

    ZVarHisto(const ZVarHisto& old)
    {
      zHisto = new TH1D(*(old.zHisto));
      zVar = old.zVar;
      zFS = old.zFS;
      zMassWindow = old.zMassWindow;
      zSign = old.zSign;
    }

    TH1* H() {return zHisto;}
    TString V() {return zVar;}
    int FS() {return zFS;}
    bool MW() {return zMassWindow;}
    int Sign() {return zSign;}
};

void FillHistos(std::vector<ZVarHisto>& VecVarHisto, double w, const int fs, TLorentzVector* vecLep, TLorentzVector* vecJ, TLorentzVector* vecD, const int sign, const bool massW, const double deltaM)
{
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    // check final state
    if(VecVarHisto[h].FS() != fs)
      continue;
    // check mass window
    if(VecVarHisto[h].MW() && !massW)
      continue;
    // check sign
    if(VecVarHisto[h].Sign() && VecVarHisto[h].Sign() != sign)
      continue;
    // variable
    TString var = VecVarHisto[h].V();
    TH1* histo = VecVarHisto[h].H();
    //printf("var = %s  histo = %s\n", var.Data(), histo->GetName());
    if(var == "etaj") 
      histo->Fill(TMath::Abs(vecJ->Eta()), w);
    else if(var == "z") 
      histo->Fill(TMath::Abs(vecD->P() / vecJ->P()), w);
    else if(var == "md") 
      histo->Fill(vecD->M(), w);
    else if(var == "dmd") 
      histo->Fill(deltaM, w);
    else if(var == "gen") 
      histo->Fill(0.0, w);
    else
    {
      //printf("Error: unknown variable %s\n", var.Data());
      //exit(1);
      continue;
    } // end variable for histo
  } // end loop over histos
}

void StoreHistos(std::vector<ZVarHisto>& VecVarHisto, const int fs)
{
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    // check final state
    if(VecVarHisto[h].FS() != fs)
      continue;
    //printf("h = %d\n", h);
    TH1* histo = VecVarHisto[h].H();
    //TString name = TString::Format("%s_t%d_ch%d+%s\n", histo->GetName(), in.Type, in.Channel, in.Name.Data());
    //TString title = TString::Format("%s-t%d-ch%d-%s\n", histo->GetTitle(), in.Type, in.Channel, in.Name.Data());
    TString name = histo->GetName();
    TString title = histo->GetTitle();
    histo->SetNameTitle(name, title);
    histo->Write();
  }
}

class ZEventRecoInput
{
  public:
    TString Name;
    std::vector<ZVarHisto> VecVarHisto;
    int Channel; // 1 mu, 2 e
    int FinalState; // 1 dstar, 2 dch, 3 mu
    int Type; // 1 data, 2 MC signal, 3 MC other, 4 MC background
    bool Gen;
    std::vector<TString> VecInFile;
    double Weight;
    long MaxNEvents;
    
    ZEventRecoInput()
    {
      Weight = 1.0;
      MaxNEvents = 100e10;
      Gen = false;
    }
    
    void AddToChain(const TString& str)
    {
      VecInFile.push_back(str);
    }
};

void eventreco(ZEventRecoInput in)
{ 
  printf("****** EVENTRECO ******\n");
  printf("input sample: %s (type %d)\n", in.Name.Data(), in.Type);
  printf("channel: %d  finalstate: %d\n", in.Channel, in.FinalState);
  
  // some steerings
  TString outDir = "/nfs/dust/cms/user/zenaiev/opendata/wcharm/outEventreco";

  // some flags
  bool flagMC = (in.Type > 1);
  
  // output file
  TFile* fout = TFile::Open(TString::Format("%s/%s-c%d-f%d.root", outDir.Data(), in.Name.Data(), in.Channel, in.FinalState), "recreate");
  
  // input tree
  TChain* chain = new TChain("tree");
  for(int f = 0; f < in.VecInFile.size(); f++)
    chain->Add(in.VecInFile[f]);
  ZTreeWcharm* preselTree = new ZTreeWcharm(flagMC);
  preselTree->Init(chain);
  // disable unneeded
  chain->SetBranchStatus("metPx", 0);
  chain->SetBranchStatus("metPy", 0);
  chain->SetBranchStatus("Nevent", 0);
  chain->SetBranchStatus("Nrun", 0);
  chain->SetBranchStatus("lumiblock", 0);
  chain->SetBranchStatus("el_convDist", 0);
  chain->SetBranchStatus("el_convDcot", 0);
  chain->SetBranchStatus("mu_stations", 0);
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
  if(in.Channel == 2 && in.FinalState != 3)
    chain->SetBranchStatus("mu_*", 0);
  if(in.Channel == 1)
    chain->SetBranchStatus("el_*", 0);
  if(in.FinalState == 1)
  {
    chain->SetBranchStatus("d_*", 0);
  }
  if(in.FinalState == 2)
  {
    chain->SetBranchStatus("ds_*", 0);
    chain->SetBranchStatus("D0*", 0);
  }
  if(in.FinalState == 3)
  {
    chain->SetBranchStatus("d_*", 0);
    chain->SetBranchStatus("D*", 0);
    chain->SetBranchStatus("ds_*", 0);
  }
  // gen level
  if(in.Gen)
  {
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("mcEventType", 1);
    chain->SetBranchStatus("CGen", 1);
  }
    
  // event count
  long nSel = 0;
  long nGen = 0;

  // event loop
  long nEvents = chain->GetEntries();
  if(nEvents > in.MaxNEvents)
    nEvents = in.MaxNEvents;
  printf("nEvents: %ld\n", nEvents);
  for(int e = 0; e < nEvents; e++)
  {
    //printf("e: %d\n", e);
    chain->GetEntry(e);
    //continue;
    //printf("%d %d %d %d %d\n", preselTree->Nel, preselTree->Nmu, preselTree->Njet, preselTree->DsCand, preselTree->DCand);
    if(flagMC)
    {
      //printf("%d %d\n", preselTree->mcEventType, preselTree->CGen);
      //preselTree->b_mcEventType->GetEntry(e);
      //printf("%d", preselTree->mcEventType);
      int ev = TMath::Abs(preselTree->mcEventType);
      int sign = (preselTree->mcEventType > 0) ? +1 : -1;
      // MC signal?
      if(in.Type == 2)
      {
        if(in.Channel == 1 && sign == -1) continue;
        if(in.Channel == 2 && sign == +1) continue;
        if (!ev || (ev % 3) != (in.FinalState % 3) || (preselTree->cGen % 2) == 0) continue;
      }
      // MC other?
      if(in.Type == 3)
      {
        if(((in.Channel == 1 && sign == +1) || (in.Channel == 2 && sign == -1)) && (ev && (ev % 3) == (in.FinalState % 3)) && (preselTree->cGen % 2) != 0) continue;
      }
      nGen++;
    }
    // gen level
    if(in.Gen)
    {
    // fill histos
      FillHistos(in.VecVarHisto, in.Weight, in.FinalState, NULL, NULL, NULL, 0, 0, 0);
      //printf("ev: %d\n", preselTree->mcEventType);
      continue;
    }
    
    // continue reco level
    // load event
    //chain->GetEntry(e);
    
    // primary vertex selection
    //if(preselTree->Npv < 1 || preselTree->pvNDOF < 4 || preselTree->pvRho > 2.0 || TMath::Abs(preselTree->pvZ) > 24.0)
    //  continue;
    // select lepton
    TLorentzVector vecLep, vecJ, vecD;
    int signLep = 0;
    // *****************************************
    // ***************** mu ********************
    // *****************************************
    if(in.Channel == 1)
    {
      // trigger: 0th bit
      bool trig = ((preselTree->Triggers >> 0) & 1);
      //trig = true;
      if(!trig) continue;
      signLep = SelectMu(preselTree, vecLep);
      if(!signLep) continue;
    }
    // *****************************************
    // ****************** e ********************
    // *****************************************
    else if(in.Channel == 2)
    {
      // trigger: 1st or 2nd bits
      bool trig = ((preselTree->Triggers >> 1) & 1) || ((preselTree->Triggers >> 2) & 1);
      //trig = true;
      if(!trig) continue;
      signLep = SelectEl(preselTree, vecLep);
      if(!signLep) continue;
    }
    //continue;
    // charm selection
    int signD = 0;
    bool massWindow = false;
    double deltaM = 0.0;
    // *****************************************
    // ***************** D* ********************
    // *****************************************
    if(in.FinalState == 1)
    {
      signD = SelectDstar(preselTree, vecJ, vecD, massWindow, deltaM);
      if(!signD) continue;
    }
    // *****************************************
    // ***************** D+ ********************
    // *****************************************
    else if(in.FinalState == 2)
    {
      signD = SelectDch(preselTree, vecJ, vecD, massWindow);
      if(!signD) continue;
    }
    // *****************************************
    // *************** c->mu *******************
    // *****************************************
    else if(in.FinalState == 3)
    {
      double bestpt = -1.0;
      for(int j = 0; j < preselTree->Njet; j++)
      {
        if(TMath::Abs(preselTree->j_pt[j]) < bestpt)
          continue;
        // check jet
        if(TMath::Abs(preselTree->j_pt[j]) < 25.0)
          continue;
        if(TMath::Abs(preselTree->j_eta[j]) > 2.5)
          continue;
        vecJ.SetPtEtaPhiM(TMath::Abs(preselTree->j_pt[j]), preselTree->j_eta[j], preselTree->j_phi[j], preselTree->j_m[j]);
        TLorentzVector vecMu;
        int sign = SelectCharmMu(preselTree, vecJ, vecLep, (in.Channel == 1), vecMu);
        if(!sign) 
          continue;
        bestpt = TMath::Abs(preselTree->j_pt[j]);
        vecD = vecMu;
        signD = sign;
      }
      //signD = SelectCMu(preselTree, vecJ, vecD) continue;
      if(bestpt < 0.0) continue;
    }

    // reco level selected
    nSel++;
    //if(nSel % 100 == 0) printf("nSel: %d\n", nSel);
    
    // fill histos
    int sign = signLep * signD;
    FillHistos(in.VecVarHisto, in.Weight, in.FinalState, &vecLep, &vecJ, &vecD, sign, massWindow, deltaM);
  }
  // print out
  printf("nSel  : %ld\n", nSel);
  if(in.Type == 2) 
  {
    printf("nGen  : %ld\n", nGen);
    printf("C = %.2f%%\n", 100. * nSel / nGen);
  }

  // store histos
  fout->cd();
  StoreHistos(in.VecVarHisto, in.FinalState);
  fout->Close();
}

#endif

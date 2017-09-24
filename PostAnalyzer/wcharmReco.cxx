#include "eventreco_wcharm.h"


TString baseDir = "/nfs/dust/cms/user/zenaiev/opendata/wcharm";
TString dataDir = baseDir + "/data";
TString mcDir = baseDir + "/mc";
TString outDir = baseDir + "/outEventreco";

// flags what to run
bool flagData = 1;
bool flagMCsig = 1;
bool flagMCdy = 1;

// flag channel
bool flagChMu = 1;
bool flagChEl = 1;

// flag final state
bool flagDstar = 1;
bool flagDch = 1;
bool flagMu = 1;

int main(int argc, char** argv)
{
  // histograms
  TH1::SetDefaultSumw2();
  std::vector<ZVarHisto> vecVH, vecVHGen;
  // control plots
  vecVH.push_back(ZVarHisto("dmd", new TH1D("h_mds_os", "m(D*) - m(D0) OS", 17, 0.136, 0.17), 1, -1, false));
  vecVH.push_back(ZVarHisto("dmd", new TH1D("h_mds_ss", "m(D*) - m(D0) SS", 17, 0.136, 0.17), 1, +1, false));
  vecVH.push_back(ZVarHisto("md", new TH1D("h_mdc_os", "mass(D+) OS", 50, 1.6, 2.2), 2, -1, false));
  vecVH.push_back(ZVarHisto("md", new TH1D("h_mdc_ss", "mass(D+) SS", 50, 1.6, 2.2), 2, +1, false));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_ds_os", "D* |eta(jet)| OS", 6, 0.0, 2.5), 1, -1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_ds_ss", "D* |eta(jet)| SS", 6, 0.0, 2.5), 1, +1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_dc_os", "D+ |eta(jet)| OS", 6, 0.0, 2.5), 2, -1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_dc_ss", "D+ |eta(jet)| SS", 6, 0.0, 2.5), 2, +1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_mu_os", "mu |eta(jet)| OS", 6, 0.0, 2.5), 3, -1, false));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_mu_ss", "mu |eta(jet)| SS", 6, 0.0, 2.5), 3, +1, false));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_ds_os", "D* z OS", 12, 0.0, 1.0), 1, -1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_ds_ss", "D* z SS", 12, 0.0, 1.0), 1, +1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_dc_os", "D+ z OS", 12, 0.0, 1.0), 2, -1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_dc_ss", "D+ z SS", 12, 0.0, 1.0), 2, +1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zmu_cp_mu_os", "mu z OS", 12, 0.0, 0.6), 3, -1, false));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zmu_cp_mu_ss", "mu z SS", 12, 0.0, 0.6), 3, +1, false));
  // generator level
  vecVHGen.push_back(ZVarHisto("gen", new TH1D("h_gen_ds", "GEN", 1, -1.0, 1.0), 1, 0, false));
  vecVHGen.push_back(ZVarHisto("gen", new TH1D("h_gen_dc", "GEN", 1, -1.0, 1.0), 2, 0, false));
  vecVHGen.push_back(ZVarHisto("gen", new TH1D("h_gen_mu", "GEN", 1, -1.0, 1.0), 3, 0, false));
  
  for(int ch = 1; ch <= 2; ch++)
  {
    if(ch == 1 && !flagChMu) continue;
    if(ch == 2 && !flagChEl) continue;
    for(int fs = 1; fs <= 3; fs++)
    {
      if(fs == 1 && !flagDstar) continue;
      if(fs == 2 && !flagDch) continue;
      if(fs == 3 && !flagMu) continue;
      // *****************************************
      // **************** DATA *******************
      // *****************************************
      if(flagData)
      {
        ZEventRecoInput in;
        //in.MaxNEvents = 100000;
        in.Name = "data";
        in.Type = 1;
        in.Channel = ch;
        in.FinalState = fs;
        in.VecVarHisto = vecVH;
        if(ch == 1) // mu
        {
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_10000_FILES-1001-1500/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_10000_FILES-1-500/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_10000_FILES-1501-1882/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_10000_FILES-501-1000/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_10001/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_20000_FILES-1-500/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_20000_FILES-501-1000/*.root");
          in.AddToChain(dataDir + "/SingleMu/wcharmSelected_SingleMu-v1_20001/*.root");
        }
        else if(ch == 2) // e
        {
          in.AddToChain(dataDir + "/SingleEl/wcharmSelected_SingleEl-v1_10000/*.root");
          in.AddToChain(dataDir + "/SingleEl/wcharmSelected_SingleEl-v1_20000_FILES-1-500/*.root");
          in.AddToChain(dataDir + "/SingleEl/wcharmSelected_SingleEl-v1_20000_FILES-501-1032/*.root");
          in.AddToChain(dataDir + "/SingleEl/wcharmSelected_SingleEl-v1_20001/*.root");
        }
        eventreco(in);
      }
      
      // *****************************************
      // ************** MC signal ****************
      // *****************************************
      // events: 72165352
      // sigma: 25431 -> 31314
      // weight: 0.5 * 5000.0 / (72165352. / (25430. * 0.32)) * (31314. / 25431.) = 0.3197
      if(flagMCsig)
      {
        ZEventRecoInput in;
        //in.MaxNEvents = 1e5;
        in.Weight = 0.3197;
        in.Name = "mcSigReco";
        in.Type = 2;
        in.Channel = ch;
        in.FinalState = fs;
        in.VecVarHisto = vecVH;
        in.AddToChain(mcDir + "/W1Jet_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_W1Jet_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00000/*.root");
        in.AddToChain(mcDir + "/W1Jet_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_W1Jet_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_30000/*.root");
        // signal reco level
        eventreco(in);
        // MC other
        in.Name = "mcSigOtherReco";
        in.Type = 3;
        eventreco(in);
        // signal gen level
        in.Name = "mcSigGen";
        in.Type = 2;
        in.VecVarHisto = vecVHGen;
        in.Gen = true;
        eventreco(in);
      }
      // *****************************************
      // **************** MC DY ******************
      // *****************************************
      // Events: 36408225
      // Sigma: 2513 -> 3048
      // weight: 0.5 * 5000. / (36408225. / (2513. * 0.1)) * (3048. / 2513.) = 0.02093
      if(flagMCdy)
      {
        ZEventRecoInput in;
        //in.MaxNEvents = 1000;
        in.Type = 4;
        in.Channel = ch;
        in.FinalState = fs;
        in.VecVarHisto = vecVH;
        in.Name = "mcDYReco";
        in.Weight = 0.02093;
        in.AddToChain(mcDir + "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/*.root");
        eventreco(in);
      }
    }
  }

  return 0;
}

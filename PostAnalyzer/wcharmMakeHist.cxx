// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT ntuples for W+c analysis (see
// Analyzer/src/Analyzer.cc) and produces histograms, which are
// further used to make final plots (see wcharmMakePlots.cxx).
// Run: ./wcharmMakeHist
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// additional files from this analysis (look there for description)
#include "wcharm_eventreco.h"
#include "wcharm_settings.h"
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv)
{
  //
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>> Settings >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // set directories to data and MC ntuples
  TString dataDir = gDataDir;
  TString mcDir = gMcDir;

  // flags what to run
  bool flagData = 1; // if 1, data will be processed
  bool flagMCsig = 1; // if 1, MC will be processed
  /*bool flagMCdy = 1;*/

  // flag channel
  bool flagWtoMu = 1; // if 1, W->mu decays will be processed
  bool flagWtoEl = 1; // if 1, W->e decays will be processed

  // flag final state
  bool flagDstar = 1; // if 1, c->D* will be processed
  bool flagDch = 1; // if 1, c->D+ will be processed
  bool flagMu = 1; // if 1, c->mu will be processed

  // histograms
  TH1::SetDefaultSumw2(); // keep histogram weights by default
  // ZVarHisto is a simple class which incorporates a histogram and a variable name.
  // This class is used to store needed input settings for control plots and cross sections.
  std::vector<ZVarHisto> vecVH, vecVHGen;
  // histograms and variables for control plots.
  // Below OS is opposite sign (W-c, W+cbar), SS is same sign (W+c, W-c).
  // mass difference M(D*) - M(D0)
  vecVH.push_back(ZVarHisto("dmd", new TH1D("h_mds_os", "m(D*) - m(D0) OS", 17, 0.136, 0.17), 1, -1, false));
  vecVH.push_back(ZVarHisto("dmd", new TH1D("h_mds_ss", "m(D*) - m(D0) SS", 17, 0.136, 0.17), 1, +1, false));
  // mass M(D+)
  vecVH.push_back(ZVarHisto("md", new TH1D("h_mdc_os", "mass(D+) OS", 50, 1.6, 2.2), 2, -1, false));
  vecVH.push_back(ZVarHisto("md", new TH1D("h_mdc_ss", "mass(D+) SS", 50, 1.6, 2.2), 2, +1, false));
  // jet eta
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_ds_os", "D* |eta(jet)| OS", 6, 0.0, 2.5), 1, -1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_ds_ss", "D* |eta(jet)| SS", 6, 0.0, 2.5), 1, +1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_dc_os", "D+ |eta(jet)| OS", 6, 0.0, 2.5), 2, -1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_dc_ss", "D+ |eta(jet)| SS", 6, 0.0, 2.5), 2, +1, true));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_mu_os", "mu |eta(jet)| OS", 6, 0.0, 2.5), 3, -1, false));
  vecVH.push_back(ZVarHisto("etaj", new TH1D("h_etaj_cp_mu_ss", "mu |eta(jet)| SS", 6, 0.0, 2.5), 3, +1, false));
  // z fraction (pT of meson or muon over pT of jet)
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_ds_os", "D* z OS", 12, 0.0, 1.0), 1, -1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_ds_ss", "D* z SS", 12, 0.0, 1.0), 1, +1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_dc_os", "D+ z OS", 12, 0.0, 1.0), 2, -1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zD_cp_dc_ss", "D+ z SS", 12, 0.0, 1.0), 2, +1, true));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zmu_cp_mu_os", "mu z OS", 12, 0.0, 0.6), 3, -1, false));
  vecVH.push_back(ZVarHisto("z", new TH1D("h_zmu_cp_mu_ss", "mu z SS", 12, 0.0, 0.6), 3, +1, false));
  // generator level numbers of D*, D+, mu (histograms with one bin)
  vecVHGen.push_back(ZVarHisto("gen", new TH1D("h_gen_ds", "GEN", 1, -1.0, 1.0), 1, 0, false));
  vecVHGen.push_back(ZVarHisto("gen", new TH1D("h_gen_dc", "GEN", 1, -1.0, 1.0), 2, 0, false));
  vecVHGen.push_back(ZVarHisto("gen", new TH1D("h_gen_mu", "GEN", 1, -1.0, 1.0), 3, 0, false));
  
  // loop over final states for W decays (wLep = 1 W->mu, wLep = 2 W->e)
  for(int wLep = 1; wLep <= 2; wLep++)
  {
    // check if we should process this final state
    if(wLep == 1 && !flagWtoMu) continue;
    if(wLep == 2 && !flagWtoEl) continue;

    // loop over final states for charm (charm = 1 D*, charm = 2 D+, charm = 3 for mu)
    for(int charm = 1; charm <= 3; charm++)
    {
      // check if we should process this final state
      if(charm == 1 && !flagDstar) continue;
      if(charm == 2 && !flagDch) continue;
      if(charm == 3 && !flagMu) continue;

      // *****************************************
      // **************** DATA *******************
      // *****************************************
      if(flagData)
      {
        // ZEventWcharmRecoInput is a class for event reconstruction, see its description in wcharm_eventreco.h
        ZEventWcharmRecoInput in;
        //in.MaxNEvents = 1000; // if you need to limit the number of processed events
        in.Name = "data"; // name pattern for output histograms
        in.Type = 1; // type = 1 for data, 2 for MC signal, 3 for MC 'ttbar other', 4 for the rest of MC background samples
        in.Channel = wLep; // W decay channel
        in.FinalState = charm; // charm final state
        in.VecVarHisto = vecVH; // need to copy it, because further will be changed
        // input ROOT ntuples
        if(wLep == 1) // W->mu
          in.AddToChain(dataDir + "/SingleMu/*.root");
        else if(wLep == 2) // W->e
          in.AddToChain(dataDir + "/SingleElectron/*.root");
        // main part: event reconstruction call
        eventreco(in);
      }
      
      // *****************************************
      // ***************** MC ********************
      // *****************************************
      //
      // MC event weights need to be changed to most precise theoretical predictions,
      // the formula is:
      // weight = lumi / (nevents / sigma_MC) * (sigma_theory / sigma_MC) = lumi / nevents * sigma_theory
      //
      // Number of events can be obtained from webpage (see http://opendata.cern.ch/collection/CMS-Simulated-Datasets),
      // but it should be checked that all events have been processed at the Analyzer step (see end of log files)
      //
      // number of events: 72165352
      // MC cross section -> theory: 25431 -> 31314
      // W->lX branching ration = 0.32
      // weight: 2500.0 / (72165352. / (25431. * 0.32)) * (31314. / 25431.) = 0.3471
      if(flagMCsig)
      {
        ZEventWcharmRecoInput in;
        //in.MaxNEvents = 1e5;
        in.Weight = 0.3197;
        in.Name = "mcSigReco";
        in.Type = 2;
        in.Channel = wLep;
        in.FinalState = charm;
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
      /*// *****************************************
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
      }*/
    }
  }

  return 0;
}

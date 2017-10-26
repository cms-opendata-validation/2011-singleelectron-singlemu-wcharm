// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT histograms for W+c analysis,
// (produced by wcharmMakeHist.cxx), and makes final plots and numbers
// (control plots and mass spectra to be compared to SMP-12-002 Fig. 2, 3 and 7,
// and W+c production cross sections to be compared to the ones quoted in SMP-12-002 Table 1).
// Run: ./wcharmMakePlots
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// additional files from this analysis
#include "wcharm_settings.h"
#include "wcharm_plots.h"
// C++ library or ROOT header files
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TMath.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv)
{
  // set user style
  SetStyle();
  
  // name patterns for W decay channels
  TString suf[3] = {"comb", "mu", "e"};
  
  // MC samples for control plots with their names, weight factors and cosmetic details for control plots
  // names are stored in MCSamples.VecMCName
  // each sample can contain multiple subsamples (container of containers)
  ZControlPlotMCSamples MCSamples;
  // W+jets other
  MCSamples.VecMCName.push_back(std::vector<TString>(1, "SigOther"));
  MCSamples.VecMCColor.push_back(68);
  MCSamples.VecMCtitle.push_back("W+jets other");
  MCSamples.VecMCFactor.push_back(1.0);
  // W+jets signal
  MCSamples.VecMCName.push_back(std::vector<TString>(1, "Sig"));
  MCSamples.VecMCColor.push_back(65);
  MCSamples.VecMCtitle.push_back("W+c");
  MCSamples.VecMCFactor.push_back(1.0);

  // produce control plots
  // vector with 2D frame histograms for each variable (6 variables as in SMP-12-002 Fig. 7)
  std::vector<TH2F*> cpHR;
  // vector with variable names (same size as cpHR)
  std::vector<TString> cpVar;
  // vector with charm final-state identifiers
  std::vector<int> cpFS;
  // jet pseudorapidity, D+ final state
  TH2F* hr_cp_etaj_dp = new TH2F("hr_cp_etaj_dc", "", 1, 0, 2.5, 1, 0, 1700);
  hr_cp_etaj_dp->GetXaxis()->SetTitle("|#eta^{jet}|");
  hr_cp_etaj_dp->GetYaxis()->SetTitle("(OS-SS) events / 0.42");
  SetCPHRange(hr_cp_etaj_dp);
  cpHR.push_back(hr_cp_etaj_dp);
  cpVar.push_back("etaj_cp_dc");
  cpFS.push_back(2);
  // momentum fraction (z), D+ final state
  TH2F* hr_cp_z_dp = new TH2F("hr_cp_z_dc", "", 1, 0, 1.0, 1, 0, 1200);
  hr_cp_z_dp->GetXaxis()->SetTitle("p^{D}/p^{jet}");
  hr_cp_z_dp->GetYaxis()->SetTitle("(OS-SS) events / 0.083");
  SetCPHRange(hr_cp_z_dp);
  cpHR.push_back(hr_cp_z_dp);
  cpVar.push_back("zD_cp_dc");
  cpFS.push_back(2);
  // jet pseudorapidity, D* final state
  TH2F* hr_cp_etaj_ds = new TH2F("hr_cp_etaj_ds", "", 1, 0, 2.5, 1, 0, 350);
  hr_cp_etaj_ds->GetXaxis()->SetTitle("|#eta^{jet}|");
  hr_cp_etaj_ds->GetYaxis()->SetTitle("(OS-SS) events / 0.42");
  SetCPHRange(hr_cp_etaj_ds);
  cpHR.push_back(hr_cp_etaj_ds);
  cpVar.push_back("etaj_cp_ds");
  cpFS.push_back(1);
  // momentum fraction (z), D* final state
  TH2F* hr_cp_z_ds = new TH2F("hr_cp_z_ds", "", 1, 0, 1.0, 1, 0, 200);
  hr_cp_z_ds->GetXaxis()->SetTitle("p^{D}/p^{jet}");
  hr_cp_z_ds->GetYaxis()->SetTitle("(OS-SS) events / 0.083");
  SetCPHRange(hr_cp_z_ds);
  cpHR.push_back(hr_cp_z_ds);
  cpVar.push_back("zD_cp_ds");
  cpFS.push_back(1);
  // jet pseudorapidity, muon from charm final state
  TH2F* hr_cp_etaj_mu = new TH2F("hr_cp_etaj_mu", "", 1, 0, 2.5, 1, 0, 7000);
  hr_cp_etaj_mu->GetXaxis()->SetTitle("|#eta^{jet}|");
  hr_cp_etaj_mu->GetYaxis()->SetTitle("(OS-SS) events / 0.42");
  SetCPHRange(hr_cp_etaj_mu);
  cpHR.push_back(hr_cp_etaj_mu);
  cpVar.push_back("etaj_cp_mu");
  cpFS.push_back(3);
  // momentum fraction (z), muon from charm final state
  TH2F* hr_cp_z_mu = new TH2F("hr_cp_z_mu", "", 1, 0, 0.6, 1, 0, 5500);
  hr_cp_z_mu->GetXaxis()->SetTitle("p^{D}/p^{jet}");
  hr_cp_z_mu->GetYaxis()->SetTitle("(OS-SS) events / 0.083");
  SetCPHRange(hr_cp_z_mu);
  cpHR.push_back(hr_cp_z_mu);
  cpVar.push_back("zmu_cp_mu");
  cpFS.push_back(3);
  // produce figure similar to Fig. 7 in SMP-12-002, see routine description in wcharm_plot.h
  FigureControlPlots(MCSamples, cpHR, cpVar, cpFS);

  // produce invariant mass distribution plots for D+ and D*
  // vector with 2D frame histograms for each variable mass distribution
  std::vector<TH2F*> mHR;
  // vector with variable names (same size as mHR)
  std::vector<TString> mVar;
  // M(D+) Fig. 2 in SMP-12-002
  // reference frame histogram
  TH2F* hFrameMassDch = new TH2F("hm_mass_dc", "", 1, 1.6, 2.2, 1, 0, 1000);
  hFrameMassDch->GetXaxis()->SetTitle("M(K#pi#pi) [GeV]");
  hFrameMassDch->GetYaxis()->SetTitle("(OS-SS) events / 0.012 GeV");
  SetCPHRange(hFrameMassDch);
  // M(D*)-M(D0) Fig. 3 in SMP-12-002
  // reference frame histogram
  TH2F* hr_mass_ds = new TH2F("hm_mass_ds", "", 1, 0.135, 0.170, 1, 0, 650);
  hr_mass_ds->GetXaxis()->SetTitle("M(K#pi#pi_{s})-M(K#pi) [GeV]");
  hr_mass_ds->GetYaxis()->SetTitle("(OS-SS) events / 0.002 GeV");
  SetCPHRange(hr_mass_ds);
  mHR.push_back(hr_mass_ds);
  mVar.push_back("mds");

  // produce plot similar to SMP-12-002 paper Fig. 2 (D+)
  // obtain event yields needed for cross-section calculation
  std::vector<ZYieldWcharm> yieldsDch = FigureMassDch(hFrameMassDch, MCSamples);
  // store final state identifier needed for cross-section calculation
  for(int i = 0; i < yieldsDch.size(); i++)
    yieldsDch[i].CharmFinalState = 2;

  // produce plot similar to SMP-12-002 paper Fig. 3 (D*+)
  // obtain event yields needed for cross-section calculation
  std::vector<ZYieldWcharm> yieldsDstar = FigureMassDstar(hr_mass_ds, MCSamples);
  // store final state identifier needed for cross-section calculation
  for(int i = 0; i < yieldsDstar.size(); i++)
    yieldsDstar[i].CharmFinalState = 1;

  // now calculate cross section to be compared to Tab. 1 SMP-12-002
  // constants
  double lumi = 2500.0; //luminosity
  double brDs = 0.00622; // branching ratio Br(D*->D0pi)*Br(D0->Kpi)
  double brDc = 0.0208; // branching ratio (D+->Kpipi)
  // loop over W decay channels
  for(int i = 0; i < 2; i++)
  {
    // D+
    CalculateCrossSection(yieldsDch[i], 1.0 / lumi / brDc);
    // D*
    CalculateCrossSection(yieldsDstar[i], 1.0 / lumi / brDs);
  }

  return 0;
}

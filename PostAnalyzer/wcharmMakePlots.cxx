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
  TH2F* hr_cp_etaj_mu = new TH2F("hr_cp_etaj_mu", "", 1, 0, 2.5, 1, 0, 6000);
  hr_cp_etaj_mu->GetXaxis()->SetTitle("|#eta^{jet}|");
  hr_cp_etaj_mu->GetYaxis()->SetTitle("(OS-SS) events / 0.42");
  SetCPHRange(hr_cp_etaj_mu);
  cpHR.push_back(hr_cp_etaj_mu);
  cpVar.push_back("etaj_cp_mu");
  cpFS.push_back(3);
  // momentum fraction (z), muon from charm final state
  TH2F* hr_cp_z_mu = new TH2F("hr_cp_z_mu", "", 1, 0, 0.6, 1, 0, 4500);
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
  // vector with charm final-state identifiers
  std::vector<int> mFS;
  // M(D+)
  TH2F* hr_mass_dc = new TH2F("hm_mass_dc", "", 1, 1.6, 2.2, 1, 0, 1000);
  hr_mass_dc->GetXaxis()->SetTitle("M(K#pi#pi) [GeV]");
  hr_mass_dc->GetYaxis()->SetTitle("(OS-SS) events / 0.012 GeV");
  SetCPHRange(hr_mass_dc);
  // M(D*)-M(D0)
  TH2F* hr_mass_ds = new TH2F("hm_mass_ds", "", 1, 0.135, 0.170, 1, 0, 650);
  hr_mass_ds->GetXaxis()->SetTitle("M(K#pi#pi_{s})-M(K#pi) [GeV]");
  hr_mass_ds->GetYaxis()->SetTitle("(OS-SS) events / 0.002 GeV");
  SetCPHRange(hr_mass_ds);
  mHR.push_back(hr_mass_ds);
  mVar.push_back("mds");
  mFS.push_back(1);
  
  // number of signal events needed for cross-section calculation
  double NData[3][2]; // in data
  double NReco[3][2]; // in MC at reco level
  double NGen[3][2]; // in MC at generator level
  
  // produce plot similar to SMP-12-002 paper Fig. 2
  // canvas: 0 for combined W->e and W-mu, 1 for W->mu, 2 for W->e
  TCanvas* c_mds[3];
  for(int ch = 0; ch < 3; ch++)
  {
    c_mds[ch] = new TCanvas(TString::Format("cmds%d", ch), "", 800, 800);
  }
  {
    TString var = "mds";
    int fs = 1;
    std::vector<TH1D*> hcp;
    hcp.resize(MCSamples.VecMCName.size() + 1);
    TLegend* leg = new TLegend(0.46, 0.52, 0.90, 0.92);
    leg->SetTextSize(0.036);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (int ch = 1; ch < 3; ch++)
    {
      c_mds[ch]->cd();
      hr_mass_ds->Draw();
      // MC
      std::vector<TH1D*> vecHMC;
      for(int mc = 0; mc < MCSamples.VecMCName.size(); mc++)
      {
        if(mc > 0)
          vecHMC.push_back(new TH1D(*vecHMC[mc - 1]));
        TString filename = TString::Format("%s/mc%sReco-c%d-f%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][0].Data(), ch, fs);
        printf("filename: %s\n", filename.Data());
        TFile* f = TFile::Open(filename);
        TH1D* hos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
        TH1D* hss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
        TH1D* h = new TH1D(*hos);
        h->Add(hss, -1.0);
        h->Scale(MCSamples.VecMCFactor[mc]);
        for(int mcf = 1; mcf < MCSamples.VecMCName[mc].size(); mcf++)
        {
          TFile* f = TFile::Open(TString::Format("%s/mc%sReco-c%d.root", gPlotsDir.Data(), MCSamples.VecMCName[mc][mcf].Data(), ch));
          TH1D* hhos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
          TH1D* hhss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
          TH1D* hh = new TH1D(*hhos);
          hh->Add(hss, -1.0);
          hh->Scale(MCSamples.VecMCFactor[mc]);
          h->Add(hh);
        }
        if(MCSamples.VecMCName[mc][0] == "Sig")
        {
          if(ch == 1)
          {
            // fit signal MC
            TF1* f_fit= new TF1("f_fit","gaus(0) + pol1(3)",0.136,0.170);
            f_fit->SetParameters(200.0,0.145,0.001);
            //f_fit->FixParameter(2,0.0005);
            h->Fit("f_fit","IR0");
            f_fit->SetLineColor(4);
            NReco[1][0] = f_fit->GetParameter(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002 / MCSamples.VecMCFactor[mc];
            NReco[1][1] = f_fit->GetParError(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002 / MCSamples.VecMCFactor[mc];
            printf("MC RECO: %.0f +- %.0f\n", NReco[1][0], NReco[1][1]);
          }
        }
        if(mc == 0)
          vecHMC.push_back(h);
        else
          vecHMC[mc]->Add(h);
      }
      // data
      TFile* fData = TFile::Open(TString::Format("%s/data-c%d-f%d.root", baseDir.Data(), ch, fs));
      TH1D* hDataos = (TH1D*)fData->Get(TString::Format("h_%s_os", var.Data()));
      TH1D* hDatass = (TH1D*)fData->Get(TString::Format("h_%s_ss", var.Data()));
      TH1D* hData = new TH1D(*hDataos);
      hData->Add(hDatass, -1.0);
      //printf("%s %s %f\n", fData->GetName(), hData->GetName(), hData->Integral());
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1.5);
      hData->SetLineColor(1);
      hData->SetMarkerColor(1);
      if(ch == 1)
        leg->AddEntry(hData, "Data", "pe");
      for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
      {
        vecHMC[mc]->SetFillColor(MCSamples.VecMCColor[mc]);
        vecHMC[mc]->SetLineColor(1);
        vecHMC[mc]->Draw("hist same");
        if(ch == 1)
          leg->AddEntry(vecHMC[mc], MCSamples.VecMCtitle[mc], "f");
        if(ch == 1)
          hcp[mc] = new TH1D(*vecHMC[mc]);
        else
          hcp[mc]->Add(vecHMC[mc]);
      }
      if(ch == 1)
      {
        // fit data
        TF1* f_fit= new TF1("f_fit","gaus(0) + pol1(3)",0.136,0.170);
        f_fit->SetParameters(200.0,0.145,0.001);
        //f_fit->FixParameter(2,0.0005);
        f_fit->SetLineColor(1);
        f_fit->SetLineWidth(2);
        f_fit->SetNpx(10000);
        hData->Fit(f_fit,"MIR0");
        leg->AddEntry(f_fit, "Gauss + polynomial ", "l");
        TString str = TString::Format("#DeltaM = %.1f #pm %.1f MeV", 1000 * f_fit->GetParameter(1), 1000 * f_fit->GetParError(1));
        leg->AddEntry((TObject*)NULL, str, "");
        NData[1][0] = f_fit->GetParameter(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002;
        NData[1][1] = f_fit->GetParError(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002;
        //NData[1][0] = f_fit->GetParameter(0);
        //NData[1][1] = f_fit->GetParError(0);
        str = TString::Format("N(D*) = %.0f #pm %.0f", NData[1][0], NData[1][1]);
        //leg->AddEntry((TObject*)NULL, str, "");
        f_fit->Draw("same");
      }
      // draw data
      hData->Draw("e0 same");
      leg->Draw();
      hr_mass_ds->Draw("axis same");
      if(ch == 1)
        hcp[hcp.size() - 1] = new TH1D(*hData);
      else
        hcp[hcp.size() - 1]->Add(hData);
    }
    // combined
    c_mds[0]->cd();
    hr_mass_ds->Draw();
    for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
      hcp[mc]->Draw("hist same");
    hcp[hcp.size() - 1]->Draw("e0 same");
    leg->Draw();
    hr_mass_ds->Draw("axis same");
  }
  // save plots
  for(int ch = 1; ch < 3; ch++)
  {
    c_mds[ch]->SaveAs(TString::Format("%s/mass_ds-c%d.eps", gPlotsDir.Data(), ch));
    c_mds[ch]->SaveAs(TString::Format("%s/mass_ds-c%d.pdf", gPlotsDir.Data(), ch));
  }
  c_mds[0]->SaveAs(TString::Format("%s/mass_ds.eps", gPlotsDir.Data()));
  c_mds[0]->SaveAs(TString::Format("%s/mass_ds.pdf", gPlotsDir.Data()));

  // *** paper Fig. 3 ***
  TCanvas* c_mdc[3];
  for(int ch = 0; ch < 3; ch++)
  {
    c_mdc[ch] = new TCanvas(TString::Format("cmdc%d", ch), "", 800, 800);
  }
  {
    TString var = "mdc";
    int fs = 2;
    std::vector<TH1D*> hcp;
    hcp.resize(MCSamples.VecMCName.size() + 1);
    TLegend* leg = new TLegend(0.40, 0.58, 0.90, 0.92);
    leg->SetTextSize(0.036);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (int ch = 1; ch < 3; ch++)
    {
      c_mdc[ch]->cd();
      hr_mass_dc->Draw();
      // MC
      std::vector<TH1D*> vecHMC;
      for(int mc = 0; mc < MCSamples.VecMCName.size(); mc++)
      {
        if(mc > 0)
          vecHMC.push_back(new TH1D(*vecHMC[mc - 1]));
        TString filename = TString::Format("%s/mc%sReco-c%d-f%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][0].Data(), ch, fs);
        printf("filename: %s\n", filename.Data());
        TFile* f = TFile::Open(filename);
        TH1D* hos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
        TH1D* hss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
        TH1D* h = new TH1D(*hos);
        h->Add(hss, -1.0);
        h->Scale(MCSamples.VecMCFactor[mc]);
        for(int mcf = 1; mcf < MCSamples.VecMCName[mc].size(); mcf++)
        {
          TFile* f = TFile::Open(TString::Format("%s/mc%sReco-c%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][mcf].Data(), ch));
          TH1D* hhos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
          TH1D* hhss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
          TH1D* hh = new TH1D(*hhos);
          hh->Add(hss, -1.0);
          hh->Scale(MCSamples.VecMCFactor[mc]);
          h->Add(hh);
        }
        if(MCSamples.VecMCName[mc][0] == "Sig")
        {
          if(ch == 1)
          {
            // fit signal MC
            TF1* f_fit= new TF1("f_fit","gaus(0) + pol1(3)",1.6,2.2);
            f_fit->SetParameters(300.0,1.8,0.02);
            //f_fit->FixParameter(2,0.015);
            h->Fit("f_fit","MIR0");
            f_fit->SetLineColor(4);
            f_fit->Draw("same");
            NReco[0][0] = f_fit->GetParameter(0) * (TMath::Sqrt(2.*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012 / MCSamples.VecMCFactor[mc];
            NReco[0][1] = f_fit->GetParError(0) * (TMath::Sqrt(2.*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012 / MCSamples.VecMCFactor[mc];
            //NReco[0][0] = f_fit->GetParameter(0) / vecMCFactorBR[fs - 1] / vecMCFactor[mc];
            //NReco[0][1] = f_fit->GetParError(0) / vecMCFactorBR[fs - 1] / vecMCFactor[mc];
            //NReco[0][0]*=TMath::Sqrt(3.);
            //NReco[0][0] = h->Integral() / vecMCFactorBR[fs - 1] / vecMCFactor[mc];
            printf("MC RECO: %.0f +- %.0f\n", NReco[0][0], NReco[0][1]);
          }
        }
        if(mc == 0)
          vecHMC.push_back(h);
        else
          vecHMC[mc]->Add(h);
      }
      // data
      TFile* fData = TFile::Open(TString::Format("%s/data-c%d-f%d.root", baseDir.Data(), ch, fs));
      TH1D* hDataos = (TH1D*)fData->Get(TString::Format("h_%s_os", var.Data()));
      TH1D* hDatass = (TH1D*)fData->Get(TString::Format("h_%s_ss", var.Data()));
      TH1D* hData = new TH1D(*hDataos);
      hData->Add(hDatass, -1.0);
      //printf("%s %s %f\n", fData->GetName(), hData->GetName(), hData->Integral());
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1.5);
      hData->SetLineColor(1);
      hData->SetMarkerColor(1);
      if(ch == 1)
        leg->AddEntry(hData, "Data", "pe");
      for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
      {
        vecHMC[mc]->SetFillColor(MCSamples.VecMCColor[mc]);
        vecHMC[mc]->SetLineColor(1);
        vecHMC[mc]->Draw("hist same");
        if(ch == 1)
          leg->AddEntry(vecHMC[mc], MCSamples.VecMCtitle[mc], "f");
        if(ch == 1)
          hcp[mc] = new TH1D(*vecHMC[mc]);
        else
          hcp[mc]->Add(vecHMC[mc]);
      }
      if(ch == 1)
      {
        // fit data
        TF1* f_fit= new TF1("f_fit","gaus(0) + pol1(3)",1.6,2.2);
        f_fit->SetParameters(500.0,1.8,0.02);
        //f_fit->FixParameter(2,0.015);
        f_fit->SetLineColor(1);
        f_fit->SetLineWidth(2);
        f_fit->SetNpx(10000);
        hData->Fit(f_fit,"MIR0");
        leg->AddEntry(f_fit, "Gauss + polynomial ", "l");
        TString str = TString::Format("M(D^{+}) = %.1f #pm %.1f MeV", 1000 * f_fit->GetParameter(1), 1000 * f_fit->GetParError(1));
        leg->AddEntry((TObject*)NULL, str, "");
        NData[0][0] = f_fit->GetParameter(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012;
        NData[0][1] = f_fit->GetParError(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012;
        //NData[0][0] = f_fit->GetParameter(0);
        //NData[0][1] = f_fit->GetParError(0);
        str = TString::Format("N(D^{+}) = %.0f #pm %.0f", NData[0][0], NData[0][1]);
        //leg->AddEntry((TObject*)NULL, str, "");
        f_fit->Draw("same");
      }
      // draw data
      hData->Draw("e0 same");
      leg->Draw();
      hr_mass_dc->Draw("axis same");
      if(ch == 1)
        hcp[hcp.size() - 1] = new TH1D(*hData);
      else
        hcp[hcp.size() - 1]->Add(hData);
    }
    // combined
    c_mdc[0]->cd();
    hr_mass_dc->Draw();
    for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
      hcp[mc]->Draw("hist same");
    hcp[hcp.size() - 1]->Draw("e0 same");
    leg->Draw();
    hr_mass_dc->Draw("axis same");
  }
  // save plots
  for(int ch = 1; ch < 3; ch++)
  {
    c_mdc[ch]->SaveAs(TString::Format("%s/mass_dc-c%d.eps", gPlotsDir.Data(), ch));
    c_mdc[ch]->SaveAs(TString::Format("%s/mass_dc-c%d.pdf", gPlotsDir.Data(), ch));
  }
  c_mdc[0]->SaveAs(TString::Format("%s/mass_dc.eps", gPlotsDir.Data()));
  c_mdc[0]->SaveAs(TString::Format("%s/mass_dc.pdf", gPlotsDir.Data()));


  // calculate cross section
  double lumi = 2500.0;
  double brDs = 0.00622;
  // D*
  {
    TFile* f = TFile::Open(baseDir + "/mcSigGen-c1-f1.root");
    TH1D* h = (TH1D*)f->Get("h_gen_ds");
    NGen[1][0] = h->Integral();
    double acc = NReco[1][0] / NGen[1][0];
    double accerr = NReco[1][1] / NGen[1][0];
    double cs = NData[1][0] / acc / lumi / brDs;
    double cserr = NData[1][1] / acc / lumi / brDs;
    printf("Data D* = %.0f +- %.0f\n", NData[1][0], NData[1][1]);
    printf("Reco D* = %.0f +- %.0f\n", NReco[1][0], NReco[1][1]);
    printf("Gen D* = %.0f\n", NGen[1][0]);
    printf("Acceptance D* = %.3f +- %.3f\n", acc, accerr);
    printf("Cross section D* = %.1f +- %.1f\n", cs, cserr);
  }

  double brDc = 0.0208;
  // D+
  {
    TFile* f = TFile::Open(baseDir + "/mcSigGen-c1-f2.root");
    TH1D* h = (TH1D*)f->Get("h_gen_dc");
    NGen[0][0] = h->Integral();
    double acc = NReco[0][0] / NGen[0][0];
    double accerr = NReco[0][1] / NGen[0][0];
    double cs = NData[0][0] / acc / lumi / brDc;
    double cserr = NData[0][1] / acc / lumi / brDc;
    printf("Data D+ = %.0f +- %.0f\n", NData[0][0], NData[0][1]);
    printf("Reco D+ = %.0f +- %.0f\n", NReco[0][0], NReco[0][1]);
    printf("Gen D+ = %.0f\n", NGen[0][0]);
    printf("Acceptance D+ = %.3f +- %.3f\n", acc, accerr);
    printf("Cross section D+ = %.1f +- %.1f\n", cs, cserr);
  }

  return 0;
}

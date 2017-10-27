// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> Helper for plotter >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifndef WCHARM_PLOTS_H
#define WCHARM_PLOTS_H

// additional files from this analysis
#include "wcharm_settings.h"
// C++ library or ROOT header files
#include <vector>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TMath.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> Prepare plot style >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Modify as you want, if needed consult
// https://root.cern.ch/doc/master/classTStyle.html
//
void SetStyle()
{
  gStyle->SetOptStat(000000000);
  gStyle->SetTitle(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(5);
  TGaxis::SetMaxDigits(4);
  gStyle->SetErrorX(0.0);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetNdivisions(206, "xyz");
}

// set style of frame histogram (font sizes etc.),
// modify as you want
//
// Arguments:
//    h: frame histogram
void SetCPHRange(TH2* h)
{
  TAxis* x = h->GetXaxis();
  x->SetTitleFont(62);
  x->SetLabelFont(62);
  x->SetTitleSize(0.045);
  x->SetLabelSize(0.045);
  x->SetTitleOffset(1.20);
  TAxis* y = h->GetYaxis();
  y->SetTitleFont(62);
  y->SetLabelFont(62);
  y->SetTitleSize(0.045);
  y->SetLabelSize(0.045);
  y->SetTitleOffset(1.70);
}

// helper class which contains MC samples which appear on control plots
class ZControlPlotMCSamples
{
  public:
    // names of MC sample for control plots:
    //   1st dimension: samples
    //   2nd dimension: subsamples
    // e.g. { signal { signalSubsampleA, signalSubsampleB}, background { bgSubsampleA, bgSubsampleB} }
    std::vector<std::vector<TString> > VecMCName;
    std::vector<int> VecMCColor; // colors, same size as the number of samples above
    std::vector<TString> VecMCtitle; // titles, same size as the number of samples above
    std::vector<double> VecMCFactor; // weight factors, same size as the number of samples above
};

// *** SMP-12-002 Fig. 7 ***
// Arguments:
//    const ZControlPlotMCSamples& MCSamples: MC samples (see description in class ZControlPlotMCSamples)
//    const std::vector<double>& vecMCFactor: vector with weights of MC samples (same size as vecMCName)
//    const std::vector<TH2F*>& cpHR: vector with 2D frame histograms for each variable
//    const std::vector<TString>& cpVar: vector with variable names (same size as cpHR)
//    const std::vector<int>& cpFS: vector with charm final-state identifiers
//
void FigureControlPlots(
    const ZControlPlotMCSamples& MCSamples,
    const std::vector<TH2F*>& cpHR,
    const std::vector<TString>& cpVar,
    const std::vector<int>& cpFS
    )
{
  // canvas: 0 for combined W->e and W-mu, 1 for W->mu, 2 for W->e
  TCanvas* c_cp[3];
  for(int ch = 0; ch < 3; ch++)
  {
    c_cp[ch] = new TCanvas(TString::Format("c%d", ch), "", 800, 1200);
    c_cp[ch]->Divide(2, 3);
  }
  // loop over 6 variables
  for(int v = 0; v < 6; v++)
  {
    // charm final state for plot labeling
    TString fsLabel;
    if(v / 2 == 0)
      fsLabel = "D^{#pm}#rightarrowK^{#pm}#pi^{#pm}#pi^{#pm}";
    if(v / 2 == 1)
      fsLabel = "D^{*#pm}#rightarrowD^{0}(K^{#pm}#pi^{#pm})#pi^{#pm}";
    if(v / 2 == 2)
      fsLabel = "semileptonic";
    // variable
    TString var = cpVar[v];
    // charm final state
    int fs = cpFS[v];
    // vector with all histograms: the size is number of MC samples + 1 data sample
    std::vector<TH1D*> hcp;
    hcp.resize(MCSamples.VecMCName.size() + 1);
    // legend
    TLegend* leg = new TLegend(0.54, 0.62, 0.90, 0.92);
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    // loop over e and mu W decay channels
    for (int ch = 1; ch < 3; ch++)
    {
      c_cp[ch]->cd(v + 1);
      // draf frame histogram
      cpHR[v]->Draw();
      // produce cumulative MC histograms
      std::vector<TH1D*> vecHMC;

      // loop over MC samples
      for(int mc = 0; mc < MCSamples.VecMCName.size(); mc++)
      {
        if(mc > 0)
          vecHMC.push_back(new TH1D(*vecHMC[mc - 1]));
        TString filename = TString::Format("%s/mc%sReco-c%d-f%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][0].Data(), ch, fs);
        printf("opening file %s\n", filename.Data());
        TFile* f = TFile::Open(filename);
        // histogram with opposite sign (OS) Wc events
        TH1D* hos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
        // histogram with same sign (SS) Wc events
        TH1D* hss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
        // histogram with OS - SS
        TH1D* h = new TH1D(*hos);
        h->Add(hss, -1.0);
        // scale by MC weight
        h->Scale(MCSamples.VecMCFactor[mc]);
        // loop over MC subsamples
        for(int mcf = 1; mcf < MCSamples.VecMCName[mc].size(); mcf++)
        {
          TString filename = TString::Format("%s/mc%sReco-c%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][mcf].Data(), ch);
          printf("opening file %s\n", filename.Data());
          TFile* f = TFile::Open(filename);
          // histogram with opposite sign (OS) Wc events
          TH1D* hhos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
          // histogram with same sign (OS) Wc events
          TH1D* hhss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
          // histogram with OS - SS
          TH1D* hh = new TH1D(*hhos);
          hh->Add(hhss, -1.0);
          // scale by MC weight
          hh->Scale(MCSamples.VecMCFactor[mc]);
          // add this MC subsample
          h->Add(hh);
        }
        // if this is first MC sample, push new histogram, otherwise add to existing one
        if(mc == 0)
          vecHMC.push_back(h);
        else
          vecHMC[mc]->Add(h);
      }

      // data
      TString filename = TString::Format("%s/data-c%d-f%d.root", gHistDir.Data(), ch, fs);
      printf("opening file %s\n", filename.Data());
      TFile* fData = TFile::Open(filename);
      // histogram with opposite sign (OS) Wc events
      TH1D* hDataos = (TH1D*)fData->Get(TString::Format("h_%s_os", var.Data()));
      // histogram with same sign (OS) Wc events
      TH1D* hDatass = (TH1D*)fData->Get(TString::Format("h_%s_ss", var.Data()));
      // histogram with OS - SS
      TH1D* hData = new TH1D(*hDataos);
      hData->Add(hDatass, -1.0);
      // graph cosmetics
      hData->SetMarkerStyle(20);
      hData->SetMarkerSize(1);
      hData->SetLineColor(1);
      hData->SetMarkerColor(1);
      // data legend entry
      if(ch == 1)
        leg->AddEntry(hData, "Data", "pe");

      // rescale MC to data
      double rescaleFactor = hData->Integral() / vecHMC[MCSamples.VecMCName.size() - 1]->Integral();
      if(rescaleFactor != rescaleFactor)
        rescaleFactor = 1.0;
      printf("rescaleFactor: %f\n", rescaleFactor);

      // draw MC
      for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
      {
        // apply rescale-to-data factor
        vecHMC[mc]->Scale(rescaleFactor);
        // histogram cosmetics
        vecHMC[mc]->SetFillColor(MCSamples.VecMCColor[mc]);
        vecHMC[mc]->SetLineColor(1);
        // draw MC hisogram
        vecHMC[mc]->Draw("hist same");
        // MC legend entry
        if(ch == 1)
          leg->AddEntry(vecHMC[mc], MCSamples.VecMCtitle[mc], "f");
        // if this is first W decay channel, push new histogram (sum of W->mu and W->e channels), otherwise add to existing one
        if(ch == 1)
          hcp[mc] = new TH1D(*vecHMC[mc]);
        else
          hcp[mc]->Add(vecHMC[mc]);
      }

      // draw data
      hData->Draw("e0 same");
      // draw legend
      leg->Draw();
      // draw text label showing W decay channel
      TPaveText* pt = new TPaveText(0.20, 0.78, 0.55, 0.90, "ndc");
      pt->AddText(TString::Format("W#rightarrow%s", (ch == 1 ? "#mu" : "e")));
      pt->AddText(fsLabel);
      pt->SetTextAlign(11);
      pt->SetFillColor(0);
      pt->SetBorderSize(0);
      pt->SetTextSize(0.05);
      pt->Draw();
      // re-draw frame
      cpHR[v]->Draw("axis same");
      // if this is first W decay channel, push new histogram (sum of W->mu and W->e channels), otherwise add to existing one
      if(ch == 1)
        hcp[hcp.size() - 1] = new TH1D(*hData);
      else
        hcp[hcp.size() - 1]->Add(hData);
    }

    // combined channel (W->mu and W->e)
    c_cp[0]->cd(v + 1);
    cpHR[v]->Draw();
    for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
      hcp[mc]->Draw("hist same");
    hcp[hcp.size() - 1]->Draw("e0 same");
    leg->Draw();
    // draw text label showing W decay channel
    TPaveText* pt = new TPaveText(0.20, 0.78, 0.55, 0.90, "ndc");
    pt->AddText(TString::Format("W#rightarrow#mu,e"));
    pt->AddText(fsLabel);
    pt->SetTextAlign(11);
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.05);
    pt->Draw();
    // re-draw frame
    cpHR[v]->Draw("axis same");
  }

  // save plots
  for(int ch = 1; ch < 3; ch++)
  {
    c_cp[ch]->SaveAs(TString::Format("%s/cp-c%d.eps", gPlotsDir.Data(), ch));
    c_cp[ch]->SaveAs(TString::Format("%s/cp-c%d.pdf", gPlotsDir.Data(), ch));
    c_cp[ch]->SaveAs(TString::Format("%s/cp-c%d.png", gPlotsDir.Data(), ch));
  }
  c_cp[0]->SaveAs(TString::Format("%s/cp.eps", gPlotsDir.Data()));
  c_cp[0]->SaveAs(TString::Format("%s/cp.pdf", gPlotsDir.Data()));
  c_cp[0]->SaveAs(TString::Format("%s/cp.png", gPlotsDir.Data()));
}


// class to store information about event numbers needed for cross-section calculation
class ZYieldWcharm
{
  public:
    double NData; // number of event in data
    double NDataUnc; // uncertainty on the number of events in data
    double NReco; // number of events in MC at reco level
    double NRecoUnc; // uncertainty on the number of events in MC at reco level
    double NGen; // number of events in MC at generator level
    double NGenUnc; // uncertainty on the number of events in MC at generator level

    int WDecayChannel; // W decay channel (1 for mu, 2 for e)
    int CharmFinalState; // charm final state (1 for D*, 2 for D+)
};

// *** SMP-12-002 Fig. 3 ***
// Arguments:
//    TH2F* hFrame: reference frame 2D histogram
//    const ZControlPlotMCSamples& MCSamples: MC samples (see description in class ZControlPlotMCSamples)
//
// Return event yields for two W decay channels (see class ZYieldWcharm for description)
//
std::vector<ZYieldWcharm> FigureMassDstar(
    TH2F* hFrame,
    const ZControlPlotMCSamples& MCSamples
    )
{
  // event yields to be returned (muon and electron W decay channels)
  // (see class ZYieldWcharm for description)
  std::vector<ZYieldWcharm> yield(2);

  TCanvas* c[3];
  // canvas: 0 for combined W->e and W-mu (not used), 1 for W->mu, 2 for W->e
  for(int ch = 0; ch < 3; ch++)
  {
    c[ch] = new TCanvas(TString::Format("cmds%d", ch), "", 800, 800);
  }

  TString var = "mds";
  int fs = 1;
  // loop over e and mu W decay channels
  for (int ch = 1; ch < 3; ch++)
  {
    // legend
    TLegend* leg = new TLegend(0.46, 0.52, 0.90, 0.92);
    leg->SetTextSize(0.036);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    c[ch]->cd();
    hFrame->Draw();
    // MC
    std::vector<TH1D*> vecHMC;
    for(int mc = 0; mc < MCSamples.VecMCName.size(); mc++)
    {
      if(mc > 0)
        vecHMC.push_back(new TH1D(*vecHMC[mc - 1]));
      TString filename = TString::Format("%s/mc%sReco-c%d-f%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][0].Data(), ch, fs);
      TFile* f = TFile::Open(filename);
      // histogram with opposite sign (OS) Wc events
      TH1D* hos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
      // histogram with same sign (OS) Wc events
      TH1D* hss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
      // histogram with OS - SS
      TH1D* h = new TH1D(*hos);
      h->Add(hss, -1.0);
      h->Scale(MCSamples.VecMCFactor[mc]);
      for(int mcf = 1; mcf < MCSamples.VecMCName[mc].size(); mcf++)
      {
        TFile* f = TFile::Open(TString::Format("%s/mc%sReco-c%d.root", gPlotsDir.Data(), MCSamples.VecMCName[mc][mcf].Data(), ch));
        TH1D* hhos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
        TH1D* hhss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
        TH1D* hh = new TH1D(*hhos);
        hh->Add(hhss, -1.0);
        hh->Scale(MCSamples.VecMCFactor[mc]);
        h->Add(hh);
      }
      if(MCSamples.VecMCName[mc][0] == "Sig")
      {
        // fit signal MC
        TF1* f_fit= new TF1("f_fit","gaus(0) + pol1(3)",0.136,0.170);
        f_fit->SetParameters(200.0,0.145,0.001);
        //f_fit->FixParameter(2,0.0005);
        h->Fit("f_fit","IR0");
        f_fit->SetLineColor(4);
        //f_fit->Draw("same");
        // number of events in the peak
        yield[ch - 1].NReco = f_fit->GetParameter(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002 / MCSamples.VecMCFactor[mc];
        // uncertainty on this number
        yield[ch - 1].NRecoUnc = f_fit->GetParError(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002 / MCSamples.VecMCFactor[mc];
        printf("MC RECO: %.0f +- %.0f\n", yield[ch - 1].NReco, yield[ch - 1].NRecoUnc);
        // store W decay channel
        yield[ch - 1].WDecayChannel = ch;
      }
      if(mc == 0)
        vecHMC.push_back(h);
      else
        vecHMC[mc]->Add(h);
    }
    // data
    TFile* fData = TFile::Open(TString::Format("%s/data-c%d-f%d.root", gHistDir.Data(), ch, fs));
    TH1D* hDataos = (TH1D*)fData->Get(TString::Format("h_%s_os", var.Data()));
    TH1D* hDatass = (TH1D*)fData->Get(TString::Format("h_%s_ss", var.Data()));
    TH1D* hData = new TH1D(*hDataos);
    hData->Add(hDatass, -1.0);
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.5);
    hData->SetLineColor(1);
    hData->SetMarkerColor(1);

    // rescale MC to data
    double rescaleFactor = hData->Integral() / vecHMC[MCSamples.VecMCName.size() - 1]->Integral();
    if(rescaleFactor != rescaleFactor)
      rescaleFactor = 1.0;
    printf("rescaleFactor: %f\n", rescaleFactor);

    // print W decay channel in lefend
    leg->AddEntry((TObject*)(NULL), (ch == 1) ? "W#rightarrowmu" : "W#rightarrowe", "");

    // add data legend entry
    leg->AddEntry(hData, "Data", "pe");

    for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
    {
      // apply rescale-to-data factor
      vecHMC[mc]->Scale(rescaleFactor);
      // plotting cosmetics
      vecHMC[mc]->SetFillColor(MCSamples.VecMCColor[mc]);
      vecHMC[mc]->SetLineColor(1);
      vecHMC[mc]->Draw("hist same");
      leg->AddEntry(vecHMC[mc], MCSamples.VecMCtitle[mc], "f");
    }

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
    // number of events in the peak
    yield[ch - 1].NData = f_fit->GetParameter(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002;
    // uncertainty on this number
    yield[ch - 1].NDataUnc = f_fit->GetParError(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.002;
    str = TString::Format("N(D*) = %.0f #pm %.0f", yield[ch - 1].NReco, yield[ch - 1].NRecoUnc);
    f_fit->Draw("same");

    // draw data
    hData->Draw("e0 same");
    leg->Draw();
    hFrame->Draw("axis same");

    // save plots
    c[ch]->SaveAs(TString::Format("%s/mass_ds-c%d.eps", gPlotsDir.Data(), ch));
    c[ch]->SaveAs(TString::Format("%s/mass_ds-c%d.pdf", gPlotsDir.Data(), ch));
    c[ch]->SaveAs(TString::Format("%s/mass_ds-c%d.png", gPlotsDir.Data(), ch));
  }

  return yield;
}

// *** SMP-12-002 Fig. 2 ***
// Arguments:
//    TH2F* hFrame: reference frame 2D histogram
//    const ZControlPlotMCSamples& MCSamples: MC samples (see description in class ZControlPlotMCSamples)
//
// Return event yields for two W decay channels (see class ZYieldWcharm for description)
//
std::vector<ZYieldWcharm> FigureMassDch(
    TH2F* hFrame,
    const ZControlPlotMCSamples& MCSamples
    )
{
  // event yields to be returned (muon and electron W decay channels)
  // (see class ZYieldWcharm for description)
  std::vector<ZYieldWcharm> yield(2);

  TCanvas* c[3];
  // canvas: 0 for combined W->e and W-mu (not used), 1 for W->mu, 2 for W->e
  for(int ch = 1; ch < 3; ch++)
  {
    c[ch] = new TCanvas(TString::Format("cmdc%d", ch), "", 800, 800);
    TString var = "mdc";
    int fs = 2;
    // legend
    TLegend* leg = new TLegend(0.40, 0.58, 0.90, 0.92);
    leg->SetTextSize(0.036);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    // loop over e and mu W decay channels
    c[ch]->cd();
    hFrame->Draw();
    // MC
    std::vector<TH1D*> vecHMC;
    for(int mc = 0; mc < MCSamples.VecMCName.size(); mc++)
    {
      if(mc > 0)
        vecHMC.push_back(new TH1D(*vecHMC[mc - 1]));
      TString filename = TString::Format("%s/mc%sReco-c%d-f%d.root", gHistDir.Data(), MCSamples.VecMCName[mc][0].Data(), ch, fs);
      printf("filename: %s\n", filename.Data());
      TFile* f = TFile::Open(filename);
      // histogram with opposite sign (OS) Wc events
      TH1D* hos = (TH1D*)f->Get(TString::Format("h_%s_os", var.Data()));
      // histogram with same sign (OS) Wc events
      TH1D* hss = (TH1D*)f->Get(TString::Format("h_%s_ss", var.Data()));
      // histogram with OS - SS
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
        // fit signal MC
        TF1* f_fit= new TF1("f_fit","gaus(0) + pol1(3)",1.6,2.2);
        f_fit->SetParameters(50.0,1.87,0.02);
        f_fit->SetParLimits(1,1.85,1.90);
        f_fit->SetParLimits(2,0.005,0.03);
        //f_fit->FixParameter(2,0.015);
        h->Fit("f_fit","MIR0");
        f_fit->SetLineColor(4);
        //f_fit->Draw("same");
        // number of events in the peak
        yield[ch - 1].NReco = f_fit->GetParameter(0) * (TMath::Sqrt(2.*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012 / MCSamples.VecMCFactor[mc];
        // uncertainty on this number
        yield[ch - 1].NRecoUnc = f_fit->GetParError(0) * (TMath::Sqrt(2.*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012 / MCSamples.VecMCFactor[mc];
        printf("MC RECO: %.0f +- %.0f\n", yield[ch - 1].NReco, yield[ch - 1].NRecoUnc);
        // store W decay channel
        yield[ch - 1].WDecayChannel = ch;
      }
      if(mc == 0)
        vecHMC.push_back(h);
      else
        vecHMC[mc]->Add(h);
    }
    // data
    TFile* fData = TFile::Open(TString::Format("%s/data-c%d-f%d.root", gHistDir.Data(), ch, fs));
    TH1D* hDataos = (TH1D*)fData->Get(TString::Format("h_%s_os", var.Data()));
    TH1D* hDatass = (TH1D*)fData->Get(TString::Format("h_%s_ss", var.Data()));
    TH1D* hData = new TH1D(*hDataos);
    hData->Add(hDatass, -1.0);
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.5);
    hData->SetLineColor(1);
    hData->SetMarkerColor(1);

    // rescale MC to data
    double rescaleFactor = hData->Integral() / vecHMC[MCSamples.VecMCName.size() - 1]->Integral();
    if(rescaleFactor != rescaleFactor)
      rescaleFactor = 1.0;
    printf("rescaleFactor: %f\n", rescaleFactor);

    // print W decay channel in lefend
    leg->AddEntry((TObject*)(NULL), (ch == 1) ? "W#rightarrowmu" : "W#rightarrowe", "");

    // add data legend entry
    leg->AddEntry(hData, "Data", "pe");
    for(int mc = MCSamples.VecMCName.size() - 1; mc >= 0; mc--)
    {
      // apply rescale-to-data factor
      vecHMC[mc]->Scale(rescaleFactor);
      // plotting cosmetics
      vecHMC[mc]->SetFillColor(MCSamples.VecMCColor[mc]);
      vecHMC[mc]->SetLineColor(1);
      vecHMC[mc]->Draw("hist same");
      leg->AddEntry(vecHMC[mc], MCSamples.VecMCtitle[mc], "f");
    }

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
    // number of events in the peak
    yield[ch - 1].NData = f_fit->GetParameter(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012;
    // uncertainty on this number
    yield[ch - 1].NDataUnc = f_fit->GetParError(0) * (TMath::Sqrt(2*TMath::Pi()) * f_fit->GetParameter(2)) / 0.012;
    str = TString::Format("N(D^{+}) = %.0f #pm %.0f", yield[ch - 1].NData, yield[ch - 1].NDataUnc);
    f_fit->Draw("same");

    // draw data
    hData->Draw("e0 same");
    leg->Draw();
    hFrame->Draw("axis same");
  }
  // save plots
  for(int ch = 1; ch < 3; ch++)
  {
    c[ch]->SaveAs(TString::Format("%s/mass_dc-c%d.eps", gPlotsDir.Data(), ch));
    c[ch]->SaveAs(TString::Format("%s/mass_dc-c%d.pdf", gPlotsDir.Data(), ch));
    c[ch]->SaveAs(TString::Format("%s/mass_dc-c%d.png", gPlotsDir.Data(), ch));
  }

  return yield;
}


// calculate cross section
// Arguments:
//    const ZYieldWcharm& yield: event yields at different levels (see class ZYieldWcharm for description)
//    const double scale: scaling factor
//
void CalculateCrossSection(const ZYieldWcharm& yield, const double scale)
{
  // get generator level histogram to get the number og events at generator level
  TString filename = TString::Format(gHistDir + "/mcSigGen-c%d-f%d.root", yield.WDecayChannel, yield.CharmFinalState);
  TFile* f = TFile::Open(filename);
  TString genHistoName = "h_gen_";
  if(yield.CharmFinalState == 1)
    genHistoName += "ds";
  else if(yield.CharmFinalState == 2)
    genHistoName += "dc";
  TH1D* h = (TH1D*)f->Get(genHistoName);
  double nGen = h->Integral();

  // calculate acceptance (also referred to as efficiency)
  double acc = yield.NReco / nGen;
  // uncertainty on acceptance (MC statistics is not unlimited)
  // (uncertainty on the number of generated events is ignored)
  double accerr = yield.NRecoUnc / nGen;

  // calculate cross section
  double cs = yield.NData / acc * scale;
  // uncertainty on corss section
  double cserr = yield.NDataUnc / acc * scale;

  // meson name for printing
  TString meson;
  if(yield.CharmFinalState == 1)
    meson = "c->D*";
  else if(yield.CharmFinalState == 2)
    meson = "c->D+";
  // W decay channel for printing
  TString decay;
  if(yield.WDecayChannel == 1)
    decay = "W->mu";
  else if(yield.WDecayChannel == 2)
    decay = "W->e";
  // print event yields and calculated quantities
  printf("%s\n", std::string(100, '*').c_str());
  printf("%s %s\n", decay.Data(), meson.Data());
  printf("Data = %.0f +- %.0f\n", yield.NData, yield.NDataUnc);
  printf("Reco = %.0f +- %.0f\n", yield.NReco, yield.NRecoUnc);
  printf("Gen = %.0f\n", meson.Data(), nGen);
  printf("Acceptance = %.3f +- %.3f\n", acc, accerr);
  printf("Cross section = %.1f +- %.1f\n", cs, cserr);
  printf("%s\n", std::string(100, '*').c_str());
}

//
// be aware: there are certainly memory leaks above (not removing
// dynamically allocated objects), although it does not matter here,
// because these routines ire called in a small program
// and all memory is free when execution finished
//

#endif // WCHARM_PLOTS_H

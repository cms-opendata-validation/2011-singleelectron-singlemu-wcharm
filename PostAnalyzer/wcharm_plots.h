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
    cpHR[v]->Draw("axis same");
  }

  // save plots
  for(int ch = 1; ch < 3; ch++)
  {
    c_cp[ch]->SaveAs(TString::Format("%s/cp-c%d.eps", gPlotsDir.Data(), ch));
    c_cp[ch]->SaveAs(TString::Format("%s/cp-c%d.pdf", gPlotsDir.Data(), ch));
  }
  c_cp[0]->SaveAs(TString::Format("%s/cp.eps", gPlotsDir.Data()));
  c_cp[0]->SaveAs(TString::Format("%s/cp.pdf", gPlotsDir.Data()));
  //
  // be aware: there are certainly memory leaks above (not removing
  // dynamically allocated objects), although it does not matter here,
  // because this routine is called in a small routine
  // and all memory is free when execution finished
  //
}
#endif // WCHARM_PLOTS_H

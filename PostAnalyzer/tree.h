//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 19 13:41:08 2016 by ROOT version 5.34/04
// from TTree tree/spam
// found on file: ../mc/W1Jet_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_W1Jet_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00000/wcharmSel_6.root
//////////////////////////////////////////////////////////

#ifndef ZTree_h
#define ZTree_h

#include <TROOT.h>
#include <TChain.h>
//#include <TFile.h>
//#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ZTree {//: public TSelector {
public :

   // MC flag
   bool _flagMC;

   // variable array max sizes
   static const int maxNel = 10;
   static const int maxNmu = 10;
   static const int maxNjet = 100;
   static const int maxDsCand = 500;
   static const int maxDCand = 50000;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          Nevent;
   UInt_t          Nrun;
   Int_t           lumiblock;
   Float_t         metPx;
   Float_t         metPy;
   Int_t           Nel;
   Float_t         el_pt[maxNel];   //[Nel]
   Float_t         el_eta[maxNel];   //[Nel]
   Float_t         el_isolat[maxNel];   //[Nel]
   Float_t         el_phi[maxNel];   //[Nel]
   Float_t         el_mt[maxNel];   //[Nel]
   Float_t         el_convDist[maxNel];   //[Nel]
   Float_t         el_convDcot[maxNel];   //[Nel]
   Int_t           el_mh[maxNel];   //[Nel]
   Int_t           Nmu;
   Float_t         mu_pt[maxNmu];   //[Nmu]
   Float_t         mu_eta[maxNmu];   //[Nmu]
   Float_t         mu_isolat[maxNmu];   //[Nmu]
   Float_t         mu_phi[maxNmu];   //[Nmu]
   Float_t         mu_mt[maxNmu];   //[Nmu]
   Int_t           mu_ValidHits[maxNmu];   //[Nmu]
   Int_t           mu_PixelHits[maxNmu];   //[Nmu]
   Float_t         mu_imp[maxNmu];   //[Nmu]
   Float_t         mu_nChi2[maxNmu];   //[Nmu]
   Int_t           mu_stations[maxNmu];   //[Nmu]
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   Int_t           DsCand;
   Float_t         dmDs[maxDsCand];   //[DsCand]
   Float_t         vcx[maxDsCand];   //[DsCand]
   Float_t         vcy[maxDsCand];   //[DsCand]
   Float_t         vcz[maxDsCand];   //[DsCand]
   Float_t         ds_pt[maxDsCand];   //[DsCand]
   Float_t         ds_eta[maxDsCand];   //[DsCand]
   Float_t         ds_m[maxDsCand];   //[DsCand]
   Float_t         ds_phi[maxDsCand];   //[DsCand]
   Float_t         D0dl[maxDsCand];   //[DsCand]
   Float_t         D0dlEr[maxDsCand];   //[DsCand]
   Float_t         D0dls[maxDsCand];   //[DsCand]
   Float_t         D0sv_x[maxDsCand];   //[DsCand]
   Float_t         D0sv_y[maxDsCand];   //[DsCand]
   Float_t         D0sv_z[maxDsCand];   //[DsCand]
   Int_t           DCand;
   Float_t         d_pt[maxDCand];   //[DCand]
   Float_t         d_eta[maxDCand];   //[DCand]
   Float_t         d_phi[maxDCand];   //[DCand]
   Float_t         d_m[maxDCand];   //[DCand]
   Float_t         Ddl[maxDCand];   //[DCand]
   Float_t         DdlEr[maxDCand];   //[DCand]
   Float_t         Ddls[maxDCand];   //[DCand]
   Float_t         Dsv_x[maxDCand];   //[DCand]
   Float_t         Dsv_y[maxDCand];   //[DCand]
   Float_t         Dsv_z[maxDCand];   //[DCand]
   Int_t           Njet;
   Float_t         j_pt[maxNjet];   //[Njet]
   Float_t         j_eta[maxNjet];   //[Njet]
   Float_t         j_m[maxNjet];   //[Njet]
   Float_t         j_phi[maxNjet];   //[Njet]
   Int_t           j_DsInJ[maxNjet];   //[Njet]
   Int_t           j_DInJ[maxNjet];   //[Njet]
   Int_t           Triggers;
   // MC
   Int_t           mcEventType;
   Int_t           CGen;

   // List of branches
   TBranch        *b_Nevent;   //!
   TBranch        *b_Nrun;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_Nel;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_isolat;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_mt;   //!
   TBranch        *b_el_convDist;   //!
   TBranch        *b_el_convDcot;   //!
   TBranch        *b_el_mh;   //!
   TBranch        *b_Nmu;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_isolat;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_mt;   //!
   TBranch        *b_mu_ValidHits;   //!
   TBranch        *b_mu_PixelHits;   //!
   TBranch        *b_mu_imp;   //!
   TBranch        *b_mu_nChi2;   //!
   TBranch        *b_mu_stations;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_DsCand;   //!
   TBranch        *b_dmDs;   //!
   TBranch        *b_vcx;   //!
   TBranch        *b_vcy;   //!
   TBranch        *b_vcz;   //!
   TBranch        *b_ds_pt;   //!
   TBranch        *b_ds_eta;   //!
   TBranch        *b_ds_m;   //!
   TBranch        *b_ds_phi;   //!
   TBranch        *b_D0dl;   //!
   TBranch        *b_D0dlEr;   //!
   TBranch        *b_D0dls;   //!
   TBranch        *b_D0sv_x;   //!
   TBranch        *b_D0sv_y;   //!
   TBranch        *b_D0sv_z;   //!
   TBranch        *b_DCand;   //!
   TBranch        *b_d_pt;   //!
   TBranch        *b_d_eta;   //!
   TBranch        *b_d_phi;   //!
   TBranch        *b_d_m;   //!
   TBranch        *b_Ddl;   //!
   TBranch        *b_DdlEr;   //!
   TBranch        *b_Ddls;   //!
   TBranch        *b_Dsv_x;   //!
   TBranch        *b_Dsv_y;   //!
   TBranch        *b_Dsv_z;   //!
   TBranch        *b_Njet;   //!
   TBranch        *b_j_pt;   //!
   TBranch        *b_j_eta;   //!
   TBranch        *b_j_m;   //!
   TBranch        *b_j_phi;   //!
   TBranch        *b_j_DsInJ;   //!
   TBranch        *b_j_DInJ;   //!
   TBranch        *b_Triggers;   //!
   // MC
   TBranch        *b_mcEventType;   //!
   TBranch        *b_CGen;   //!

   ZTree(bool flagMC = false) : fChain(0), _flagMC(flagMC) { }
   virtual ~ZTree() { }
   /*virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);*/
   virtual void    Init(TTree *tree);
   /*virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();*/

   //ClassDef(tree,0);
};

//#endif

//#ifdef tree_cxx
void ZTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
   fChain->SetBranchAddress("Nrun", &Nrun, &b_Nrun);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("Nel", &Nel, &b_Nel);
   fChain->SetBranchAddress("el_pt", el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_isolat", el_isolat, &b_el_isolat);
   fChain->SetBranchAddress("el_phi", el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_mt", el_mt, &b_el_mt);
   fChain->SetBranchAddress("el_convDist", el_convDist, &b_el_convDist);
   fChain->SetBranchAddress("el_convDcot", el_convDcot, &b_el_convDcot);
   fChain->SetBranchAddress("el_mh", el_mh, &b_el_mh);
   fChain->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
   fChain->SetBranchAddress("mu_pt", mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_isolat", mu_isolat, &b_mu_isolat);
   fChain->SetBranchAddress("mu_phi", mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_mt", mu_mt, &b_mu_mt);
   fChain->SetBranchAddress("mu_ValidHits", mu_ValidHits, &b_mu_ValidHits);
   fChain->SetBranchAddress("mu_PixelHits", mu_PixelHits, &b_mu_PixelHits);
   fChain->SetBranchAddress("mu_imp", mu_imp, &b_mu_imp);
   fChain->SetBranchAddress("mu_nChi2", mu_nChi2, &b_mu_nChi2);
   fChain->SetBranchAddress("mu_stations", mu_stations, &b_mu_stations);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("DsCand", &DsCand, &b_DsCand);
   fChain->SetBranchAddress("dmDs", dmDs, &b_dmDs);
   fChain->SetBranchAddress("vcx", vcx, &b_vcx);
   fChain->SetBranchAddress("vcy", vcy, &b_vcy);
   fChain->SetBranchAddress("vcz", vcz, &b_vcz);
   fChain->SetBranchAddress("ds_pt", ds_pt, &b_ds_pt);
   fChain->SetBranchAddress("ds_eta", ds_eta, &b_ds_eta);
   fChain->SetBranchAddress("ds_m", ds_m, &b_ds_m);
   fChain->SetBranchAddress("ds_phi", ds_phi, &b_ds_phi);
   fChain->SetBranchAddress("D0dl", D0dl, &b_D0dl);
   fChain->SetBranchAddress("D0dlEr", D0dlEr, &b_D0dlEr);
   fChain->SetBranchAddress("D0dls", D0dls, &b_D0dls);
   fChain->SetBranchAddress("D0sv_x", D0sv_x, &b_D0sv_x);
   fChain->SetBranchAddress("D0sv_y", D0sv_y, &b_D0sv_y);
   fChain->SetBranchAddress("D0sv_z", D0sv_z, &b_D0sv_z);
   fChain->SetBranchAddress("DCand", &DCand, &b_DCand);
   fChain->SetBranchAddress("d_pt", d_pt, &b_d_pt);
   fChain->SetBranchAddress("d_eta", d_eta, &b_d_eta);
   fChain->SetBranchAddress("d_phi", d_phi, &b_d_phi);
   fChain->SetBranchAddress("d_m", d_m, &b_d_m);
   fChain->SetBranchAddress("Ddl", Ddl, &b_Ddl);
   fChain->SetBranchAddress("DdlEr", DdlEr, &b_DdlEr);
   fChain->SetBranchAddress("Ddls", Ddls, &b_Ddls);
   fChain->SetBranchAddress("Dsv_x", Dsv_x, &b_Dsv_x);
   fChain->SetBranchAddress("Dsv_y", Dsv_y, &b_Dsv_y);
   fChain->SetBranchAddress("Dsv_z", Dsv_z, &b_Dsv_z);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
   fChain->SetBranchAddress("j_pt", j_pt, &b_j_pt);
   fChain->SetBranchAddress("j_eta", j_eta, &b_j_eta);
   fChain->SetBranchAddress("j_m", j_m, &b_j_m);
   fChain->SetBranchAddress("j_phi", j_phi, &b_j_phi);
   fChain->SetBranchAddress("j_DsInJ", j_DsInJ, &b_j_DsInJ);
   fChain->SetBranchAddress("j_DInJ", j_DInJ, &b_j_DInJ);
   fChain->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
   // MC
   if(_flagMC) fChain->SetBranchAddress("mcEventType", &mcEventType, &b_mcEventType);
   if(_flagMC) fChain->SetBranchAddress("CGen", &CGen, &b_CGen);
}

/*Bool_t tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}*/

#endif // #ifdef ZTree_h

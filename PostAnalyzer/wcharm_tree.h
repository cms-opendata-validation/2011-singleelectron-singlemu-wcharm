// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This class contains all needed variables to read ROOT ntuples for W+c analysis
// (automaticlly produced by ROOT, then slightly tuned manually)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifndef WCHARM_TREE_H
#define WCHARM_TREE_H

// ROOT header files
#include <TROOT.h>
#include <TChain.h>

// Class which gives access to all information in each event stored in ntuples
class ZTree
{
  public :
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // public members (for direct access outside the class)
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // pointer to the analyzed TTree or TChain
    TTree *fChain;

    // MC flag (true for MC, false for data)
    bool _flagMC;

    // variable array max sizes
    static const int maxNel = 10; // electrons
    static const int maxNmu = 10; // muons
    static const int maxNjet = 100; // jets
    static const int maxDsCand = 500; // D* candidates
    static const int maxDCand = 50000; // D+ candidates

    // Ntuple variables (their description can be found also in Analyzer/src/Analyzer.cc)
    //[N] means that this is fixed size array with N elements
    Int_t           evRunNumber; // run number
    Int_t           evEventNumber; // event number

    Float_t         metPx; // missing transverse energy x component
    Float_t         metPy; // missing transverse energy y component

    Int_t           Nel; // number of electrons
    Float_t         el_pt[maxNel];   //[Nel] electron pT
    Float_t         el_eta[maxNel];   //[Nel] electron eta
    Float_t         el_isolat[maxNel];   //[Nel] electron isolation
    Float_t         el_phi[maxNel];   //[Nel] electron phi
    Float_t         el_mt[maxNel];   //[Nel] missing transverse mass built from electron pT and missing transverse energy
    Float_t         el_convDist[maxNel];   //[Nel] (not used) electron conversion distance
    Float_t         el_convDcot[maxNel];   //[Nel] (not used) electron conversion cotangent
    Int_t           el_mh[maxNel];   //[Nel] electron missing hits number

    Int_t           Nmu; // number of muons
    Float_t         mu_pt[maxNmu];   //[Nmu] muon pT
    Float_t         mu_eta[maxNmu];   //[Nmu] muon eta
    Float_t         mu_isolat[maxNmu];   //[Nmu] muon isolation
    Float_t         mu_phi[maxNmu];   //[Nmu] muon phi
    Float_t         mu_mt[maxNmu];   //[Nmu] missing transverse mass built from muon pT and missing transverse energy
    Int_t           mu_ValidHits[maxNmu];   //[Nmu] muon valid hits number
    Int_t           mu_PixelHits[maxNmu];   //[Nmu] muon pixel hits number
    Float_t         mu_imp[maxNmu];   //[Nmu] transverse impact parameter of the muon with respect to the beam spot
    Float_t         mu_nChi2[maxNmu];   //[Nmu] muon track number of degrees of freedom

    Float_t         pv_x; // (not used) x component of the primary vertex
    Float_t         pv_y; // (not used) y component of the primary vertex
    Float_t         pv_z; // (not used) z component of the primary vertex

    Int_t           DsCand; // number of D* candidates
    Float_t         dmDs[maxDsCand];   //[DsCand] mass difference M(D*) - M(D0)
    Float_t         vcx[maxDsCand];   //[DsCand] (not used) x component of distance between D0 vertex and slow pion track
    Float_t         vcy[maxDsCand];   //[DsCand] (not used) y component of distance between D0 vertex and slow pion track
    Float_t         vcz[maxDsCand];   //[DsCand] (not used) z component of distance between D0 vertex and slow pion track
    Float_t         ds_pt[maxDsCand];   //[DsCand] D* pT
    Float_t         ds_eta[maxDsCand];   //[DsCand] D* eta
    Float_t         ds_m[maxDsCand];   //[DsCand] D* mass
    Float_t         ds_phi[maxDsCand];   //[DsCand] D* phi
    Float_t         D0dl[maxDsCand];   //[DsCand] D0 decay length
    Float_t         D0dlEr[maxDsCand];   //[DsCand] (not used) D0 decay length uncertainty
    Float_t         D0dls[maxDsCand];   //[DsCand] D0 decay length significance
    Float_t         D0sv_x[maxDsCand];   //[DsCand] (not used) x component of D0 vertex
    Float_t         D0sv_y[maxDsCand];   //[DsCand] (not used) y component of D0 vertex
    Float_t         D0sv_z[maxDsCand];   //[DsCand] z component of D0 vertex

    Int_t           DCand; // number of D+ candidates
    Float_t         d_pt[maxDCand];   //[DCand] D+ pT
    Float_t         d_eta[maxDCand];   //[DCand] D+ eta
    Float_t         d_phi[maxDCand];   //[DCand] D+ phi
    Float_t         d_m[maxDCand];   //[DCand] D+ mass
    Float_t         Ddl[maxDCand];   //[DCand] D+ decay length
    Float_t         DdlEr[maxDCand];   //[DCand] (not used) D+ decay length uncertainty
    Float_t         Ddls[maxDCand];   //[DCand] D+ decay length significance
    Float_t         Dsv_x[maxDCand];   //[DCand] (not used) x component of D+ vertex
    Float_t         Dsv_y[maxDCand];   //[DCand] (not used) y component of D+ vertex
    Float_t         Dsv_z[maxDCand];   //[DCand] x component of D+ vertex

    Int_t           Njet; // number of jets
    Float_t         j_pt[maxNjet];   //[Njet] jet pT
    Float_t         j_eta[maxNjet];   //[Njet] jet eta
    Float_t         j_m[maxNjet];   //[Njet] jet mass
    Float_t         j_phi[maxNjet];   //[Njet] jet phi
    Int_t           j_DsInJ[maxNjet];   //[Njet] D* candidate associated with jet
    Int_t           j_DInJ[maxNjet];   //[Njet] D+ candidate associated with jet
    Int_t           Triggers; // trigger bits

    // variables for MC only (generator level)
    Int_t           mcEventType; // detailed description below
    Int_t           cGen; // number of charm quarks
    // mcEventType is type of event:
    // mcEventType = 1: D*+ W-
    // mcEventType = 2: D+ W-
    // mcEventType = 3: mu+ W-
    // mcEventType = 4: D*- W+
    // mcEventType = 5: D- W+
    // mcEventType = 6: mu- W+
    // mcEventType = 0: anything else
    // D mesons and muons are from charm
    // positive values for W decaying to muons, negative values for W decaying to electrons

    // List of branches (their names follow variable names with prefix b_)
    TBranch        *b_evRunNumber;   //!
    TBranch        *b_evEventNumber;   //!
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
    // for MC only only
    TBranch        *b_mcEventType;   //!
    TBranch        *b_cGen;   //!

    // constructor
    // argument: true for MC, false (default) for data
    ZTree(bool flagMC = false) : fChain(0), _flagMC(flagMC) { }

    // destructor
    virtual ~ZTree() { }

    // initialise with provided tree pointer
    virtual void    Init(TTree *tree);
};

// initialise with provided tree pointer
void ZTree::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

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
  if(_flagMC) fChain->SetBranchAddress("cGen", &cGen, &b_cGen);
}

#endif // #ifdef ZTree_wcharm_h

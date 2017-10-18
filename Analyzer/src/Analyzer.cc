// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc demo/Analyzer/src/Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue Jul 19 12:43:42 CEST 2016
// $Id$
//
//


// system include files
#include <memory>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
//added 
//
//------ ADD EXTRA HEADER FILES--------------------//
#include "math.h"
#include <vector>
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for vertices refit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"

// for muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

//for electron informaton
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

// MET
#include "DataFormats/METReco/interface/PFMET.h"

// jets
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/BTauReco/interface/JetTag.h>
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//MC
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//ROOT
#include "TH1.h"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>

using namespace std;
using namespace edm;
using namespace reco;
//
// class declaration
//


class Analyzer : public edm::EDAnalyzer {
  public:
    explicit Analyzer(const edm::ParameterSet&);
    ~Analyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    //
    void SelectEvt(const edm::Event& iEvent);
    void SelectEl(const edm::Event& iEvent);
    void SelectMu(const edm::Event& iEvent,const reco::BeamSpot &beamSpot);
    void SelectJet(const edm::Event& iEvent,const edm::EventSetup& iSetup,const reco::BeamSpot &beamSpot);\
    float CalcMt(float Phi_l,float Pt);
    void SelectMET(const edm::Event& iEvent);
    void ReconstrDstar(const Handle<reco::TrackCollection> &tracks,const ESHandle<TransientTrackBuilder> &theB, const reco::BeamSpot &beamSpot,const TLorentzVector &J);
    void ReconstrD(const Handle<reco::TrackCollection> &tracks,const ESHandle<TransientTrackBuilder> &theB,const reco::BeamSpot &beamSpot, const TLorentzVector &J);
    vector<float> CalcDLS(vector<TransientTrack> tracksVector,  TransientVertex CMSFittedVtx,const reco::BeamSpot &beamSpot);
    float GetSimilarity(float *v1, float *v2, ROOT::Math::SMatrix<float,3> M);
    void FindTriggerBits(const HLTConfigProvider& trigConf);
    void SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR);
    void PrintTriggerBits();
    // mc
    const reco::Candidate* GetFinalState(const reco::Candidate* particle, const int id);
    void SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles);
    bool CheckParents(const reco::Candidate *c1,const reco::Candidate *c2);

    TFile* file;
    TTree* tree;
    //General info
    unsigned int _evRunNumber;
    unsigned int _evEventNumber;
    int Nevt;
    int Nsel;
    int _flagGEN;
    // Jet correction label
    std::string mJetCorr;
    // Missing momentum
    float metPx;
    float metPy;
    float metPt;
    float metPhi;

    //Electrons
    static const int Maxel=100;
    int Nel;
    float el_pt[Maxel];
    float el_eta[Maxel];
    float el_isolat[Maxel];
    float el_theta[Maxel];
    float el_phi[Maxel];
    float el_mt[Maxel];
    ///
    float el_convDist[Maxel];
    float el_convDcot[Maxel];
    int el_mh[Maxel];

    //Muons
    static const int Maxmu=100;
    int Nmu;
    float mu_pt[Maxmu];
    float mu_eta[Maxmu];
    float mu_isolat[Maxmu];
    float mu_theta[Maxmu];
    float mu_phi[Maxmu];
    float mu_mt[Maxmu];
    ///extras
    int mu_ValidHits[Maxmu];
    int mu_PixelHits[Maxmu];
    float mu_imp[Maxmu];
    float mu_nChi2[Maxmu];

    // primary vertex
    float pv_x, pv_y, pv_z;

    // jets
    static const int Maxjet=100;
    int Njet;
    float j_pt[Maxjet];
    float j_eta[Maxjet];
    float j_phi[Maxjet];
    float j_m[Maxjet];
    int j_DsInJ[Maxjet];
    int j_DInJ[Maxjet];
    //D*
    int DsCand;
    static const int MaxDs=1000;
    float ds_pt[MaxDs];
    float ds_eta[MaxDs];
    float ds_phi[MaxDs];
    float ds_m[MaxDs];
    float vcx[MaxDs];
    float vcy[MaxDs];
    float vcz[MaxDs];
    float dmDs[MaxDs];
    //d0
    float D0dl[MaxDs];
    float D0dlEr[MaxDs];
    float D0dls[MaxDs];
    // secondary vertex
    float D0sv_x[MaxDs];
    float D0sv_y[MaxDs];
    float D0sv_z[MaxDs];
    // D+-
    int DCand;
    static const int MaxD=20000;
    float d_pt[MaxD];
    float d_eta[MaxD];
    float d_phi[MaxD];
    float d_m[MaxD];
    float Ddl[MaxD];
    float DdlEr[MaxD];
    float Ddls[MaxD];
    // secondary vertex
    float Dsv_x[MaxD];
    float Dsv_y[MaxD];
    float Dsv_z[MaxD];

    int _triggers;
    edm::InputTag _inputTagTriggerResults;
    std::vector<std::vector<int> > _vecTriggerBits;
    std::vector<std::string> _vecTriggerNames;
    edm::InputTag _inputTagMCgen;


    //mc
    int mcDspWm=0;
    int mcDsmWp=0;
    int mcDpWm=0;
    int mcDmWp=0;
    int mcCmupWm=0;
    int mcCmumWp=0;
    int mcEventType=0;
    int cGen;

    int CWm=0;
    int CbWp=0;


    //
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig)

{
  setbuf(stdout, NULL);
  std::string fileout = iConfig.getParameter<std::string>("outFile");
  file = new TFile(fileout.c_str(), "recreate");
  tree = new TTree("tree", "spam");

  //h1  = new TH1F("h1","{p_{t}(#mu)}",300,0.,100.);
  //h2  = new TH1F("h2","{p_{t}(e)}",300,0.,100.);
  //now do what ever initialization is needed
  Nevt=0;
  Nsel=0;
  //general
  _inputTagTriggerResults = edm::InputTag("TriggerResults", "", "HLT");
  _inputTagMCgen = edm::InputTag("genParticles");
  _flagGEN = iConfig.getParameter<int>("gen");
  // jet correction label
  mJetCorr = "ak5PFL1FastL2L3Residual";


  tree->Branch("evRunNumber",&_evRunNumber,"evRunNumber/i");
  tree->Branch("evEventNumber",&_evEventNumber,"evEventNumber/i");

  // missing nrg
  tree->Branch("metPx",&metPx,"metPx/F");
  tree->Branch("metPy",&metPy,"metPy/F");

  //el
  tree->Branch("Nel",&Nel,"Nel/I");
  tree->Branch("el_pt",&el_pt,"el_pt[Nel]/F");
  tree->Branch("el_eta",&el_eta,"el_eta[Nel]/F");
  tree->Branch("el_isolat",&el_isolat,"el_isolat[Nel]/F");
  tree->Branch("el_phi",&el_phi,"el_phi[Nel]/F");
  tree->Branch("el_mt",&el_mt,"el_mt[Nel]/F");
  ///extras
  tree->Branch("el_convDist",&el_convDist,"el_convDist[Nel]/F");
  tree->Branch("el_convDcot",&el_convDcot,"el_convDcot[Nel]/F");
  tree->Branch("el_mh",&el_mh,"el_mh[Nel]/I");


  //mu
  tree->Branch("Nmu",&Nmu,"Nmu/I");
  tree->Branch("mu_pt",&mu_pt,"mu_pt[Nmu]/F");
  tree->Branch("mu_eta",&mu_eta,"mu_eta[Nmu]/F");
  tree->Branch("mu_isolat",&mu_isolat,"mu_isolat[Nmu]/F");
  tree->Branch("mu_phi",&mu_phi,"mu_phi[Nmu]/F");
  tree->Branch("mu_mt",&mu_mt,"mu_mt[Nmu]/F");
  ///extras
  tree->Branch("mu_ValidHits",&mu_ValidHits,"mu_ValidHits[Nmu]/I");
  tree->Branch("mu_PixelHits",&mu_PixelHits,"mu_PixelHits[Nmu]/I");
  tree->Branch("mu_imp",&mu_imp,"mu_imp[Nmu]/F");
  tree->Branch("mu_nChi2",&mu_nChi2,"mu_nChi2[Nmu]/F");

  // primaru vertex
  tree->Branch("pv_x",&pv_x,"pv_x/F");
  tree->Branch("pv_y",&pv_y,"pv_y/F");
  tree->Branch("pv_z",&pv_z,"pv_z/F");

  // D*
  tree->Branch("DsCand",&DsCand,"DsCand/I");
  tree->Branch("dmDs",&dmDs,"dmDs[DsCand]/F");
  tree->Branch("vcx",&vcx,"vcx[DsCand]/F");
  tree->Branch("vcy",&vcy,"vcy[DsCand]/F");
  tree->Branch("vcz",&vcz,"vcz[DsCand]/F");
  tree->Branch("ds_pt",&ds_pt,"ds_pt[DsCand]/F");
  tree->Branch("ds_eta",&ds_eta,"ds_eta[DsCand]/F");
  tree->Branch("ds_m",&ds_m,"ds_m[DsCand]/F");
  tree->Branch("ds_phi",&ds_phi,"ds_phi[DsCand]/F");
  // D0
  //	tree->Branch("D0Cand",&D0Cand,"D0Cand/I");
  tree->Branch("D0dl",&D0dl,"D0dl[DsCand]/F");
  tree->Branch("D0dlEr",&D0dlEr,"D0dlEr[DsCand]/F");
  tree->Branch("D0dls",&D0dls,"D0dls[DsCand]/F");
  // secondary vertex
  tree->Branch("D0sv_x",D0sv_x,"D0sv_x[DsCand]/F");
  tree->Branch("D0sv_y",D0sv_y,"D0sv_y[DsCand]/F");
  tree->Branch("D0sv_z",D0sv_z,"D0sv_z[DsCand]/F");
  // D+-
  tree->Branch("DCand",&DCand,"DCand/I");
  tree->Branch("d_pt",&d_pt,"d_pt[DCand]/F");
  tree->Branch("d_eta",&d_eta,"d_eta[DCand]/F");
  tree->Branch("d_phi",&d_phi,"d_phi[DCand]/F");
  tree->Branch("d_m",&d_m,"d_m[DCand]/F");
  tree->Branch("Ddl",&Ddl,"Ddl[DCand]/F");
  tree->Branch("DdlEr",&DdlEr,"DdlEr[DCand]/F");
  tree->Branch("Ddls",&Ddls,"Ddls[DCand]/F");
  // secondary vertex
  tree->Branch("Dsv_x",Dsv_x,"Dsv_x[DCand]/F");
  tree->Branch("Dsv_y",Dsv_y,"Dsv_y[DCand]/F");
  tree->Branch("Dsv_z",Dsv_z,"Dsv_z[DCand]/F");


  //jets
  tree->Branch("Njet",&Njet,"Njet/I");
  tree->Branch("j_pt",&j_pt,"j_pt[Njet]/F");
  tree->Branch("j_eta",&j_eta,"j_eta[Njet]/F");
  tree->Branch("j_m",&j_m,"j_m[Njet]/F");
  tree->Branch("j_phi",&j_phi,"j_phi[Njet]/F");
  tree->Branch("j_DsInJ",&j_DsInJ,"j_DsInJ[Njet]/I");
  tree->Branch("j_DInJ",&j_DInJ,"j_DInJ[Njet]/I");
  //triggers
  tree->Branch("Triggers", &_triggers, "Triggers/I");
  _vecTriggerNames.push_back("HLT_IsoMu24");
  _vecTriggerNames.push_back("HLT_Ele27");
  _vecTriggerNames.push_back("Ele32");

  // MC generated info
  if(_flagGEN)
  {
  /*	tree->Branch("mcDsmWp",&mcDsmWp,"mcDsmWp/I");
        tree->Branch("mcDpWm",&mcDpWm,"mcDpWm/I");
        tree->Branch("mcDmWp",&mcDmWp,"mcDmWp/I");
        tree->Branch("mcCmupWm",&mcCmupWm,"mcCmupWm/I");
        tree->Branch("mcCmumWp",&mcCmumWp,"mcCmumWp/I");*/
  tree->Branch("mcEventType",&mcEventType,"mcEventType/I");
  tree->Branch("cGen",&cGen,"cGen/I");
  /*tree->Branch("mcDspWm",&mcDspWm,"mcDspWm/I");
        tree->Branch("mcDspWm",&mcDspWm,"mcDspWm/I");
        tree->Branch("mcDspWm",&mcDspWm,"mcDspWm/I");
        tree->Branch("mcDspWm",&mcDspWm,"mcDspWm/I"); */
  }
}


Analyzer::~Analyzer()
{


  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  file->cd();
  tree->Write();
  file->Close();
  if(this->_flagGEN)
  {
    printf("D*+ W-: %f%%\n",100.*mcDspWm/CWm);
    printf("D*- W+: %f%%\n",100.* mcDsmWp/CbWp);
    printf("D+ W-: %f%%\n",100.*mcDpWm/CWm);
    printf("D- W+: %f%%\n",100.*mcDmWp/CbWp);
    printf("mu+ W-: %f%%\n",100.*mcCmupWm/CWm);
    printf("mu- W+: %f%%\n",100.*mcCmumWp/CbWp);
    printf("c W-   %d, cbar W+     %d \n",CWm,CbWp);
  }


  printf("Processed %d events, selected %d\n", Nevt, Nsel);
}


//
// member functions
//


void Analyzer::SelectEvt(const edm::Event& iEvent){
  _evRunNumber=iEvent.id().event();
  //cout << "Event: " << evt << endl;
  _evEventNumber=iEvent.id().run();
}

void Analyzer::SelectMET(const edm::Event& iEvent){
  edm::Handle<edm::View<reco::PFMET> > pfmets; // Handle<reco::PFMET> mets;  didn't work for some reason
  iEvent.getByLabel("pfMet",pfmets);
  metPx = (pfmets->front()).px();
  metPy = (pfmets->front()).py();
  metPt=(pfmets->front()).pt();
  metPhi=(pfmets->front()).phi();
}

const reco::Candidate* Analyzer::GetFinalState(const reco::Candidate* particle, const int id){
  for(unsigned int i = 0; i < particle->numberOfDaughters(); i++)
  {
    const reco::Candidate* daughter = particle->daughter(i);
    if(daughter->pdgId() == id )
      return daughter;
    const reco::Candidate* result = GetFinalState(daughter, id);
    if(result)
      return result;
  }
  return NULL;
}

bool Analyzer::CheckParents(const reco::Candidate *c1,const reco::Candidate *c2){
  if (c1->numberOfMothers()!=c2->numberOfMothers())
    return false;
  else{
    for(unsigned int x=0;x<c2->numberOfMothers();x++){
      if (c2->mother(x)!=c1->mother(x))
        return false;
    }
  }
  return true;
}

void Analyzer::SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles){
  const reco::Candidate* Wp = NULL;
  const reco::Candidate* Wm = NULL;
  const reco::Candidate* Lp = NULL;
  const reco::Candidate* Lm = NULL;
  const reco::Candidate* Nu = NULL;
  const reco::Candidate* Nub = NULL;
  const reco::Candidate* C = NULL;
  const reco::Candidate* Cb = NULL;

  const reco::Candidate* Dp = NULL;
  const reco::Candidate* Dm = NULL;
  const reco::Candidate* Dsp = NULL;
  const reco::Candidate* Dsm = NULL;
  const reco::Candidate* CMup = NULL;
  const reco::Candidate* CMum = NULL;

  const reco::Candidate* Kp = NULL;
  const reco::Candidate* Km = NULL;
  const reco::Candidate* Pi1p = NULL;
  const reco::Candidate* Pi1m = NULL;
  const reco::Candidate* Pi2p = NULL;
  const reco::Candidate* Pi2m = NULL;
  const reco::Candidate* D0 = NULL;
  const reco::Candidate* D0b = NULL;
  mcEventType=0;
  cGen=0;
  int WGen=0;
  //const vector <reco::Candidate*> Cwcand;
  bool sign=false;
  const reco::Candidate* W = NULL;

  for(unsigned int p = 0; p < genParticles->size(); p++) {
    const reco::Candidate* prt = &genParticles->at(p);
    if (abs(prt->pdgId())==24 && prt->status()==3){
      WGen+=1;
      sign=((prt->pdgId())<0);
      W= (sign ? Wm : Wp);
      W = prt;
    }
  }
  //	printf("1 prnt\n");

  if (WGen!=1 ){
    printf("Achtung! Number of W in event isn't ONE! This event will be skipped.\n");
    return;
  }
  const reco::Candidate*& c = (sign ? C : Cb);

  // check if the numbers of c(cbar) is odd
  for(unsigned int p = 0; p < genParticles->size(); p++) {
    const reco::Candidate* prt = &genParticles->at(p);
    if ((prt->pdgId()==(sign? 4: -4)) && prt->status()==3){
      if (CheckParents(W,prt) && prt->pt()>25. && abs(prt->eta())<2.5){
        c=prt;
        cGen++;
      }
    }
  }
  //	printf("charm %p \n",c);
  //	printf("2 prnt\n");
  if (cGen%2==0 ){
    //printf("Warning! Number of Charm quarks in event is even! This event will be skipped.\n");
    return;
  }

  // in case the numbers of c(cbar) is odd and single W are present in event
  if ((cGen%2)!=0 && WGen==1 ){

    const reco::Candidate*& l = (sign ? Lm : Lp);
    const reco::Candidate*& nu = (sign ? Nub : Nu);
    const reco::Candidate*& D = (sign ? Dp : Dm);
    const reco::Candidate*& Ds = (sign ? Dsp : Dsm);
    const reco::Candidate*& cmu = (sign ? CMup : CMum);
    const reco::Candidate*& K = (sign ? Km : Kp);
    const reco::Candidate*& Pi1 = (sign ? Pi1p : Pi1m);
    const reco::Candidate*& Pi2 = (sign ? Pi2p : Pi2m);
    const reco::Candidate*& d0 = (sign ? D0 : D0b);


    // find W -> lnu decay
    for(unsigned int d = 0; d < W->numberOfDaughters(); d++) {
      const reco::Candidate* daughter = W->daughter(d);
      if(daughter->status() != 3)
        continue;
      // electron
      if(daughter->pdgId() == (sign ? 11 : -11) && daughter->pt()>35.){
        if(l) printf("Error: multiple hard-scattering l\n");
        l = daughter;
      }
      if(daughter->pdgId() == (sign ? -12 : 12)){
        if(nu) printf("Error: multiple hard-scattering nu\n");
        nu = daughter;
      }
      // muon
      if(daughter->pdgId() == (sign ? 13 : -13) && daughter->pt()>25.){
        if(l) printf("Error: multiple hard-scattering l\n");
        l = daughter;
      }
      if(daughter->pdgId() == (sign ? -14 : 14)){
        if(nu) printf("Error: multiple hard-scattering nu\n");
        nu = daughter;
      }

    }

    if(!l || !nu)
      return;
    if (abs(l->eta())>2.1)
      return;
    //	printf("prnt 3: l %p\n",l);

    (sign? CWm++ : CbWp++);

    // get final states

    D=(sign ? GetFinalState( c, 411): GetFinalState( c, -411)); //411 D+
    Ds=(sign ? GetFinalState( c, 413): GetFinalState( c, -413)); //413 D*+
    cmu=(sign ? GetFinalState( c, -13): GetFinalState( c, 13)); // 13 mu-
    //	printf("prnt 4 D: %p Ds: %p cmu: %p\n",D,Ds,cmu);
    if (D){
      //	printf("Found D\n");
      if (D->numberOfDaughters()==3){

        //		printf("3 daughters: %d %d %d\n", D->daughter(0)->pdgId(), D->daughter(1)->pdgId(), D->daughter(2)->pdgId());
        for(unsigned int d = 0; d < D->numberOfDaughters(); d++){
          const reco::Candidate* prt = D->daughter(d);
          if (prt->pdgId()==(sign ? -321 : 321)) // 321  K+
            K=prt;

          if (prt->pdgId()==(sign ? 211 : -211) && !Pi1) // 211 Pi+
            Pi1=prt;
          if (prt->pdgId()==(sign ? 211 : -211) && (*&prt!= *&Pi1))
            Pi2=prt;
        }
        //			printf("3 daughters: %p %p %p\n", K, Pi1, Pi2);
        if(K && Pi1 && Pi2){ //C W-
          //printf("3 dauhters found\n");

          if(sign){
            //	printf("Found D+ W-\n");
            if(Km->pdgId()==-321 && Pi1p->pdgId()==211 && Pi2p->pdgId()==211){
              mcEventType=2;
              mcDpWm++;
            }
          }
          // Cbar W+
          else{
            //		printf("Found D W+\n");
            if(Kp->pdgId()==321 && Pi1m->pdgId()==-211 && Pi2m->pdgId()==-211){
              mcEventType=5;
              mcDmWp++;
            }
          }// W+
        }// all tracks are ok
      }//num of candiadate == 3
    }

    //	printf("c: %p\n", c);
    if (Ds){
      //		printf("Found Ds\n");
      if(Ds->numberOfDaughters()==2){

        for(unsigned int d = 0; d < Ds->numberOfDaughters(); d++){
          const reco::Candidate* prt = Ds->daughter(d);
          if (prt->pdgId()==(sign ? 421 : -421))  //421  D0
            d0=prt;
          if (prt->pdgId()==(sign ? 211 : -211))  // 211 Pi+
            Pi1=prt;
        }
        if(d0 && Pi1 && d0->numberOfDaughters()==2){
          for(unsigned int d2 = 0; d2 < d0->numberOfDaughters(); d2++){
            const reco::Candidate* prt2 = d0->daughter(d2);
            if ((prt2->pdgId())==(sign? 211 : -211))
              Pi2=prt2;
            if (prt2->pdgId()==(sign ? -321 : 321))
              K=prt2; // 321  K+
          }

          if (K && Pi1 && Pi2){
            if (sign){ // c W-
              //printf("Found D*+ W-\n");
              if(Km->pdgId()==-321 && Pi1p->pdgId()==211 && Pi2p->pdgId()==211 && D0->pdgId()==421){
                mcEventType=1;
                mcDspWm++;
              }
            }

            else{ // cbar W+
              //printf("Found D*- W+\n");
              if(Kp->pdgId()==321 && Pi1m->pdgId()==-211 && Pi2m->pdgId()==-211 && D0b->pdgId()==-421){
                mcEventType=4;
                mcDsmWp++;
              }
            }
          }// final candidates are ok
        }// number of daughter of D0 is 2
      } // number of daugthers of D* is 2
    }//D* end

    // c->mu
    if (cmu){
      // printf("Found c->mu\n");
      if (sign)
      { // c W-
        if (cmu->pdgId()==-13)
        {
          // printf("Found c->mu+ W-\n");
          mcEventType=3;
          mcCmupWm++;
        }
      }
      else
      { // cbar W+
        // printf("Found cbar->mu- W+\n");
        if (cmu->pdgId()==13)
        {
          mcEventType=6;
          mcCmumWp++;
        }
      }
    }
    if (abs(l->pdgId())==11)
      mcEventType*=-1;
  }/// W and charm have same parent
}



void Analyzer::SelectEl(const edm::Event& iEvent){

  Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons",electrons);
  Nel=0;
  float pt,eta, phi_l,iso,mt;

  // loop through electrons in the event
  for( reco::GsfElectronCollection::const_iterator it= electrons->begin(); it !=electrons->end(); it++) {
    pt=it->pt();
    eta=it->eta();
    phi_l=it->phi();
    iso=(it->dr03TkSumPt()+it->dr03EcalRecHitSumEt() +it->dr03HcalTowerSumEt())/(it->pt());

    if (Nel<Maxel){
      if (pt>20. && (abs(eta)<1.44 || (abs(eta)>1.57 && abs(eta)<2.5) ) && iso < 0.05){

        // calc mt
        mt=sqrt(2.0*pt*metPt*(1.0-cos(phi_l-metPhi)));

        el_pt[Nel]=pt*it->charge();
        el_eta[Nel]=eta;
        el_theta[Nel]=it->theta();
        el_phi[Nel]=phi_l;
        el_isolat[Nel]=iso;
        el_mt[Nel]=mt;
        el_convDist[Nel]=it->convDist();
        el_convDcot[Nel]=it->convDcot();
        el_mh[Nel]=((it->gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
        Nel++;
      }
    }
    else {
      printf("Maximum number of electrons %d reached, skipping the rest\n", Maxel);
      return;
    }

  }	// end the loop through electrons in the event

}

void Analyzer::SelectMu(const edm::Event& iEvent, const reco::BeamSpot &beamSpot){

  //Handle<reco::TrackCollection> gmuons;
  //iEvent.getByLabel("globalMuons", gmuons);
  edm::Handle<reco::MuonCollection> gmuons;
  iEvent.getByLabel("muons", gmuons);
  Nmu=0;
  float pt,eta, phi_l,iso,mt;//,phi_Em, mt, Em;

  math::XYZPoint point(beamSpot.position());

  for( reco::MuonCollection::const_iterator it= gmuons->begin(); it !=gmuons->end(); it++) {
    if((it->globalTrack()).isNull())
      continue;
    pt=it->pt();
    eta=it->eta();
    phi_l=it->phi();


    if (Nmu<Maxmu){
      if ((abs(eta)<2.4)){
        const reco::MuonPFIsolation& iso04 = it->pfIsolationR04();
        iso = (iso04.sumChargedParticlePt + iso04.sumNeutralHadronEt)/it->pt();
        
        mt=sqrt(2.0*pt*metPt*(1.0-cos(phi_l-metPhi)));

        mu_pt[Nmu]=pt*it->charge();
        mu_eta[Nmu]=eta;
        mu_theta[Nmu]=it->theta();
        mu_isolat[Nmu]=iso;
        mu_phi[Nmu]=phi_l;
        mu_mt[Nmu]=mt;

        /// extras
        mu_ValidHits[Nmu]=0;
        mu_PixelHits[Nmu]=0;
        const reco::TrackRef& track=it->globalTrack();
        const reco::HitPattern& p = (track)->hitPattern();
        /// if the hit is valid and in pixel
        for (int n=0; n<p.numberOfHits(); n++) {
          uint32_t hit = p.getHitPattern(n);

          if (p.validHitFilter(hit) && p.pixelHitFilter(hit)) mu_PixelHits[Nmu]++;
          if (p.validHitFilter(hit)) mu_ValidHits[Nmu]++;
        }
        ///

        mu_imp[Nmu]=track->dxy(point);
        mu_nChi2[Nmu]=(track->chi2())/(track->ndof());
        ///extras

        Nmu++;
      }

    }
    else {
      printf("Maximum number of muons %d reached, skipping the rest\n", Maxmu);
      return;
    }



  } //  end the loop through muons in the event


} // ends SelectMu func


vector<float> Analyzer::CalcDLS(vector<TransientTrack> tracksVector,TransientVertex CMSFittedVtx,const reco::BeamSpot &beamSpot)     {

  vector <float> DLValues;
  float Lx = 0.;
  float Ly = 0.;
  float Lz = CMSFittedVtx.position().z() - beamSpot.z0();

  float PxSum = 0.;
  float PySum = 0.;
  float PtSum = 0.;

  for (vector<TransientTrack>::iterator itSumP = tracksVector.begin(); itSumP != tracksVector.end(); itSumP++)
  {
    PxSum += (itSumP->track()).px();
    PySum += (itSumP->track()).py();
  }

  PtSum = sqrt(PxSum*PxSum + PySum*PySum);

  Lx = (CMSFittedVtx.position().x() - (beamSpot.x0() + beamSpot.dxdz()*Lz) );
  Ly = (CMSFittedVtx.position().y() - (beamSpot.y0() + beamSpot.dydz()*Lz) );

  float ProjectedDecayLength  = 0.;
  float ProjectedDecayLengthError = 0.;
  /// calculate DL xy :
  ProjectedDecayLength = (Lx*PxSum + Ly*PySum)/PtSum;
  // cout<<"DL xy proj = "<<ProjectedDecayLength<<endl;
  if (ProjectedDecayLength<0){
    DLValues.push_back(-10.);
    return DLValues;}
  typedef ROOT::Math::SMatrix<float,3>       SMatrix33;

  SMatrix33 DataVertexM;

  DataVertexM(0,0) = beamSpot.covariance(0, 0);
  DataVertexM(0,1) = beamSpot.covariance(0, 1);
  DataVertexM(0,2) = 0;

  DataVertexM(0,2) = 0;
  DataVertexM(1,0) = beamSpot.covariance(1, 0);
  DataVertexM(1,1) = beamSpot.covariance(1, 1);
  DataVertexM(1,2) = 0;
  DataVertexM(2,0) = 0;
  DataVertexM(2,1) = 0;
  DataVertexM(2,2) = 0;
  // }


  /// error DL calculation in matrices:
  /// *********************************  *********************************
  SMatrix33 SumVertexCovM;

  float MyVertexError[3][3];

  MyVertexError[0][0] = CMSFittedVtx.positionError().cxx();
  MyVertexError[0][1] = CMSFittedVtx.positionError().cyx();
  MyVertexError[0][2] = 0;
  MyVertexError[1][0] = CMSFittedVtx.positionError().cyx();
  MyVertexError[1][1] = CMSFittedVtx.positionError().cyy();
  MyVertexError[1][2] = 0;
  MyVertexError[2][0] = 0;
  MyVertexError[2][1] = 0;
  MyVertexError[2][2] = 0;


  for (int i =0; i<3; i++)        {
    for (int j=0; j<3; j++)         {

      SumVertexCovM[i][j] = MyVertexError[i][j] + DataVertexM[i][j];
    }
  }


  /// Extract DLxy sigma;
  float DLSigma = 0.;
  float PClusterV[3] = {PxSum, PySum, 0.};                                       /// vector components of the DL 3D

  DLSigma = GetSimilarity(PClusterV, PClusterV, SumVertexCovM);
  DLSigma = DLSigma/(PtSum*PtSum);
  DLSigma = sqrt(DLSigma);
  ProjectedDecayLengthError = DLSigma;

  /// *********************************  *********************************

  DLValues.push_back(ProjectedDecayLength);
  DLValues.push_back(ProjectedDecayLengthError);
  DLValues.push_back(ProjectedDecayLength/ProjectedDecayLengthError);

  /// couts:

  //  cout<<"Refit Vtx = "<<CMSFittedVtx.position().x()<<"    VtxY = "<<CMSFittedVtx.position().y()<<"        VtxZ = "<<CMSFittedVtx.position().z()<<endl;

  return DLValues;

}


float Analyzer::GetSimilarity(float *v1, float *v2, ROOT::Math::SMatrix<float,3> M)       {
  float results = 0.;
  float v_tmp[3] = {0.,0.,0.};

  for (int i =0 ; i<3; i++)       {
    for (int j=0;j<3;j++)           {
      v_tmp[i] += M[i][j]*v2[j];
    }
  }

  for (int z=0;z<3; z++)  {
    results +=      v1[z]*v_tmp[z];
  }

  return results;
}

void Analyzer::ReconstrDstar(const Handle<reco::TrackCollection> &tracks,const ESHandle<TransientTrackBuilder> &theB,
                             const reco::BeamSpot &beamSpot, const TLorentzVector &J){ // float RecMass,float deltaM, float DL, float PK, float Ppi,float Ps){
  vector<reco::TransientTrack> genralTracks = theB->build(tracks); //      declare new track builder for my new Transient track collection  ;
  float dR,MD0,MDstar,LXY;//,lxyx,lxyy,lxy,x0,y0;
  float vcD0[3]={0.,0.,0.};
  vector <float> DLvar;
  TLorentzVector vK;
  TLorentzVector vPi;
  //	DL=0.001;
  float mK=0.49367;
  float mPi=0.13957;
  int D0Cand=0;

  // D0 secondary vertex coordinates
  float d0sv[3] = {0., 0., 0.};

  int c1,c2;
  //x0 = beamSpot.x0();
  //  y0 = beamSpot.y0();
  for( reco::TrackCollection::const_iterator it1 = tracks->begin(); it1 != tracks->end(); it1++) {
    //cout << "it1 :" << (it1-tracks->begin()) << endl;
    if(it1->pt()<1. )
      continue;
    vK.SetPtEtaPhiM(it1->pt(),it1->eta(),it1->phi(),mK);
    if (vK.DeltaR(J)>0.3)
      continue;
    c1=it1->charge();

    for( reco::TrackCollection::const_iterator it2 = tracks->begin(); it2 != tracks->end(); it2++){

      if( it2->pt()<1. || *&it1==*&it2)
        continue;
      c2=it2->charge();
      if (c2*c1==1)
        continue;
      vPi.SetPtEtaPhiM(it2->pt(),it2->eta(),it2->phi(),mPi);
      if (vPi.DeltaR(J)>0.3)
        continue;
      TLorentzVector vD0;
      vD0=vPi+vK;
      MD0=(vD0).M();
      if (abs(MD0-1.864)>0.07)
        continue;
      TLorentzVector vPs;
      TLorentzVector Ds;

      LXY=-777.;

      for( reco::TrackCollection::const_iterator it3 = tracks->begin(); it3 != tracks->end(); it3++){
        // check charges of K and Pis. They must b same
        if(c1==it3->charge() || it3->pt()<0.3 || *&it3==*&it1 || *&it3==*&it2)
          continue;

        vPs.SetPtEtaPhiM(it3->pt(),it3->eta(),it3->phi(),mPi);
        if (vPs.DeltaR(J)>0.3)
          continue;
        Ds=vD0+vPs;
        MDstar=Ds.M();
        float dM=abs(MDstar-MD0);
        if (dM>0.17 || dM<0.135)
          continue;
        dR=vPs.DeltaR(vD0);


        if (dR>0.1)
          continue;
        // Check the slow pion originates from the same region as the D0

        if(LXY==-777.){
          // init vector for tracks of D0's products
          vector<TransientTrack> myD0Tracks;
          // result TransientTrack vector filling
          for(vector<TransientTrack>::iterator gt_trans = genralTracks.begin(); gt_trans !=genralTracks.end(); gt_trans++){
            const reco::TrackRef trackRef1 = (gt_trans->trackBaseRef()).castTo<reco::TrackRef>();
            /// check if candidates from loops over TrackCollaction
            /// exist in TransientTrack
            if ( &*it1 == trackRef1.get() || &*it2 == trackRef1.get() )	{
              TransientTrack  transientTrack1 = theB->build(trackRef1);
              myD0Tracks.push_back(transientTrack1);
            }
          }//ends filling vector TransientTrack  for D0

          TransientVertex D0Vertex;
          /// check if it's enough tracks to reconstruct vertex D0
          //cout << "MD0: " << MD0 << endl;
          if (myD0Tracks.size()<2 ){
            cout<<"Error! Not enough tracks to reconstruct vertex"<<endl;
            break;}
          AdaptiveVertexFitter theFitter;
          /// TransientVertex myVertexBS = theFitter.vertex(mytracks, vertexBeamSpot);                             // if you want the beam constraint
          D0Vertex = theFitter.vertex(myD0Tracks);                                                             // if you don't want the beam constraint
          // check validity of D0 vertex
          if (!D0Vertex.isValid()){
            //cout<<"Error! Vertex is invalid"<<endl;
            break;}

          DLvar=CalcDLS(myD0Tracks,D0Vertex,beamSpot);
          LXY=DLvar[0];

          /*	lxyx = (D0Vertex.position()).x()-x0;
                                                lxyy = (D0Vertex.position()).y()-y0;
                                                                /// direction of candidates regarding candidate
                                                lxy = sqrt(lxyx*lxyx + lxyy*lxyy);

//								float px,py,p;
//								px=it1->px()+it2->px(); 
//								py=it1->py()+it2->py();
//								p=sqrt(px*px+py*py);

                                                                /// calc decay length projected on momentum direction
                                                                /// LXY = LXY * Cos(LXY^p)
                                                lxy*=(vD0.Px()*lxyx+vD0.Py()*lxyy)/(lxy*vD0.Pt());
                                                        */
          // if (decay length of D0 condition)
          // Achtung! DEFINE DL!
          if (LXY<0.)
            break; // which exactly limit? DEFINE DL!
          ///	cout<<"Nazar's func "<<LXY<<"	 direct calc "<<lxy<<endl;
          D0Cand++;
          //locate D0 'vertex' using average coordinate of tracks
          vcD0[0] = 0.5*(it1->vx()+it2->vx()); vcD0[1] = 0.5*(it1->vy()+it2->vy()); vcD0[2] = 0.5*(it1->vz()+it2->vz());

          // secondary vertex
          d0sv[0] = D0Vertex.position().x();
          d0sv[1] = D0Vertex.position().y();
          d0sv[2] = D0Vertex.position().z();
        }// dl CHECK
        float vc2[3];
        vc2[0] = abs(vcD0[0]-it3->vx()); vc2[1]= abs(vcD0[1]- it3->vy()) ;vc2[2] = abs(vcD0[2]-it3->vz());
        vcx[DsCand]=vc2[0]; vcy[DsCand]=vc2[1];vcz[DsCand]=vc2[2];
        D0dl[DsCand]=LXY;
        D0dlEr[DsCand]=DLvar[1];
        D0dls[DsCand]=DLvar[2];

        //	if (sqrt(vc2[0]*vc2[0] + vc2[1]*vc2[1])>10. || vc2[2] >10.) // check if tracks start close
        //		continue;

        dmDs[DsCand]=MDstar-MD0;
        ds_pt[DsCand]=Ds.Pt()*it3->charge();
        ds_eta[DsCand]=Ds.Eta();
        ds_m[DsCand]=MDstar;
        ds_phi[DsCand]=Ds.Phi();

        // secondary vertex
        D0sv_x[DsCand] = d0sv[0];
        D0sv_y[DsCand] = d0sv[1];
        D0sv_z[DsCand] = d0sv[2];

        DsCand++;

        if (DsCand>=MaxDs){
          cout<<"Max number of Ds candidates exceeded! The next ones will be ignored!"<<endl;
          break;
        }


      } //TRACK 3 LOOP
      if (DsCand>=MaxDs){
        cout<<"Max number of Ds candidates exceeded! The next ones will be ignored!"<<endl;
        break;
      }
    } //track 2 loop ends
    if (DsCand>=MaxDs){
      cout<<"Max number of Ds candidates exceeded! The next ones will be skipped!"<<endl;
      break;
    }
  }	//track 1 loop ends
} // func ends

void Analyzer::ReconstrD(const Handle<reco::TrackCollection> &tracks,const ESHandle<TransientTrackBuilder> &theB,
                         const reco::BeamSpot &beamSpot, const TLorentzVector &J){ // float RecMass,float deltaM, float DL, float PK, float Ppi,float Ps){
  vector<reco::TransientTrack> genralTracks = theB->build(tracks); //      declare new track builder for my new Transient track collection  ;
  float MD,LXY;//,x0,y0,lxy,lxyx,lxyy;
  //	float vcD0[3]={0.,0.,0.};

  TLorentzVector vK;
  TLorentzVector vPi1;
  TLorentzVector vPi2;
  //	DL=0.001;
  float mK=0.49367;
  float mPi=0.13957;
  //D0Cand=0;

  int c1,c2;
  //	x0 = beamSpot.x0();
  //	y0 = beamSpot.y0();
  for( reco::TrackCollection::const_iterator it1 = tracks->begin(); it1 != tracks->end(); it1++) {

    if(it1->pt()<1.)
      continue;
    vK.SetPtEtaPhiM(it1->pt(),it1->eta(),it1->phi(),mK);
    if (vK.DeltaR(J)>0.3)
      continue;
    c1=it1->charge();

    for( reco::TrackCollection::const_iterator it2 = tracks->begin(); it2 != tracks->end(); it2++){
      c2=it2->charge();
      if (c2*c1==1  || it2->pt()<1.)
        continue;
      if(*&it1==*&it2)
        continue;
      vPi1.SetPtEtaPhiM(it2->pt(),it2->eta(),it2->phi(),mPi);
      if (vPi1.DeltaR(J)>0.3)
        continue;


      TLorentzVector vD;

      //for( reco::TrackCollection::const_iterator it3 = tracks->begin(); it3 != tracks->end(); it3++){
      for( reco::TrackCollection::const_iterator it3 = it2 + 1; it3 != tracks->end(); it3++){
        // check charges of K and Pis. They must b same
        if(c1==it3->charge() || it3->pt()<1. )
          continue;
        if (*&it3==*&it1 || *&it3==*&it2)
        {
          //printf("skipping\n");
          continue;
        }
        vPi2.SetPtEtaPhiM(it3->pt(),it3->eta(),it3->phi(),mPi);
        if (vPi2.DeltaR(J)>0.3)
          continue;
        vD=vK+vPi1+vPi2;
        MD=vD.M();
        //	float dM=abs(MD-1.8695);
        if (MD>2.2 || MD<1.6)
          continue;

        // init vector for tracks of D0's products
        vector<TransientTrack> myDTracks;
        // result TransientTrack vector filling
        for(vector<TransientTrack>::iterator gt_trans = genralTracks.begin(); gt_trans !=genralTracks.end(); gt_trans++){
          const reco::TrackRef trackRef1 = (gt_trans->trackBaseRef()).castTo<reco::TrackRef>();
          /// check if candidates from loops over TrackCollaction
          /// exist in TransientTrack
          if ( &*it1 == trackRef1.get() || &*it2 == trackRef1.get() || &*it3 == trackRef1.get())	{
            TransientTrack  transientTrack1 = theB->build(trackRef1);
            myDTracks.push_back(transientTrack1);
          }
        }//ends filling vector TransientTrack  for D0

        /// check if it's enough tracks to reconstruct vertex D0
        LXY=-1.;
        if (myDTracks.size()<3 ){
          cout<<"ERROR! Not enough tracks to reconstruct vertex."<<endl;
          continue;}
        AdaptiveVertexFitter theFitter;
        /// TransientVertex myVertexBS = theFitter.vertex(mytracks, vertexBeamSpot);                             // if you want the beam constraint
        auto DVertex = theFitter.vertex(myDTracks);                                                             // if you don't want the beam constraint
        // check validity of D0 vertex
        if (!DVertex.isValid()){
          //cout<<endl<<"ERROR! Vertex isn't valid."<<endl;
          continue;}
        const std::vector<AdaptiveVertexFitter::RefCountedVertexTrack> fittedTracks = DVertex.tracks();
        if(fittedTracks.size() != 3)
          continue;
        if((*fittedTracks[0]).weight() < 0.001 || (*fittedTracks[1]).weight() < 0.001 || (*fittedTracks[2]).weight() < 0.001)
          continue;
        vector <float> DLvar=CalcDLS(myDTracks,DVertex,beamSpot);
        LXY=DLvar[0];

        /*
                                        lxyx = (DVertex.position()).x()-x0;
                                        lxyy = (DVertex.position()).y()-y0;
                                                /// direction of candidates regarding candidate
                                        lxy = sqrt(lxyx*lxyx + lxyy*lxyy);

                                                /// calc decay length projected on momentum direction
                                                /// LXY = LXY * Cos(LXY^p)
                                        lxy*=(vD.Px()*lxyx+vD.Py()*lxyy)/(lxy*vD.Pt());
                                                                */
        // if (decay length of D0 condition)
        // Achtung! DEFINE DL!

        if (LXY<0.)
          continue; // which exactly limit? DEFINE DL!
        ///		Ddl[DCand]=LXY;
        ///		DCand++;
        //		cout<<"Nazar's func "<<LXY<<"	 direct calc "<<lxy<<endl;
        Ddl[DCand]=LXY;
        DdlEr[DCand]=DLvar[1];
        Ddls[DCand]=DLvar[2];
        d_m[DCand]=MD;

        d_pt[DCand]=vD.Pt()*it3->charge();
        d_eta[DCand]=vD.Eta();
        d_phi[DCand]=vD.Phi();

        // secondary vertex
        Dsv_x[DCand] = DVertex.position().x();
        Dsv_y[DCand] = DVertex.position().y();
        Dsv_z[DCand] = DVertex.position().z();

        DCand++;

        if (DCand>=MaxD){
          cout<<"Max number of D candidates exceeded! "<< DCand<<" The next ones will be ignored!"<<endl;
          break;
        }

      } //TRACK 3 LOOP

      if (DCand>=MaxD){
        cout<<"Max number of D candidates exceeded! "<< DCand<<"  The next ones will be ignored!"<<endl;
        break;
      }
    } //track 2 loop ends

    if (DCand>=MaxD){
      cout<<"Max number of D candidates exceeded! "<< DCand<<"  The next ones will be ignored!"<<endl;
      break;
    }
  }	//track 1 loop ends
} 

void Analyzer::SelectJet(const edm::Event& iEvent,const edm::EventSetup& iSetup,const reco::BeamSpot &beamSpot)
{
  // jets
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel("ak5PFJets", jets);
  Njet=0;
  // track 2 associate with jet
  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);

  /// Load the tools to work with vertices:

  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  TLorentzVector Jet;
  float pt,eta;
  DCand=0;
  DsCand=0;
  float dtmp=0;
  float dstmp=0;
  // Load jet energy correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector(mJetCorr, iSetup);
  // Loop over jets
  for(reco::PFJetCollection::const_iterator i_pfjet = jets->begin(); i_pfjet != jets->end(); i_pfjet++)
  {
    // Apply jet energy correction (JEC)
    double jec = corrector->correction(*i_pfjet, iEvent, iSetup);
    // copy original (uncorrected) jet;
    reco::PFJet corjet = *i_pfjet;
    // apply JEC
    corjet.scaleEnergy(jec);
    // check jet
    pt=corjet.pt();eta=corjet.eta();
    if(pt<25. || abs(eta)>2.5)
      continue;
    j_pt[Njet]=pt;
    j_eta[Njet]=eta;
    j_phi[Njet]=corjet.phi();
    j_m[Njet]=corjet.mass();

    Jet.SetPtEtaPhiE(pt,eta,corjet.phi(),corjet.energy());
    if (DsCand<MaxDs){
      ReconstrDstar(tracks,theB,beamSpot,Jet);
    }


    j_DsInJ[Njet]=DsCand-dstmp;
    dstmp=DsCand;
    if (DCand<MaxD){
      ReconstrD(tracks,theB,beamSpot,Jet);
    }

    j_DInJ[Njet]=DCand-dtmp;
    dtmp=DCand;
    //	cout<<"DsCand ="<<DsCand<<"  "<<" j_DsInJ="<<j_DsInJ[Njet]<<" DCand ="<<DCand<<"  j_DInJ="<<j_DInJ[Njet]<<" Njet="<<Njet<<endl<<endl;
    Njet++;
    if (Njet>=Maxjet){
      cout<<endl<<"Warning! Number of jets exceed maximum expected! Next ones will be skipped"<<endl;
      break;
    }
  }// jets loop
}

// returns vector of integers which are needed trigger bits
// (called in the beginning of each run)
void Analyzer::FindTriggerBits(const HLTConfigProvider& trigConf)
{
  _vecTriggerBits.clear();
  _vecTriggerBits.resize(_vecTriggerNames.size());
  std::vector<std::string> trigNames;
  trigNames = trigConf.triggerNames();
  for(unsigned int i = 0; i < trigNames.size(); i++)
  {
    std::string currentName = trigConf.triggerNames()[i];
    //printf("%5d  %s\n", i, currentName.c_str());
    for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
    {
      if(currentName.find(_vecTriggerNames[n]) != std::string::npos)
      {
        _vecTriggerBits[n].push_back(i);
      }
    }
  }
}

void Analyzer::PrintTriggerBits()
{
  printf("********* Trigger Bits: **********\n");
  for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
  {
    printf("%s: ", _vecTriggerNames[n].c_str());
    for(unsigned int i = 0; i < _vecTriggerBits[n].size(); i++)
      printf(" %d ", _vecTriggerBits[n][i]);
    printf("\n");
  }
}

// fill trigger bits
void Analyzer::SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR)
{
  for(unsigned int i = 0; i < _vecTriggerBits.size(); i++)
  {
    int status = 0;
    for(unsigned int j = 0; j < _vecTriggerBits[i].size(); j++)
    {
      status = status || HLTR->accept(_vecTriggerBits[i][j]);
    }
    //printf("TRIGGER:  %d %s\n", status, _vecTriggerNames[i].c_str());
    _triggers ^= (-status ^ _triggers) & (1 << i);
  }
  //printf("*************\n");
}



void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace reco;
  using namespace std;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
  // beamspot 2 calc decay length
  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  reco::BeamSpot beamSpot;
  if ( beamSpotHandle.isValid() )	{
    beamSpot = *beamSpotHandle;
  } else	{
    edm::LogInfo("Demo")
        << "No beam spot available from EventSetup \n";}

  // primary vertex (beamspot for x, y)
  Handle<reco::VertexCollection> primVertex;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", primVertex);
  reco::VertexCollection::const_iterator pv = primVertex->begin();
  pv_z = primVertex->begin()->z();
  pv_x = beamSpot.x0() + beamSpot.dxdz() * pv_z;
  pv_y = beamSpot.y0() + beamSpot.dydz() * pv_z;

  edm::Handle<reco::GenParticleCollection> genParticles;
  //edm::Handle<reco::MuonCollection> gmuons;
  //iEvent.getByLabel("muons", gmuons);

  SelectMET(iEvent);
  SelectEvt(iEvent);
  Nel=0;
  SelectEl(iEvent);
  SelectMu(iEvent,beamSpot);
  if(_flagGEN)
  {
    iEvent.getByLabel(_inputTagMCgen, genParticles);
    SelectMCGen(genParticles);
  }
  if (Nel>0 || Nmu>0 || mcEventType>0){
    SelectJet(iEvent, iSetup,beamSpot);
    if ((((DCand>0) || (Nmu>0) || (DsCand>0))&& Njet>0) || mcEventType>0)
    {
      //mc


      // triggers
      Handle<TriggerResults> HLTR;
      iEvent.getByLabel(_inputTagTriggerResults, HLTR);
      SelectTriggerBits(HLTR);
      // fill
      tree->Fill();
      Nsel++;
    }
  }//Nel = Nmu = Njet = DsCand = D0Cand = 0;
  //printf("before Fill()\n");
  //if (Nevt%1000==0)
  //cout<<endl<<"Nevt= "<<Nevt;

  Nevt++;
  if( (Nevt % 1000) == 0)
  {
    //printf("*****************************************************************\n");
    printf("************* NEVENTS = %d K, selected = %d *************\n", Nevt / 1000, Nsel);
    //printf("*****************************************************************\n");
  }
}



// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob()
{

}

// ------------ method called when starting to processes a run  ------------
void Analyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){
  HLTConfigProvider triggerConfig;
  bool changed = true;
  triggerConfig.init(iRun, iSetup, _inputTagTriggerResults.process(), changed);
  FindTriggerBits(triggerConfig);
  PrintTriggerBits();
}

// ------------ method called when ending the processing of a run  ------------
void 
Analyzer::endRun(edm::Run const& , edm::EventSetup const& )
{

}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);

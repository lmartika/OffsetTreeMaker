// -*- C++ -*-
//
// Package:    treemaker/OffsetTreeMaker
// Class:      OffsetTreeMaker
//
/**\class OffsetTreeMaker OffsetTreeMaker.cc treemaker/OffsetTreeMaker/plugins/OffsetTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  charles harrington
//         Created:  Mon, 09 Nov 2015 17:09:43 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "DataFormats/JetReco/interface/PFJet.h"
#include "parsePileUpJSON2.h"
#include "jetVetoMap.h"
#include <vector>
#include "TMath.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//root files
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

using namespace std;

const int ETA_BINS = 82;
const int ETA_BINS_GME = 18;
const int PHI_BINS_GME = 12;
const int MAXNPV = 100;
const int MAXJETS = 4;


float etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

float phibins[PHI_BINS_GME+1] = 
  {-M_PI, -M_PI*5./6., -M_PI*2./3., -M_PI/2., -M_PI/3., -M_PI/6., 0., M_PI/6., M_PI/3., M_PI/2., M_PI*2./3., M_PI*5./6., M_PI};

bool sortJetPt(const pair<int, float>& jet1, const pair<int, float>& jet2){ return jet1.second > jet2.second; }

class OffsetTreeMaker : public edm::one::EDAnalyzer<> {
  public:
    explicit OffsetTreeMaker(const edm::ParameterSet&);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    int getEtaIndex(float eta);
    enum Flavor{
      chm = 0, chu, nh, ne, hfh, hfe, lep, untrk, numFlavors, X //undefined
    };
    Flavor getFlavor(reco::PFCandidate::ParticleType id);

    int counter;
    TFile* root_file;
    TTree* tree;

    TH1F* h;
    TH1F* h1;
    TH2F* h2_GME; // 2d ET_eta_phi histo 
    TH2F* et_etaphi;
    TH2F* et_etaphi_chs;
    TRandom3* rand;

    TH1F* h_nFoundVertex;
    TH1F* h_bestweight;
    TH1F* h_bestweight1;

    int nEta;
    float energy[ETA_BINS], eRMS[ETA_BINS], et[ETA_BINS], etMED[ETA_BINS], etMEAN[ETA_BINS], etMEDchs[ETA_BINS], etMEANchs[ETA_BINS];
    float et_gme[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t ch_et_gme[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t f[numFlavors][ETA_BINS];  //energy fraction by flavor


    ULong64_t event;
    int run, lumi, bx;
//    int PileUpbx = 1;
    float mu;
//    float mubx;  // uncomment for lumibx
    float rho, rhoC0, rhoCC, rhoC;

    int mua[16];
    float puz[50];

    int nPVall, nPV;
    float pv_ndof[MAXNPV], pv_z[MAXNPV], pv_rho[MAXNPV];

    float ht;
    int nJets;
    float jet_eta[MAXJETS], jet_phi[MAXJETS], jet_pt[MAXJETS], jet_area[MAXJETS], jet_mass[MAXJETS];
    float jet_ch[MAXJETS], jet_nh[MAXJETS], jet_ne[MAXJETS], jet_hfh[MAXJETS], jet_hfe[MAXJETS], jet_lep[MAXJETS];

    vector<int> pf_type;
    vector<float> pf_pt, pf_eta, pf_phi, pf_et;
    //vector<string> eras;

//----------------------------Modification for hcal depth analysis----------------------------------

    vector<double>  nh_eta, nh_energy;
    vector<double>  nh_ECAL, nh_HCAL, nh_HO;
    vector<double>  nh_rawECAL, nh_rawHCAL, nh_rawHO;
    vector<double> nh_HCAL_depth1_fraction, nh_HCAL_depth2_fraction, nh_HCAL_depth3_fraction, nh_HCAL_depth4_fraction,
                   nh_HCAL_depth5_fraction, nh_HCAL_depth6_fraction, nh_HCAL_depth7_fraction;

    vector<double>  ne_eta, ne_energy;
    vector<double>  ne_ECAL, ne_HCAL, ne_HO;
    vector<double>  ne_rawECAL, ne_rawHCAL, ne_rawHO;

    float nh_sum_HCAL[ETA_BINS], nh_sum_ECAL[ETA_BINS], nh_sum_HO[ETA_BINS];
    float nh_sum_rawHCAL[ETA_BINS], nh_sum_rawECAL[ETA_BINS], nh_sum_rawHO[ETA_BINS];

    float nh_HCAL_depth1_fraction_SUM[ETA_BINS], nh_HCAL_depth2_fraction_SUM[ETA_BINS], nh_HCAL_depth3_fraction_SUM[ETA_BINS], nh_HCAL_depth4_fraction_SUM[ETA_BINS], nh_HCAL_depth5_fraction_SUM[ETA_BINS], nh_HCAL_depth6_fraction_SUM[ETA_BINS], nh_HCAL_depth7_fraction_SUM[ETA_BINS];

    float ne_sum_HCAL[ETA_BINS], ne_sum_ECAL[ETA_BINS], ne_sum_HO[ETA_BINS];
    float ne_sum_rawHCAL[ETA_BINS], ne_sum_rawECAL[ETA_BINS], ne_sum_rawHO[ETA_BINS];

    bool writeEnergyDeposition_;

//--------------------------------------------------------------------------------------


    TString RootFileName_, jetVetoMapFileName_, mapName2_; // mapName6_;
    vector<vector<double>> JetVetoMap ;
 
    string puFileName_, pubxFileName_, jet_type_; //era_,
    int numSkip_;
    bool isMC_, writeCands_, doL1L2L3Res_, dojetVetoMap_;

    //string jet_type;

    edm::EDGetTokenT< vector<reco::Vertex> > pvTag_;
    edm::EDGetTokenT< vector<reco::Track> > trackTag_;
    edm::EDGetTokenT< vector<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT< vector<reco::PFCandidate> > pfTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<double> rhoC0Tag_;
    edm::EDGetTokenT<double> rhoCCTag_;
    edm::EDGetTokenT<double> rhoCTag_;
    edm::EDGetTokenT< vector<reco::PFJet> > pfJetTag_;

};

OffsetTreeMaker::OffsetTreeMaker(const edm::ParameterSet& iConfig)
{
  numSkip_ = iConfig.getParameter<int> ("numSkip");
  RootFileName_ = iConfig.getParameter<string>("RootFileName");
  puFileName_ = iConfig.getParameter<string>("puFileName");
//  pubxFileName_ = iConfig.getParameter<string>("pubxFileName");  // uncomment for lumibx
  jetVetoMapFileName_ = iConfig.getParameter<string>("jetVetoMapFileName");
  mapName2_ = iConfig.getParameter<string>("mapName2");
//  mapName6_ = iConfig.getParameter<string>("mapName6");
  isMC_ = iConfig.getParameter<bool>("isMC");
  writeCands_ = iConfig.getParameter<bool>("writeCands");
  pvTag_ = consumes< vector<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
  trackTag_ = consumes< vector<reco::Track> >( iConfig.getParameter<edm::InputTag>("trackTag") );
  muTag_ = consumes< vector<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  pfTag_ = consumes< vector<reco::PFCandidate> >( iConfig.getParameter<edm::InputTag>("pfTag") );
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  rhoC0Tag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoC0Tag") );
  rhoCCTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCCTag") );
  rhoCTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCTag") );
  pfJetTag_ = consumes< vector<reco::PFJet> >( iConfig.getParameter<edm::InputTag>("pfJetTag") );
//  era_ = iConfig.getParameter<string>("era");
  jet_type_ = iConfig.getParameter<string>("jet_type");
  doL1L2L3Res_ = iConfig.getParameter<bool>("doL1L2L3Res");
  dojetVetoMap_ = iConfig.getParameter<bool>("dojetVetoMap");

  writeEnergyDeposition_ = iConfig.getParameter<bool>("writeEnergyDeposition");  // Modification for hcal depth

}

// ------------ method called once each job just before starting event loop  ------------
void  OffsetTreeMaker::beginJob() {

//  if (!isMC_){  // uncomment for lumibx
//    parsePileUpbxJSON2( pubxFileName_ );
//    cout << "ParsePileupbx Passed" <<  endl;
//  }
  root_file = new TFile(RootFileName_,"RECREATE");
  tree = new TTree("T","Offset Tree");

  counter = -1;
  h = new TH1F("mu", "mu", 200, 0, 100);
  rand = new TRandom3;
  h2_GME = new TH2F("GME_2D", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , -1*M_PI , M_PI);
  et_etaphi = new TH2F("et_etaphi", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);
  et_etaphi_chs = new TH2F("et_etaphi_chs", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);
  //Histos for checking CHS
  h_nFoundVertex = new TH1F ("nFoundVertex", "nFoundVertex", 3, -0.5, 2.5);
  h_bestweight = new TH1F ("bestweight", "bestweight", 50, 0., 1.);
  h_bestweight1 = new TH1F ("bestweight1", "bestweight1", 50, 0., 1.);
  
  if (dojetVetoMap_){
    jetVetoMap( jetVetoMapFileName_,mapName2_ ); //mapName6_
    JetVetoMap = getVetoMap();
  }

  if (!isMC_){
    //parsePileUpbxJSON2( pubxFileName_ ); // uncomment for lumibx
    //cout << "ParsePileupbx Passed" <<  endl;  // uncomment for lumibx
    parsePileUpJSON2( puFileName_ );
    cout << "ParsePU Passed" << endl;

    tree->Branch("run", &run, "run/I");
    tree->Branch("lumi", &lumi, "lumi/I");
    tree->Branch("bx", &bx, "bx/I");
    tree->Branch("event", &event, "event/l");
  }

  if (writeCands_) {
    tree->Branch("pf_type", "std::vector<int>",   &pf_type);
    tree->Branch("pf_pt",   "std::vector<float>", &pf_pt);
    tree->Branch("pf_eta",  "std::vector<float>", &pf_eta);
    tree->Branch("pf_phi",  "std::vector<float>", &pf_phi);
    tree->Branch("pf_et",  "std::vector<float>", &pf_et);
  }

  tree->Branch("mu", &mu, "mu/F");
//  tree->Branch("mubx", &mubx, "mubx/F");  // uncomment for lumibx
  tree->Branch("mua", mua, "mua[16]/I");
  tree->Branch("puz", puz, "puz[50]/F");

  tree->Branch("rho",   &rho,   "rho/F");
  tree->Branch("rhoC0", &rhoC0, "rhoC0/F");
  tree->Branch("rhoCC", &rhoCC, "rhoCC/F");

  tree->Branch("nPVall",  &nPVall, "nPVall/I");
  tree->Branch("nPV",     &nPV,    "nPV/I");
  tree->Branch("pv_ndof", pv_ndof, "pv_ndof[nPVall]/F");
  tree->Branch("pv_z",    pv_z,    "pv_z[nPVall]/F");
  tree->Branch("pv_rho",  pv_rho,  "pv_rho[nPVall]/F");

  tree->Branch("nEta",   &nEta,  "nEta/I");
  tree->Branch("energy", energy, "energy[nEta]/F");
  tree->Branch("et",     et,     "et[nEta]/F");
  tree->Branch("eRMS",   eRMS,   "eRMS[nEta]/F");

  tree->Branch("etMED",     etMED,     "etMED[nEta]/F");
  tree->Branch("etMEAN",    etMEAN,    "etMEAN[nEta]/F");
  tree->Branch("etMEDchs",     etMEDchs,     "etMEDchs[nEta]/F");
  tree->Branch("etMEANchs",    etMEANchs,    "etMEANchs[nEta]/F");

  tree->Branch("et_gme",       et_gme,     "et_gme[18][11]/F");
  tree->Branch("ch_et_gme",       ch_et_gme,     "ch_et_gme[18][11]/b");

  tree->Branch("fchm",   f[chm],   "fchm[nEta]/b");
  tree->Branch("fchu",   f[chu],   "fchu[nEta]/b");
  tree->Branch("fnh",    f[nh],    "fnh[nEta]/b");
  tree->Branch("fne",    f[ne],    "fne[nEta]/b");
  tree->Branch("fhfh",   f[hfh],   "fhfh[nEta]/b");
  tree->Branch("fhfe",   f[hfe],   "fhfe[nEta]/b");
  tree->Branch("flep",   f[lep],   "flep[nEta]/b");
  tree->Branch("funtrk", f[untrk], "funtrk[nEta]/b");


  tree->Branch("ht", &ht, "ht/F");
  tree->Branch("nJets",    &nJets,   "nJets/I");
  tree->Branch("jet_eta",  jet_eta,  "jet_eta[nJets]/F");
  tree->Branch("jet_phi",  jet_phi,  "jet_phi[nJets]/F");
  tree->Branch("jet_pt",   jet_pt,   "jet_pt[nJets]/F");
  tree->Branch("jet_area", jet_area, "jet_area[nJets]/F");

  tree->Branch("jet_ch",  jet_ch,  "jet_ch[nJets]/F");
  tree->Branch("jet_nh",  jet_nh,  "jet_nh[nJets]/F");
  tree->Branch("jet_ne",  jet_ne,  "jet_ne[nJets]/F");
  tree->Branch("jet_hfh", jet_hfh, "jet_hfh[nJets]/F");
  tree->Branch("jet_hfe", jet_hfe, "jet_hfe[nJets]/F");
  tree->Branch("jet_lep", jet_lep, "jet_lep[nJets]/F");

//----------------------------Modification for hcal depth analysis----------------------------------

  if (writeEnergyDeposition_){
    tree->Branch("nh_eta", "std::vector<double>", &nh_eta);
    tree->Branch("nh_energy", "std::vector<double>", &nh_energy);
    tree->Branch("nh_ECAL", "std::vector<double>", &nh_ECAL);
    tree->Branch("nh_HCAL", "std::vector<double>", &nh_HCAL);
    tree->Branch("nh_rawECAL", "std::vector<double>", &nh_rawECAL);
    tree->Branch("nh_rawHCAL", "std::vector<double>", &nh_rawHCAL);
    tree->Branch("nh_HO", "std::vector<double>", &nh_HO);
    tree->Branch("nh_rawHO", "std::vector<double>", &nh_rawHO);

    tree->Branch("nh_HCAL_depth1_fraction", "std::vector<double>", &nh_HCAL_depth1_fraction);
    tree->Branch("nh_HCAL_depth2_fraction", "std::vector<double>", &nh_HCAL_depth2_fraction);
    tree->Branch("nh_HCAL_depth3_fraction", "std::vector<double>", &nh_HCAL_depth3_fraction);
    tree->Branch("nh_HCAL_depth4_fraction", "std::vector<double>", &nh_HCAL_depth4_fraction);
    tree->Branch("nh_HCAL_depth5_fraction", "std::vector<double>", &nh_HCAL_depth5_fraction);
    tree->Branch("nh_HCAL_depth6_fraction", "std::vector<double>", &nh_HCAL_depth6_fraction);
    tree->Branch("nh_HCAL_depth7_fraction", "std::vector<double>", &nh_HCAL_depth7_fraction);

    tree->Branch("nh_sum_HCAL", nh_sum_HCAL, "nh_sum_HCAL[nEta]/F");
    tree->Branch("nh_sum_ECAL", nh_sum_ECAL, "nh_sum_ECAL[nEta]/F");
    tree->Branch("nh_sum_HO", nh_sum_HO, "nh_sum_HO[nEta]/F");
    tree->Branch("nh_sum_rawHCAL", nh_sum_rawHCAL, "nh_sum_rawHCAL[nEta]/F");
    tree->Branch("nh_sum_rawECAL", nh_sum_rawECAL, "nh_sum_rawECAL[nEta]/F");
    tree->Branch("nh_sum_rawHO", nh_sum_rawHO, "nh_sum_rawHO[nEta]/F");

    tree->Branch("nh_HCAL_depth1_fraction_SUM", nh_HCAL_depth1_fraction_SUM, "nh_HCAL_depth1_fraction_SUM[nEta]/F");
    tree->Branch("nh_HCAL_depth2_fraction_SUM", nh_HCAL_depth2_fraction_SUM, "nh_HCAL_depth2_fraction_SUM[nEta]/F");
    tree->Branch("nh_HCAL_depth3_fraction_SUM", nh_HCAL_depth3_fraction_SUM, "nh_HCAL_depth3_fraction_SUM[nEta]/F");
    tree->Branch("nh_HCAL_depth4_fraction_SUM", nh_HCAL_depth4_fraction_SUM, "nh_HCAL_depth4_fraction_SUM[nEta]/F");
    tree->Branch("nh_HCAL_depth5_fraction_SUM", nh_HCAL_depth5_fraction_SUM, "nh_HCAL_depth5_fraction_SUM[nEta]/F");
    tree->Branch("nh_HCAL_depth6_fraction_SUM", nh_HCAL_depth6_fraction_SUM, "nh_HCAL_depth6_fraction_SUM[nEta]/F");
    tree->Branch("nh_HCAL_depth7_fraction_SUM", nh_HCAL_depth7_fraction_SUM, "nh_HCAL_depth7_fraction_SUM[nEta]/F");

    tree->Branch("ne_eta", "std::vector<double>", &ne_eta);
    tree->Branch("ne_energy", "std::vector<double>", &ne_energy);
    tree->Branch("ne_ECAL", "std::vector<double>", &ne_ECAL);
    tree->Branch("ne_HCAL", "std::vector<double>", &ne_HCAL);
    tree->Branch("ne_rawECAL", "std::vector<double>", &ne_rawECAL);
    tree->Branch("ne_rawHCAL", "std::vector<double>", &ne_rawHCAL);
    tree->Branch("ne_HO", "std::vector<double>", &ne_HO);
    tree->Branch("ne_rawHO", "std::vector<double>", &ne_rawHO);

    tree->Branch("ne_sum_HCAL", ne_sum_HCAL, "ne_sum_HCAL[nEta]/F");
    tree->Branch("ne_sum_ECAL", ne_sum_ECAL, "ne_sum_ECAL[nEta]/F");
    tree->Branch("ne_sum_HO", ne_sum_HO, "ne_sum_HO[nEta]/F");
    tree->Branch("ne_sum_rawHCAL", ne_sum_rawHCAL, "ne_sum_rawHCAL[nEta]/F");
    tree->Branch("ne_sum_rawECAL", ne_sum_rawECAL, "ne_sum_rawECAL[nEta]/F");
    tree->Branch("ne_sum_rawHO", ne_sum_rawHO, "ne_sum_rawHO[nEta]/F");
  }

//--------------------------------------------------------------------------------------


}

// ------------ method called for each event  ------------
void OffsetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  counter++;
  if (counter%numSkip_ != 0) return;

//------------ Pileup ------------//
  //cout <<"\n\nEvent # :" << counter << endl;
  if (isMC_){
    edm::Handle< vector<PileupSummaryInfo> > pileups;
    iEvent.getByToken(muTag_, pileups);

    mu = pileups->at(1).getTrueNumInteractions();
    // mubx = mu;  // uncomment for lumibx
    int j = 0;
    std::vector<PileupSummaryInfo>::const_iterator pileupinfo;
    for(pileupinfo = pileups->begin(); pileupinfo != pileups->end(); ++pileupinfo){
      mua[j] = pileupinfo->getPU_NumInteractions();
      /*cout << "index = "<< j << endl;
      cout << "BX = "<< pileupinfo->getBunchCrossing()<<endl;
      cout << "mu = "<<pileupinfo->getTrueNumInteractions() <<endl;
      cout << "muactual = "<< pileupinfo->getPU_NumInteractions()<<endl;*/

      if (pileupinfo->getBunchCrossing() == 0){
        for(int i = 0; i != pileupinfo->getPU_NumInteractions() && i < 50 ; ++i){
          puz[i] = (pileupinfo->getPU_zpositions())[i]; 
        }
      } 
      ++j;
    }
  }
  else{
    run = int(iEvent.id().run());
    lumi = int(iEvent.getLuminosityBlock().luminosityBlock());
    bx = iEvent.bunchCrossing();
    event = iEvent.id().event();

    //mubx = getAvgPUbx( run, lumi, bx );  // uncomment for lumibx
    //cout << "Fetched PUbx = " << mubx << endl;
    mu = getAvgPU( run, lumi );
    //cout << "Fetched PU = " << mu << endl;
    if (mu==0) return;
  }  

//------------ Primary Vertices ------------//

  edm::Handle< vector<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPVall = primaryVertices->size();
  if(nPVall > MAXNPV) nPVall = MAXNPV ;
  nPV = 0;

  for (int i = 0; i != nPVall; ++i){
    reco::Vertex pv = primaryVertices->at(i);

    pv_ndof[i] = pv.ndof();
    pv_z[i] = pv.z();
    pv_rho[i] = pv.position().rho();

    if( !pv.isFake() && pv_ndof[i] >= 4.0 && fabs(pv_z[i]) <= 24.0 && fabs(pv_rho[i]) <= 2.0 )
      nPV++;
  }

//------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;

  edm::Handle<double> rhoC0Handle;
  iEvent.getByToken(rhoC0Tag_, rhoC0Handle);
  rhoC0 = *rhoC0Handle;

  edm::Handle<double> rhoCCHandle;
  iEvent.getByToken(rhoCCTag_, rhoCCHandle);
  rhoCC = *rhoCCHandle;

  edm::Handle<double> rhoCHandle;
  iEvent.getByToken(rhoCTag_, rhoCHandle);
  rhoC = *rhoCHandle;

//------------ PF Particles ------------//

  edm::Handle< vector<reco::PFCandidate> > pfCandidates;
  iEvent.getByToken(pfTag_, pfCandidates);

//----------------------------Modification for hcal depth analysis----------------------------------


  memset(nh_sum_HCAL, 0, sizeof(nh_sum_HCAL));
  memset(nh_sum_ECAL, 0, sizeof(nh_sum_ECAL));
  memset(nh_sum_HO, 0, sizeof(nh_sum_HO));
  memset(nh_sum_rawHCAL, 0, sizeof(nh_sum_rawHCAL));
  memset(nh_sum_rawECAL, 0, sizeof(nh_sum_rawECAL));
  memset(nh_sum_rawHO, 0, sizeof(nh_sum_rawHO));


  memset(nh_HCAL_depth1_fraction_SUM, 0, sizeof(nh_HCAL_depth1_fraction_SUM));
  memset(nh_HCAL_depth2_fraction_SUM, 0, sizeof(nh_HCAL_depth2_fraction_SUM));
  memset(nh_HCAL_depth3_fraction_SUM, 0, sizeof(nh_HCAL_depth3_fraction_SUM));
  memset(nh_HCAL_depth4_fraction_SUM, 0, sizeof(nh_HCAL_depth4_fraction_SUM));
  memset(nh_HCAL_depth5_fraction_SUM, 0, sizeof(nh_HCAL_depth5_fraction_SUM));
  memset(nh_HCAL_depth6_fraction_SUM, 0, sizeof(nh_HCAL_depth6_fraction_SUM));
  memset(nh_HCAL_depth7_fraction_SUM, 0, sizeof(nh_HCAL_depth7_fraction_SUM));

  memset(ne_sum_HCAL, 0, sizeof(ne_sum_HCAL));
  memset(ne_sum_ECAL, 0, sizeof(ne_sum_ECAL));
  memset(ne_sum_HO, 0, sizeof(ne_sum_HO));
  memset(ne_sum_rawHCAL, 0, sizeof(ne_sum_rawHCAL));
  memset(ne_sum_rawECAL, 0, sizeof(ne_sum_rawECAL));
  memset(ne_sum_rawHO, 0, sizeof(ne_sum_rawHO));
 
//------------------------------------------------------------------------------------------------------

  memset(energy, 0, sizeof(energy));    //reset arrays to zero
  memset(et, 0, sizeof(et));
  memset(etMED, 0, sizeof(etMED));
  memset(etMEAN, 0, sizeof(etMEAN));
  memset(etMEDchs, 0, sizeof(etMEDchs));
  memset(etMEANchs, 0, sizeof(etMEANchs));

  nEta = ETA_BINS;

  float eFlavor[numFlavors][ETA_BINS] = {};
  float e2[ETA_BINS] = {};  //energy squared
  int nPart[ETA_BINS] = {}; //number of particles per eta bin
  memset(et_gme,        0, sizeof(et_gme));
  memset(ch_et_gme,     0, sizeof(ch_et_gme));

  pf_type.clear(); pf_pt.clear(); pf_eta.clear(); pf_phi.clear(); pf_et.clear();

//----------------------------Modification for hcal depth analysis----------------------------------


  nh_eta.clear(); nh_energy.clear(); nh_ECAL.clear(); nh_HCAL.clear();
  nh_rawECAL.clear(); nh_rawHCAL.clear();
  nh_HO.clear(); nh_rawHO.clear();
  nh_HCAL_depth1_fraction.clear();  nh_HCAL_depth2_fraction.clear();  nh_HCAL_depth3_fraction.clear(); 
  nh_HCAL_depth4_fraction.clear();
  nh_HCAL_depth5_fraction.clear();  nh_HCAL_depth6_fraction.clear();  nh_HCAL_depth7_fraction.clear(); 

  ne_eta.clear(); ne_energy.clear(); ne_ECAL.clear(); ne_HCAL.clear();
  ne_rawECAL.clear(); ne_rawHCAL.clear();
  ne_HO.clear(); ne_rawHO.clear();
 
//------------------------------------------------------------------------------------------------------

  h2_GME->Reset();
  et_etaphi->Reset();
  et_etaphi_chs->Reset();

  vector<reco::PFCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
  for (i_pf = pfCandidates->begin(); i_pf != endpf; ++i_pf) {

    bool pfVeto = false;
    if (dojetVetoMap_){
      for(int j=0; j != int((JetVetoMap[0]).size()); ++j){
        if ((i_pf->eta()>= JetVetoMap[0][j]) && (i_pf->eta()<= JetVetoMap[1][j]) && (i_pf->phi()>= JetVetoMap[2][j]) && (i_pf->phi()<= JetVetoMap[3][j]))
          pfVeto=true;
      } 
    }
    if (pfVeto) continue;

    int etaIndex = getEtaIndex( i_pf->eta() );
    Flavor flavor = getFlavor( i_pf->particleId() );

    if (etaIndex == -1 || flavor == X) continue;

    bool attached1 = false;
    bool attached = false;

    reco::TrackRef pftrack( i_pf->trackRef() );

    if (flavor == chm && !pftrack.isNull() ) { //check charged hadrons ONLY
     // unsigned int iVertex = 0;
      //unsigned int index=0;
      unsigned int nFoundVertex = 0;
      float bestweight = 0;

      vector<reco::Vertex>::const_iterator i_pv, endpv = primaryVertices->end();
      for (i_pv = primaryVertices->begin(); i_pv != endpv ; ++i_pv) {
        
        if ( !i_pv->isFake() && i_pv->ndof() >= 4.0 && fabs(i_pv->z()) <= 24.0 && fabs(i_pv->position().rho())<=2.0 ) {

          float w = i_pv->trackWeight(pftrack);
          //cout << " weight = " << w <<endl;

          if (w > bestweight){
	    bestweight=w;
	//    iVertex=index;
	    nFoundVertex++;
          }
        //  ++index;

          reco::Vertex::trackRef_iterator i_vtxTrk, endvtxTrk = i_pv->tracks_end();
          for(i_vtxTrk = i_pv->tracks_begin(); i_vtxTrk != endvtxTrk && !attached1; ++i_vtxTrk) {
              
            reco::TrackRef vtxTrk(i_vtxTrk->castTo<reco::TrackRef>());
            if (vtxTrk == pftrack)
              attached1 = true;
          } 
        }
      }
      if (nFoundVertex > 0) attached = true;
      
      //cout << endl;
      //cout << "Best Weight, nFoundVertex, attached = "<< bestweight << "\t" << nFoundVertex << "\t" << attached1 << endl;
      //if(nFoundVertex != int(attached1)) cout << "  =======> Different decision " << endl ;
      if (!attached) flavor = chu; //unmatched charged hadron

      h_nFoundVertex->Fill(nFoundVertex);
      h_bestweight->Fill(bestweight);
      if (attached1) h_bestweight1->Fill(bestweight);
    }
    float e = i_pf->energy();

    energy[etaIndex] += e;
    et[etaIndex] += i_pf->et();
    eFlavor[flavor][etaIndex] += e;

    h2_GME->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
    et_etaphi->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
    if (flavor != chm) et_etaphi_chs->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
    e2[etaIndex] += (e*e);
    nPart[etaIndex] ++;

    if (writeCands_) {
      pf_type.push_back( static_cast<int>(flavor) );
      pf_pt.push_back( i_pf->pt() );
      pf_eta.push_back( i_pf->eta() );
      pf_phi.push_back( i_pf->phi() );
      pf_et.push_back( i_pf->et() );
    }

//----------------------------Modification for hcal depth analysis----------------------------------

    if (writeEnergyDeposition_){
      if (flavor == nh){ 
        // per PF candidate quantities
        nh_eta.push_back(i_pf->eta());
        nh_energy.push_back(i_pf->energy());
        nh_ECAL.push_back(i_pf->ecalEnergy());
        nh_HCAL.push_back(i_pf->hcalEnergy());
        nh_HO.push_back(i_pf->hoEnergy());
        nh_rawECAL.push_back(i_pf->rawEcalEnergy());
        nh_rawHCAL.push_back(i_pf->rawHcalEnergy());
        nh_rawHO.push_back(i_pf->rawHoEnergy());

        nh_HCAL_depth1_fraction.push_back(i_pf->hcalDepthEnergyFraction(1)); 
        nh_HCAL_depth2_fraction.push_back(i_pf->hcalDepthEnergyFraction(2)); 
        nh_HCAL_depth3_fraction.push_back(i_pf->hcalDepthEnergyFraction(3)); 
        nh_HCAL_depth4_fraction.push_back(i_pf->hcalDepthEnergyFraction(4)); 
        nh_HCAL_depth5_fraction.push_back(i_pf->hcalDepthEnergyFraction(5)); 
        nh_HCAL_depth6_fraction.push_back(i_pf->hcalDepthEnergyFraction(6)); 
        nh_HCAL_depth7_fraction.push_back(i_pf->hcalDepthEnergyFraction(7)); 

        //per eta slice quantities
        int etaIndex = getEtaIndex( i_pf->eta() );
        if (etaIndex == -1) continue;
        nh_sum_HCAL[etaIndex] += i_pf->hcalEnergy();
        nh_sum_ECAL[etaIndex] += i_pf->ecalEnergy();
        nh_sum_HO[etaIndex] += i_pf->hoEnergy();
        nh_sum_rawHCAL[etaIndex] += i_pf->rawHcalEnergy();
        nh_sum_rawECAL[etaIndex] += i_pf->rawEcalEnergy();
        nh_sum_rawHO[etaIndex] += i_pf->rawHoEnergy();

        nh_HCAL_depth1_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(1); 
        nh_HCAL_depth2_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(2); 
        nh_HCAL_depth3_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(3); 
        nh_HCAL_depth4_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(4); 
        nh_HCAL_depth5_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(5); 
        nh_HCAL_depth6_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(6); 
        nh_HCAL_depth7_fraction_SUM[etaIndex] += i_pf->hcalEnergy() * i_pf->hcalDepthEnergyFraction(7); 
      }

      if (flavor == ne){
        // per PF candidate quantities
        ne_eta.push_back(i_pf->eta());
        ne_energy.push_back(i_pf->energy());
        ne_ECAL.push_back(i_pf->ecalEnergy());
        ne_HCAL.push_back(i_pf->hcalEnergy());
        ne_HO.push_back(i_pf->hoEnergy());
        ne_rawECAL.push_back(i_pf->rawEcalEnergy());
        ne_rawHCAL.push_back(i_pf->rawHcalEnergy());
        ne_rawHO.push_back(i_pf->rawHoEnergy());

        //per eta slice quantities
        int etaIndex = getEtaIndex( i_pf->eta() );
        if (etaIndex == -1) continue;
        ne_sum_HCAL[etaIndex] += i_pf->hcalEnergy();
        ne_sum_ECAL[etaIndex] += i_pf->ecalEnergy();
        ne_sum_HO[etaIndex] += i_pf->hoEnergy();
        ne_sum_rawHCAL[etaIndex] += i_pf->rawHcalEnergy();
        ne_sum_rawECAL[etaIndex] += i_pf->rawEcalEnergy();
        ne_sum_rawHO[etaIndex] += i_pf->rawHoEnergy();

      }
    }
//------------------------------------------------------------------------------------------------------

  }

  for (int ieta = 1; ieta != (ETA_BINS_GME+1); ++ieta){
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      et_gme[ieta-1][iphi-1] =  h2_GME->GetBinContent(ieta, iphi);
      ch_et_gme[ieta-1][iphi-1] =  char(min(255, int( h2_GME->GetBinContent(ieta, iphi)/0.1)));
    }
  } 
  for (int ieta = 1; ieta != (ETA_BINS+1); ++ieta){
    vector<double> x, x_chs;
    double et_sum = 0., et_sum_chs = 0.;
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      //if (ieta == ETA_BINS/2) cout << " iphi " << et_etaphi->GetBinContent(ieta, iphi) <<endl;
      x.push_back(et_etaphi->GetBinContent(ieta, iphi));
      et_sum += et_etaphi->GetBinContent(ieta, iphi);
      x_chs.push_back(et_etaphi_chs->GetBinContent(ieta, iphi));
      et_sum_chs += et_etaphi_chs->GetBinContent(ieta, iphi);
    }
    sort(x.begin(),x.end());
    sort(x_chs.begin(),x_chs.end()); 
    etMED[ieta-1]  = 0.5 * (x[5]+x[6]);
    etMEAN[ieta-1] = float(et_sum)/float(PHI_BINS_GME);
    etMEDchs[ieta-1]  = 0.5 * (x_chs[5]+x_chs[6]);
    etMEANchs[ieta-1] = float(et_sum_chs)/float(PHI_BINS_GME);
  }


//------------ Tracks ------------//

  edm::Handle< vector<reco::Track> > tracks;
  iEvent.getByToken(trackTag_, tracks);

  vector<reco::Track>::const_iterator i_trk, endtrk = tracks->end();
  for (i_trk = tracks->begin(); i_trk != endtrk; ++i_trk) {

    if ( !i_trk->quality(reco::Track::tight) ) continue;
    bool matched = false;

    vector<reco::PFCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
    for (i_pf = pfCandidates->begin();  i_pf != endpf && !matched; ++i_pf) {

      if ( &(*i_trk) == i_pf->trackRef().get() )
        matched = true;      
    }
    if (matched) continue;

    int etaIndex = getEtaIndex( i_trk->eta() );
    if (etaIndex == -1) continue;

    float e = i_trk->p();

    energy[etaIndex] += e;
    et[etaIndex] += i_trk->pt();
    eFlavor[untrk][etaIndex] += e;
    e2[etaIndex] += (e*e);
    nPart[etaIndex] ++;
  }

  for (int i=0; i != nEta; ++i){

    for (int flav = 0; flav != numFlavors; ++flav){
      UChar_t f_value; float eFlav = eFlavor[flav][i]; float E = energy[i];

      if (eFlav == 0)      f_value = 0;
      else if (eFlav == E) f_value = 255;
      else                 f_value = int( eFlav * 256 / E );

      f[flav][i] = f_value;
    }

    nPart[i] == 0 ? eRMS[i] = 0 : eRMS[i] = sqrt( e2[i]/nPart[i] );
  }


//------------ PF Jets ------------//

  edm::Handle< vector<reco::PFJet> > pfJets;
  iEvent.getByToken(pfJetTag_, pfJets);

  ht = 0;
  vector<reco::PFJet>::const_iterator i_jet, endjet = pfJets->end();
  for (i_jet = pfJets->begin(); i_jet != endjet; ++i_jet) {
    float pt = i_jet->pt();
    if (pt > 10) ht += pt; // jet pt cut of 10 GeV
  }

  // Correcteing jets for pileup and sorting wrt pT
//  vector<JetCorrectorParameters> jetPars;
//  FactorizedJetCorrector* jetCorrectors;

//  JetCorrectorParameters* L1JetPars  = new JetCorrectorParameters(era_ + "/" + era_ + "_L1FastJet_"    +  jet_type_ + ".txt");
//  JetCorrectorParameters* L2JetPars  = new JetCorrectorParameters(era_ + "/" + era_ + "_L2Relative_"   +  jet_type_ + ".txt");
//  JetCorrectorParameters* L3JetPars  = new JetCorrectorParameters(era_ + "/" + era_ + "_L3Absolute_"   +  jet_type_ + ".txt");
//  JetCorrectorParameters* ResJetPars = new JetCorrectorParameters(era_ + "/" + era_ + "_L2L3Residual_" +  jet_type_ + ".txt");
//  jetPars.push_back( *L1JetPars );
//  if (doL1L2L3Res_){
//    jetPars.push_back( *L2JetPars );
//    jetPars.push_back( *L3JetPars );
//    jetPars.push_back( *ResJetPars );
//  }
//  jetCorrectors = new FactorizedJetCorrector( jetPars );

  vector<pair<int, double> > jet_index_corrpt; // Pair of jet index and corrected jet pT

  nJets = pfJets->size();
  for (int i=0; i != nJets; ++i){ // Looping through ALL the jets
    reco::PFJet jet = pfJets->at(i);


/*  This part breaks the code. Temporarily commented
    bool jetVeto = false;
    if (dojetVetoMap_){
      for(int j=0; j != int((JetVetoMap[0]).size()); ++j){
        if ((jet.eta()>= JetVetoMap[0][j]) && (jet.eta()<= JetVetoMap[1][j]) && (jet.phi()>= JetVetoMap[2][j]) && (jet.phi()<= JetVetoMap[3][j]))
          jetVeto=true;
      }
    }
    if (jetVeto) continue;
*/

    // Applying L1 corrections
//    jetCorrectors->setJetEta( jet.eta() );
//    jetCorrectors->setJetPt( jet.pt() );
//    jetCorrectors->setJetA( jet.jetArea() );
//    jetCorrectors->setRho(rho);

//    jet_index_corrpt.push_back( make_pair(i, jet.pt()*jetCorrectors->getCorrection() ) );
    jet_index_corrpt.push_back( make_pair(i, jet.pt() ) );
  }

  sort(jet_index_corrpt.begin(), jet_index_corrpt.end(), sortJetPt); // Sorting jets based on pT

  pfJets->size()<MAXJETS ? nJets = pfJets->size() : nJets = MAXJETS; // Only storing FOUR jets
  for (int i=0; i != nJets; ++i){ // Looping through the jets
    
    int is = jet_index_corrpt[i].first ; 
    reco::PFJet jet = pfJets->at(is);
 
    jet_eta[i] = jet.eta();
    jet_phi[i] = jet.phi();
    jet_pt[i] = jet.pt();
    jet_mass[i] = jet.mass();
    jet_area[i] = jet.jetArea();

    jet_ch[i] = jet.chargedHadronEnergyFraction();
    jet_nh[i] = jet.neutralHadronEnergyFraction();
    jet_ne[i] = jet.photonEnergyFraction();
    jet_hfh[i] = jet.HFHadronEnergyFraction();
    jet_hfe[i] = jet.HFEMEnergyFraction();
    jet_lep[i] = jet.electronEnergyFraction() + jet.muonEnergyFraction();
  }  

//------------ Fill Tree ------------//

  tree->Fill();
}


// ------------ method called once each job just after ending the event loop  ------------
void OffsetTreeMaker::endJob() {
  
  if (root_file !=0) {
    cout << "About to write histograms" << endl;
    h_nFoundVertex->Write();
    cout << "nFoundVertex written..." << endl;
    h_bestweight->Write();
    h_bestweight1->Write();
    
    root_file->Write();
    delete root_file;
    root_file = 0;
  }
  else {cout << "No files to write" << endl;}
}

int OffsetTreeMaker::getEtaIndex(float eta){

  for (int i=0; i != ETA_BINS; ++i){
    if (etabins[i] <= eta && eta < etabins[i+1]) return i;
  }
  if (eta == etabins[ETA_BINS]) return ETA_BINS-1;
  else return -1;
}


OffsetTreeMaker::Flavor OffsetTreeMaker::getFlavor(reco::PFCandidate::ParticleType id)
{
    if (id == reco::PFCandidate::h)
        return chm;     //initially matched charged hadron
    else if (id == reco::PFCandidate::e)
        return lep;
    else if (id == reco::PFCandidate::mu)
        return lep;
    else if (id == reco::PFCandidate::gamma)
        return ne;
    else if (id == reco::PFCandidate::h0)
        return nh;
    else if (id == reco::PFCandidate::h_HF)
        return hfh;
    else if (id == reco::PFCandidate::egamma_HF)
        return hfe;
    else
        return X;
}

//define this as a plug-in
DEFINE_FWK_MODULE(OffsetTreeMaker);

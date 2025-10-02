// -*- C++ -*-
//
// Package:    test/OffsetTreeMaker_mAOD
// Class:      OffsetTreeMaker_mAOD
//
/**\class OffsetTreeMaker_mAOD OffsetTreeMaker_mAOD.cc test/OffsetTreeMaker_mAOD/plugins/OffsetTreeMaker_mAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Garvita Agarwal
//         Created:  Sat, 04 Jan 2020 05:00:48 GMT
//
// 


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Ref.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"

#include "parsePileUpJSON2.h"
#include "jetVetoMap.h"

#include <vector>
#include "TMath.h"
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
const int MAXNPV = 50;
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

class OffsetTreeMaker_mAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit OffsetTreeMaker_mAOD(const edm::ParameterSet&);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    int getEtaIndex(float eta);
    enum Flavor{
      chm = 0, chu, nh, ne, hfh, hfe, lep, untrk, numFlavors, X //undefined
    };
    Flavor getFlavor(int id);

    int counter;
    TFile* root_file;
    TTree* tree;

    TH1F* h;
    TH2F* h2_GME; // 2d ET_eta_phi histo 
    TH2F* h2_finnereta;
    TH2F* et_etaphi;
    TH2F* et_etaphi_chs;
    TRandom3* rand;

    int nEta;
    float energy[ETA_BINS], eRMS[ETA_BINS], et[ETA_BINS], etMED[ETA_BINS], etMEAN[ETA_BINS], etMEDchs[ETA_BINS], etMEANchs[ETA_BINS];
    float et_gme[ETA_BINS_GME][PHI_BINS_GME];
    //float et_gme_gen[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t ch_et_gme[ETA_BINS_GME][PHI_BINS_GME];
    //UChar_t ch_et_gme_gen[ETA_BINS_GME][PHI_BINS_GME];
    UChar_t f[numFlavors][ETA_BINS];  //energy fraction by flavor


    ULong64_t event;
    int run, lumi, bx;
    float mu;
    float rho, rhoC0, rhoCC;
//-------------------------------------->>>>>>>>>>
    int mua[16];
    float puz[50];

    unsigned int nPVall, nPV;
    float pv_ndof[MAXNPV], pv_z[MAXNPV], pv_rho[MAXNPV];

    float ht;
    int nJets;
    float jet_eta[MAXJETS], jet_phi[MAXJETS], jet_pt[MAXJETS], jet_area[MAXJETS];
    float jet_ch[MAXJETS], jet_nh[MAXJETS], jet_ne[MAXJETS], jet_hfh[MAXJETS], jet_hfe[MAXJETS], jet_lep[MAXJETS];

    vector<int> pf_type;
    vector<float> pf_pt, pf_eta, pf_phi, pf_et;

    TString RootFileName_, jetVetoMapFileName_, mapName2_;
    vector<vector<double>> JetVetoMap;
    string puFileName_;
    int numSkip_;
  bool isMC_, writeCands_, dojetVetoMap_;

    edm::EDGetTokenT< edm::View<reco::Vertex> > pvTag_;
    edm::EDGetTokenT< edm::View<pat::PackedCandidate> > trackTag_;
    edm::EDGetTokenT< edm::View<PileupSummaryInfo> > muTag_;
    edm::EDGetTokenT< edm::View<pat::PackedCandidate> > pfTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<double> rhoC0Tag_;
    edm::EDGetTokenT<double> rhoCCTag_;
    edm::EDGetTokenT< edm::View<pat::Jet> > pfJetTag_;
};

OffsetTreeMaker_mAOD::OffsetTreeMaker_mAOD(const edm::ParameterSet& iConfig)
{
  pvTag_ = consumes< edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag"));
  trackTag_ = consumes< edm::View<pat::PackedCandidate> >( iConfig.getParameter<edm::InputTag>("trackTag"));
  muTag_ = consumes< edm::View<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("muTag") );
  pfTag_ = consumes< edm::View<pat::PackedCandidate> >( iConfig.getParameter<edm::InputTag>("pfTag"));
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag"));
  rhoC0Tag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoC0Tag"));
  rhoCCTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoCCTag"));
  pfJetTag_ = consumes< edm::View<pat::Jet> >( iConfig.getParameter<edm::InputTag>("pfJetTag"));

  numSkip_ = iConfig.getParameter<int> ("numSkip");
  RootFileName_ = iConfig.getParameter<string>("RootFileName");
  puFileName_ = iConfig.getParameter<string>("puFileName");
  isMC_ = iConfig.getParameter<bool>("isMC");
  writeCands_ = iConfig.getParameter<bool>("writeCands");

  jetVetoMapFileName_ = iConfig.getParameter<string>("jetVetoMapFileName");
  mapName2_ = iConfig.getParameter<string>("mapName2");
  dojetVetoMap_ = iConfig.getParameter<bool>("dojetVetoMap");
}
// ------------ method called for each event  ------------
void
OffsetTreeMaker_mAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  counter++;
  if (counter%numSkip_ != 0) return;

//------------ Pileup ------------//
  // cout <<"\n\nEvent # :" << counter << endl;
  if (isMC_){
    edm::Handle< edm::View<PileupSummaryInfo> > pileups;
    iEvent.getByToken(muTag_, pileups);

    mu = pileups->at(1).getTrueNumInteractions();

    int j = 0;
    edm::View<PileupSummaryInfo>::const_iterator pileupinfo;
    for(pileupinfo = pileups->begin(); pileupinfo != pileups->end(); ++pileupinfo){
      mua[j] = pileupinfo->getPU_NumInteractions();
      //cout << "index = "<< j << endl;
      //cout << "BX = "<< pileupinfo->getBunchCrossing()<<endl;
      //cout << "mu = "<<pileupinfo->getTrueNumInteractions() <<endl;
      //cout << "muactual = "<< pileupinfo->getPU_NumInteractions()<<endl;

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

    mu = getAvgPU( run, lumi );

    if (mu==0) return;

  }

//------------ Primary Vertices ------------//

  edm::Handle< edm::View<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  nPVall = primaryVertices->size();
  bool isGoodPV[nPVall];
  
  for (unsigned int i = 0; i < nPVall; ++i){
    isGoodPV[i] = false;
  }

  nPV = 0;
  for (unsigned int i = 0; i < nPVall; ++i){
    reco::Vertex pv = primaryVertices->at(i);

    pv_ndof[i] = pv.ndof();
    pv_z[i] = pv.z();
    pv_rho[i] = pv.position().rho();

    if( !pv.isFake() && pv_ndof[i] >= 4.0 && fabs(pv_z[i]) <= 24.0 && fabs(pv_rho[i]) <= 2.0 ){
      nPV++;
      isGoodPV[i] = true;
    }
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


//------------ PF Particles ------------//

  edm::Handle< edm::View<pat::PackedCandidate> > pfCandidates;
  iEvent.getByToken(pfTag_, pfCandidates);

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
  //float et_gme_gen[ETA_BINS_GME][PHI_BINS_GME] = {};
  memset(et_gme,        0, sizeof(et_gme));
  memset(ch_et_gme,     0, sizeof(ch_et_gme));
  //memset(et_gme_gen,    0, sizeof(et_gme_gen));
  //memset(ch_et_gme_gen, 0, sizeof(ch_et_gme_gen));

  pf_type.clear(); pf_pt.clear(); pf_eta.clear(); pf_phi.clear(); pf_et.clear();
  h2_GME->Reset();
  et_etaphi->Reset();
  et_etaphi_chs->Reset();
  h2_finnereta->Reset();

  edm::View<pat::PackedCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
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
    Flavor flavor = getFlavor( i_pf->pdgId() );

    if (etaIndex == -1 || flavor == X) continue;
    if (flavor == chm) { //check charged hadrons ONLY
      bool attached = false;
      for (unsigned int ipv = 0; ipv < nPVall && !attached; ipv++) {
        if (isGoodPV[ipv] && i_pf->fromPV(ipv) == 3)
          attached = true;  //pv used in fit
      }
      if (!attached) flavor = chu;
    }

    float e = i_pf->energy();

    energy[etaIndex] += e;
    et[etaIndex] += i_pf->et();
    eFlavor[flavor][etaIndex] += e;

    h2_GME->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
    h2_finnereta->Fill(i_pf->eta(),i_pf->phi(),i_pf->et());
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

  /*for (int ieta = 1; ieta != (ETA_BINS_GME+1); ++ieta){
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      et_gme[ieta-1][iphi-1] =  h2_GME->GetBinContent(ieta, iphi);
      ch_et_gme[ieta-1][iphi-1] =  char(min(255, int( h2_GME->GetBinContent(ieta, iphi)/0.1)));
    }
  } 
  for (int ieta = 1; ieta != (ETA_BINS+1); ++ieta){
    vector<double> x;
    double et_sum = 0;
    for (int iphi = 1; iphi != PHI_BINS_GME+1; ++iphi){
      x.push_back(h2_finnereta->GetBinContent(ieta, iphi));
      et_sum += h2_finnereta->GetBinContent(ieta, iphi);
    }
    sort(x.begin(),x.end()); 
    float median = x[5]; // for 11 phi bins
    etMED[ieta-1] = median;
    etMEAN[ieta-1] = float(et_sum)/(11);
    
    } */

//------------ Tracks ------------//

  // Difference to AOD: the lowest pT tracks are not saved, losttracks contains only tracks not associated to a PFcandidate
  
  edm::Handle< edm::View<pat::PackedCandidate> > tracks;
  iEvent.getByToken(trackTag_, tracks);

  edm::View<pat::PackedCandidate>::const_iterator i_trk, endtrk = tracks->end();
  for (i_trk = tracks->begin(); i_trk != endtrk; ++i_trk) {

    //    const pat::PackedCandidate& track = (*tracks)[i_trk];

    if (!i_trk->hasTrackDetails()) continue;
    
    reco::Track const& pseudotrack = i_trk->pseudoTrack();

    bool goodtrack = pseudotrack.quality(reco::TrackBase::qualityByName("highPurity"));
    if (!goodtrack)            continue;

    /*
    bool matched = false;

    edm::View<pat::PackedCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
    for (i_pf = pfCandidates->begin();  i_pf != endpf && !matched; ++i_pf) {

      //      if ( &(*i_trk) == i_pf->trackRef().get() )
      const reco::Track *constitTrack = i_pf->bestTrack();
      // Want to check if constit track matches with??
      if (constitTrack and abs(constitTrack->pt()-i_trk->pt()) < 1e-3 and abs(constitTrack->eta()-i_trk->eta()) < 1e-3 and abs(constitTrack->phi()-i_trk->phi()) < 1e-3) {
	std::cout << "lost track matched" << std::endl;
	matched = true;
      }

      //      if ( &(*i_trk) == i_pf->trackRef().get() )  matched = true;
      //      if ( &(*i_trk) == i_pf->bestTrack() )  matched = true;      
    }

    if (matched) continue ; */

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

  //Are JEC applied?
  
  edm::Handle< edm::View<pat::Jet> > pfJets;
  iEvent.getByToken(pfJetTag_, pfJets);

  ht = 0;
  edm::View<pat::Jet>::const_iterator i_jet, endjet = pfJets->end();
  for (i_jet = pfJets->begin(); i_jet != endjet; ++i_jet) {

    float pt = i_jet->pt();
    if (pt > 10) ht += pt;
  }

  pfJets->size()<MAXJETS ? nJets = pfJets->size() : nJets = MAXJETS;
  for (int i=0; i != nJets; ++i){
    pat::Jet jet = pfJets->at(i);

    jet_eta[i] = jet.eta();
    jet_phi[i] = jet.phi();
    jet_pt[i] = jet.pt();
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


// ------------ method called once each job just before starting event loop  ------------
void
OffsetTreeMaker_mAOD::beginJob()
{
  root_file = new TFile(RootFileName_,"RECREATE");
  tree = new TTree("T","Offset Tree");

  counter = -1;
  h = new TH1F("mu", "mu", 100, 0, 50);
  rand = new TRandom3;
  h2_GME = new TH2F("GME_2D", "GME_2D_yPhi_xEta", ETA_BINS_GME, -5.0, 5.0, PHI_BINS_GME , -1*M_PI , M_PI);
  et_etaphi = new TH2F("et_etaphi", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);
  et_etaphi_chs = new TH2F("et_etaphi_chs", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);
  h2_finnereta = new TH2F("finnereta_ET", "yPhi_xEta", ETA_BINS, etabins, PHI_BINS_GME, phibins);

  if (dojetVetoMap_){
    jetVetoMap( jetVetoMapFileName_,mapName2_ );
    JetVetoMap = getVetoMap();
  }
  
  if (!isMC_){
    parsePileUpJSON2( puFileName_ );

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
  //tree->Branch("et_gme_gen",   et_gme_gen, "et_gme_gen[18][11]/F");
  tree->Branch("ch_et_gme",       ch_et_gme,     "ch_et_gme[18][11]/b");
  //tree->Branch("ch_et_gme_gen",   ch_et_gme_gen, "ch_et_gme_gen[18][11]/b");

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
}

// ------------ method called once each job just after ending the event loop  ------------
void
OffsetTreeMaker_mAOD::endJob()
{
  if (root_file !=0) {
    root_file->Write();
    delete root_file;
    root_file = 0;
  }
}
int OffsetTreeMaker_mAOD::getEtaIndex(float eta){

  for (int i=0; i != ETA_BINS; ++i){
    if (etabins[i] <= eta && eta < etabins[i+1]) return i;
  }
  if (eta == etabins[ETA_BINS]) return ETA_BINS-1;
  else return -1;
}


OffsetTreeMaker_mAOD::Flavor OffsetTreeMaker_mAOD::getFlavor(int id)
{
    if (id == 211)
        return chm;     //initially matched charged hadron
    else if (id == 11)
        return lep;
    else if (id == 13)
        return lep;
    else if (id == 22)
        return ne;
    else if (id == 130)
        return nh;
    else if (id == 1)
        return hfh;
    else if (id == 2)
        return hfe;
    else
        return X;
}

DEFINE_FWK_MODULE(OffsetTreeMaker_mAOD);

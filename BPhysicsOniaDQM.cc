/*
 *  See header file for a description of this class.
 *
 *  $Date: 2010/01/04 11:50:13 $
 *  $Revision: 1.5 $
 *  \author S. Bolognesi, Erik - CERN
 */

#include "DQM/Physics/src/BPhysicsOniaDQM.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;
using namespace edm;
using namespace reco;

BPhysicsOniaDQM::BPhysicsOniaDQM(const ParameterSet& parameters) {
  cout << "constructor\n";  
  // Muon Collection Label
//  theMuonCollectionLabel = parameters.getParameter<InputTag>("MuonCollection");
  cout << "after muon collection label\n";
  /*glbSigCut = NULL;
  glbSigNoCut = NULL;
  glbBkgNoCut = NULL;
  staSigCut = NULL;
  staSigNoCut = NULL;
  staBkgNoCut = NULL;
  trkSigCut = NULL;
  trkSigNoCut = NULL;
  trkBkgNoCut = NULL;

  JPsiGlbYdLumi = NULL;
  JPsiStaYdLumi = NULL;
  JPsiTrkYdLumi = NULL;*/
}

BPhysicsOniaDQM::~BPhysicsOniaDQM() { 
}

void BPhysicsOniaDQM::beginJob() {
  cout << "BeginJob \n";
  // the services
  theDbe = Service<DQMStore>().operator->();

  metname = "oniaAnalyzer";
  LogTrace(metname)<<"[BPhysicsOniaDQM] Parameters initialization";

  if(theDbe!=NULL){
    theDbe->setCurrentFolder("Physics/BPhysics");  // Use folder with name of PAG
    glbSigCut = theDbe->book1D("glbSigCut", "Opposite-sign glb-glb dimuon mass", 750, 0, 15);
    glbSigNoCut = theDbe->book1D("glbSigNoCut", "Opposite-sign glb-glb dimuon mass (no cut)", 750, 0, 15);
    glbBkgNoCut = theDbe->book1D("glbBkgNoCut", "Same-sign glb-glb dimuon mass (no cut)", 750, 0, 15);
    staSigCut = theDbe->book1D("staSigCut", "Opposite-sign sta-sta dimuon mass", 500, 0, 15);
    staSigNoCut = theDbe->book1D("staSigNoCut", "Opposite-sign sta-sta dimuon mass (no cut)", 500, 0, 15);
    staBkgNoCut  = theDbe->book1D("staBkgNoCut", "Same-sign sta-sta dimuon mass (no cut)", 500, 0, 15);
    trkSigCut = theDbe->book1D("trkSigCut", "Opposite-sign trk-trk dimuon mass", 750, 0, 15);
    trkSigNoCut = theDbe->book1D("trkSigNoCut", "Opposite-sign trk-trk dimuon mass (no cut)", 750, 0, 15);
    trkBkgNoCut = theDbe->book1D("trkBkgNoCutt", "Same-sign trk-trk dimuon mass (no cut)", 750, 0, 15);
  }

}

void BPhysicsOniaDQM::analyze(const Event& iEvent, const EventSetup& iSetup) {

  LogTrace(metname)<<"[BPhysicsOniaDQM] Analysis of event # ";
  cout << "analyze\n";
  
  // Take the STA muon container
  Handle<MuonCollection> muons;
  iEvent.getByLabel(theMuonCollectionLabel,muons);

  Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel("offlinePrimaryVertices",privtxs);
  VertexCollection::const_iterator privtx;

  if(privtxs->begin() != privtxs->end()){
    privtx = privtxs->begin();
    RefVtx = privtx->position();
  } else {
    RefVtx.SetXYZ(0.,0.,0.);
  }

  if(muons.isValid()){
    for (MuonCollection::const_iterator recoMu1 = muons->begin(); recoMu1!=muons->end(); ++recoMu1){

      // only loop over the remaining muons if recoMu1 is one of the following
      if(recoMu1->isGlobalMuon() || recoMu1->isTrackerMuon() || recoMu1->isStandAloneMuon()){
        for (MuonCollection::const_iterator recoMu2 = recoMu1+1; recoMu2!=muons->end(); ++recoMu2){

          // fill the relevant histograms if recoMu2 satisfies one of the following
          if (recoMu1->isGlobalMuon() && recoMu2->isGlobalMuon()){
            math::XYZVector vec1 = recoMu1->globalTrack()->momentum();
            math::XYZVector vec2 = recoMu2->globalTrack()->momentum();
            float massJPsi = computeMass(vec1,vec2);

            // if opposite charges, fill glbSig, else fill glbBkg
            if (((*recoMu1).charge()*(*recoMu2).charge())<0) {
              if(glbSigNoCut!=NULL){
                glbSigNoCut->Fill(massJPsi);
                if (selGlobalMuon(*recoMu1) && selGlobalMuon(*recoMu2)) {
                  if (glbSigCut!=NULL) glbSigCut->Fill(massJPsi);
                  if (massJPsi >= 3.0 && massJPsi <= 3.2) jpsiGlbSigPerLS++;
                }
              }
            } else {
              if(glbBkgNoCut!=NULL){
                glbBkgNoCut->Fill(massJPsi);
              }
            }
          }
          
          if(recoMu1->isStandAloneMuon() && recoMu2->isStandAloneMuon() &&
            fabs(recoMu1->outerTrack()->d0()) < 5 && fabs(recoMu1->outerTrack()->dz()) < 30 &&
            fabs(recoMu2->outerTrack()->d0()) < 5 && fabs(recoMu2->outerTrack()->dz()) < 30){
            math::XYZVector vec1 = recoMu1->outerTrack()->momentum();
            math::XYZVector vec2 = recoMu2->outerTrack()->momentum();
            float massJPsi = computeMass(vec1,vec2);

            // if opposite charges, fill staSig, else fill staBkg
            if (((*recoMu1).charge()*(*recoMu2).charge())<0) {
              if(staSigNoCut!=NULL){
                staSigNoCut->Fill(massJPsi);
                /*if (selStandaloneMuon(*recoMu1) && selStandaloneMuon(*recoMu2)) {
                  if (staSigCut!=NULL) staSigCut->Fill(massJPsi);
                  if (massJPsi >= 3.0 && massJPsi <= 3.2) jpsiStaSigPerLS++;
                }*/
              }
            } else {
              if(staBkgNoCut!=NULL){
                staBkgNoCut->Fill(massJPsi);
              }
            }
          }

          if(recoMu1->isTrackerMuon() && recoMu2->isTrackerMuon() &&
            muon::isGoodMuon(*recoMu1, muon::TrackerMuonArbitrated) &&
            muon::isGoodMuon(*recoMu2, muon::TrackerMuonArbitrated)){
            math::XYZVector vec1 = recoMu1->innerTrack()->momentum();
            math::XYZVector vec2 = recoMu2->innerTrack()->momentum();
            float massJPsi = computeMass(vec1,vec2);

            // if opposite charges, fill trkSig, else fill trkBkg
            if (((*recoMu1).charge()*(*recoMu2).charge())<0) {
              if(trkSigNoCut!=NULL){
                trkSigNoCut->Fill(massJPsi);
                if (selTrackerMuon(*recoMu1) && selTrackerMuon(*recoMu2)) {
                  if (trkSigCut!=NULL) trkSigCut->Fill(massJPsi);
                  if(massJPsi >= 3.0 && massJPsi <= 3.2) jpsiTrkSigPerLS++;
                }
              }
            } else {
              if(trkBkgNoCut!=NULL){
                trkBkgNoCut->Fill(massJPsi);
              }
            }
          }

        }//end of 2nd MuonCollection
      }//end of GLB,STA,TRK muon check
    }//end of 1st MuonCollection
  }//Is this MuonCollection vaild?

}

void BPhysicsOniaDQM::endJob(void) {
  LogTrace(metname)<<"[BPhysicsOniaDQM] EndJob";
}

void BPhysicsOniaDQM::beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &iSetup)
{
  LogTrace(metname)<<"[BPhysicsOniaDQM] Start of a LuminosityBlock";
  cout << "Begin luminosity block\n";
  jpsiGlbSigPerLS = 0;
  jpsiStaSigPerLS = 0;
  jpsiTrkSigPerLS = 0;
}

void BPhysicsOniaDQM::endLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &iSetup)
{
  LogTrace(metname)<<"[BPhysicsOniaDQM] Start of a LuminosityBlock";
  cout << "End luminosity block\n";

  edm::Handle<LumiSummary> lumiSummary;
  lumiBlock.getByLabel("lumiProducer",lumiSummary);

  int LBlockNum = lumiBlock.id().luminosityBlock();
  
  jpsiGlbSig.insert( pair<int,int>(LBlockNum, jpsiGlbSigPerLS) );
  jpsiStaSig.insert( pair<int,int>(LBlockNum, jpsiStaSigPerLS) );
  jpsiTrkSig.insert( pair<int,int>(LBlockNum, jpsiTrkSigPerLS) );
  cout << "lumi: " << LBlockNum << "\t" << jpsiGlbSig[LBlockNum] << "\t" << jpsiStaSig[LBlockNum] << "\t" << jpsiTrkSig[LBlockNum] << endl;
  
  theDbe->setCurrentFolder("Physics/BPhysics");
  if(JPsiGlbYdLumi!=NULL) {
    theDbe->removeElement("JPsiGlbYdLumi");   // Remove histograms from previous run
    theDbe->removeElement("JPsiStaYdLumi");
    theDbe->removeElement("JPsiTrkYdLumi");
  }

  int xmin = (*jpsiGlbSig.begin()).first;
  int xmax = (*jpsiGlbSig.rbegin()).first;
  int nx   = xmax - xmin + 1;
  cout << "x-axis " << xmin << " " << xmax << endl;

  JPsiGlbYdLumi = theDbe->book1D("JPsiGlbYdLumi", "JPsi yield from global-global dimuon", nx, xmin, xmax);
  JPsiStaYdLumi = theDbe->book1D("JPsiStaYdLumi", "JPsi yield from standalone-standalone dimuon", nx, xmin, xmax);
  JPsiTrkYdLumi = theDbe->book1D("JPsiTrkYdLumi", "JPsi yield from tracker-tracker dimuon", nx, xmin, xmax);

  map<int,int>::iterator glb;
  map<int,int>::iterator sta;
  map<int,int>::iterator trk;
  for (glb = jpsiGlbSig.begin(); glb != jpsiGlbSig.end(); ++glb)
  {
    int bin = (*glb).first - xmin + 1;  //X-axis bin #
    sta = jpsiStaSig.find((*glb).first);
    trk = jpsiTrkSig.find((*glb).first);
    JPsiGlbYdLumi->setBinContent(bin,(*glb).second);
    JPsiStaYdLumi->setBinContent(bin,(*sta).second);
    JPsiTrkYdLumi->setBinContent(bin,(*trk).second);
    cout << "glb: " << bin << "\t" << (*glb).first << "\t" << (*glb).second << endl;
    cout << "sta: " << bin << "\t" << (*sta).first << "\t" << (*sta).second << endl;
    cout << "trk: " << bin << "\t" << (*trk).first << "\t" << (*trk).second << endl;
  }
}

void BPhysicsOniaDQM::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  LogTrace(metname)<<"[BPhysicsOniaDQM] Start of a Run";
  cout << "Start of a Run\n";
}

void BPhysicsOniaDQM::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  LogTrace(metname)<<"[BPhysicsOniaDQM] End of a Run";
  cout << "End of a Run\n";
  if (!jpsiGlbSig.empty()) {
    jpsiGlbSig.clear();
    jpsiStaSig.clear();
    jpsiTrkSig.clear();
  }
}

float BPhysicsOniaDQM::computeMass(const math::XYZVector &vec1,const math::XYZVector &vec2){
  // mass of muon
  float massMu = 0.10566;
  float eMu1 = -999;
  if(massMu*massMu + vec1.Mag2()>0)
    eMu1 = sqrt(massMu*massMu + vec1.Mag2());
  float eMu2 = -999;
  if(massMu*massMu + vec2.Mag2()>0)
    eMu2 = sqrt(massMu*massMu + vec2.Mag2());

  float pJPsi = -999;
  if((vec1+vec2).Mag2()>0)
    pJPsi = sqrt((vec1+vec2).Mag2());
  float eJPsi = eMu1 + eMu2;

  float massJPsi = -999;
  if((eJPsi*eJPsi - pJPsi*pJPsi) > 0)
    massJPsi = sqrt(eJPsi*eJPsi - pJPsi*pJPsi);
 
 return massJPsi;
}

bool BPhysicsOniaDQM::isMuonInAccept(const reco::Muon &recoMu)
{
  return (fabs(recoMu.eta()) < 2.4 &&
         ((fabs(recoMu.eta()) < 1.3 && recoMu.pt() > 3.3) ||
          (fabs(recoMu.eta()) > 1.3 && fabs(recoMu.eta()) < 2.2 && recoMu.p() > 2.9) ||
          (fabs(recoMu.eta()) > 2.2 && recoMu.pt() > 0.8)));
}

bool BPhysicsOniaDQM::selGlobalMuon(const reco::Muon &recoMu)
{
  TrackRef iTrack = recoMu.innerTrack();
  const reco::HitPattern &p = iTrack->hitPattern();
  
  TrackRef gTrack = recoMu.globalTrack();
  const reco::HitPattern &q = gTrack->hitPattern();

  return (isMuonInAccept(recoMu) &&
          iTrack->found() > 11 &&
          gTrack->chi2()/gTrack->ndof() < 20.0 &&
          q.numberOfValidMuonHits() > 0 &&
          iTrack->chi2()/iTrack->ndof() < 4.0 &&
          //recoMu.muonID("TrackerMuonArbitrated") &&
          //recoMu.muonID("TMLastStationAngTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
          fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool BPhysicsOniaDQM::selTrackerMuon(const reco::Muon &recoMu)
{
  TrackRef iTrack = recoMu.innerTrack();
  const reco::HitPattern &p = iTrack->hitPattern();

  return (isMuonInAccept(recoMu) &&
          iTrack->found() > 11 &&
          iTrack->chi2()/iTrack->ndof() < 4.0 &&
          //recoMu.muonID("TrackerMuonArbitrated") &&
          //recoMu.muonID("TMLastStationAngTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
          fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}
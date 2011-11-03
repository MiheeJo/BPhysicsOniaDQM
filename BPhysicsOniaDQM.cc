/*
 *  See header file for a description of this class.
 *
 *  $Date: 2010/11/11 17:33:03 $
 *  $Revision: 1.6 $
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
  // Muon Collection Label
  theMuonCollectionLabel = parameters.getParameter<InputTag>("MuonCollection");
  vertex = parameters.getParameter<InputTag>("vertex");

  global_background = NULL;
  diMuonMass_global = NULL;
  tracker_background = NULL;
  diMuonMass_tracker = NULL;
  standalone_background = NULL;
  diMuonMass_standalone = NULL;
  
  glbSigCut = NULL;
  glbSigNoCut = NULL;
  glbBkgNoCut = NULL;
  staSigNoCut = NULL;
  staBkgNoCut = NULL;
  trkSigCut = NULL;
  trkSigNoCut = NULL;
  trkBkgNoCut = NULL;
}

BPhysicsOniaDQM::~BPhysicsOniaDQM() { 
}

void BPhysicsOniaDQM::beginJob() {
  // the services
  theDbe = Service<DQMStore>().operator->();

  metname = "oniaAnalyzer";
  LogTrace(metname)<<"[BPhysicsOniaDQM] Parameters initialization";

  if(theDbe!=NULL){
    theDbe->setCurrentFolder("Physics/BPhysics");  // Use folder with name of PAG
    global_background = theDbe->book1D("global_background", "Same-sign global-global dimuon mass", 750, 0, 15);
    diMuonMass_global = theDbe->book1D("diMuonMass_global", "Opposite-sign global-global dimuon mass", 750, 0, 15);
    tracker_background = theDbe->book1D("tracker_background", "Same-sign tracker-tracker (arbitrated) dimuon mass", 750, 0, 15);
    diMuonMass_tracker = theDbe->book1D("diMuonMass_tracker", "Opposite-sign tracker-tracker (arbitrated) dimuon mass", 750, 0, 15);
    standalone_background = theDbe->book1D("standalone_background", "Same-sign standalone-standalone dimuon mass", 500, 0, 15);
    diMuonMass_standalone = theDbe->book1D("diMuonMass_standalone", "Opposite-sign standalone-standalone dimuon mass", 500, 0, 15);

    glbSigCut = theDbe->book1D("glbSigCut", "Opposite-sign glb-glb dimuon mass", 650, 0, 130);
    glbSigNoCut = theDbe->book1D("glbSigNoCut", "Opposite-sign glb-glb dimuon mass (no cut)", 650, 0, 130);
    glbBkgNoCut = theDbe->book1D("glbBkgNoCut", "Same-sign glb-glb dimuon mass (no cut)", 650, 0, 130);
    staSigNoCut = theDbe->book1D("staSigNoCut", "Opposite-sign sta-sta dimuon mass (no cut)", 430, 0, 129);
    staBkgNoCut  = theDbe->book1D("staBkgNoCut", "Same-sign sta-sta dimuon mass (no cut)", 430, 0, 129);
    trkSigCut = theDbe->book1D("trkSigCut", "Opposite-sign trk-trk dimuon mass", 650, 0, 130);
    trkSigNoCut = theDbe->book1D("trkSigNoCut", "Opposite-sign trk-trk dimuon mass (no cut)", 650, 0, 130);
    trkBkgNoCut = theDbe->book1D("trkBkgNoCutt", "Same-sign trk-trk dimuon mass (no cut)", 650, 0, 130);
  }

}

void BPhysicsOniaDQM::analyze(const Event& iEvent, const EventSetup& iSetup) {

  LogTrace(metname)<<"[BPhysicsOniaDQM] Analysis of event # ";
  
  // Take the STA muon container
  Handle<MuonCollection> muons;
  iEvent.getByLabel(theMuonCollectionLabel,muons);

  Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(vertex,privtxs);
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
              if(diMuonMass_global!=NULL){  // BPhysicsOniaDQM original one
                diMuonMass_global->Fill(massJPsi);
              }

              if(glbSigNoCut!=NULL){
                glbSigNoCut->Fill(massJPsi);
                if (selGlobalMuon(*recoMu1) && selGlobalMuon(*recoMu2)) {
                  if (glbSigCut!=NULL) glbSigCut->Fill(massJPsi);
                }
              }
            } else {
              if(global_background!=NULL){  // BPhysicsOniaDQM original one
                global_background->Fill(massJPsi);
              }

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
              if(diMuonMass_standalone!=NULL){
                diMuonMass_standalone->Fill(massJPsi);
              }

              if(staSigNoCut!=NULL){
                staSigNoCut->Fill(massJPsi);
              }
            } else {
              if(standalone_background!=NULL){
                standalone_background->Fill(massJPsi);
              }

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
              if(diMuonMass_tracker!=NULL){
                diMuonMass_tracker->Fill(massJPsi);
              }

              if(trkSigNoCut!=NULL){
                trkSigNoCut->Fill(massJPsi);
                if (selTrackerMuon(*recoMu1) && selTrackerMuon(*recoMu2)) {
                  if (trkSigCut!=NULL) trkSigCut->Fill(massJPsi);
                }
              }
            } else {
              if(tracker_background!=NULL){
                tracker_background->Fill (massJPsi);
              }

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

void BPhysicsOniaDQM::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  LogTrace(metname)<<"[BPhysicsOniaDQM] Start of a Run";
}

void BPhysicsOniaDQM::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  LogTrace(metname)<<"[BPhysicsOniaDQM] End of a Run";
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
        ((fabs(recoMu.eta()) < 1.0 && recoMu.pt() >= 3.4) ||
         (1.0 <= fabs(recoMu.eta()) && fabs(recoMu.eta()) < 1.5 && recoMu.pt() >= 5.8-2.4*fabs(recoMu.eta())) ||
         (1.5 <= fabs(recoMu.eta()) && recoMu.pt() >= 3.3667-7.0/9.0*fabs(recoMu.eta()))));
}


bool BPhysicsOniaDQM::selGlobalMuon(const reco::Muon &recoMu)
{
  TrackRef iTrack = recoMu.innerTrack();
  const reco::HitPattern &p = iTrack->hitPattern();
  
  TrackRef gTrack = recoMu.globalTrack();
  const reco::HitPattern &q = gTrack->hitPattern();

  return (isMuonInAccept(recoMu) &&
          iTrack->found() > 10 &&
          gTrack->chi2()/gTrack->ndof() < 20.0 &&
          iTrack->chi2()/iTrack->ndof() < 4.0 &&
          p.pixelLayersWithMeasurement() > 0 &&
          fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool BPhysicsOniaDQM::selTrackerMuon(const reco::Muon &recoMu)
{
  TrackRef iTrack = recoMu.innerTrack();
  const reco::HitPattern &p = iTrack->hitPattern();

  return (isMuonInAccept(recoMu) &&
          iTrack->found() > 10 &&
          iTrack->chi2()/iTrack->ndof() < 4.0 &&
          p.pixelLayersWithMeasurement() > 0 &&
          fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}


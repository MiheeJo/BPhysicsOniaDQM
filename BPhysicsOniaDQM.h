#ifndef BPhysicsOniaDQM_H
#define BPhysicsOniaDQM_H


/** \class BPhysicsOniaDQM
 *
 *  DQM offline for quarkonia
 *
 *  $Date: 2010/11/11 17:33:03 $
 *  $Revision: 1.6 $
 *  \author S. Bolognesi, Eric - CERN
 */

#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include <string>
#include <cmath>
#include <map>

class DQMStore;
class MonitorElement;

class BPhysicsOniaDQM : public edm::EDAnalyzer {
 public:

  /// Constructor
  BPhysicsOniaDQM(const edm::ParameterSet&);

  /// Destructor
  virtual ~BPhysicsOniaDQM();

  /// Inizialize parameters for histo binning
  void beginJob();

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&);
  void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  void endJob(void);

 private:

  float computeMass(const math::XYZVector &vec1,const math::XYZVector &vec2);
  bool isMuonInAccept(const reco::Muon &recoMu);
  bool selGlobalMuon(const reco::Muon &recoMu);
  bool selTrackerMuon(const reco::Muon &recoMu);

  // ----------member data ---------------------------

  DQMStore* theDbe;
  
  edm::InputTag vertex;
  // Switch for verbosity
  std::string metname;

  // Muon Label
  edm::InputTag theMuonCollectionLabel;

  //The histos
  MonitorElement* diMuonMass_global;
  MonitorElement* diMuonMass_tracker;
  MonitorElement* diMuonMass_standalone;
  MonitorElement* global_background;
  MonitorElement* tracker_background;
  MonitorElement* standalone_background;

  MonitorElement* glbSigCut;
  MonitorElement* glbSigNoCut;
  MonitorElement* glbBkgNoCut;
  MonitorElement* staSigNoCut;
  MonitorElement* staBkgNoCut;
  MonitorElement* trkSigCut;
  MonitorElement* trkSigNoCut;
  MonitorElement* trkBkgNoCut;

  math::XYZPoint RefVtx;
};
#endif


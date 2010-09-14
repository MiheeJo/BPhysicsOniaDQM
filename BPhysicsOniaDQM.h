#ifndef BPhysicsOniaDQM_H
#define BPhysicsOniaDQM_H


/** \class BPhysicsOniaDQM
 *
 *  DQM offline for quarkonia
 *
 *  $Date: 2010/07/19 19:54:13 $
 *  $Revision: 1.5 $
 *  \author S. Bolognesi, Eric - CERN
 */

#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include <vector>
#include <algorithm>
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
  void beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &iSetup);
  void endLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &iSetup);
  void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup);

  /// Save the histos
  void endJob(void);

 private:

  float computeMass(const math::XYZVector &vec1,const math::XYZVector &vec2);
  bool isMuonInAccept(const reco::Muon &recoMu);
  bool selGlobalMuon(const reco::Muon &recoMu);
  bool selTrackerMuon(const reco::Muon &recoMu);

  // ----------member data ---------------------------
  
  DQMStore* theDbe;
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
  MonitorElement* yieldJPsi_lumi;

  //Yield of dimuon objects
  int jpsiGlbSigPerLS;
  int jpsiTrkSigPerLS;
  std::map<int,int> jpsiGlbSig;
  std::map<int,int> jpsiTrkSig;

  math::XYZPoint RefVtx;
};
#endif  

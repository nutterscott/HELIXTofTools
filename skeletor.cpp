#include <algorithm>
#include <iostream>

/// ROOTy bits
#include <TFile.h>
#include <TTree.h>
#include <TEnv.h>
#include <TProfile.h>  //SN
#include <TApplication.h>

///////////////////////////////////////////////////////////
#include "HLX_TimeStamp.h"
#include "HLX_CLI.h"
#include "HLX_Logger.h"
#include "HLX_Stopwatch.h"
#include "HLX_RunHeader.h"

#include "HLX_Geometry.h"
#include "HLX_Stg1_Event.h"
#include "HLX_Stg2_Event.h"
#include "HLX_Stg3_Event.h"
#include "HLX_Stg3_Hit.h"

// Analysis stuff
#include "HLX_Getch.h"
#include "HLX_HitViewer.h"
#include "HXA_Stg2_HitBuilder.h"
#include "HXA_HoughLineFinder.h"
#include "HXA_HoughCircleFinder.h"
#include "HXA_TrackFitterLinear.h"
#include "HXA_TrackFitterGenFit.h"
#include "HXA_TrackFitterRK.h"
#include "HXA_CKFTrackFinder.h"
#include "HXA_StaggerTrackFinder.h"

#include "HXA_ResidualSummary.h"
#include "HXA_CutSummary.h"
#include "HXA_CalManager_HodoHit.h"
#include "HXA_CalManager_HSK.h"
#include "HXA_IsochroneManager.h"
#include "HXA_RICHRecon.h"
#include "HXA_Stg2_DCTHit_AlgoTag.h"
#include "HXA_RICHReconMinimizer.h"
#include "HXA_FullTrackRICHMinimization.h"
#include "HXA_SegmentFinder.h"

#include "HXA_ResultHeader.h"
#include "HLX_ExitRootApp.h"

using namespace std;

void splash();

/// This is skeleton code designed to be a starting point for analysis code, starting with a stg1 file.
/// It
int main(int argc, char *argv[])
{
  hlog::noWarts();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// handle any command lineFunc bois
  auto cli = HLX_CLI::get();

  cli->add_flag("--splash", [](size_t t){splash();exit(0);}, "Show the splash screen");

  bool kForceBON{false};
  cli->add_flag("--bon", kForceBON, "Forces B-ON track fitting, even on B-OFF runs (default: auto-select)");

  string kIPFName;
  cli->add_option("ipfname", kIPFName, "Input stg1 file name or run_list.toml name (if --rl is passed)")->required();

  string kOPFName{"skeletor_out.root"};
  cli->add_option("-o, --out", kOPFName, "Output file name");

  long kMaxCount{-1};
  cli->add_option("-N", kMaxCount, "Only process N events");

  long kStartIndex{-1};
  cli->add_option("-S", kStartIndex, "Start at event x");

  bool kClobber{false};
  cli->add_flag("-c", kClobber, "clobber output file");

  bool kTreatIPFNameAsRunList{false};
  cli->add_flag("--rl, --run-list", kTreatIPFNameAsRunList, "Treats input file name as a run list");

  string kIPFDirectory{"data/stg1/stg1.merged.fixed"};
  cli->add_option("--dd, --data-dir", kIPFDirectory, "Search path for stg1 data files (requires --rl)")->needs("--rl");

  int kRunIdx{0};
  cli->add_option("--ri, --run-idx", kRunIdx, "Run index to look at within the provided run list (requires --rl)")->needs("--rl");

  bool kDoIsochrone2Pass{false};
  cli->add_flag("--iso2", kDoIsochrone2Pass, "Do a 2-pass isochrone fit (slope-matching)");

  bool kCutByAlgoTagClamp{false};
  cli->add_flag("--tag", kCutByAlgoTagClamp, "Get rid of the hits that were clamped when reconstructed using isochrone fits");

  bool kShowPlots{false};
  cli->add_flag("--plots", kShowPlots, "Show plots at the end (blocking)");

  bool kInteractive{false};
  cli->add_flag("-I", kInteractive, "Interactive - event by event mode")->excludes("--plots");

  bool kUseRK{false};
  cli->add_flag("--rk", kUseRK, "Use Runge-Kutta for B-ON track fitting (default: GenFit)");

  bool kFiducialCut{false};
  cli->add_flag("--fid", kFiducialCut, "Use fiducial cuts for track fitting (default: no)");

  bool kDoRICHTrackMinimization{false};
  cli->add_flag("--rich-min", kDoRICHTrackMinimization, "Do a combined track /rich minimization in RICH reconstruction (default: no)");

  bool kWriteOutTree{false};
  cli->add_flag("--tree", kWriteOutTree, "Write an output tree (big outputs!)");

  bool kNoCuts{false};
  cli->add_flag("--no-cuts", kNoCuts, "Do not perform any cuts in the event loop");

  bool kUseSingleCol{false};
  cli->add_flag("--scol", kUseSingleCol, "Use Tracks with only single columns or don't cross into another column");

  bool kUseSingleSide{false};
  cli->add_flag("--sside", kUseSingleSide, "Use Single Side tracks only (no crossing sense planes or cathodes)");

  bool kUseTwoCol{false};
  cli->add_flag("--twocol", kUseTwoCol, "Use Tracks with only 2 columns")->excludes("--scol");

  //=============================================================================================================
  // track finder/hit selection options
  auto selGroup = cli->add_option_group("Tracking finding Algorithms");
  bool kUseSegmentFinder{false};
  selGroup->add_flag("--seg", kUseSegmentFinder, "Use SegmentFinder for track finding/hit selection");

  bool kUseCKF{false};
  selGroup->add_flag("--ckf", kUseCKF, "Use CKF for track finding/hit selection");

  bool kUseHough{false};
  selGroup->add_flag("--hough", kUseHough, "Use Hough for track finding/hit selection");
  selGroup->require_option(1);

  bool kUseStagger{false};
  selGroup->add_flag("--stagger", kUseStagger, "Use stagger scoring for track finding/hit selection");

  double kCKFEnvelope{-1};
  cli->add_option("--ckf-env", kCKFEnvelope, "CKF envelope size (sigma)")->needs("--ckf");

  /// allow this to override everything
  bool kCountAndExit{false};
  cli->add_flag("--count",
                [&](size_t n) {if(n>0) {kCountAndExit = true; selGroup->require_option(0);}},
                "Just count the entries in the file and exit")->trigger_on_parse();

  CLI11_PARSE(*cli, argc, argv)
  if(!kCountAndExit) hlog::info(hlog::kOrange, "Input file name: '{}' | Mode: {}", kIPFName, (kTreatIPFNameAsRunList) ? "Run List" : "Single File");

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// spawn the any key watcher
  if(kInteractive) HLX_Getch::spawnWatchThread();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// TApplication not always needed
  TApplication *theApp{};
  if(kInteractive || kShowPlots) {
    hlog::warn("Interactive or plot mode requested - make sure your DISPLAY var is set and X server is running");
    theApp = new TApplication("App", &argc, argv);
    theApp->SetReturnFromRun(true); // don't exit() on end
  }

  if (kTreatIPFNameAsRunList){
    HLX_TOML_Table runListTable(kIPFName);
    vector<int> runlistNumbers;
    runListTable.getValueAtKeyPath("Runlist", runlistNumbers);
    if (kRunIdx >= runlistNumbers.size()) {
      hlog::error("Run index {} too big for subset size {}", kRunIdx, runlistNumbers.size());
      if(theApp) theApp->Terminate();
      return 1;
    }
    kIPFName = format("{}/run_{}.stg1.root", kIPFDirectory, runlistNumbers[kRunIdx]);
    hlog::info(hlog::kOrange, "Run index {} from list has ID {} and path {}", kRunIdx, runlistNumbers[kRunIdx], kIPFName);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Load a stg1 File
  auto *_file0 = TFile::Open(kIPFName.c_str(), "+read");
  if(!_file0 || !_file0->IsOpen() || _file0->IsZombie()){
    hlog::error("Could not open file {}", kIPFName);
    if(theApp) theApp->Terminate();
    return 1;
  }

  /// Stg1 Tree
  auto t = _file0->Get<TTree>("t");
  if(!t){
    hlog::error("Could not get the tree at namecycle 't' !");
    if(theApp) theApp->Terminate();
    return 1;
  } else {
    if(kCountAndExit) {
      printf("%lld\n", t->GetEntries());
      if(theApp) theApp->Terminate();
      return 0;
    }
    hlog::info("Tree loaded with entry count: {}!", t->GetEntries());
  }

  /// get the run summary - needed because we need the pedestal file names
  unique_ptr<HLX_Stg1_RunSummary> runSummary(_file0->Get<HLX_Stg1_RunSummary>("rs"));
  if(!runSummary){
    hlog::error("Could not get the HLX_Stg1_RunSummary from file {}!", kIPFName);
    if(theApp) theApp->Terminate();
    return 1;
  } else {
    hlog::info("Run Summary loaded!");
  }

  /// Input branch setup
  auto *evtS1 = new HLX_Stg1_Event();
  int rval = t->SetBranchAddress("evtS1", &evtS1);
  if(rval != 0){
    hlog::error("Could not set the branch address for evtS1");
    if(theApp) theApp->Terminate();
    return 1;
  }

  //get the first event to check magnetic field
  t->GetEntry(0);
  HXA_CalManager_HSK* calManHSK = HXA_CalManager_HSK::loadFromRoot();
  
  bool kMagFieldIsOn = static_cast<bool>(calManHSK->fMagnetChargeState.getValueAtTime(evtS1->getTrigger()->getTimeStamp()));
//  bool kMagFieldIsOn = static_cast<bool>(calManHSK->fMagnetChargeState.getValueAtTime(evtS1->fTimeStamp));
  hlog::info("Magnetic Field for this run is reported by CalMan to be {}", kMagFieldIsOn ? "ON" : "OFF");
  if (kMagFieldIsOn) {
  } else {
    if(kForceBON) {
      kMagFieldIsOn = true;
      hlog::warn("FORCING magnetic field to be ON, though CalMan reports it OFF");
    }
  }

  /// Get the geometry singleton - it will autoload the files.
  auto theGeom = HLX_Geometry::getInstance();

  // You MUST configure a magnetic field model if you want tracking to work.  This is new behavior as of July 2026.
  // previously, no configuration would default to a zero-field model
  auto fieldManager = HXA_MagneticFieldManager::getInstance();
  if(kMagFieldIsOn){
    fieldManager->setModel(FieldType::kHEAT, true);
  } else {
    fieldManager->setModel(FieldType::kZeroField);
  }
  fieldManager->getModel()->dump(false, true);

  /// Report on hit selection choice
  if(kUseCKF){
    hlog::info(hlog::kBurlyWood, "Using CKF for DCT hit selection");
  } else if (kUseSegmentFinder) {
    hlog::info(hlog::kBurlyWood, "Using Segment Finder for DCT hit selection");
  } else if (kUseStagger) {
    hlog::info(hlog::kBurlyWood, "Using Stagger Score for DCT hit selection");
  } else if (kUseHough) {
    string algoType = kMagFieldIsOn ? "Circular Hough Transform" : "Linear Hough Transform";
    hlog::info(hlog::kBurlyWood, "Using {} for DCT hit selection", algoType);
  }

  if(kUseTwoCol) hlog::info(hlog::kViolet, "WARNING! TWO COLUMN TRKS ONLY!");
  if(kUseSingleSide) hlog::info(hlog::kViolet, "WARNING! SINGLE SIDE TRKS ONLY!");
  if(kUseSingleCol) hlog::info(hlog::kViolet, "WARNING! SINGLE COLUMN TRKS ONLY!");

  //Figure out whether isochrones are needed and load them if so
  HLX_TOML_Table steeringInfo;
  steeringInfo.loadTableFromFile("config/steeringS1S2.toml");
  string algoName = steeringInfo.getValueOrDefaultAtKeyPath("DCT.ProcessingAlgorithm", "", true);
  string hybridAlgoName = steeringInfo.getValueOrDefaultAtKeyPath("DCT.kHybrid.Bending", "", false);

  hlog::info("DCT algo name: {}, DCT hybrid algo name: {}", algoName, hybridAlgoName);
  const bool kUseIsochrones = (algoName == "kIsochrone") || (hybridAlgoName == "kIsochrone") || (hybridAlgoName == "kIsochroneHEATFits") || (hybridAlgoName == "kIsochroneHEATFitsScaled") || (hybridAlgoName == "kIsochroneHEATFitsRefined");

  // Track the cuts
  HXA_CutSummary cutSummary;

  if (kUseIsochrones) {
    auto isoMan = HXA_IsochroneManager::getInstance();
    if(kMagFieldIsOn){
      hlog::info("Loading isochrone fit files for HEAT field with along wire adjustments...");
      //isoMan->loadAllIsochroneFitsFromFile();
      isoMan->loadAllIsochroneFitsFromSingleFile();
    }
    else{
      hlog::info("Loading isochrones for no field...");
      isoMan->loadIsochronesByState(HXA_IsochroneManager::eMagnetState::kBOFF, 299.0); // we can just load the no-field isochrones for now since we are only doing a quick check of the code, but eventually we will want to load all the isochrones for all the states and temps and have the manager pick the right one based on the conditions at the time of the event
    }
    hlog::info(hlog::kOrange, "After loading isochrone data -> Memory usage: {:6.3f} MB", HLX_Util::getCurrentMemoryUsageMB());
  }

  /// Load hodo mapping info
  auto hodoCalMan = HXA_CalManager_HodoHit::getInstance();

  /// A lambda to build DCT hits from a Hough line
  auto buildDCTStage3HitsFromHoughLine = [](const HLX_Stg2_Event *evtS2, const HXA_HoughLine &inLine, bool kAreWeUsingIsochrones, bool invertSign=false) -> vector<HLX_Stg3_Hit> {
    auto lineHits = kAreWeUsingIsochrones ? inLine.getContributingPoints(true, 1) : inLine.getContributingPoints();
    vector<HLX_Stg3_Hit> stage3HitsDCT;
    for(auto pointInfo: lineHits) {
      if (pointInfo.type != HXA_HoughPoint::Type::kDCT) continue;

      int idx = pointInfo.idx;
      int sign = pointInfo.sign;
      if(invertSign) sign *= -1; // invert the sign if requested
      if (idx < 0 || idx >= evtS2->fDCTHits.size()) {
        hlog::error("Index out of range - this should not be possible");
        continue;
      }
      auto hit = evtS2->fDCTHits[idx];
      hit.setDriftSign(sign); // optionally set the drift sign
      stage3HitsDCT.push_back({{hit.getX(), hit.getYWithDriftSign(sign), hit.getZWithDriftSign(sign)}, {hit.getXErr(), hit.getYErr(), 0.0}, HLX_Stg3_Hit::Type::kDCT});
      stage3HitsDCT.back().fHitID = hit.getUniqueID();
      stage3HitsDCT.back().fDriftTime = hit.getDriftTime(); // copy the drift time
      stage3HitsDCT.back().fDriftSign = hit.getDriftSign(); // copy the drift sign
    }
    return stage3HitsDCT;
  };

  /// A lambda to build DCT hits from a Hough circle
  auto buildDCTStage3HitsFromHoughCircle = [=](const HLX_Stg2_Event *evtS2, const HXA_HoughCircle &inCircle, bool kAreWeUsingIsochrones, bool invertSign=false) -> vector<HLX_Stg3_Hit> {
    auto circleHits = kAreWeUsingIsochrones ? inCircle.getContributingPoints(true, 1) : inCircle.getContributingPoints();
    vector<HLX_Stg3_Hit> stage3HitsDCT;
    for(auto pointInfo: circleHits) {
      if (pointInfo.type != HXA_HoughPoint::Type::kDCT) continue;

      int idx = pointInfo.idx;
      int sign = pointInfo.sign;
      if(invertSign) sign *= -1; // invert the sign if requested
      if (idx < 0 || idx >= evtS2->fDCTHits.size()) {
        hlog::error("Index out of range - this should not be possible");
        continue;
      }
      auto hit = evtS2->fDCTHits[idx];
      hit.setDriftSign(sign); // optionally set the drift sign
      stage3HitsDCT.push_back({{hit.getX(), hit.getYWithDriftSign(sign), hit.getZWithDriftSign(sign)}, {hit.getXErr(), hit.getYErr(), 0.0}, HLX_Stg3_Hit::Type::kDCT});
      stage3HitsDCT.back().fHitID = hit.getUniqueID();
      stage3HitsDCT.back().fDriftTime = hit.getDriftTime(); // copy the drift time
      stage3HitsDCT.back().fDriftSign = hit.getDriftSign(); // copy the drift sign

    }
    return stage3HitsDCT;
  };

  /// a lambda that uses segment finder
  auto buildDCTStage3HitsWithSegmentFinder = [](HLX_Stg2_Event *evtS2, bool flipSign=false) {
    /// use segment finder to select only the hits from the chosen segment,
    /// but do not mutate evtS2->fDCTHits so downstream code still sees the full event.
    auto dctSegments  = HXA_SegmentFinder::findSegment(evtS2->fDCTHits);
    auto hitsInSegments = HXA_SegmentFinder::findDCTHitsFromSegment(dctSegments, 0);
    if (flipSign) {
      for (auto &hit : hitsInSegments) {
        hit.fDriftSign *= -1;
      }
    }

    // Ken always seems to replace the stg2 hits here, but that seems wrong
    //evtS2->fDCTHits.swap(hitsInSegments);

    // Use the shared helper to convert stage2 DCT hits to stage3 hits.
    auto stage3DCTHits = HXA_SegmentFinder::convertDCTHitsToStage3Hits(hitsInSegments);

    return stage3DCTHits;
  };

  /// Set the ROOT error ignore level to suppress warnings about missing branches, etc.
  if(cli->count("-V") < 1) gErrorIgnoreLevel = kSysError;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///=================================================================================================================
  /// Open the output file before any histograms/trees are created below and they will get magically associated
  /// with it in that uniquely awful root way
  TFile opf(kOPFName.c_str(), kClobber ? "RECREATE" : "CREATE");
  if(!opf.IsOpen()){
    hlog::error("Could not open output file '{}' for writing (pass -c to clobber an existing file)", kOPFName);
    if(theApp) theApp->Terminate();
    return 1;
  }

  HXA_Track *outTrack{};
  TTree *outTree{};
  if(kWriteOutTree){
    outTree = new TTree("t", "Output Tree");
    outTree->SetDirectory(&opf);
    outTree->Branch("trk", &outTrack, 32000, 0);
  }

  /// TOF info -- SLN
//  vector<HLX_Stg2_TOFHit>* TOFHits = new vector<HLX_Stg2_TOFHit>(); //top, bore, bottom
  HLX_Stg2_TOFHit* topTOF = new HLX_Stg2_TOFHit();
  HLX_Stg2_TOFHit* boreTOF = new HLX_Stg2_TOFHit();
  HLX_Stg2_TOFHit* botTOF = new HLX_Stg2_TOFHit(); //top, bore, bottom
//Hits:
  Int_t Zhi, Zlo; //trig bits
  struct TOFStruct {
    Double_t ttofX, ttofY, boreX, boreY, btofX, btofY,richX, richY;
    Double_t pathLength, topPLC, borePLC, botPLC;
  };
  double evttime;
  int runnum;
  Double_t richZEst, DCTChi2, hitCount;
  Double_t rig, mdr;  //rigidity, mdr
  TOFStruct tvars;
  double topdircos[3], boredircos[3], botdircos[3];

  auto *toft = new TTree("toft", "TOF tree");
//  toft->Branch("TOFHits", &TOFHits);
  toft->Branch("topTOF", &topTOF);
  toft->Branch("boreTOF", &boreTOF);
  toft->Branch("botTOF", &botTOF);
  toft->Branch("Zhi",&Zhi,"Zhi/I");
  toft->Branch("Zlo",&Zlo,"Zlo/I");
  toft->Branch("tvars",&tvars,"ttofX/D:ttofY/D:boreX/D:boreY/D:btofX/D:btofY/D:richX/D:richY/D:pathLength/D:topPLC/D:borePLC/D:botPLC/D");
//  toft->Branch("tvars",&tvars);
  toft->Branch("runnum",&runnum,"runnum/I");
  toft->Branch("evttime",&evttime,"evttime/D");
  toft->Branch("richZEst",&richZEst,"richZEst/D");
  toft->Branch("DCTChi2",&DCTChi2,"DCTChi2/D");
  toft->Branch("rig",&rig,"rig/D");
  toft->Branch("mdr",&mdr,"mdr/D");
  toft->Branch("hitCount",&hitCount,"hitCount/D");
  toft->Branch("topdircos", topdircos,   "topdircos[3]/D");
  toft->Branch("boredircos",boredircos, "boredircos[3]/D");
  toft->Branch("botdircos", botdircos,   "botdircos[3]/D");



  /// Some Histos
  auto ttofMap = new TH2D("h_TTOFMap", "ttofMap",     200, -800, 800, 200, -800, 800);
  auto boreMap = new TH2D("h_BoreMap", "boreMap",     200, -800, 800, 200, -800, 800);
  auto btofMap = new TH2D("h_BTOFMap", "btofMap",     200, -800, 800, 200, -800, 800);
  auto hHodoC  = new TH2D("h_HodoC"  , "Hodo Coarse", 200, -500, 500, 640, -320, 320 );
  auto hHodoF  = new TH2D("h_HodoF"  , "Hodo Fine",   200, -500, 500, 640, -320, 320 );
  auto hTOFdT  = new TH1D("h_TOFdT"  , "TOF dT"  ,    200, -50, 50);
  auto h_Chi2Pre   = new TH1D("h_Chi2Pre"   , "Fit chi2",    100, 0, 200);
  hHodoC->SetTitle("Coarse Precision;Fiber Hit Position (mm);Residual (trk-hit) (mm);");
  hHodoF->SetTitle("Fine Precision;Fiber Hit Position (mm);Residual (trk-hit) (mm);"  );

  auto hTTOF = new TH1D("h_TTOF", "TTOF Z;Z Est;Counts", 200, 0, 12);
  auto hBTOF = new TH1D("h_BTOF", "BTOF Z;Z Est;Counts", 200, 0, 12);
  auto hHitCount = new TH1D("h_HitCount", "Track Hit Count;Hits/Track;Counts", 100, 0, 100);

  gStyle->SetOptFit(1101);
  auto h_ResVPosZ   = new TH2D("h_ResVPosZ",  "Y Residuals vs Z; Hit Z Coord [mm];Y Resid [mm]", 75,  -300, 300, 200, -5, 5);
  auto h_ResVDrift  = new TH2D("h_ResVDrift", "Y Residuals vs Drift Time; Drift Time [us];Y Resid [mm]", 240,  -15,  15, 200, -5, 5);
  auto h_ResVDrift_filtered = new TH2D("h_ResVDrift_filtered", "Y Residuals vs Large Drift Time; Drift Time [us];Y Resid [mm]", 240,  -15,  15, 200, -5, 5);
  auto h_ResAll     = new TH1D("h_ResAll",     "", 200,  -5,  5);
  auto h_ResPos     = new TH1D("h_ResPos",     "Y Residuals Pos Stag; Y Resid [mm];Counts", 200,  -5,  5);
  auto h_ResNeg     = new TH1D("h_ResNeg",     "Y Residuals Neg Stag; Y Resid [mm];Counts", 200,  -5,  5);
  auto h_FOM        = new TH1D("h_FOM",        "FOM; FOM     ;Counts", 505,  -5,  500);
  auto h_Rig        = new TH1D("h_Rig",        "Rig; Rig [GV];Counts", 1005,  -505,  500);
  auto h_Eta        = new TH1D("h_Eta",        "Eta; Eta [1/GV];Counts", 200,   -0.5,  1.5);
  auto h_RigErr     = new TH1D("h_RigErr",     "RigErr; RigErr [?];Counts", 1005,  -505,  500);
  auto h_MDR        = new TH1D("h_MDR",        "MDR; MDR(GV);Counts", 200,  0,  1000);
  auto h_Chi2Post   = new TH1D("h_Chi2Post",   "Reduced #chi^2; #chi^2/ndf;Counts", 200,  0,  10);
  auto h_YErr       = new TH1D("h_YErr",       "Y Errors; Y Err Val(mm);Counts", 200,  0,  5.00);
  auto h_XErr       = new TH1D("h_XErr",       "X Errors; X Err Val(mm);Counts", 200,  0,  100.00);
  auto h_RICHBeta   = new TH1D("h_RICHBeta",   "RICH Beta; RICH Beta;Counts", 200,  0.5,  1.5);
  auto h_RICHRes    = new TH1D("h_RICHRes",    "RICH Residuals; RICH Resid (mm);Counts", 200,  -0.01,  0.01);
  auto h_RICHCharge = new TH1D("h_RICHCharge", "RICH Charge; RICH Charge (e);Counts", 200,  0.0,  1500.0);

  /// For tracking resolution
  HXA_ResidualSummary resSum("Residual Summary");
  HXA_ResidualSummary pullSum("Pull Summary", true);

  /// RICH Recon
  HXA_RICHReconMinimizer richRecon;
  HXA_FullTrackRICHMinimization fullTrackRichMin;

  richRecon.setRICHGoodBadPixels("../Analysis/Selector/selector_user/reference_data/rich_pixel_trims_20240105_50.0ns_06.7C.toml");

  if (kDoRICHTrackMinimization) {
    fullTrackRichMin.setXStepSize(17.3); //mm
    if (kMagFieldIsOn) fullTrackRichMin.setYStepSize(2.144); //mm
    else fullTrackRichMin.setYStepSize(0.924); //mm
    fullTrackRichMin.setAzimuthalStepSize(0.0568); //rad
    fullTrackRichMin.setIncidenceStepSize(0.0234); //rad

    richRecon.setTrackXError(21.6); //mm
    if (kMagFieldIsOn) richRecon.setTrackYError(0.984); //mm
    else richRecon.setTrackYError(0.564); //mm
    richRecon.setTrackAzimuthError(0.0524); //rad
    richRecon.setTrackIncidenceError(0.0152); //rad
  }
  else{
    richRecon.setTrackXError(17.3); //mm
    if (kMagFieldIsOn) richRecon.setTrackYError(2.144); //mm
    else richRecon.setTrackYError(0.924); //mm
    richRecon.setTrackAzimuthError(0.0568); //rad
    richRecon.setTrackIncidenceError(0.0234); //rad
  }
  
  /// Make and configure the hit builder
  HXA_Stg2_HitBuilder hitBuilder2;
  hitBuilder2.configure("config/steeringS1S2.toml", runSummary.get());

  /// calculate loop range
  long entryCount = t->GetEntries();
  long startEntry = std::max(kStartIndex, 0L);
  long endEntry = (kMaxCount < 0) ? entryCount : std::min(entryCount, startEntry + kMaxCount);
  long totalEntryCount = endEntry - startEntry;
  if (totalEntryCount <= 0) {
    hlog::error("No entries to process!");
    if(theApp) theApp->Terminate();
    return 1;
  }
  hlog::info("Processing {} entries [{} to {}]", totalEntryCount, startEntry, endEntry-1);

  /// Start Loop
  HLX_Stopwatch rateSW, totalSW;
  const int kReportInterval = std::max(totalEntryCount/10,1L);
  long selectedEventCount{};

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇👇
  /// Here is where the action starts - main event loop
  for(int eventIdx=startEntry; eventIdx < endEntry; eventIdx++){
    t->GetEntry(eventIdx);
    if(eventIdx%((kReportInterval))==0){
      hlog::info("Event: {} | Parse Rate: {:6.2f} evts/sec", eventIdx, rateSW.getRate(kReportInterval));
      rateSW.restart();
    }

    /// Do something.
    if(!evtS1) {
      hlog::error("No event found at index {}", eventIdx);
      continue;
    }

    ///=================================================================================================================
    /// Do the actual stage1 => stage2 conversion
    /// the stage2 hit-builder will call all DCT/TOF/RICH hit builders for you
    HLX_Stg2_Event evtS2_orig;
    int rval2 = hitBuilder2.processStage1(eventIdx, evtS1, &evtS2_orig);
    if(rval2 != 0){
      hlog::error("Bad stage1 => stage2 conversion!!!");
      continue;
    }
    cutSummary.registerCutPass("entry");

    ///=================================================================================================================
    /// Optionally kill the canonical bad wires - add your own if you want!
    vector<int> badWireList = {11, 143, 155, 168, 190};
    auto killBadWires = [](HLX_Stg2_Event *evtS2, const vector<int> &inBadWireList) {
      auto hitVec = &evtS2->fDCTHits;
      hitVec->erase(std::remove_if(hitVec->begin(), hitVec->end(), [&](const HLX_Stg2_DCTHit& dctHit) {
        return std::find(inBadWireList.begin(), inBadWireList.end(), dctHit.fWireID) != inBadWireList.end();
      }), hitVec->end());
    };
    killBadWires(&evtS2_orig, badWireList);

    ///=================================================================================================================
    /// Report some event info
    double eventTimeDays = HLX_TimeStamp::dT(evtS2_orig.fTimeStamp, true);
    auto dctHitCount = evtS2_orig.fDCTHits.size();

    for(auto const &hit: evtS2_orig.fDCTHits) {
      h_XErr->Fill(hit.fXErr);
      h_YErr->Fill(hit.fYErr);
    }

//    hlog::sep(hlog::kBurlyWood);
//    hlog::info(hlog::kBurlyWood, "Evt idx: {} | Run/Id: {}/{} | DCT Hits: {:3d} | Evt time: {}" , eventIdx, evtS2_orig.fRunNumber, evtS2_orig.fEventNumber, dctHitCount, eventTimeDays);
    auto evtS2 = evtS2_orig; // changed such that things can be cut before any track finding or track recon occurs

    if(kCutByAlgoTagClamp){
      hlog::info(hlog::kViolet, "WARNING! Removing stage 2 hits that were reconstructed from out of bounds of the isochroneHEATFits!");
       auto it = std::remove_if(evtS2.fDCTHits.begin(), evtS2.fDCTHits.end(),
                       [](const HLX_Stg2_DCTHit& h) { return (HXA_Stg2_DCTHit_AlgoTag(h.fAlgoTag).getOptions() & 0xFF00) != 0; });
                       evtS2.fDCTHits.erase(it, evtS2.fDCTHits.end());
       dctHitCount = evtS2.fDCTHits.size();
    }


    /// do some simple filtering to find clean events
    if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "dctHitCount < 36", dctHitCount < 36,[](){});

    if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "dctHitCount > 300", dctHitCount > 300, [](){});

    ///================================================================================================================
    /// Select events with clean TOF hits
    vector<HLX_Stg2_TOFHit> ttofHits, boreHits, btofHits;
    vector<HLX_Stg3_Hit> stage3TOFHits; // save these too, for residual analysis later
    for(auto tofHit2: evtS2.fTOFHits){
      if (!tofHit2.isHiZHit() && !tofHit2.isLoZHit()) continue;
      auto layerName = theGeom->fTOF.getLayerNameIDByZ(tofHit2.getZ());
      if(layerName == HLX_Geometry_TOF::LayerName::kTop){
        ttofHits.push_back(tofHit2);
        stage3TOFHits.push_back({{tofHit2.getX(), tofHit2.getY(), tofHit2.getZ()}, {tofHit2.getXErr(), tofHit2.getYErr(), tofHit2.getZErr()}, HLX_Stg3_Hit::Type::kTTOF});
      }
      if(layerName == HLX_Geometry_TOF::LayerName::kBore){
        boreHits.push_back(tofHit2);
        stage3TOFHits.push_back({{tofHit2.getX(), tofHit2.getY(), tofHit2.getZ()}, {tofHit2.getXErr(), tofHit2.getYErr(), tofHit2.getZErr()}, HLX_Stg3_Hit::Type::kBore});
      }
      if(layerName == HLX_Geometry_TOF::LayerName::kBottom){
        btofHits.push_back(tofHit2);
        stage3TOFHits.push_back({{tofHit2.getX(), tofHit2.getY(), tofHit2.getZ()}, {tofHit2.getXErr(), tofHit2.getYErr(), tofHit2.getZErr()}, HLX_Stg3_Hit::Type::kBTOF});
      }
    }

    /// single hits per layer only
    CutSummary_APPLY_CUT(cutSummary, "Single TTOF hits", ttofHits.size() != 1 || boreHits.size() != 1 || btofHits.size() != 1, [](){});

    double dTOF = btofHits[0].fTime - ttofHits[0].fTime;
    hTOFdT->Fill(dTOF);

    ///=================================================================================================================
    /// Find the SiPM ID of hodoscope hit
    double   hodoMax   = 0;
    uint32_t hodoID    = 0xFFFFFFFF;
    for (auto &hit: evtS2.fHodoHits){ 
      auto signal = hit.fSignal;
      auto sipmID = hit.fChanID; // begin with 0xA

      if( hodoMax < signal ) {
        hodoMax = signal;
        hodoID  = sipmID;
      }
    }
    if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "Bad HodoID", hodoID == 0xFFFFFFFF, [](){});

    ///================================================================================================================
    /// build some stg3 hits to use for tracking
    vector<HLX_Stg3_Hit> stage3DCTHits;

    /// Pick your favorite hit selection method
    if(kUseCKF) {
      HXA_CKFTrackFinder ckf;
      if(kCKFEnvelope>0) ckf.fCfg.fEnvelopeCutSigma = kCKFEnvelope;
      ckf.fVerbose = cli->count("-V");
//      hlog::info(hlog::kBlue, "Skeletor|CKF => Using CKF on {} DCT Hits:", evtS2.fDCTHits.size());
      ckf.setHits(evtS2.fDCTHits); //CKF will build the CKF hits and layers internally

      int hitsPerLayerCount = 5;
      int maxSeedCount = 200;
      int candCount = ckf.findCandidateTracks(hitsPerLayerCount, maxSeedCount);
//      hlog::info(hlog::kBlue, "Skeletor|CKF => Found {} candidate tracks", candCount);
      if (candCount > 0) {
        auto candTrack0 = ckf.getCandidateTrack(0, false); //get the best candidate track
        stage3DCTHits = candTrack0.getHits();
//        hlog::info(hlog::kBlue, "Skeletor|CKF => Selected Track: {:3d} | Charge: {:f} | NHits: {:3d}", 0, candTrack0.getCharge(), stage3DCTHits.size());
      } else {
//        hlog::warn("CKF => No Candidate Tracks found");
      }
    } else
    if(kUseSegmentFinder) {
      /// works for both BON and BOFF
      stage3DCTHits = buildDCTStage3HitsWithSegmentFinder(&evtS2, false);
    } else if(kUseStagger) {
      HXA_StaggerTrackFinder stagger;
      stage3DCTHits = stagger.buildStage3Hits(evtS2);
    }else {
      /// Use Hough transform to select hits, deciding between lines and circles based on Magnetic field
      if(kMagFieldIsOn){
        HXA_HoughCircleFinder houghCircleFinder;
        houghCircleFinder.addDCTHits(evtS2.fDCTHits, true); // add hits, filtering bad ones
        /// Perform the Hough transform for isochrones
        vector<HXA_HoughCircle> houghCircles = houghCircleFinder.findCircles();

        vector<HXA_HoughPoint> hitsFromBestHoughCircle;
        if (houghCircles.size() > 0){
          hitsFromBestHoughCircle = kUseIsochrones ? houghCircles[0].getContributingPoints(1) : houghCircles[0].getContributingPoints();
          hlog::info("DCT Hough Circles:");
          int houghIndex{0};
          for(auto const &circle: houghCircles) {
            hlog::info("Circle {:2d} => Score: {:4d}, Radius: {:8.3f}, CenterX: {:8.3f} CenterY: {:8.3f}", houghIndex++, circle.fScore, circle.fRadiusXY,circle.fCenterX,circle.fCenterY);
          }
        } else {
          hlog::info("No Hough circles found");
          continue;
        }
        stage3DCTHits = buildDCTStage3HitsFromHoughCircle(&evtS2, houghCircles[0], kUseIsochrones);
      }
      else {
        HXA_HoughLineFinder houghLineFinder;
        houghLineFinder.addDCTHits(evtS2.fDCTHits, true); // add hits, filtering bad ones.
        hlog::info(hlog::kOrange, "Hough Line Finder has {} points", houghLineFinder.fPoints.size());

        /// Perform the Hough transform
        vector<HXA_HoughLine> houghLines = houghLineFinder.findLines();

        vector<HXA_HoughPoint> hitsFromBestHoughLine;
        if (houghLines.size() > 0){
          hitsFromBestHoughLine = kUseIsochrones ? houghLines[0].getContributingPoints(1) : houghLines[0].getContributingPoints();

          hlog::info("DCT Hough Lines:");
          int houghIndex{0};
          for(auto const &line: houghLines) {
            hlog::info("  Line {:2d} => Score: {}", houghIndex++, line.fScore);
          }
        } else {
          hlog::info("No Hough lines found");
          continue;
        }
        stage3DCTHits = buildDCTStage3HitsFromHoughLine(&evtS2, houghLines[0], kUseIsochrones);
      }
    }

    /// do a final filter if you like
    if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "Stg3 HitCount < 36", stage3DCTHits.size() < 36,
      [&](){ hlog::info("Not enough Stg3 DCT hits ({}) for tracking", stage3DCTHits.size());}
    );
    selectedEventCount++;

    ///***********************************************************************************************//
    /// do the actual track fit
    ///***********************************************************************************************//
    /// build a track from the hits.
    HXA_Track track;
    track.setHits(stage3DCTHits);
    track.setCharge(1.0);

    /// build a track fitter with the track
    std::unique_ptr<HXA_TrackFitterBase> myFitter;
    if(kMagFieldIsOn){
      if(kUseRK) {
        myFitter = std::make_unique<HXA_TrackFitterRK>(&track);
      } else {
        myFitter = std::make_unique<HXA_TrackFitterGenFit>(&track);
      }
    } else {
      myFitter = std::make_unique<HXA_TrackFitterLinear>(&track);
    }

    /// do the fit and get the result
    HXA_TrackFitResult fit = myFitter->doFit();
    fit.dump(kMagFieldIsOn ? (kUseRK ? "RK Fit:" : "GF Fit:") : "Linear Fit:");

    // need to come up with a column count for the track object?
    if(kUseSingleCol && !kNoCuts){
      int numcolumns=track.getColumnCount();
      CutSummary_APPLY_CUT(cutSummary, "Num Cols > 1",
                         numcolumns > 1,
                         [&](){ hlog::info("THIS TRACK Crossed: {} columns", numcolumns); }
      );
      //hlog::info("THIS TRACK Crossed: {} columns", numcolumns);
    }
    if(kUseTwoCol && !kNoCuts){
      int numcolumns=track.getColumnCount();
      CutSummary_APPLY_CUT(cutSummary, "Num Cols != 2",
                         numcolumns !=2,
                         [&](){ hlog::info("THIS TRACK Crossed: {} columns", numcolumns); }
      );
      //hlog::info("THIS TRACK Crossed: {} columns", numcolumns);
    }

    if(kUseSingleSide && !kNoCuts){
      hlog::info("CHECKING IF THIS TRACK is on a single side");
      int numcolumns=track.getColumnCount();
      double driftsign_n=0;
      double driftsign_p=0;
      auto thesehits = track.getHits();
      for(auto hit: thesehits){
        if(hit.fDriftSign>0) driftsign_p+=1.0;
        else driftsign_n+=1.0;
      }
      //hlog::info("CHECKING IF THIS TRACK is on a single side {:.3f},{:.3f}",driftsign_p,driftsign_n);
      CutSummary_APPLY_CUT(cutSummary, "More than 1 side",
                         (driftsign_p>0 && driftsign_n>0) || numcolumns > 1,
                         [&](){ hlog::info("THIS TRACK crossed planes"); }
      );
    }


    ///=================================================================================================================
    /// we can optionally move around the stage3 hits based on the isochrones and the original fit.
    /// better would be to do this on stg2 hits, then re-do the Hough, but that is for another time.
    if(kDoIsochrone2Pass && !kUseIsochrones) {
      hlog::warn("You have requested a 2-pass isochrone fit, but isochrones were not used in the first pass. Ignoring...");
    } else if (kDoIsochrone2Pass && kUseIsochrones) {
      hlog::info(hlog::kYellow, "Doing 2nd pass isochrone adjustment of hit positions based on initial fit");
      
      auto trackHits = track.getHits(); //grab copy of the hits to work on
      auto iso2stepper = kMagFieldIsOn ? (kUseRK ? HXA_StepperBase::Type::kNystrom : HXA_StepperBase::Type::kGenFit) : HXA_StepperBase::Type::kLinear;
      for(auto &hit : trackHits){
        if(hit.fType != HLX_Stg3_Hit::Type::kDCT) continue;
        double slope = fit.ps().stepToZPlane(hit.fZ, iso2stepper).getMy();
        auto hitID = hit.fHitID;
        auto wireID = theGeom->fDCT.getWireIDFromHitID(hitID);
        /// destagger has already occurred, so there should only be a single Y position left
        double signedDriftTime = hit.fDriftTime*hit.fDriftSign;
        HXA_IsochroneCurve isoCurve;
        if(kMagFieldIsOn){
          isoCurve = HXA_IsochroneManager::getInstance()->getIsochroneFittedCurve(signedDriftTime, wireID, hit.fX);
          //isoCurve.dump();
        }
        else{
          isoCurve = HXA_IsochroneManager::getInstance()->getIsochroneCurve(signedDriftTime, wireID);
        }

        if(!isoCurve.fIsOK) {
          if(cli->count("-V")>1) hlog::error("Process => Isochrone curve is empty for wireID: {:3d} and drift time: {:8.3f}", wireID, signedDriftTime);
          continue; // skip this hit
        }
        auto [newY, newZ] = isoCurve.getYZAtClosestSlope(slope, 0.25,0.75);
        //hlog::info("HitID: {:8d} | WireID: {:3d} | Orig Y: {:8.3f} | Z: {:8.3f} | DriftTime: {:8.3f} | Sign: {:2d} | SlopeYZ: {:8.5f} | New Y:{:8.3f} | New Z: {:8.3f}", hitID, wireID, hit.fY, hit.fZ, hit.fDriftTime, hit.fDriftSign, slope, newY, newZ);
        hit.fY = newY; // set the new Y position
        hit.fZ = newZ; // set the new Z position
      }// hit loop
      ///==============================================================================================================
      /// replace the hits in the track and redo the fit
      track.setHits(trackHits);
      myFitter->setTrack(&track); // make sure the fitter has the updated track
      fit = myFitter->doFit();
      fit.dump(kMagFieldIsOn ? "GenFit 2nd Pass" : "Linear 2nd Pass");
    }// kDoIsochrone2Pass

    if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "Fit Status < 0",
                         fit.fStatus < 0 || std::isinf(fit.fStatus) || std::isnan(fit.fStatus),
                         [&](){ hlog::info("Bad fit status: {}", fit.fStatus); }
    );
    track.addFit(fit);
    track.dumpHits();

    /// fill the chi2 histogram before cutting
    double thisRedChi2 = fit.getReducedChi2();
    h_Chi2Pre->Fill(thisRedChi2);

    if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "Bad Chi2", thisRedChi2 < 0.0 || thisRedChi2 > 200.0 || !std::isfinite(thisRedChi2),
                         [&](){ hlog::info("Bad reduced chi2: {}", thisRedChi2); }
    );

    /// Save hitcount
    hHitCount->Fill(track.getHits().size());
    h_MDR->Fill(fit.ps().getMDR());
    h_Chi2Post->Fill(fit.getReducedChi2());

    /// use the fit to get hit locations at the TOF layers - this is the linear extrapolation
    auto stepper = kMagFieldIsOn ? HXA_StepperBase::Type::kGenFit : HXA_StepperBase::Type::kLinear;
    auto [ttofX, ttofY] = fit.ps().getXYAtZ(theGeom->fTOF.fTopLocZ,         stepper);
    auto [boreX, boreY] = fit.ps().getXYAtZ(theGeom->fTOF.fBorePaddleLocZ,  stepper);
    auto [btofX, btofY] = fit.ps().getXYAtZ(theGeom->fTOF.fBottomLocZ,      stepper);
    auto [richX, richY] = fit.ps().getXYAtZ(theGeom->fRICH.fFocalPlaneLocZ, stepper);

    /// For DCT Ficucial
    pair<double, double> dctZs = theGeom->fDCT.getTopAndBotZ();
    double tdctZ = dctZs.first;
    double bdctZ = dctZs.second;
    pair<double, double> tdct = fit.ps().getXYAtZ(tdctZ, stepper);
    double tdctX = tdct.first;
    double tdctY = tdct.second;
    pair<double, double> bdct = fit.ps().getXYAtZ(bdctZ, stepper);
    double bdctX = bdct.first;
    double bdctY = bdct.second;
    bool isInsideDCTTop = theGeom->fDCT.isInside(tdctX, tdctY, tdctZ, 10.0); // top hit is within 10mm of vol
    bool isInsideDCTBot = theGeom->fDCT.isInside(bdctX, bdctY, bdctZ, 10.0);
    if(kFiducialCut) {
      if(!kNoCuts)CutSummary_APPLY_CUT(cutSummary, "DCT Fiduc.", !isInsideDCTBot || !isInsideDCTTop,
         [=]() {hlog::info("DCT Fid Fail! Top:{:8.3f}/{:8.3f}/{:8.3f} | Bot:{:8.3f}/{:8.3f}/{:8.3f}", tdctX, tdctY, tdctZ, bdctX, bdctY, bdctZ);}
      );
    }

    ///================================================================================================================
    ///TOF Beta and Charge
    
    double pathLength = fit.ps().getPathlengthBetweenZPlanes(ttofHits[0].getZ(), btofHits[0].getZ(), stepper);
    double timeDiff = btofHits[0].fTime - ttofHits[0].fTime;
    double topSig = ttofHits[0].getSignal();
    double botSig = btofHits[0].getSignal();
    double boreSig = boreHits[0].getSignal();
    double topPLC = fit.ps().stepToZPlane(ttofHits[0].getZ()).getPLC();
    double botPLC = fit.ps().stepToZPlane(btofHits[0].getZ()).getPLC();
    double borePLC = fit.ps().stepToZPlane(boreHits[0].getZ()).getPLC();

    double ttofZ = pow(topSig / topPLC, 1.0/1.7);
    double btofZ = pow(botSig / botPLC, 1.0/1.7);
    double boreZ = pow(boreSig / borePLC, 1.0/1.7);
    hTTOF->Fill(ttofZ);
    hBTOF->Fill(btofZ);
    
    /// record the residuals, using the stepper of your choice (also used below for non-contributing hits)
    //HXA_StepperBase::Type stepperType = kMagFieldIsOn ? HXA_StepperBase::Type::kGenFit : HXA_StepperBase::Type::kLinear;
    HXA_StepperBase::Type stepperType = HXA_StepperBase::Type::kNystrom;  // another opt: HXA_StepperBase::Type::kNystrom
    auto fitResiduals = track.getResiduals(fit, stepperType);
    auto thisFOM=track.getFigureOfMerit(fitResiduals);
    h_FOM->Fill(thisFOM);
    auto thisRig = fit.ps().getRig();
    auto thisRigErr = fit.ps().getRigErr();
    h_Rig->Fill(thisRig);
    h_RigErr->Fill(thisRigErr);
    h_Eta->Fill(fit.ps().getEta());

    resSum.fillHistos(fitResiduals);
    pullSum.fillHistos(fitResiduals);

    // loop over residuals to fill histos
    for (auto res : fitResiduals) {
      h_ResVPosZ->Fill(res.fPosZ, res.fResY);
      h_ResVDrift->Fill(res.fDriftTime, res.fResY);
      h_ResAll->Fill(res.fResY);
      if(res.fDriftTime>10.0 || res.fDriftTime<-10.0){
        h_ResVDrift_filtered->Fill(res.fDriftTime, res.fResY);
      }
      if (res.fStaggerDist > 0) {
        h_ResPos->Fill(res.fResY);
      } else {
        h_ResNeg->Fill(res.fResY);
      }
    }/// residual loop

    /// extrapolate the fit to the TOF layers and fill the maps
    resSum.fillNonContributingHistos(fit.ps(), stage3TOFHits, stepperType);
    pullSum.fillNonContributingHistos(fit.ps(), stage3TOFHits, stepperType);

    //hlog::info(hlog::kBurlyWood, "TOF Top: ({:8.3f}, {:8.3f}) | Bore: ({:8.3f}, {:8.3f}) | Bottom: ({:8.3f}, {:8.3f})", ttofX, ttofY, boreX, boreY, btofX, btofY);

    ttofMap->Fill(ttofX, ttofY);
    boreMap->Fill(boreX, boreY);
    btofMap->Fill(btofX, btofY);
    
    double tofY[3] = {ttofY, boreY, btofY}; //Order: top, bore, bottom 

    ///================================================================================================================
    /// hodoscope hit
    /// use the track fit to get a hit location at the hodo
    auto [trkX_hodo, trkY_hodo] = fit.ps().getXYAtZ(theGeom->fHodo.fFiberPlaneLocZ, HXA_StepperBase::Type::kLinear);

    /// Hodoscope - coarse precision
    auto [ hitX_hodoC, hitXerr_hodoC ] = hodoCalMan->getBestFiber_CoarseX(hodoID, trkY_hodo );
    Double_t resX_hodoC = hitX_hodoC - trkX_hodo;
    hHodoC->Fill( hitX_hodoC, resX_hodoC );
    
    /// Hodoscope - fine precision. 
    /// [2024/10/02] Xerr is fixed as -1 temporally. Will be updated later
    auto [ fiberIDF, hitX_hodoF, hitXerr_hodoF ] = hodoCalMan->getBestFiber_FineX(hodoID, trkX_hodo, trkY_hodo );
    Double_t resX_hodoF = hitX_hodoF - trkX_hodo;
    hHodoF->Fill( hitX_hodoF, resX_hodoF );

    ///================================================================================================================
    /// RICH Direct hits
    //  RICH reconstruction requires an estimate of particle charge
    //  Use TOF to get this however you please
    int elementalChargeEstimate = 1; 
    map<int, double> richDirectHitThresholdMap = {
      {1, 20.0}, // Hydrogen
      {2, 30.0}, // Helium
      {3, 40.0}, // Lithium
      {4, 50.0}, // Beryllium
      {5, 60.0}, // Boron
      {6, 70.0}, // Carbon
      {7, 80.0}, // Nitrogen
      {8, 90.0}, // Oxygen
    };
    auto richDirectHits = richRecon.getRICHDirectHits(&evtS2, richDirectHitThresholdMap.at(elementalChargeEstimate), {richX, richY}, 200.0);
    /// extrapolate the fit to the RICH layer and fill the maps
    resSum.fillNonContributingHistos(fit.ps(), richDirectHits, stepperType);
    pullSum.fillNonContributingHistos(fit.ps(), richDirectHits, stepperType);

    ///=================================================================================================================
    /// RICH Beta Reconstruction
    HLX_Stg3_Event evtS3;
    evtS3.fTrackingHits = stage3DCTHits;
    for (auto &hit : evtS2.fRICHHits) {
      HLX_Stg3_RICHHit h(hit.getX(), hit.getY(), hit.getZ());
      h.fSignal = hit.fSignal;
      h.fTDC = static_cast<int8_t>(hit.fTDC.to_ulong());
      h.fChanID = hit.fChanID;
      h.fStatus = HXA_RICHReconMinimizer::fHitStatusMap.at(HXA_RICHHitStatusType::kGood); //initialize all hits as good
      evtS3.fRICHHits.push_back(h);
    }

    richRecon.filterBadPixels(&evtS3, richDirectHitThresholdMap.at(elementalChargeEstimate));
    richRecon.filterCherenkovRing(&evtS3);
    richRecon.setParticleState(&fit.ps(), stepperType);

    if (kDoRICHTrackMinimization){
            fullTrackRichMin.setEvent(&evtS3);
            fullTrackRichMin.setParticleState(&richRecon, stepperType);
            fullTrackRichMin.minimizeParticleState(elementalChargeEstimate);
            fullTrackRichMin.updateParticleState(&richRecon);
    }

    richRecon.reconstructPixelHits(&evtS3);
    double medianPhi = HXA_RICHReconMinimizer::getMedianPhi(&evtS3);
    double pixelSignalThreshold = 0.0; // do not reject any pixels based on minimum signal size
    pair<double, double> richBetaEst = richRecon.getBetaEstimate(&evtS3, elementalChargeEstimate, pixelSignalThreshold);
    h_RICHBeta->Fill(richBetaEst.first);

    //Get two additional estimates for beta, considering half the pixels only each time for self correlation
    pair<double, double> richBetaEstLeft = richRecon.getBetaEstimate(&evtS3, elementalChargeEstimate, pixelSignalThreshold, 0.0, medianPhi);
    pair<double, double> richBetaEstRight = richRecon.getBetaEstimate(&evtS3, elementalChargeEstimate, pixelSignalThreshold, medianPhi, 2.0 * hlx::units::pi);
    if (richBetaEst.first > 0 && isfinite(richBetaEst.first)) {
        h_RICHRes->Fill((richBetaEstLeft.first - richBetaEstRight.first) / richBetaEst.first);
    }


    ///=================================================================================================================
    /// RICH Charge Reconstruction
    double richChargeEst = -1.0;
    if (richRecon.fIndexUsed * richBetaEst.first > 1.0){
      double cherenkovAngleReconstructed = acos(1 / (richRecon.fIndexUsed * richBetaEst.first));
      richChargeEst = richRecon.getChargeEstimate(&evtS3, cherenkovAngleReconstructed);
      h_RICHCharge->Fill(richChargeEst);
    }

    // TOF TTree - SN
    ///---///---///-
    //Tree variables
    topTOF  = &ttofHits[0];
    boreTOF = &boreHits[0];
    botTOF  = &btofHits[0];
    
    Zhi = Zlo = -1;
    if (evtS1->fTrigger.size() == 1)
    {
      Zhi = (int)evtS1->getTrigger()->getTriggerInfo().fBits.fZHi;
      Zlo = (int)evtS1->getTrigger()->getTriggerInfo().fBits.fZLo;      
    }
//    cout << endl << " **************** ZHi = " << Zhi << " ***************" << endl << endl;
//    cout << endl << " **************** ZLo = " << Zlo << " ***************" << endl << endl;    
    tvars.ttofX = ttofX; 
    tvars.ttofY = ttofY;
    tvars.boreX = boreX;
    tvars.boreY = boreY;
    tvars.btofX = btofX;
    tvars.btofY = btofY;
    tvars.richX = richX;
    tvars.richY = richY;
    tvars.pathLength = pathLength;
    tvars.topPLC = topPLC;
    tvars.borePLC = borePLC;
    tvars.botPLC = botPLC;
    
    runnum = evtS2.fRunNumber;
    evttime = HLX_TimeStamp::dTD(evtS2_orig.fTimeStamp); //in days since official launch

    richZEst = richChargeEst;
    DCTChi2 = fit.getReducedChi2();
    rig = fit.ps().getRig();
    rig = ( std::isfinite(rig) )? rig : 0; //save reasonable values to ttree
    mdr = fit.ps().getMDR();
    hitCount = track.getHits().size();
    for (int ij=0; ij<3; ij++)
    {
      topdircos[ij]  = (fit.ps().stepToZPlane(ttofHits[0].getZ()).getDirection())(ij);
      boredircos[ij] = (fit.ps().stepToZPlane(boreHits[0].getZ()).getDirection())(ij);
      botdircos[ij]  = (fit.ps().stepToZPlane(btofHits[0].getZ()).getDirection())(ij);
    }
    
    toft->Fill();
    
    //Watch this trick: set up arrays to avoid lots of repeated code.
    HLX_Stg2_TOFHit* tofHits[3] = {&ttofHits[0], &boreHits[0], &btofHits[0]};  //not copies, so use pointer ->
    
    ///---///---///


    ///=================================================================================================================
    if(kInteractive) {
      auto dctViewer = new HLX_HitViewer("DCT Viewer"); //DCT Zoom
      dctViewer->drawDetectorDCT(evtS2.fEventNumber);
      //dctViewer->fViewerDCT->drawStage2HitsLR(&evtS2Original); //draws all the hits with both LR options
      dctViewer->fViewerDCT->drawStage2Event(&evtS2, false); //draws all LR hits
      dctViewer->fViewerDCT->drawStage3Hits(stage3DCTHits, 28, 47, 1.25); //draws the hits that made it to stg3 after selection
      dctViewer->fViewerDCT->drawFitLine(fit, kGreen, 1, 3.0);
      dctViewer->fCanv->Update();
      HLX_Getch::pressAnyKey(format(hlog::kOrange, "Press any key to for next event!"));
      delete dctViewer; // recreated every event, must not leak
    }
    /// if writing a tree, do it here
    if(kWriteOutTree) {
      outTrack = &track;
      outTree->Fill();
    }

  }//main event loop
  totalSW.stop();

  hlog::info(hlog::kOrange, "Final Memory usage: {:6.3f} MB", HLX_Util::getCurrentMemoryUsageMB());
  hlog::info(hlog::kOrange, "Total Processing Time: {} | Rate: {:6.2f} Hz", totalSW.getBestStartStopTime(), double(totalEntryCount)/totalSW.getStartStopTime());
  hlog::info(hlog::kOrange, "Selected {}/{} events (efficiency: {:5.3f})", selectedEventCount, totalEntryCount, (double)selectedEventCount/(double)totalEntryCount);

  /// Dump cut summary
  cutSummary.dump();
  /// Dump cut summary
    /// Save the residual summary to a file

  if(kShowPlots) {
    TString pdfout = kOPFName;
    pdfout.ReplaceAll(".root",".pdf");
    TString pdfoutopen = pdfout+"(";

    auto c0 = new TCanvas("c0", "c0", 0, 0, 800, 800);
    c0->cd(1);
    h_Chi2Pre->SetTitle("Chi2 Distribution;#chi^{2};Entries");
    h_Chi2Pre->Draw();
    c0->Print(pdfoutopen.Data());

    auto c1 = new TCanvas("c1", "c1", 0,0, 1600, 500);
    c1->Divide(3,1);
    c1->cd(1);
    ttofMap->Draw("colz");
    c1->cd(2);
    boreMap->Draw("colz");
    c1->cd(3);
    btofMap->Draw("colz");
    c1->Modified();
    c1->Update();
    c1->Print(pdfout.Data());

    auto c2 = new TCanvas("c2", "c2", 0, 600, 1600, 500);
    c2->Divide(3,1);
    c2->cd(1); gPad->SetGridy();
    gPad->SetLeftMargin(0.13);
    hHodoC->Draw("colz");
    c2->cd(2); gPad->SetGridy();
    gPad->SetLeftMargin(0.13);
    hHodoF->Draw("colz");
    c2->cd(3); gPad->SetGridx();
    gPad->SetLeftMargin(0.13);
    hHodoF->ProjectionY()->Draw();
    c2->Print(pdfout.Data());

    /// TOF delta-T
    auto c3 = new TCanvas("c3", "c3", 0, 0, 800, 600);
    hTOFdT->Draw();
    c3->Print(pdfout.Data());

    /// TOF Zeds
    auto c6 = new TCanvas("c6", "c6", 1200, 0, 1150, 600);
    c6->Divide(2,1);
    c6->cd(1);
    hTTOF->Draw();
    c6->cd(2);
    hBTOF->Draw();
    c6->Print(pdfout.Data());

    /// Residuals
    resSum.drawCanvas(true);

    double hiRes = 1.0;
    gStyle->SetOptStat(1101);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1101);
    auto cc = new TCanvas("cc", "Stagger Plots", 0, 0, hiRes*1200, hiRes*800);
    cc->Divide(2, 2);
    cc->cd(1);
    h_ResVPosZ->Draw("COLZ2");
    cc->cd(2);
    h_ResAll->Draw();
    h_ResPos->SetFillColorAlpha(8, 0.5);
    h_ResPos->Draw("sames");
    h_ResNeg->SetFillColorAlpha(9, 0.5);
    h_ResNeg->Draw("sames");
    cc->Modified();
    cc->Update();//call to cause the stats to be drawn
    

    auto fit = new TF1("fit", "gaus", -0.5, 0.5);
    fit->SetLineWidth(2);
    h_ResAll->Fit("fit", "R", "", -0.5, 0.5);
    cc->cd(3);
    h_ResVDrift->Draw("COLZ2");
    gPad->Modified();
    gPad->Update();
    TString pdfoutclose = pdfout+")";
    cc->Print(pdfoutclose.Data());

  }

  ///=================================================================================================================
  /// Save all the shit to the output file
  opf.cd();
  resSum.Write("residualSummary");
  pullSum.Write("pullSummary");
  cutSummary.Write("cutSummary");

  ///---///---///
  //tof tree -SN
  toft->Write();
  ///---///---///
  
  unique_ptr<HXA_ResultHeader> resHdr(HXA_ResultHeader::build());
  resHdr->Write();

  opf.Write(); //this will write any objects connected to the file
  opf.Close();

  if(theApp) delete theApp;

  return 0;
}


















































void splash(){
  string splash = R"(%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@@@@@@@@@@@@@@@@@@@ Skeletor has come to analyze your data! @@@@@@@@@@@@@@@@@@@@
&&%%%##%#(@@@@@@@@@@@@@@@@@@@@@@@@@&&@&&&&%%%&@@@@@@@&&&&&&&&%%%%%%%%%%%%%%%%%&%
&&&%%%##%#&@@@@@@@@@@@@@@@@@@@@@@@@&&&@&&%%%#&@@@@@@@@&%%%%%&&&&%%&%%%%%%%%%%##&
&%%#(#(((%@@@@@@@@@@@@@@@@@@@@@@@@@@&&&@&%%(&@@@@@@@@@@@@@@&&&&%%%&%%%%%%%%#%@@@
&%###((((&@@@@@@@@@@@@@@@@@@@@@@@@@@&&&&%&%##@@@@@@@@@@@@@@@%&&%%%%%%##%&@@@@@@@
%#%#((((&@@@@@@@@@@@@@@@@@@@&%@&%#%%@@@@&%#&&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
&#%(((&@@@@@@@@@@@@@@@@@&%&@@&#/*/#%###%%%####%%&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#(##&@@@@@@@@@@@@@@@@@&%@@/.........,#%####%#########%%@@@@@@@@@@@@@@@@@@@@@@@@@
%%@@@@@@@@@@@@@@@@@@@%@@(..............,%#################%&@@@@@@@@@@@@@@@@@@@@
&@@@@@@@@@@@@@@@@@@&&@&,..................####################%@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@%&@%......................#################%##%@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@%&@@,.......................,%###############&%##%@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@%%@%,../#*........,(/....,**,..#%##############&%###%@@@@@@@@@@@
@@@@@@@@@@@@@@@&%@#,*......*....,,..../&@@@%..../%##############&%###%@@@@@@@@@@
@@@@@@@@@&&&&&@#@%..#@@@%/..,/(/..*%@@@@@@@@#....#%#############%@%%##%@@@@@@@@@
@@@@@@@&&%%%%%&%@@/,&@@@@@@&*(....(&@@@@@@@@#.....&##############%@@@&%%@@@@@@@@
@@@@@&%%%%%%#%%%@@(./%@@@@&*&@@%...,/#&@&(,.......(%##############@@@@@@@@@@@@@@
@@@@&&%%%%##%%%&@#........,&&*#&#...,............,#@##############%@@@@@@@@@@@@@
@&&&%%&%####(%%&@(.../,................,(/,*(&@@@@@@##############%@@@@@@@@@@@@@
&&&%#%####(#&@#@@@@@%,....................,&@#,,#&@@###############&@@@@@@@@@@@@
####%(/%&&%%%@#&@@@%%#///(,*,.(,/**/,/,..#@@&..../@&#########%%####%@@@@@@@@@@@@
(###&@@@&%&%&@%%@@@@@*#,((.*,,%,,*./.(&@@@@@*...,&@&##########%&####@@@@@@@@@@@@
(&&&@@@@@&%&%@%%@@@@@@/(@@@@@@@@@@@@@@@@@@@&....*&@#############&&##&@@@@@@@@@@@
%%%%&@@@@&%%&&%#&@@@@@@#/,**,*,.*.,,.%@@@@@*....,&&###############&&%@@@@@@@@@@&
%%%%&%@@&@&%&&%#&@@@@@@#...,...,.,..,...........,&##################%@&(&@@@@&&&
%%%%%%&@&&&&#&&##@@@@@@#.......................(&#####################%@&@&&&&@&
%%%&%%&@@&&&&%&##%@@@@@&,............*#&@@@@@@@&########################%@&&&&&&
%%%%%%%&@@&@&%&%##%@@@@@@@@@@@@@@@@@@@@@@@@@@@%###########################&@&&&@
%%%%%%%&@@@@&%%&####%@@@@@@@@@@@@@@@@@@@@@@@&######################%%%#####&&&&@
%%%%%%%%@@@@@&&@&######&@@@@@@@@@@@@@@@@@@%################%%&&%%%%&&@@&%%%#@&&&
%%%%%%%%@@@@%%%%%%&%######%&@@@@@@@@@@@&##############%%&%%%%%%%%%%%%%%%%&@@%@&&
&@@@@@&%%&%%%%%%%%%%%%&%########%%&%#############%&&%%%%%%%%%%%%%%%%%%%%%%%%%&@&
%%%%%@%%%%%%%%%%%%%%%%%%%%&%###############%&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%@%%%%%%%%%%%%%%%%%%%%#%%%&@&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
/////@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//#//&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%)";

  hlog::info(hlog::kGreen, splash);
}

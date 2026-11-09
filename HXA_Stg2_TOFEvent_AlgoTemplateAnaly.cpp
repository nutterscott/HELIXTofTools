////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Filename: HXA_Stg2_TOFEvent_AlgoTemplateAnaly.cpp
/// Created:  2026/04/10 by SLN (modified Sakai code)
/// Created:  2025/08/05 by sakai (refactored from original TOF hitbuilder)
/// PSU/Carl's original TOF processing algorithm implementation
/// but updated with Coutu's time dependent calibrations and ramp fit using Gerard Visser formula.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "HXA_Stg2_TOFEvent_AlgoTemplateAnaly.h"
#include "HLX_Stg1_Event.h"
#include "HLX_Stg2_Event.h"

//#include <TMath.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Configure the PSU algorithm with root file and steering parameters
/// \param steeringInfo - TOML table with configuration parameters
/// \returns 0 on success, negative on error
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::configure(const string& pedFileName, const HLX_TOML_Table& steeringInfo)
{
  fProcessMode = eProcessingMode::kProcess; 

  // Load calibration data
  loadSCTOFCalibs();
  
  //define TGRaph and ramp fitting functions:
  defineGVRampFitFuncs();
  
  //use dynamic peds as calculated in this routine
  fDynamicPeds = true;
  
  hlog::info("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::configure => TemplateAnaly algorithm configured successfully");
  return kSuccess;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Configure the PSU algorithm while ignoring pedestal file and steering parameters 
/// \returns 0 on success, negative on error
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
int HXA_Stg2_TOFEvent_AlgoTemplateDirect::configure()
{
  fProcessMode = eProcessingMode::kReprocess;

  // Load calibration data
  if (!fTOFcal) loadSCTOFCalibs();
  //define TGRaph and ramp fitting functions:
  if (!fGVup) defineGVRampFitFuncs();

  hlog::info("HXA_Stg2_TOFEvent_AlgoTemplateDirect::configure => TemplateDirect Reprocess algorithm configured successfully");
  return kSuccess;
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Process a stage1 event to produce TOF hits in the stage2 event
/// \param evtS1 - input stage1 event
/// \param evtS2 - output stage2 event
/// \returns 0 on success, negative on error
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::processEvent(HLX_Stg1_Event* evtS1, HLX_Stg2_Event* evtS2)
{
  if (!evtS1 || !evtS2) {
    hlog::error("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::processEvent => Bad s1 or s2 event!");
    return kErrorBadEvent;
  }
  if (evtS1->fTOF.empty()) {
    hlog::error("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::processEvent => No TOF subevents found in input S1 event!");
    return kErrorBadEvent;
  }
  
  ///DEBUG:
  if (DebugMode > 0) cout << endl << endl <<  "******* HXA_Stg2_TOFEvent_AlgoTemplateAnaly PROCESSING NEW EVENT " << fEvtCnt << " ********" << endl;
  fEvtCnt++;
  clearDebug();

  //update calibs for current run if needed
  int runnum = evtS2->fRunNumber;
  getSCTOFCalibsForRun(runnum);
  
  // Step 1: Build FEE hits from stage1 data
  vector<HLX_Stg2_TOFFEEHit> feeHitVec;
  int result = assembleFEEHits(evtS1, evtS2, feeHitVec);
  if (result != kSuccess) return result;

  // Step 2: Build paddle-end hits from FEE hits
  map<HLX_Geometry_TOF_PaddleEndTag, HLX_Stg2_TOFPaddleEndHit> padEndHitMap;
  result = assemblePaddleEndHits(feeHitVec, padEndHitMap);
  if (result != kSuccess) return result;

  // Step 3: Build paddle hits from paddle-end hits and add to event
  result = assemblePaddleHits(padEndHitMap, evtS2);
  
  ///DEBUG output summary
  //only output for events with one top and one bottom hit
  if (DebugMode > 0)
  {
    auto theGeom = HLX_Geometry::getInstance();
    
    vector<HLX_Stg2_TOFHit> ttofHits, boreHits, btofHits;
    //collect top bottom and bore hits
    for(auto tofHit2: evtS2->fTOFHits){
      if (!tofHit2.isHiZHit() && !tofHit2.isLoZHit()) continue;
      auto layerName = theGeom->fTOF.getLayerNameIDByZ(tofHit2.getZ());
      if(layerName == HLX_Geometry_TOF::LayerName::kTop){
        ttofHits.push_back(tofHit2);
      }
      if(layerName == HLX_Geometry_TOF::LayerName::kBore){
        boreHits.push_back(tofHit2);
      }
      if(layerName == HLX_Geometry_TOF::LayerName::kBottom){
        btofHits.push_back(tofHit2);
      }
    }
    //skip if too many paddles hit
    if (ttofHits.size() != 1 || boreHits.size() != 1 || btofHits.size() != 1) return result; //skip DEBUG output for this event
    
    //output relevant info. Assume proper end and fee ordering in paddle hit object
    HLX_Stg2_TOFPaddleEndHit tE = ttofHits[0].getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit tW = ttofHits[0].getPaddleEndHit(HLX_Geometry_TOF::kWest);
    HLX_Stg2_TOFPaddleEndHit bE = btofHits[0].getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit bW = btofHits[0].getPaddleEndHit(HLX_Geometry_TOF::kWest);
    HLX_Stg2_TOFPaddleEndHit boreE = boreHits[0].getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit boreW = boreHits[0].getPaddleEndHit(HLX_Geometry_TOF::kWest);

    HLX_Stg2_TOFPaddleEndHit endHit[6] = {tE, tW, bE, bW, boreE, boreW};
  
    cout << "pad, endID, feeID, SCfeeID, Rtype, Zped, Rped, Status, GVfitT, udCorr, t" << endl;
    for (int i=0; i<6; i++)
    {
      int nFees     = endHit[i].getFEEHitCount();
      int paddleID  = endHit[i].fTag.getPaddleID();
      int endID     =  endHit[i].fTag.getPaddleEndID();
      
//      if (paddleID == 16) continue; //i don't know how to deal with dT for bore paddle
      if ( nFees < 2 ) continue;
//      if ( endHit[i].getFEEHit(0).fStatusFlag != HLX_Stg2_TOFFEEHit::kGood ||  endHit[i].getFEEHit(1).fStatusFlag != HLX_Stg2_TOFFEEHit::kGood) continue;
      double dT[4] = {0};
      for (int j=0; j<nFees; j++)
      {
        HLX_Stg2_TOFFEEHit feeHit = endHit[i].getFEEHit(j);
        string feeStat = "GOOD";
        if (feeHit.fStatusFlag == HLX_Stg2_TOFFEEHit::kBad) feeStat="BAD";
        int feeID = feeHit.fTag.getFEEID();
        auto [laneIdx, chanIdx] = getGeom()->fTOF.fGeoMap[{paddleID, endID,  feeID}];
        double feeTime = endHit[i].getFEEHit(j).fTime;
        dT[j] = feeTime;
        int globFeeID = laneIdx*8+chanIdx;
        int scFeeID = globFeeID; //tofcal12
//        int scFeeID = fSCFeeIDMap[globFeeID]; //tofcal11
//        std::streamsize original_precision = std::cout.precision();
        cout << paddleID << "  " << endID << "  " << feeID << "  " << scFeeID << "  " << fDBRType[globFeeID] << "  " 
          << feeHit.fPulsePedestal << " " << feeHit.fRampPedestal << "  " << feeStat << "  " << fDBFT[globFeeID] << "  " << fDBRampdT[globFeeID] << std::fixed << std::setprecision(3) << "  " << feeHit.fTime << endl;
//        std::cout.precision(original_precision);
        // Important: Reset 'fixed' or 'scientific' if they were used
//        std::cout << std::defaultfloat; 
      }
      cout << paddleID << "  " << endID << " dT (1-0)  (2-1) (3-2): ";
      for (int j=1; j<nFees; j++)
      {
        cout << dT[j] - dT[j-1] << "    ";
      }
      cout << endl;
    }
    
    cout << "TimeDiff, Pos: " << ttofHits[0].fTimeDiff << "  " << ttofHits[0].fTimeDiffPos << endl;
    cout << "ChargeDiv, Pos: " << ttofHits[0].fChargeDivFrac << "  " << ttofHits[0].fChargeDivPos << endl;
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Assemble FEE hits from stage1 TOF data
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::assembleFEEHits(HLX_Stg1_Event* evtS1, HLX_Stg2_Event* evtS2, vector<HLX_Stg2_TOFFEEHit>& feeHitVec)
{
  for(const auto &tofBoard: evtS1->fTOF[0].fTOFBoards) {
    auto laneIdx  = tofBoard.fFEELaneID;
    auto tofclock = tofBoard.fTOFClock;

    // Should always be 8 channels for working boards
    for(int chanIdx = 0; chanIdx < tofBoard.fChan.size(); chanIdx++) {
      auto theChan = tofBoard.fChan[chanIdx];

      // Get physical/geo info
      auto [paddleID, endID, feeID] = getGeom()->fTOF.getGeoInfo(laneIdx, chanIdx);
      auto tag = HLX_Geometry_TOF_FEETag(paddleID, endID, feeID);
      auto [loZHitCount, hiZHitCount] = theChan.countPulseTrigFlags();

      // Add it all to a hit
      HLX_Stg2_TOFFEEHit feeHit;
      feeHit.fTag = HLX_Geometry_TOF_FEETag(paddleID, endID, feeID);
      feeHit.fTOFClock     = tofclock;
      feeHit.fPulseData    = theChan.getPulseData();
      feeHit.fRampData     = theChan.getRampData();
      feeHit.fHitCountLoZ  = loZHitCount;
      feeHit.fHitCountHiZ  = hiZHitCount;
      feeHit.fStatusFlag   = HLX_Stg2_TOFFEEHit::kGood;

      // This is where the ped subtraction, gain correction, etc. happens
      if(processTOFFEEHit_SC(feeHit) < 0)
        feeHit.fStatusFlag = HLX_Stg2_TOFFEEHit::kBad;

      feeHitVec.push_back(feeHit);
    }
  }
  return kSuccess;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Assemble paddle-end hits from FEE hits
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::assemblePaddleEndHits(const vector<HLX_Stg2_TOFFEEHit>& feeHitVec, map<HLX_Geometry_TOF_PaddleEndTag, HLX_Stg2_TOFPaddleEndHit>& padEndHitMap)
{
  for(const auto &feeHit: feeHitVec) {
    auto paddleEndTag = feeHit.fTag.getPaddleEndTag();
    // See if this is the first hit for this paddle end
    if (padEndHitMap.count(paddleEndTag) == 0) {
      padEndHitMap[paddleEndTag] = HLX_Stg2_TOFPaddleEndHit{};
      padEndHitMap[paddleEndTag].fTag = paddleEndTag;
    }
    // Add the hit to the vector
    padEndHitMap[paddleEndTag].fFEEHits.push_back(feeHit);
  }

  // Sort the FEE hits in each paddle-end by geo order for simple access
  for(auto &[tag, endHit]: padEndHitMap) {
    std::sort(endHit.fFEEHits.begin(), endHit.fFEEHits.end(), [](const auto &a, const auto &b) {
      return a.fTag.getFEEID() < b.fTag.getFEEID();
    });
  }
  
  if (DebugMode > 4)
  {
    cout << "HXA_Stg2_TOFEvent_AlgoTemplateAnaly::assemblePaddleEndHits" << endl;
    for (const auto &[endTag, endHit]: padEndHitMap) 
    {
      int paddleID    = endTag.getPaddleID();
      int paddleEndID = endTag.getPaddleEndID();
    
      std::cout << "Paddle, end: " << paddleID << ", " << paddleEndID << std::endl;
      for (int i=0; i<endHit.getFEEHitCount(); i++) 
      {
        int FEEID       = endHit.fFEEHits[i].fTag.getFEEID();
        cout << "Fee ID" << FEEID << endl;
      }
    }
  }
      
  //// DO UP/DOWN RAMP TIMING CORRECTIONS HERE IN SAMPLE UNITS!!
  //SC corrects for up/down timing in FEE pairs at ends. 
  // Here goes:
  // Top and bottom FEEs:
  //  First FEE: up -  no correction
  //      Second FEE: up -  subtract peakuu[] in SAMPLE units, not ns!!
  //                  down- subtract peakud[]
  //  First FEE: down- subtract (peakud[]-peakdd[])
  //      Second FEE: up -  subtract peakud[]-peakdd[]+peakdu[]
  //                  down- subtract peakud[]
  // Bore paddle : A whole different animal, with lots of 
  
  for(auto &[endTag, endHit]: padEndHitMap) 
  {
    int paddleID    = endTag.getPaddleID();
    int paddleEndID = endTag.getPaddleEndID();
    
    //identify ramp type and modify times accordingly
    bool feeUp[4] = {true}; 
    int scID[4] = {0}; //based on SC's counting scheme identifying paddle ends (0-40)

    //loop over FEEs on this end, set ramptype
    for (int i=0; i<endHit.getFEEHitCount(); i++) 
    {
      HLX_Stg2_TOFFEEHit feeHit = endHit.fFEEHits[i];
      int FEEID  = feeHit.fTag.getFEEID();
      auto [laneIdx, chanIdx] = getGeom()->fTOF.fGeoMap[{paddleID, paddleEndID, FEEID}];

      if (paddleID<=15) {scID[i] = paddleID*2 + paddleEndID;}   // so simple now
      else if (paddleID==16) {scID[i] = 32 + paddleEndID*4 + FEEID;}        // BP, still very simple

      ///DEBUG::
      if (DebugMode > 3) {    
        cout << "-----------" << endl;
        cout << "DEBUG: NEW END FEE padID, endID, feeID, lane, chan: " << paddleID << ", " << paddleEndID << ", " << FEEID << ", " << laneIdx << ", " << chanIdx << endl; 
        cout << "DEBUG: i, scID, scID after: " << i << ", " << scID[i];
        cout << ", " << scID[i] << endl;
      }

      feeUp[i] = true;  //default up ramp
      if (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ) feeUp[i] = false; //down ramp
      if (DebugMode > 3) cout << "Setting feeUp= " << std::boolalpha << feeUp[i]  << std::boolalpha << endl;

      ///DEBUG::
      if (DebugMode > 3) {    
        if (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kLowToHigh ) cout << "UP UP UP UP UP UP UP UP UP UP ramp" << endl;
        if (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ) cout << "DOWN DOWN DOWN DOWN DOWN DOWN DOWN ramp" << endl;
      }
      int globalFeeID = laneIdx*8+chanIdx;
      fDBRType[globalFeeID] = 1; //up
      if (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ) fDBRType[globalFeeID] = -1; //down ramp

    }
    
    ///DEBUG::
    if (DebugMode > 3) {    
      cout << "DEBUG: feeUp[0-4]: " << std::boolalpha << feeUp[0]  << ", " << feeUp[1]  << ", " << feeUp[2]  << ", " << feeUp[3] << std::noboolalpha << endl;
      cout << "DEBUG: scID[0-4]: " << scID[0]  << ", " << scID[1]  << ", " << scID[2]  << ", " << scID[3] << endl;
    }

    // Now correct the times based on ramp order.
    for (int i=0; i<endHit.getFEEHitCount(); i++) 
    {
      HLX_Stg2_TOFFEEHit &feeHit = endHit.fFEEHits[i];
      int FEEID  = feeHit.fTag.getFEEID();
      
//      double beforeTime = feeHit.fTime;
      double corr = 0;
      
      if (FEEID == 0 && !feeUp[0]) corr = ((fCal.peakud[scID[1]] - fCal.peakdd[scID[1]])*fSampleToNanoSecond);
      if (FEEID == 1)
      {
        if ( feeUp[0] &&  feeUp[1]) corr = ((fCal.peakuu[scID[1]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1]) corr = ((fCal.peakud[scID[1]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1]) corr = ((fCal.peakud[scID[1]])*fSampleToNanoSecond);
      }
      if (FEEID == 2) //bore paddle only. It gets interesting here....
      {
        if ( feeUp[0] &&  feeUp[1] &&  feeUp[2]) corr = ((fCal.peakuu[scID[1]]+fCal.peakuu[scID[2]])*fSampleToNanoSecond);
        if ( feeUp[0] &&  feeUp[1] && !feeUp[2]) corr = ((fCal.peakuu[scID[1]]+fCal.peakud[scID[2]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1] &&  feeUp[2]) corr = ((fCal.peakud[scID[1]]+fCal.peakdu[scID[2]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1] && !feeUp[2]) corr = ((fCal.peakud[scID[1]]+fCal.peakdd[scID[2]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1] &&  feeUp[2]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]]+fCal.peakuu[scID[2]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1] && !feeUp[2]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]]+fCal.peakud[scID[2]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1] &&  feeUp[2]) corr = ((fCal.peakud[scID[1]]+fCal.peakdu[scID[2]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1] && !feeUp[2]) corr = ((fCal.peakud[scID[1]]+fCal.peakdd[scID[2]])*fSampleToNanoSecond);
      }
      if (FEEID == 3) //bore paddle only. It gets interesting here....
      {
        if ( feeUp[0] &&  feeUp[1] &&  feeUp[2] &&  feeUp[3]) corr = ((fCal.peakuu[scID[1]]+fCal.peakuu[scID[2]]+fCal.peakuu[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] &&  feeUp[1] &&  feeUp[2] && !feeUp[3]) corr = ((fCal.peakuu[scID[1]]+fCal.peakuu[scID[2]]+fCal.peakud[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] &&  feeUp[1] && !feeUp[2] &&  feeUp[3]) corr = ((fCal.peakuu[scID[1]]+fCal.peakud[scID[2]]+fCal.peakdu[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] &&  feeUp[1] && !feeUp[2] && !feeUp[3]) corr = ((fCal.peakuu[scID[1]]+fCal.peakud[scID[2]]+fCal.peakdd[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1] &&  feeUp[2] &&  feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdu[scID[2]]+fCal.peakuu[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1] &&  feeUp[2] && !feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdu[scID[2]]+fCal.peakud[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1] && !feeUp[2] &&  feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdd[scID[2]]+fCal.peakdu[scID[3]])*fSampleToNanoSecond);
        if ( feeUp[0] && !feeUp[1] && !feeUp[2] && !feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdd[scID[2]]+fCal.peakdd[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1] &&  feeUp[2] &&  feeUp[3]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]]+fCal.peakuu[scID[2]]+fCal.peakuu[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1] &&  feeUp[2] && !feeUp[3]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]]+fCal.peakuu[scID[2]]+fCal.peakud[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1] && !feeUp[2] &&  feeUp[3]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]]+fCal.peakud[scID[2]]+fCal.peakdu[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] &&  feeUp[1] && !feeUp[2] && !feeUp[3]) corr = ((fCal.peakud[scID[1]]-fCal.peakdd[scID[1]]+fCal.peakdu[scID[1]]+fCal.peakud[scID[2]]+fCal.peakdd[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1] &&  feeUp[2] &&  feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdu[scID[2]]+fCal.peakuu[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1] &&  feeUp[2] && !feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdu[scID[2]]+fCal.peakud[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1] && !feeUp[2] &&  feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdd[scID[2]]+fCal.peakdu[scID[3]])*fSampleToNanoSecond);
        if (!feeUp[0] && !feeUp[1] && !feeUp[2] && !feeUp[3]) corr = ((fCal.peakud[scID[1]]+fCal.peakdd[scID[2]]+fCal.peakdd[scID[3]])*fSampleToNanoSecond);
      }
      
      feeHit.fTime -= corr;
      
      ////DEBUG:
      auto [laneIdx, chanIdx] = getGeom()->fTOF.fGeoMap[{paddleID, paddleEndID, FEEID}];
      fDBRampdT[laneIdx*8+chanIdx] = corr;
      if (DebugMode > 1) {
        string ramp = (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kLowToHigh)? "UP": "DOWN";
//        cout << "scID[1 - 3]: " << scID[1] << "  " << scID[2] << "  " << scID[3] << endl;
//        std::streamsize original_precision = std::cout.precision();
        cout << paddleID << "  " << paddleEndID << "  " << FEEID << " " << ramp << " ramp, dTime: " << corr << std::fixed << std::setprecision(3) << "  " << feeHit.fTime << endl;
//        std::cout.precision(original_precision);
        // Important: Reset 'fixed' or 'scientific' if they were used
//        std::cout << std::defaultfloat;
//        cout << "Corrections scID[1]: peakuu, peakdd, peakud: " << fCal.peakuu[scID[1]] << ", " << fCal.peakdd[scID[1]] << ", " << fCal.peakud[scID[1]] << endl;
      }
    }
    
  }
  
    
  return kSuccess;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Assemble paddle hits from paddle-end hits and add to stage2 event
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::assemblePaddleHits(const map<HLX_Geometry_TOF_PaddleEndTag, HLX_Stg2_TOFPaddleEndHit>& padEndHitMap,
                                                     HLX_Stg2_Event* evtS2)
{
  map<HLX_Geometry_TOF_PaddleTag, HLX_Stg2_TOFHit> padHitMap;

  for(const auto &[tag, endHit]: padEndHitMap) {
    auto paddleTag = tag.getPaddleTag();
    if (padHitMap.count(paddleTag) == 0) {
      // This is the first hit for this paddle
      padHitMap[paddleTag] = HLX_Stg2_TOFHit{};
      padHitMap[paddleTag].fTag = paddleTag;
    }
    // Add the hit to the appropriate end
    if (tag.getPaddleEndID() == HLX_Geometry_TOF::kEast) {
      padHitMap[paddleTag].fHitEast = endHit;
    } else {
      padHitMap[paddleTag].fHitWest = endHit;
    }

    ///DEBUG:
    if (DebugMode > 1) {
      int paddleID = paddleTag.getPaddleID();
      int endID    = tag.getPaddleEndID();
      for (int j=0; j<endHit.getFEEHitCount(); j++)
      {
        HLX_Stg2_TOFFEEHit feeHit = endHit.getFEEHit(j);
        if (feeHit.fStatusFlag == HLX_Stg2_TOFFEEHit::kBad) continue;
        int FEEID = feeHit.fTag.getFEEID();
        auto [laneIdx, chanIdx] = getGeom()->fTOF.fGeoMap[{paddleID, endID,  FEEID}];
        double pretime = feeHit.fTime - feeHit.fTOFClock*fSampleToNanoSecond;
        double nobboff = pretime+fBboff[laneIdx];
//        cout << endl << endl << "Final: Paddle " << paddleID << ", End " << endID << ", FEE " << FEEID << ", SCFeeID " << fSCFeeIDMap[laneIdx*8+chanIdx] << endl; //tofcal11
        cout << endl << endl << "Final: Paddle " << paddleID << ", End " << endID << ", FEE " << FEEID << ", SCFeeID " << laneIdx*8+chanIdx << endl; //tofcal12
        if (tag.getPaddleEndID() == HLX_Geometry_TOF::kEast) cout << "EAST" << endl;  else  cout << "WEST" << endl; 
        cout << "time, w/UD, noBBoff: " << feeHit.fTime << "  " << pretime << "  " << nobboff << endl;
        cout << "X,Y,Z= " << feeHit.getX() << "  " << feeHit.getY() << "  " << feeHit.getZ() << endl;
//        cout << "   Signal ped, ramp ped: " << feeHit.fPulsePedestal << "   " << feeHit.fRampPedestal << endl;
//        cout << "   Signalsum: " << feeHit.fSignalSum << endl;
        int globFeeID = laneIdx*8+chanIdx;
        cout << "Ramptype, pulse ped, ramp ped : " << fDBRType[globFeeID] << ", " << feeHit.fPulsePedestal << ", " << feeHit.fRampPedestal << ", Signalsum " << feeHit.fSignalSum << endl;
        cout << "fittime, up/dn dT, : " << fDBFT[globFeeID] << ", " << fDBRampdT[globFeeID] << endl;
      }
    } //end DEBUG

  }

  // Finally, add the paddle hits to the S2 event
  for(auto &[tag, padHit]: padHitMap) {
    padHit.fTime          = padHit.getTime(fMaximumScaler*fSampleToNanoSecond);     //average time on end between FEEs
    padHit.fChargeDivFrac = padHit.getChargeDivFrac();                              //(E-W)/(E+W) for W,E = geometric mean signal of W,E end FEEs
    padHit.fTimeDiff      = padHit.getTimeDiff(fMaximumScaler*fSampleToNanoSecond); //E-W of FEE avg times --- put in SC paddle time shifts here
    
    double west = padHit.fHitWest.getY();
    double east = padHit.fHitEast.getY();
    double center = (west + east) / 2;
    double length = (west - east);
    
    int paddleID = padHit.fTag.getPaddleID(); //HLX numbering scheme


    //ChargeDivFrac position from fit to DCTY vs fChargeDivFrac
//    padHit.fChargeDivPos = fPolZ[paddleID].second * padHit.fChargeDivFrac + fPolZ[paddleID].first; //=DCTY - center
    //tofcal12:
    double Zdiff = padHit.fChargeDivFrac;
    padHit.fChargeDivPos = fCal.ZDm[paddleID] * Zdiff + fCal.ZDb[paddleID]; // = DCTY - center
    padHit.fChargeDivPos += center;

    //TimeDiff postion from fit to DCTY vs fTimeDiff
//    padHit.fTimeDiffPos = fPolT[paddleID].second * padHit.fTimeDiff + fPolT[paddleID].first; //=DCTY - center
    //tofcal12:
    padHit.fTimeDiffPos = fCal.TDm[paddleID] * padHit.fTimeDiff + fCal.TDb[paddleID]; // = DCTY - center
    padHit.fTimeDiffPos += center;

    evtS2->fTOFHits.push_back(padHit);
    
    ///DEBUG:
    if (DebugMode > 1) {
      cout << endl << "assemblePaddleHits: " << paddleID << "  " << west << "  " << east << "  " << length << "  " << center << endl;
      cout << "TimeDiffPos  : " << padHit.fTimeDiffPos << " = " << fCal.TDm[paddleID]  << " * " << padHit.fTimeDiff  << " + " << fCal.TDb[paddleID] << endl;
      cout << "fChargeDivPos: " << padHit.fChargeDivPos << " = " << fCal.ZDm[paddleID]  << " * " << Zdiff  << " + " << fCal.ZDb[paddleID] << endl << endl;
    }

  }
  
  return kSuccess;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Process individual FEE hit using SC algorithm
/// This is Carl's TOF analysis code, extracted from the TOFanalysis library developed at PSU,
/// implemented here by Dr. Kenichi Sakai.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::processTOFFEEHit_SC(HLX_Stg2_TOFFEEHit& feeHit)
{
  //[0] Pedestal calculation
  //[1] Pedestal subtraction
  //[2] Search Pulse
  //[3] Charge calculation
  //[4] Search Ramp
  //[5] Fit GV ramp
  //[6] Time calculation
  
  //What SC uses:
  //const int fNumFEE = 72;  // 32 FEEs for top, 8 FEEs for bore, 32 FEEs for bottom
  // FEE ordering in arrays here:  0-3 top0,  4-7 top1,  8-11 top2,  12-15 top3,  16-19 top4,  20-23 top5,
  // 24-27 top6,  28-31 top7,  32-39 bore0-bore7,  40-43 bot0,  44-47 bot1,  48-51 bot2,  52-55 bot3,  56-59 bot4,
  // 60-63 bot5,  64-67 bot6,  68-71 bot7
  // and eventually 40-71 bot0-31
  // For top/bottom: sequence of 4 consecutive array entries track the 0-y, 1-y, 0+y, 1+y FEEs on the same paddle
  // For bore, sequence is 0-y, 1-y, 2-y, 3-y, 0+y, 1+y, 2+y, 3+y
  //Paddles: 0-7 (top), 8-15 (bot), 16 (bore) from SC
  
  //what everyone else uses:
  int paddleID    = feeHit.fTag.getPaddleID();
  int FEEID       = feeHit.fTag.getFEEID();
  int paddleEndID = feeHit.fTag.getPaddleEndID();
  auto [laneIdx, chanIdx] = getGeom()->fTOF.fGeoMap[{paddleID, paddleEndID, FEEID}];
//  int SCfeeID = fSCFeeIDMap[laneIdx*8 + chanIdx]; //SC labels!! for use with current fCal only! tofcal11
  int SCfeeID = laneIdx*8 + chanIdx; //HLX labels. tofcal12

  //[0] Pedestal calculation: dynamic calculation or load pedestal file
//  auto [chanPedPulse, chanPedRamp] = calculatePedestal(feeHit);
  auto [chanPedPulse, chanPedRamp] = getPedSubSignalAndRamp(feeHit, SCfeeID);
  if (DebugMode > 1) {
    cout << endl << endl << "-----  NEW FEE pad, end, fee, SCfeeID: " << paddleID << ", " << paddleEndID << ", " << FEEID << ", " << SCfeeID << endl;
    cout << "Pedestals: pulse, ramp: " << chanPedPulse << "\t" << chanPedRamp << endl;
  }

  //following ped checks no longer meaningful when using SC version of calculatePedestalSC()
  if (chanPedPulse < 0 || chanPedRamp < 0 || feeHit.fStatusFlag == HLX_Stg2_TOFFEEHit::kBad) 
  {
    if (DebugMode > 2) {
      cout << "KILLED: [0] " << " chanPedPulse < 0 || chanPedRamp < 0 || feeHit.fStatusFlag == HLX_Stg2_TOFFEEHit::kBad" << endl;
      if (DebugMode > 4) {
        cout << "i, Signal, Ramp DATA:" << endl;
        for (int idx = 0; idx < feeHit.fRampData.size(); idx++) cout << idx << "     " << feeHit.fPulseData[idx] << "     " << feeHit.fRampData[idx] << endl;
      }
    }
    return -1;
  }

  if (DebugMode > 1) {
    if(feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kLowToHigh ) cout << "UP UP UP UP UP UP UP UP UP UP ramp" << endl;
    if(feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ) cout << "DOWN DOWN DOWN DOWN DOWN DOWN DOWN ramp" << endl;
  }

  //[1] Pedestal subtraction
  for (int idx = 0; idx < feeHit.fPulseData.size(); idx++)
    feeHit.fPulseSignal.push_back(feeHit.fPulseData[idx]-chanPedPulse);

  for (int idx = 0; idx < feeHit.fRampData.size(); idx++)
    feeHit.fRampSignal.push_back(feeHit.fRampData[idx]-chanPedRamp);

  if (DebugMode > 3) {
//  cout << "Signal pulse:" << endl;
//  for (int idx = 0; idx < feeHit.fPulseSignal.size(); idx++) cout << idx << "\t" << feeHit.fPulseSignal[idx] << endl;
    cout << "i, Signal, Ramp (ped sub):" << endl;
    for (int idx = 0; idx < feeHit.fRampSignal.size(); idx++) cout << idx << "     " << feeHit.fPulseSignal[idx] << "     " << feeHit.fRampSignal[idx] << endl;
  }
  
  //[2] Search Pulse - fixed end points in SC method
  int startPulseIdx{-1};
  int endPulseIdx{-1};
  //set threshold depending on channel id  
  int SCThresholdPulse = getPulseThreshold(feeHit.fTag);
  for (int idx = 0; idx < feeHit.fPulseSignal.size(); idx++)
    if (feeHit.fPulseSignal[idx] > SCThresholdPulse) {
      if (startPulseIdx == -1) startPulseIdx = idx;
      endPulseIdx = idx;
    }

  /// following cut removed after talking to Stephane 5/21/26 SN. SC wants to still calculate the timing pulse even if a signal pulse is not found.
//  if (startPulseIdx == -1 || (endPulseIdx-startPulseIdx) == 0) 
//  {
//    if (DebugMode > 2) 
//    {
//      cout << "KILLED: [2] startPulseIdx == -1 || (endPulseIdx-startPulseIdx) == 0" << endl;
//      cout << "  Threshold, start, end, end-start:" << SCThresholdPulse << "\t" << startPulseIdx << "\t"  << endPulseIdx << "\t" << (endPulseIdx-startPulseIdx) << endl;
//      cout << "i, Signal, Ramp (ped sub):" << endl;
//      for (int idx = 0; idx < feeHit.fRampSignal.size(); idx++) cout << idx << "     " << feeHit.fPulseSignal[idx] << "     " << feeHit.fRampSignal[idx] << endl;
//    }
//    return -1;
//  }

  //[3] Charge calculation
  feeHit.fSignalSumRaw = getSCPulseSum(feeHit.fPulseSignal, fPulseSumStartIdx, fPulseSumEndIdx, SCThresholdPulse);   //custom version to use SC thresholds for pulses
  
  //REMOVE Oct 11 2026 SLN -- handle in pass2 or skeletor using measured thicknesses
  //scale by scintillator thickness:
    //  double thick = fThickness[paddleID];
    //  feeHit.fSignalSum = feeHit.fSignalSumRaw/thick;
  //scale by channel gain from calibration file
  double chanGain = fCal.mpv[SCfeeID]/fCal.mpv[1];   //normalize to FEE1 (lane 0, ch 1: pad_id=16, end_id=1, fee_id=3) as of tofcal12.2.root
  feeHit.fSignalSum = feeHit.fSignalSumRaw/chanGain;
  feeHit.fSignalSum /= fCal.mippk;  //scale by global mippk to Z=1
//  feeHit.fSignalMax = HLX_Stg1_FeeEvent_TOFChan::getPulseMax(feeHit.fPulseSignal)/chanGain/thick/fCal.mippk;  //commented out because of 'thick' variable
  feeHit.fSignalMax = HLX_Stg1_FeeEvent_TOFChan::getPulseMax(feeHit.fPulseSignal)/chanGain/fCal.mippk;

  /// following cut removed after talking to Stephane 5/21/26 SN. SC wants to still calculate the timing pulse even if a signal pulse is not found.
//  if (feeHit.fSignalSum < 0 || feeHit.fSignalMax < 0) 
//  {
//    if (DebugMode > 2) {
//      cout << "KILLED: [3] feeHit.fSignalSum < 0 || feeHit.fSignalMax < 0" << endl;
//      cout << "            " << feeHit.fSignalSum << "     " << feeHit.fSignalMax << endl;
//    }
//    return -1;
//  }
  if ( feeHit.fSignalSum > 0 ) feeHit.fSignalSumErr = sqrt(feeHit.fSignalSum);

  //[4] Search Ramp to find starting location
  //Debug: remove following check (comes from KS)
//  int rampSlopeStep = feeHit.fRampSignal[startPulseIdx]-feeHit.fRampSignal[startPulseIdx-1];
//  if (fabs(rampSlopeStep) < fRampSlopeMinStep)  //check that the ramp is associated with signal pulse
//  {
  ///Debug:
////    cout << "KILLED: [4]" << endl;
//    return -1;
//  }

  feeHit.fStartRampIdx = -1;
  feeHit.fEndRampIdx = -1;

  //set upper limit to ramp search based on top vs bot/bore paddle:
  int rampSearchEndIdx = (paddleID >= 8 && paddleID <= 15)? fRampSearchEndIdxTop : fRampSearchEndIdxBotBore;
  //the following should search both up and down ramps, combined into one with a C++ style ternary operator:
  for (int idx = fRampSearchStartIdx; idx <= rampSearchEndIdx; idx++) {
    bool foundRamp = (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow)? fabs(feeHit.fRampSignal[idx]) < -fCal.pdn[SCfeeID][0]-fRampSearchMinVal : fabs(feeHit.fRampSignal[idx]) >= fRampSearchMinVal;
    if (foundRamp) {
      if (DebugMode > 4)
      {
        cout << "--Debug: Found ramp, idx = " << idx << endl;
        cout << "--idx, feeHit.fRampSignal[idx]: " << idx << ", " << feeHit.fRampSignal[idx] << endl;
      }
      feeHit.fStartRampIdx = idx;
      break;
    }
  }
  
  feeHit.fEndRampIdx = feeHit.fStartRampIdx + fRampFitLength;

  //check that reasonable end location was found
  if (feeHit.fEndRampIdx > fRampLength || feeHit.fStartRampIdx == -1) //no proper ramp found
  {
    if (DebugMode > 2) {
      cout << "KILLED: [4b] feeHit.fStartRampIdx == -1 || feeHit.fEndRampIdx > fRampLength" << endl;
      cout << "             " << feeHit.fStartRampIdx << "     " << feeHit.fEndRampIdx << endl;
    }
    return -1;
  }
  
  //[5] Fit ramp with GV template
    //init TGraph with ramp
    //init ramp fit
    //fit Tgraph
    //extract timing

  //init TGraph with ramp
  double sample[fRampFitLength], ramp[fRampFitLength];
  int knt = 0;
  for (int idx= feeHit.fStartRampIdx; idx < feeHit.fEndRampIdx; idx++)
  {
    sample[knt] = idx;
    ramp[knt]   = feeHit.fRampSignal[idx];
    knt++;
  }
  auto gRamp = TGraph(fRampFitLength, sample, ramp); //NB: out of scope outside this function
  if (DebugMode > 4) 
  {
    cout << "Graph: i, sample, Ramp (ped sub):" << endl;
    for (int idx = 0; idx < fRampFitLength; idx++) cout << idx << "     " << sample[idx] << "     " << ramp[idx] << endl;
  }

  //init ramp fit
  (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow) ? initGVdnRampFit(SCfeeID) : initGVupRampFit(SCfeeID);

  if (DebugMode > 1) {
    if(feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kLowToHigh ) cout << "Fitting UP ramp" << endl;
    if(feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ) cout << "Fitting DOWN ramp" << endl;
  }
  //fit Tgraph
  feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ? 
    gRamp.Fit("fGVdn", "QR", "", (double)feeHit.fStartRampIdx, (double)feeHit.fEndRampIdx) : 
    gRamp.Fit("fGVup", "QR", "", (double)feeHit.fStartRampIdx, (double)feeHit.fEndRampIdx); 
    
    //SC extracts the chisq here. I don't know that he uses it
  double chisq = (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow) ? fGVdn->GetChisquare() : fGVup->GetChisquare();

  //check that reasonable chi2 was found
  if (chisq <= 0 || log10(chisq) > 5) 
  {
    if (DebugMode > 2) {
      cout << "KILLED: [5] chisq > 0 && log10(chisq) > 5" << endl;
      cout << "     chi2   " << chisq << endl;
    }
    return -1;
  }
  

  //[6] Time calculation
  //extract timing
  double fittime = (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow) ? fGVdn->GetParameter(6) : fGVup->GetParameter(6);
  feeHit.fTime = (fittime+feeHit.fTOFClock)*fSampleToNanoSecond;     // SC has additional check: "check rare condition at the limit of a 16 bit integer" What is it?
  
  ///DEBUG:
  fDBFT[laneIdx*8+chanIdx] = fittime*fSampleToNanoSecond;
  if (DebugMode > 1) {
    TF1 *gv = (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow)? fGVdn : fGVup;
    cout << "Post fit: fit time, fTime " << fittime*fSampleToNanoSecond << "\t" << feeHit.fTime << endl;
    cout << "Fit parameters: par, name, min, val, max" << endl;
    for (int j = 0; j<gv->GetNpar(); j++)
    {
      double pmin, pmax;
      gv->GetParLimits(j,pmin,pmax);
      cout << j << "\t" << gv->GetParName(j) << "\t" << pmin << "\t" << gv->GetParameter(j) << "\t" << pmax << endl;
    }
  }
  

  feeHit.fTime -= fBboff[laneIdx]; //subtract board offset
  
  //following correction made by SC but not part of TOF calculation, so HLX_Stg2_TOFHit::getTimeDiff will be offset.
//  feeHit.fTime -= fCal.etoeshift[SCfeeID]; //subtract (West - East) offset  
  
  if (DebugMode > 1) {
    if(feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kLowToHigh ) cout << "SUCCESS: UP ramp  " << fittime*fSampleToNanoSecond << "\t" << std::fixed << std::setprecision(3) << feeHit.fTime << endl;
    if(feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kHighToLow ) cout << "SUCCESS: DOWN ramp " << fittime*fSampleToNanoSecond << "\t" << std::fixed << std::setprecision(3) << feeHit.fTime << endl;
  }

  return 0;
}

int HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getPulseThreshold(HLX_Geometry_TOF_FEETag tag)
{
  if (tag.getPaddleID() == 16) return 5;  //bore FEEs
  if (tag.getPaddleID() == 0 && tag.getPaddleEndID() == 1) return 10; //paddle 0 (top), end 1, both FEEs
  return 20; //all other FEEs
}

double HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getSCPulseSum(const vector<float> in, int startIdx, int endIdx, int summingThreshold)   //custom version to use SC thresholds for pulses
{
    if(startIdx >= endIdx) {
//      hlog::error("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getSCPulseSum => start index {} >= end index {}!",  startIdx, endIdx);
      return 0;
    }
    if(endIdx >= (int)in.size()){
//      hlog::error("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getSCPulseSum => end index {} out of range. Max: {}!",  endIdx, in.size());
      return 0;
    }

    double pulse{0}; //sum
    for(int n=startIdx;n<=endIdx;n++){
      if (in[n] >= summingThreshold ) pulse += (double)in[n];
    }
    return pulse;
}

// Open SCTOF parameter file (root TTree organized by run number)
bool HXA_Stg2_TOFEvent_AlgoTemplateAnaly::loadSCTOFCalibs(string ipfName) //only called once
{
  ////////  NB: put checks here to make sure file is not already open?  /////////  
  ////////  NB: put checks here to make sure correct version is opened?  /////////  

  if (!std::filesystem::exists(ipfName.c_str())) {
    hlog::error("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::loadSCTOFCalibs => file not found: {}", ipfName);
    exit(-1);
  }

  fFile = new TFile(ipfName.c_str());
  if (!fFile->IsOpen()) {
    hlog::error("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::loadSCTOFCalibs => could not open file: {}", ipfName);
    exit(-1);
  }

  hlog::info("HXA_Stg2_TOFEvent_AlgoTemplateAnaly::loadSCTOFCalibs : load from {}", ipfName);

  fTOFcal = (TTree*)fFile->Get("tofcal");

  ///////////////////
  //   VALID FOR tofcalv13 AND AFTER:
  ///////////////////
  // global ToF parameters, same for any run number
  fTOFcal->SetBranchAddress("pup",  fRtCal.pup);  // for each of 72 FEEs, up ramp parameters (7 of them) for GV function
  //  up ramp:   [0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fTOFcal->SetBranchAddress("epup", fRtCal.epup);  // for each of 72 FEEs, error on up ramp parameters (7 of them) for GV function
  fTOFcal->SetBranchAddress("pdn",  fRtCal.pdn);  // for each of 72 FEEs, down ramp parameters (7 of them) for GV function
  //  down ramp:   -[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fTOFcal->SetBranchAddress("epdn", fRtCal.epdn);  // for each of 72 FEEs, error on down ramp parameters (7 of them) for GV function
  fTOFcal->SetBranchAddress("mippk",&fRtCal.mippk);  // MIP peak, average most probable Landau fit value, used for determining the charge Z from the pulse integral between samples 11 and 25 (inclusive)
  fTOFcal->SetBranchAddress("interTW",&fRtCal.interTW);  // intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
  fTOFcal->SetBranchAddress("slopeTW",&fRtCal.slopeTW);  // slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
  fTOFcal->SetBranchAddress("einterTW",&fRtCal.einterTW);  // error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
  fTOFcal->SetBranchAddress("eslopeTW",&fRtCal.eslopeTW);  // error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)	
//  fTOFcal->SetBranchAddress("interT",&interT);  // top paddle intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
//  fTOFcal->SetBranchAddress("slopeT",&slopeT);  // top paddle slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
//  fTOFcal->SetBranchAddress("einterT",&einterT);  // top paddle error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
//  fTOFcal->SetBranchAddress("eslopeT",&eslopeT);  // top paddle error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)	
//  fTOFcal->SetBranchAddress("interB",&interB);  // bottom paddle intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
//  fTOFcal->SetBranchAddress("slopeB",&slopeB);  // bottom paddle slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
//  fTOFcal->SetBranchAddress("einterB",&einterB);  // bottom paddle error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
//  fTOFcal->SetBranchAddress("eslopeB",&eslopeB);  // bottom paddle error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)

  // ToF calibration parameters that change with run number

  fTOFcal->SetBranchAddress("runnum", &fRtCal.runnum);    // run number
  fTOFcal->SetBranchAddress("nument", &fRtCal.nument);    // number of entries in that run
  fTOFcal->SetBranchAddress("runtime",&fRtCal.runtime);   // time of the first event in the run, in days since 1st event in run 14753
  fTOFcal->SetBranchAddress("pedm",     fRtCal.pedm);    // for each of 72 FEEs, mean value of the charge pedestal for this run
  fTOFcal->SetBranchAddress("rms",      fRtCal.rms);     // for each of 72 FEEs, rms of the charge pedestal peak for this run
  fTOFcal->SetBranchAddress("mpv",      fRtCal.mpv);     // for each of 72 FEEs, most probable charge value of the pedestal-subtracted signal for this run
  fTOFcal->SetBranchAddress("empv",     fRtCal.empv);    // for each of 72 FEEs, uncertainty on the mpv value of the pedestal-subtracted charge signal for this run
  fTOFcal->SetBranchAddress("etoeshift",fRtCal.etoeshift);  // for each of 17 paddles, average time shift (in samples) between the two ends for this run (top: 0-7, bore: 8, bottom: 9-16)
  // for the following arrays of size 40: top are 0-15, bore are 16-23, bottom are 24-39;
  // top order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  // bore order is 00E, 01E, 12E, 23E, 00W, 01W, 12W, 23W   where 00E and 00W can be ignored 
  // bottom order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  fTOFcal->SetBranchAddress("peakuu", fRtCal.peakuu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+up
  fTOFcal->SetBranchAddress("widthuu",fRtCal.widthuu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+up
  fTOFcal->SetBranchAddress("peakud", fRtCal.peakud);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+dn
  fTOFcal->SetBranchAddress("widthud",fRtCal.widthud);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+dn
  fTOFcal->SetBranchAddress("peakdu", fRtCal.peakdu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+up
  fTOFcal->SetBranchAddress("widthdu",fRtCal.widthdu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+up
  fTOFcal->SetBranchAddress("peakdd", fRtCal.peakdd);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+dn
  fTOFcal->SetBranchAddress("widthdd",fRtCal.widthdd);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+dn
  
  //new branches from tofcal12:
  //Position from fChargeDivFrac by run:
  fTOFcal->SetBranchAddress("ZDm",  fRtCal.ZDm );  //paddle slope, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  fTOFcal->SetBranchAddress("ZDb",  fRtCal.ZDb );  //paddle intercept, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  //Position from fTimeDiff by run:
  fTOFcal->SetBranchAddress("TDm",  fRtCal.TDm );  //paddle slope, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b
  fTOFcal->SetBranchAddress("TDb",  fRtCal.TDb );  //paddle intercept, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b

  //better scheme for top and bottom paddle thickness measurements. NB: EVENTUALLY MOVE THIS TO SINGLE READ IN SEPARATE OBJECT IN tofcal root FILE, NOT IN TREE!
  fTOFcal->SetBranchAddress("padx", fRtCal.padx);  // x values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fTOFcal->SetBranchAddress("pady", fRtCal.pady);  // y values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fTOFcal->SetBranchAddress("padt", fRtCal.padt);  // thickness (mm) of top paddle at padx,pady  
  //keep bore stuff:
  fTOFcal->SetBranchAddress("borex", fRtCal.borex);  // x values where we have thickness measurements for bore paddle
  fTOFcal->SetBranchAddress("borey", fRtCal.borey);  // y values where we have thickness measurements for bore paddle
  fTOFcal->SetBranchAddress("boret", fRtCal.boret);  // thickness (mm) of bore paddle at topx,topy
  
  //make copy in case of mods:
  fTOFcal->GetEntry( 0 ); // load up the calibration parameters structure for first run
  fCal = fRtCal;

  getSCTOFCalibsForRun( 14765 ); //set default values
  
  return kTRUE;

}

////////////////////////////////////////////////////////
void HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getSCTOFCalibsForRun(int runnum, bool doAdjustment) //Getting fCal for run runnum
{
  static int runcurrent = 0;
//  cout << "getSCTOFCalibsForRun: Getting fCal for run " << runnum << ", runcurrent =  " << runcurrent << endl;
  if (runnum == runcurrent) return;

  //copy previous values
  PSUTOFCalibs oldCal = fCal;

//  cout << "getSCTOFCalibsForRun: fTOFcal->GetEntries() " << fTOFcal->GetEntries() << endl;

  // Now that we know which run we're dealing with, scroll down fTOFcal.root to find the relevant calibration parameters
  for (int ical = 0; ical < fTOFcal->GetEntries(); ical++) 
  {	
    fTOFcal->GetEntry( ical ); // load up the calibration parameters structure for this run
//    cout << "getSCTOFCalibsForRun: Checking entry: " << ical << ", fRtCal.runnum: " << fRtCal.runnum << ", runnum: " << runnum << endl;
    if (fRtCal.runnum >= runnum) {
//      cout << "getSCTOFCalibsForRun: Setting fTOFcal entry ical " << ical << " with fRtCal.runnum = " << fRtCal.runnum << endl;
      runcurrent = fRtCal.runnum;
      fCal = fRtCal;  //copy INPUT to OUTPUT
      break; // now we have the parameters for this run, or at least for the next higher run number in the database
    }
  }
  
//  cout << "getSCTOFCalibsForRun: found entry " << runcurrent << endl;

  ///DEBUG:
  if (DebugMode > 3) {
    cout << "Calibration parameters, run " << runnum << endl;
    cout << fCal.runnum << "  " << fCal.peakuu[20] << "  " << fCal.widthuu[14] << endl;
    cout << fCal.padx[7][3] << "  " << fCal.borex[0][5] << "  " << fCal.padt[10][3][12] << endl;
    for (int i=0; i<17; i++) cout << fCal.ZDm[i] << "  " << fCal.ZDb[i] << endl;
  }

//adjusting particular values inside fCal for run
  if (doAdjustment) adjustSCTOFCalibsForRun( oldCal ); //adjusting particular values inside fCal for run
  
}

///////////////////////////////////////////////////////////////
//do adjustment of bad parameters inside fCal for run
void HXA_Stg2_TOFEvent_AlgoTemplateAnaly::adjustSCTOFCalibsForRun(PSUTOFCalibs oldCal) 
{

  ////////  NOTE: IS THE FOLLOWING CHECK STILL NEEDED? /////////////
  for( Int_t j = 0; j < 40; j++) { // check for pathological peak?? corrections, and reset to last good known if needed
    if (fCal.widthuu[j]<0 || fCal.widthuu[j]>0.04) {fCal.peakuu[j] = oldCal.peakuu[j];}
    if (fCal.widthud[j]<0 || fCal.widthud[j]>0.04) {fCal.peakud[j] = oldCal.peakud[j];}
    if (fCal.widthdu[j]<0 || fCal.widthdu[j]>0.04) {fCal.peakdu[j] = oldCal.peakdu[j];}
    if (fCal.widthdd[j]<0 || fCal.widthdd[j]>0.04) {fCal.peakdd[j] = oldCal.peakdd[j];}
  }

  // Check mpv (fee by fee fSignalSumRaw peak value after PLC and thickness corrections)
  for( Int_t j = 0; j < 72; j++ ) {
    if ( fCal.empv[j] <= 0 || fCal.mpv[j] <= 0 || fCal.empv[j]/fCal.mpv[j] > 0.02 ) { // Need Landau peak determined to better than 2% accuracy. Assumes previous values were good. 
      fCal.mpv[j]  = oldCal.mpv[j];
      fCal.empv[j] = oldCal.empv[j];
    }
  }

  //do iffy check for reasonable fit values in ZDm, TDm etc
  for (int i=0; i<17; i++)
  {
    if (fCal.ZDm[i] == -10) {fCal.ZDm[i] = oldCal.ZDm[i]; fCal.ZDb[i] = oldCal.ZDb[i];}
    if (fCal.TDm[i] == -10) {fCal.TDm[i] = oldCal.TDm[i]; fCal.TDb[i] = oldCal.TDb[i];}
  }

}

void HXA_Stg2_TOFEvent_AlgoTemplateAnaly::defineGVRampFitFuncs()
{
  
//  fGRamp = new TGraph(fRampLength);
  //Debug:
//  cout << "HXA_Stg2_TOFEvent_AlgoTemplateAnaly::defineGVRampFitFuncs: Init" << endl;
  
  // GV up ramp function to fit the timing ramp. Note fixed start and end
  fGVup = new TF1("fGVup","[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))", 9., 22.); //SC original
//  fGVup = new TF1("fGVup","[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))", 4., 22.); //SN expanded region
  // Set parameter names (optional)
  fGVup->SetParName(0,"a");
  fGVup->SetParName(1,"b");
  fGVup->SetParName(2,"c");
  fGVup->SetParName(3,"d");
  fGVup->SetParName(4,"e");
  fGVup->SetParName(5,"f");
  fGVup->SetParName(6,"dx");
  // Set parameter start values (mandatory).
  fGVup->SetParLimits(0, 13500., 16000.);
  fGVup->SetParLimits(1, 0.5,4.);
  fGVup->SetParLimits(2, 0.,2.);
  fGVup->SetParLimits(3, 0.,2.);
  fGVup->SetParLimits(4, 0.,2.);
  fGVup->SetParLimits(5, 0.,2.);
  fGVup->SetParLimits(6, 9.,15.);
  
  // GV down ramp function to fit the timing ramp
  fGVdn = new TF1("fGVdn","-[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))", 9., 22.); //SC original
//  fGVdn = new TF1("fGVdn","-[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))", 4., 22.); //SN expanded region
  // Set parameter names (optional)
  fGVdn->SetParName(0,"a");
  fGVdn->SetParName(1,"b");
  fGVdn->SetParName(2,"c");
  fGVdn->SetParName(3,"d");
  fGVdn->SetParName(4,"e");
  fGVdn->SetParName(5,"f");
  fGVdn->SetParName(6,"dx");
  // Set parameter start values (mandatory).
  fGVdn->SetParLimits(0, -16000., -13500.);
  fGVdn->SetParLimits(1, 0.5,4.);
  fGVdn->SetParLimits(2, 0.,2.);
  fGVdn->SetParLimits(3, 0.,2.);
  fGVdn->SetParLimits(4, 0.,2.);
  fGVdn->SetParLimits(5, 0.,2.);
  fGVdn->SetParLimits(6, 9.,15.);
}

void HXA_Stg2_TOFEvent_AlgoTemplateAnaly::initGVupRampFit(int FEEID)
{
  //Debug:
  if (DebugMode > 3) cout << "initGVupRampFit 0: " << fCal.pup[FEEID][0] << endl;
  fGVup->FixParameter(0,fCal.pup[FEEID][0]);
  if (DebugMode > 3) cout << "initGVupRampFit 1: " << fCal.pup[FEEID][1] << endl;
  fGVup->FixParameter(1,fCal.pup[FEEID][1]);
  if (DebugMode > 3) cout << "initGVupRampFit 2: " << fCal.pup[FEEID][2] << endl;
  fGVup->FixParameter(2,fCal.pup[FEEID][2]);
  if (DebugMode > 3) cout << "initGVupRampFit 3: " << fCal.pup[FEEID][3] << endl;
  fGVup->FixParameter(3,fCal.pup[FEEID][3]);
  if (DebugMode > 3) cout << "initGVupRampFit 4: " << fCal.pup[FEEID][4] << endl;
  fGVup->FixParameter(4,fCal.pup[FEEID][4]);
  if (DebugMode > 3) cout << "initGVupRampFit 5: " << fCal.pup[FEEID][5] << endl;
  fGVup->FixParameter(5,fCal.pup[FEEID][5]);
}
  
void HXA_Stg2_TOFEvent_AlgoTemplateAnaly::initGVdnRampFit(int FEEID)
{
  //Debug:
  if (DebugMode > 3) cout << "initGVdnRampFit 0: " << fCal.pdn[FEEID][0] << endl;
  fGVdn->FixParameter(0,fCal.pdn[FEEID][0]);
  if (DebugMode > 3) cout << "initGVdnRampFit 1: " << fCal.pdn[FEEID][1] <<  endl;
  fGVdn->FixParameter(1,fCal.pdn[FEEID][1]);
  if (DebugMode > 3) cout << "initGVdnRampFit 2: " << fCal.pdn[FEEID][2] << endl;
  fGVdn->FixParameter(2,fCal.pdn[FEEID][2]);
  if (DebugMode > 3) cout << "initGVdnRampFit 3: " << fCal.pdn[FEEID][3] << endl;
  fGVdn->FixParameter(3,fCal.pdn[FEEID][3]);
  if (DebugMode > 3) cout << "initGVdnRampFit 4: " << fCal.pdn[FEEID][4] << endl;
  fGVdn->FixParameter(4,fCal.pdn[FEEID][4]);
  if (DebugMode > 3) cout << "initGVdnRampFit 5: " << fCal.pdn[FEEID][5] << endl;
  fGVdn->FixParameter(5,fCal.pdn[FEEID][5]);
}

pair<double, double> HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getPedSubSignalAndRamp(HLX_Stg2_TOFFEEHit& feeHit, int FEEID)
{
  if (fProcessMode == eProcessingMode::kProcess) {
    double chanPedPulse           = getSignalPedestal(feeHit.fPulseData, FEEID);
    auto [chanPedRamp, rampType]  = getRampPedestal(feeHit.fRampData); //rampType = -1: down, 0= no ramp; +1= up ramp

    if      (rampType == 1)  feeHit.fRampType = HLX_Stg2_TOFFEEHit::RampType::kLowToHigh;
    else if (rampType == -1) feeHit.fRampType = HLX_Stg2_TOFFEEHit::RampType::kHighToLow;
    else                     feeHit.fStatusFlag = HLX_Stg2_TOFFEEHit::kBad;
    
    feeHit.fTime = 0;
    if (feeHit.fStatusFlag == HLX_Stg2_TOFFEEHit::kBad) return {-1, -1};

    if (!fSubtractPeds) {
      chanPedPulse = 0;
      chanPedRamp  = 0;
    } else if (!fDynamicPeds) {
      chanPedPulse = fChanPedPulseMap[feeHit.fTag];
      chanPedRamp  = (feeHit.fRampType == HLX_Stg2_TOFFEEHit::RampType::kLowToHigh)? fChanPedRampLowMap[feeHit.fTag] : fChanPedRampHighMap[feeHit.fTag];
    }
    feeHit.fPulsePedestal = chanPedPulse;
    feeHit.fRampPedestal  = chanPedRamp;
  }

  return {feeHit.fPulsePedestal, feeHit.fRampPedestal};
}

double HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getSignalPedestal(const vector<int16_t> vec, int feeID)
{

  //avg of first fPedSamples channels:
  double ped = 0;
  for (int k = 0; k < fPedSamples; k++) {
      ped += vec[k];
  }
  
  //look for weirdness in channels from len-2 to end
  double pedlo = ped;
  for (Int_t k = fPedSamples; k < fRampLength; k++) {    // determine the running sum for the pedestal, keeping the minimum value
    pedlo = pedlo - vec[k-fPedSamples] + vec[k];
    ped = min(ped,pedlo);
  }
  ped /= fPedSamples;

  if (abs(ped-fCal.pedm[feeID])>2.*fCal.rms[feeID]) ped=fCal.pedm[feeID];    // guard against a pathological pedestal value; if more than 2.*sigma away from mean pedestal, use mean pedestal

  return ped;
}

pair<double, int> HXA_Stg2_TOFEvent_AlgoTemplateAnaly::getRampPedestal(const vector<uint16_t> vec)
{

  //find up or down configuration:
  double pedlo = 0, pedhi = 0;
  int kstart  = 1;
  int kend    = fPedSamples-1;
  int nk      = kend - kstart;
  
  for (int k = kstart; k < kend; k++) {
    pedlo += vec[kend+k];  // use samples 7 to 10 for the low part
    pedhi += vec[21-k];  // use samples 17 to 20 for the high part
  }
  
  //Use pedlo and pedhi to determine ramp direction. Check that they are at least a certain ADC distance apart.
  int upRamp = 0; //failure
  if (pedhi - pedlo >  fRampSlopeMinStep) upRamp = 1;   // up ramp
  if (pedhi - pedlo < -fRampSlopeMinStep) upRamp = -1;  // down ramp
  
  if (upRamp == 0 ) return {-1, 0};

// Need to fix the pedestal determination; some rare events don't have a real pedestal in the first or last numpedsamp,
// but instead somewhere in between (e.g., starts high, has a down ramp, then some weird late trace distortion)

  //sneak up on true pedestal through successive iterations:
  const int nIter = 5;
  int pedIter[nIter] = {1500, 100, 40, 20, 10};
  double ped = 0.;
  
  for (int itr=0; itr<nIter; itr++)
  {
    int numave = 0;
    double ave = 0.;
    for (int i=1; i < fRampLength; i++) {
      if (abs(vec[i]-ped) < pedIter[itr]) { // look at the low average baseline, for which the raw adcr values are usually around 400 or 500; looking at numbers less than 1500 should be at the bottom range of things
        numave++;
        ave += vec[i];
      }
    }
    if (numave > 0) 
    {
      ped = ave/numave;     // this low average should be pretty close to the pedestal, modulo the squiggle artifact at the bottom of the ramp
    } 
    else return {-1, upRamp};
  }
  
  return {ped, upRamp};
}

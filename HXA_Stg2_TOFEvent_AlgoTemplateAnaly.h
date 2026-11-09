////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Filename: HXA_Stg2_TOFEvent_AlgoTemplateAnaly.h
/// Created:  2026/04/10 by SLN (modified Sakai code)
/// Created:  2025/08/05 by sakai (refactored from original TOF hitbuilder)
/// PSU/Carl's original TOF processing algorithm implementation
/// but updated with Coutu's time dependent calibrations and ramp fit using Gerard Visser formula.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <stdexcept>
#include "HXA_TimeSeries.h"
#include "HXA_Stg2_TOFEvent_AlgoBase.h"

#include <TTree.h>
#include <TMath.h>
#include <TF1.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// PSU TOF processing algorithm - this contains Carl's original TOF analysis code
/// extracted from the TOFanalysis library developed at PSU, implemented by Dr. Kenichi Sakai.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HXA_Stg2_TOFEvent_AlgoTemplateAnaly : public HXA_Stg2_TOFEvent_AlgoBase {
  public:
    HXA_Stg2_TOFEvent_AlgoTemplateAnaly() { fAlgoType = eEventProcessingAlgorithm::kTemplateAnaly; }
    ~HXA_Stg2_TOFEvent_AlgoTemplateAnaly() override = default;

  public:
    /// Processing parameters
    const int    fPedSamples              {6};    //number of samples to use in calculating the pedestals
    const int    fPulseSumStartIdx        {11};
    const int    fPulseSumEndIdx          {25};
    const int    fRampSearchStartIdx      {9};    //limited range of allowed start ramps in SC world
    const int    fRampSearchEndIdxTop     {13};
    const int    fRampSearchEndIdxBotBore {14};
    const int    fRampSearchMinVal        {500};
    const int    fRampFitLength           {8};
    const int    fRampLength              {32};
    const int    fRampSlopeMinStep        {2000};
    const double fSampleToNanoSecond      {25};
    const double fMaximumScaler           {65536};

  public:
    //fee    : HLX_Geometry_TOF_FEETag
    //paddle : HLX_Geometry_TOF_PaddleTag

    static const int fNumFEE = 72;
    //Below commented out 11 Oct 2026 SLN -- detailed thickness maps now available in post-processing. Add to pass2 when available here.
    // scintillator thickness, average for now (will have to use detailed maps later when tracking is available) by paddleID. 0.956 means scintillator is 0.956cm thick.
//    Double_t fThickness[17] = {0.956,0.967,0.961,0.988,0.984,0.956,0.997,0.962,0.979,0.964,1.019,0.956,1.009,0.972,0.984,0.987,0.997};
    //Paddles: 0-7 (top), 8-15 (bot), 16 (bore) from SC
    Double_t fBboff[9] = {4.2285, -7.8465, -7.8365, -7.6625, -7.9085, 0, 0.0625, -0.0505, -0.0395};    // Timing offsets by board, in ns
    

  public:
    /// Configure the algorithm with pedestal file and steering parameters
    /// Mode: Process
    int configure(const string& pedFileName, const HLX_TOML_Table& steeringInfo) override;
    /// Mode: Reprocess
    int configure() override {
      fProcessMode = eProcessingMode::kReprocess;
      throw std::runtime_error(
        "HXA_Stg2_TOFEvent_AlgoTemplateAnaly::configure() is not implemented yet"
      );
      //return 0; //Once implemented, return 0 for success or -1 for error
    }

    /// Process a stage1 event to produce TOF hits in the stage2 event
    int processEvent(HLX_Stg1_Event* evtS1, HLX_Stg2_Event* evtS2) override;

 private:
   
   ///DEBUGGING:
   int DebugMode = 0; //higher means more. 5=all debug comments, 0=no comments
   double fDBFT[fNumFEE]={0};
   double fDBRampdT[fNumFEE]={0};
   int fDBRType[fNumFEE]={0};
   int fEvtCnt = 0;
   
   void clearDebug()
   {
     for (int ij=0; ij<fNumFEE; ij++)
     {
       fDBFT[ij] = 0;
       fDBRampdT[ij] = 0;
       fDBRType[ij] = 0;
     }
   }
   /*
   /// insane translation array from lane/chan to SC global FEEID:
   const int fSCFeeIDMap[fNumFEE] = {
     38,39,36,37,35,34,33,32, //[0,0] -> [0,7]
     0,1,4,5,13,12,9,8,       //[1,0]
     10,11,14,15,7,6,3,2,     //[2,0]
     31,30,27,26,18,19,22,23, //[3,0]
     21,20,17,16,24,25,28,29, //[4,0]
     44,45,48,49,53,52,40,41, //[5,0]
     54,55,50,51,42,43,47,46, //[6,0]
     66,67,70,71,63,62,59,58, //[7,0]
     61,60,57,56,64,65,68,69  //[8,0]
   };
   */
   //paddle numbering:
   // HLX   SC
   //  0    9 bottom
   //  1   10
   //  2   11
   //  3   12
   //  4   13
   //  5   14
   //  6   15
   //  7   16
   //  8    0 top
   //  9    1
   //  10   2
   //  11   3
   //  12   4
   //  13   5
   //  14   6
   //  15   7
   //  16   8 bore
   //test slope and intercept to compute paddle tDiff
   //Run 14801 GV: paddles 0-16 (inc. bore)
   // position = slope*var + intercept
    
    //CALIBRATION TREE VARIABLES: Saved in structure for eases of copying
    struct PSUTOFCalibs {
      Double_t pup[72][7],epup[72][7],pdn[72][7],epdn[72][7],mippk;
      Int_t runnum,nument;
      Double_t runtime,pedm[72],rms[72],mpv[72],empv[72];
      Double_t etoeshift[17];
      Double_t peakuu[40], widthuu[40], peakud[40], widthud[40], peakdu[40], widthdu[40], peakdd[40], widthdd[40];
//      Double_t interT[64],einterT[64],slopeT[64],eslopeT[64],interB[64],einterB[64],slopeB[64],eslopeB[64];
      Double_t interTW[64],einterTW[64],slopeTW[64],eslopeTW[64];
//      Int_t padnum[17];
//      Double_t topx[8][7],topy[8][17],topt[8][7][17],botx[8][7],boty[8][17],bott[8][7][17]; //,borex[1][11],borey[1][11],boret[1][11][11];
  
      //new branches: tofcal12
      Double_t ZDm[17], ZDb[17];  //paddle slope, intercept, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
      Double_t TDm[17], TDb[17];  //paddle slope, intercept, from DCTY vs fTimeDiff:      fTimeDiffPos      = m*fTimeDiff + b
  
      //tofcal13 with HELIX numbering scheme:
      Double_t padx[16][7], pady[16][17], padt[16][7][17]; //combine top and bottom paddles
      Double_t borex[1][11],borey[1][11],boret[1][11][11]; //bore paddle as before
    };
    
    
    TFile*        fFile;      //TFile with Calibs
    TTree*        fTOFcal = {};    //calibration parameters
    PSUTOFCalibs  fCal;       //global struct for use in calib getting
    PSUTOFCalibs  fRtCal;     //global struct with variables that ROOT points to

    /// Process individual FEE hit using SC algorithm
    int processTOFFEEHit_SC(HLX_Stg2_TOFFEEHit& feeHit);

    /// Assembly methods for building higher-level hits
    int     assembleFEEHits(HLX_Stg1_Event* evtS1, HLX_Stg2_Event* evtS2, vector<HLX_Stg2_TOFFEEHit>& feeHitVec);
    int     assemblePaddleEndHits(const vector<HLX_Stg2_TOFFEEHit>& feeHitVec, map<HLX_Geometry_TOF_PaddleEndTag, HLX_Stg2_TOFPaddleEndHit>& padEndHitMap);
    int     assemblePaddleHits(const map<HLX_Geometry_TOF_PaddleEndTag, HLX_Stg2_TOFPaddleEndHit>& padEndHitMap, HLX_Stg2_Event* evtS2);
    
    //local funcs:
    int     getPulseThreshold(HLX_Geometry_TOF_FEETag tag);   //get SC pulse summing thresholds
    double  getSCPulseSum(const vector<float> in, int startIdx, int endIdx, int summingThreshold);   //custom version to use SC thresholds for pulses
    void    defineGVRampFitFuncs(); //define GV template ramp fitting functions
    void    initGVupRampFit(int feeID);   //initialize GV fit template for up ramps
    void    initGVdnRampFit(int feeID);   //initialize GV fit template for down ramps 

    // TOF calibration parameters from TTree
    inline static const string kFileName {"config/calibration_files/caldata_tofhit_AlgoTemplateAnaly.root"};
    bool                  loadSCTOFCalibs(string ipfName = kFileName);                  //only called once
    void                  getSCTOFCalibsForRun(int runnum, bool doAdjustment = true);   //called for each event to update if new run encountered    
    void                  adjustSCTOFCalibsForRun(PSUTOFCalibs oldCal);                //do adjustment of bad parameters inside fCal for run
    pair<double, double>  getPedSubSignalAndRamp(HLX_Stg2_TOFFEEHit& feeHit, int feeID);// Calculate pedestal
  
    /// Pedestal processing for this technique
    double                getSignalPedestal(const vector<int16_t> vec, int feeID);
    pair<double, int>     getRampPedestal(const vector<uint16_t> vec);
    
    
//    TGraph*       fGRamp;     //ramp TGraph for fitting
    TF1*          fGVdn;      //GV fitting template, down ramp
    TF1*          fGVup = {};      //GV fitting template, up ramp
    
  };


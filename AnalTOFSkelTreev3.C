//void AnalTOFSkel(TString infile) //infile is tof output created by skeletor
//infile: skeletor.15550.root
//
// SC 10/12/2026   repackage AnalTOFSkel.C into AnalTOFSkelTree.C which will be invoked within AnalTOFSkelChain.C to
// process multiple skeletor files (in a TChain) with a list of input trees read from a text file

////////////////////////////////////////////////

//include PSU tof calib tools:
#include "PSUTofTools12.C"
/////  SET THIS IN PSUTofTools12.C!! ////
//Set up input tofcal ttree:
//string kInputTofcalName = "tofcal12.1.root";

//include PSU thickness tools for paddle thicknesses
#include "PSUThicknessTools.C"

//include signal response correction functions
#include "ResponseCorrTools.C"

//include signal calculation for arithmetic and geometric means of paddles and ends
#include "GetSignal.C"

void AnalTOFSkelTree(TTree *toft, const TString &outfile, bool recalcChargeDivPos = true, bool recalcTimeDiffPos = true )
{
  /// Get the geometry singleton - it will autoload the files.
  auto theGeom = HLX_Geometry::getInstance();
  
  ////////////////////////////////////////
  /// INPUT skeletor TOF info
//  vector<HLX_Stg2_TOFHit>* TOFHits = nullptr; //top, bore, bottom
  HLX_Stg2_TOFHit* topTOF = new HLX_Stg2_TOFHit();
  HLX_Stg2_TOFHit* borTOF = new HLX_Stg2_TOFHit();
  HLX_Stg2_TOFHit* botTOF = new HLX_Stg2_TOFHit(); //top, bore, bottom
//Hits:
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

//  toft->SetBranchAddress("TOFHits", &TOFHits);
  toft->SetBranchAddress("topTOF", &topTOF);
  toft->SetBranchAddress("boreTOF", &borTOF);
  toft->SetBranchAddress("botTOF", &botTOF);
  toft->SetBranchAddress("tvars",&tvars);
  toft->SetBranchAddress("runnum",&runnum);
  toft->SetBranchAddress("evttime",&evttime);
  toft->SetBranchAddress("richZEst",&richZEst);
  toft->SetBranchAddress("DCTChi2",&DCTChi2);
  toft->SetBranchAddress("rig",&rig);
  toft->SetBranchAddress("mdr",&mdr);
  toft->SetBranchAddress("hitCount",&hitCount);
  toft->SetBranchAddress("topdircos", topdircos);
  toft->SetBranchAddress("boredircos",boredircos);
  toft->SetBranchAddress("botdircos", botdircos);
  
  
///////////////////////////////////////////////
  // Load INPUT tofcal calibration data
  loadSCTOFCalibsv12();
  
  //set up OUTPUT tofcal calibration  
//  createSCTOFCalibs();
  
  //load response correction parameters:
  //for sqrt(signal) based corrections
//  string kInputResponseCorrFileName = "combined.Bon.sqrt.PadSigTGraphFits.txt";
//  if (outfile.Contains("Boff")) kInputResponseCorrFileName = "combined.Boff.sqrt.PadSigTGraphFits.txt";
  //for sqrt(signal) based corrections
  string kInputResponseCorrFileName = "combined.Bon.PadSigTGraphFits.txt";
  if (outfile.Contains("Boff")) kInputResponseCorrFileName = "combined.Boff.PadSigTGraphFits.txt";
  readResponseCorr(kInputResponseCorrFileName);
  
  ///////////////////////////////////////////
  //output root file for histograms
  TFile *opf = TFile::Open(outfile, "RECREATE");

  if(!opf || !opf->IsOpen() || opf->IsZombie())
  {
    cout << "Could not open file " << outfile << endl;
    exit(0);
  }
  
  //////////////////////////////////////////
  //BOOK SOME HISTOGRAMS
  //////////////////////////////////////////
  //These 3D histograms don't seem to work as TProfiles, which I prefer to fit over fitting the 2D histo.
//  TH3D* hZdct = new TH3D("hZdct","padID vs DCTY vs fChargeDivFrac;paddle;fChargeDivFrac;dctY (mm)",17,0,17, 50,-1,1, 200,-1000,1000);
//  TH3D* hTdct = new TH3D("hTdct","padID vs DCTY vs fTimeDiff;paddle;fTimeDiff;dctY (mm)",17,0,17, 50,-12,12, 200,-1000,1000);

  //pos from timing, chargeDivFrac, and difference
  TH2D** hChargeDivPos;                 //position from charge division 
  hChargeDivPos     = new TH2D*[17];
  TH2D** hSignal;                       // Charge (getSignal() = geometric mean of end fees and geometric mean of ends) vs DCTY for all paddles
  hSignal           = new TH2D*[17];
  TH2D** hSignalA;                       // Charge (arithmetic mean of end fees and arith mean of ends) vs DCTY for all paddles
  hSignalA           = new TH2D*[17];
  TH2D** hSignalG;                       // Charge (arithmetic mean of end fees and geometric mean of ends) vs DCTY for all paddles
  hSignalG           = new TH2D*[17];
  TH2D** hChg;                       // Charge (sqrt(hSignal)) vs DCTY for all paddles
  hChg           = new TH2D*[17];
  TH2D** hChgA;                       // Charge (sqrt(hSignalA)) vs DCTY for all paddles
  hChgA           = new TH2D*[17];
  TH2D** hChgG;                       // Charge (sqrt(hSignalG)) vs DCTY for all paddles
  hChgG           = new TH2D*[17];
  TH2D** hTimeDiffPos;                  //position from end to end time difference 
  hTimeDiffPos      = new TH2D*[17];
  TH2D** hPosDiff;                      //comparison of position from charge division and time diff
  hPosDiff          = new TH2D*[17];
  TH2D** hXY;                           // SC  check (X,Y) for each paddle
  hXY               = new TH2D*[17];
  TH2D** hXYshift;                      // SC  check (X,Y) for each paddle; (X,Y) shifted to a paddle position (0-20 cm, 0-180 cm)
  hXYshift          = new TH2D*[17];
  TH2D** hYXshift;                      // SC  check (Y,X) for each paddle; (Y,X) shifted to a paddle position (0-20 cm, 0-180 cm);  YX is interesting to look at instead of XY to look for possible paddle rotations
  hYXshift          = new TH2D*[17];
  TProfile** pYXshift;
  pYXshift = new TProfile*[17];
  TH2D** hXYChi;                        // SC  check (X,Y) for each paddle, after a DCTChi2 cut
  hXYChi            = new TH2D*[17];
  TH2D** hXYHitC;                       // SC  check (X,Y) for each paddle, after a DCT hitCount cut
  hXYHitC           = new TH2D*[17];
  TH1D** hthck;                         // SC  thckness at (X,Y) for each paddle
  hthck             = new TH1D*[17];
  TH2D** hZdT;                          //fChargeDivPos vs fTimeDiffPos timing differences
  hZdT = new TH2D*[17];
  TProfile** pZdT;
  pZdT = new TProfile*[17];
  TH2D** hZdct;                         //position differences: dcty vs fChargeDivFrac
  hZdct = new TH2D*[17];
  TProfile** pZdct;
  pZdct = new TProfile*[17];
  TH2D** hTdct;                         //position differences: fTimeDiff vs dcty
  hTdct = new TH2D*[17];
  TProfile** pTdct;
  pTdct = new TProfile*[17];
  

  for (int i=0;i<17;i++)
  {
    TString hname   = TString::Format("%s%d","hChargeDivPos",i);
    TString htitle  = TString::Format("%s%d%s","fChargeDivPos vs DCTY, paddle",i,";DCTY (mm);Pos (mm)");
    if (i<16) hChargeDivPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-1000,1000, 200,-1000,1000);
    else      hChargeDivPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-400,400,  200,-400,400);
    
    hname   = TString::Format("%s%d","hTimeDiffPos",i);
    htitle  = TString::Format("%s%d%s","fTimeDiffPos vs DCTY, paddle",i,";DCTY (mm);Pos (mm)");
    if (i<16) hTimeDiffPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-1000,1000, 200,-1000,1000);
    else      hTimeDiffPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-400,400,   200,-400,400);
    
    TString hname1   = TString::Format("%s%d","hSignal",i);
    TString htitle1  = TString::Format("%s%d%s","fSignal vs DCTY, paddle",i,";DCTY (mm);Signal (MIPs)");
    if (i<16) hSignal[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hSignal[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hSignalA",i);
    htitle1  = TString::Format("%s%d%s","Arithmetic mean vs DCTY, paddle",i,";DCTY (mm);Signal (MIPs)");
    if (i<16) hSignalA[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hSignalA[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hSignalG",i);
    htitle1  = TString::Format("%s%d%s","Geometric mean vs DCTY, paddle",i,";DCTY (mm);Signal (MIPs)");
    if (i<16) hSignalG[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hSignalG[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hChg",i);
    htitle1  = TString::Format("%s%d%s","Chg vs DCTY, paddle",i,";DCTY (mm);Chg (Z))");
    if (i<16) hChg[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hChg[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hChgA",i);
    htitle1  = TString::Format("%s%d%s","ChgA vs DCTY, paddle",i,";DCTY (mm);Chg (Z))");
    if (i<16) hChgA[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hChgA[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hChgG",i);
    htitle1  = TString::Format("%s%d%s","ChgG vs DCTY, paddle",i,";DCTY (mm);Chg (Z))");
    if (i<16) hChgG[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hChgG[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    TString hname2  = TString::Format("%s%d","hPosDiff",i);
    TString htitle2 = TString::Format("%s%d%s","Position diff (fTimeDiffPos-fChargeDivPos) vs DCTY, paddle",i,";DCTY (mm);Pos diff (mm)");
    if (i<16) hPosDiff[i]     = new TH2D(hname2.Data(),htitle2.Data(),200,-500,500, 200,-1000,1000);
    else      hPosDiff[i]     = new TH2D(hname2.Data(),htitle2.Data(),200,-400,400, 200,-400,400);
    
    // SC create histograms of extrapolated (X,Y) value for each paddle, see if we see the scintillator shapes
    TString hname3  = TString::Format("%s%d","hXY",i);
    TString htitle3 = TString::Format("%s%d%s","DCTY vs DCTX, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXY[i]          = new TH2D(hname3.Data(),htitle3.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXY[i]          = new TH2D(hname3.Data(),htitle3.Data(),200,-400,400, 200,-400,400);
    
    TString hname6  = TString::Format("%s%d","hXYshift",i);
    TString htitle6 = TString::Format("%s%d%s","shifted DCTY vs DCTX, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXYshift[i]     = new TH2D(hname6.Data(),htitle6.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXYshift[i]     = new TH2D(hname6.Data(),htitle6.Data(),200,-400,400, 200,-400,400);
    
    TString hname8  = TString::Format("%s%d","hYXshift",i);
    TString htitle8 = TString::Format("%s%d%s","shifted DCTX vs DCTY, paddle",i,";DCTY (mm);DCTY (mm)");
    if (i<16) hYXshift[i]     = new TH2D(hname8.Data(),htitle8.Data(),200,-1000,1000, 200,-200,200);
    else      hYXshift[i]     = new TH2D(hname8.Data(),htitle8.Data(),200,-400,400, 200,-400,400);
    hname8.ReplaceAll("hYXshift","pYXshift");
    if (i<16) pYXshift[i] = new TProfile(hname8.Data(),htitle8.Data(),200,-1000,1000, -200,200);
    else      pYXshift[i] = new TProfile(hname8.Data(),htitle8.Data(),200,-400,400,-400,400); //bore paddle

    TString hname4  = TString::Format("%s%d","hXYChi",i);
    TString htitle4 = TString::Format("%s%d%s","DCTY vs DCTX with DCTChi2 cut, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXYChi[i]       = new TH2D(hname4.Data(),htitle4.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXYChi[i]       = new TH2D(hname4.Data(),htitle4.Data(),200,-400,400, 200,-400,400);
    
    TString hname5  = TString::Format("%s%d","hXYHitC",i);
    TString htitle5 = TString::Format("%s%d%s","DCTY vs DCTX with DCT hitCount cut, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXYHitC[i]      = new TH2D(hname5.Data(),htitle5.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXYHitC[i]      = new TH2D(hname5.Data(),htitle5.Data(),200,-400,400, 200,-400,400);
    
    TString hname7  = TString::Format("%s%d","hthck",i);
    TString htitle7 = TString::Format("%s%d","Thickness (mm) for paddle",i);
    hthck[i]                  = new TH1D(hname7.Data(),htitle7.Data(),200,8.,14.);
    
    hname   = TString::Format("%s%d","hZdT",i);
    htitle  = TString::Format("%s%d%s","fChargeDivPos vs fTimeDiffPos, paddle",i,";fTimeDiffPos (mm);fChargeDivPos (mm)");
    if (i<16) hZdT[i] = new TH2D(hname.Data(),htitle.Data(),200,-1000,1000, 200,-1000,1000);
    else      hZdT[i] = new TH2D(hname.Data(),htitle.Data(),200,-400,400,   200,-1000,1000);
    hname.ReplaceAll("hZdT","pZdT");
    if (i<16) pZdT[i] = new TProfile(hname.Data(),htitle.Data(),200,-1000,1000, -1000,1000);
    else      pZdT[i] = new TProfile(hname.Data(),htitle.Data(),200,-400,400,   -1000,1000);
    
    hname   = TString::Format("%s%d","hZdct",i);
    htitle  = TString::Format("%s%d%s","DCTY vs fChargeDivFrac, paddle",i,";fChargeDivFrac;dctY (mm)");
    if (i<16) hZdct[i] = new TH2D(hname.Data(),htitle.Data(),50,-1,1,    200,-1000,1000);
    else      hZdct[i] = new TH2D(hname.Data(),htitle.Data(),50,-0.4,0.4,200,-1000,1000);
    hname.ReplaceAll("hZdct","pZdct");
    if (i<16) pZdct[i] = new TProfile(hname.Data(),htitle.Data(),50,-1,1,    -1000,1000);
    else      pZdct[i] = new TProfile(hname.Data(),htitle.Data(),50,-0.4,0.4,-1000,1000); //bore paddle
    
    hname              = TString::Format("%s%d","hTdct",i);
    htitle             = TString::Format("%s%d%s","DCTY vs fTimeDiff, paddle",i,";fTimeDiff;dctY (mm)");
    if (i<16) hTdct[i] = new TH2D(hname.Data(),htitle.Data(),60,-12,12,200,-1000,1000);
    else      hTdct[i] = new TH2D(hname.Data(),htitle.Data(),60,-4,4,  200,-1000,1000);
    hname.ReplaceAll("hTdct","pTdct");
    if (i<16) pTdct[i] = new TProfile(hname.Data(),htitle.Data(),60,-12,12,-1000,1000);
    else      pTdct[i] = new TProfile(hname.Data(),htitle.Data(),60,-4,4,  -1000,1000);
    
  }

  //end-based hists:

  TH2D** hTW;               //time walk explorations - one for each end
  hTW      = new TH2D*[72];  
  TH2D** heSig;             //ln(end signal) vs dcty, SPW getSignal()
  
  for (int i=0;i<72;i++)
  {
    TString hname   = TString::Format("%s%d","hTW",i);
    TString htitle  = TString::Format("%s%d%s","fTimeDiffPos - DCTy vs FEE signalSumRaw, fee ",i,";signalSumRaw;fTimeDiffPos - DCTY (mm)");
    hTW[i] = new TH2D(hname.Data(),htitle.Data(),400,0,2000, 200,-400,400);
  }

  
  //fee to fee end timing differences
  TH1D** hEnddT;
  hEnddT = new TH1D*[17*2];
  heSig      = new TH2D*[17*2];  
  TH2D** heSigA;             //ln(end signal) vs dcty, arithmetic mean
  heSigA      = new TH2D*[17*2];  
  TH2D** heSigG;             //ln(end signal) vs dcty, geometric mean
  heSigG      = new TH2D*[17*2];  
  
  for (int i=0;i<17;i++)
  {
    for (int j=0;j<2;j++)
    {
      TString hname = TString::Format("%s%d%s%d","hEnddT",i,"_",j);
      TString htitle = TString::Format("%s%d%s%d%s","Side by side time diff, paddle ",i,", end ",j,";dT (ns)");
//      cout << "i*2+j = " << i*2+j << "  " << hname << "   " << htitle << endl;
      hEnddT[i*2+j] = new TH1D(hname.Data(),htitle.Data(),1000,-20,20);

      TString hname1   = TString::Format("%s%d%s%d","heSig",i,"_",j);
      TString htitle1  = TString::Format("%s%d%s%d%s","loge(signal) vs (DCTY-padCenter), paddle ", i, " end ",j ,";(DCTY-padCenter) (mm);loge(signal) (ln(MIPs))");
      heSig[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(), 50,-800,800, 100,-2.,8.);

      hname1   = TString::Format("%s%d%s%d","heSigA",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signalA) vs (DCTY-padCenter),  paddle ", i, " end ",j, ";(DCTY-padCenter) (mm);loge(signalA) (ln(MIPs))");
      heSigA[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 100,-2.,8.);

      hname1   = TString::Format("%s%d%s%d","heSigG",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signalG) vs (DCTY-padCenter),  paddle ", i, " end ",j ,";(DCTY-padCenter) (mm);loge(signalG) (ln(MIPs))");
      heSigG[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 100,-2.,8.);
    }
  }
  
  //charge hists by correction level and signal construction type
  //All use sqrt(padSig)
  // Level 0: avg zenith and cos correction
  // L1: PLC and avg thickness correction
  // L2: thickness correction by position
  // L3: reponse correction
  
  const int nLZ = 4;

  TH1D** hZall;
  hZall = new TH1D*[nLZ]; //L0-3
  TH1D** hZAall;
  hZAall = new TH1D*[nLZ]; //L0-3
  TH1D** hZGall;
  hZGall = new TH1D*[nLZ]; //L0-3

  TH1D** hZtb;
  hZtb = new TH1D*[nLZ]; //L0-3
  TH1D** hZAtb;
  hZAtb = new TH1D*[nLZ]; //L0-3
  TH1D** hZGtb;
  hZGtb = new TH1D*[nLZ]; //L0-3
  
  for (int i=0; i<nLZ; i++)
  {
    TString hname   = TString::Format("%s%d","hZall",i);
    TString htitle  = TString::Format("%s%d","SPW Charge (t+b+bor), Level ", i);
    hZall[i] = new TH1D(hname.Data(),htitle.Data(),200,0.,20.);

    hname   = TString::Format("%s%d","hZAall",i);
    htitle  = TString::Format("%s%d","Arithmetic Charge (t+b+bor), Level ", i);
    hZAall[i] = new TH1D(hname.Data(),htitle.Data(),200,0.,20.);

    hname   = TString::Format("%s%d","hZGall",i);
    htitle  = TString::Format("%s%d","Geometric charge (t+b+bor), Level ", i);
    hZGall[i] = new TH1D(hname.Data(),htitle.Data(),200,0.,20.);

    hname   = TString::Format("%s%d","hZtb",i);
    htitle  = TString::Format("%s%d","SPW Charge (t+b), Level ", i);
    hZtb[i] = new TH1D(hname.Data(),htitle.Data(),200,0.,20.);

    hname   = TString::Format("%s%d","hZAtb",i);
    htitle  = TString::Format("%s%d","Arithmetic Charge (t+b), Level ", i);
    hZAtb[i] = new TH1D(hname.Data(),htitle.Data(),200,0.,20.);

    hname   = TString::Format("%s%d","hZGtb",i);
    htitle  = TString::Format("%s%d","Geometric charge (t+b), Level ", i);
    hZGtb[i] = new TH1D(hname.Data(),htitle.Data(),200,0.,20.);
  }
  
  
  // SC charge distributions, with either average thickness PLC or true thickness PLC correction
  /*
  // original charge determination scheme using ToF only: zenith and thickness corr. in an average way
  auto hallZo  = new TH1D("hallZo" , "Z t+b+bp orig." , 240, 0, 10);
  auto htbZo  = new TH1D("htbZo" , "Z t+b orig." , 240, 0, 10);
  auto hallZSCo  = new TH1D("hallZSCo" , "Z t+b+bp (SC) orig." , 240, 0, 10);
  auto hallZSCcorro  = new TH1D("hallZSCcorro" , "Z t+b+bp (SC+fudge) orig." , 240, 0, 10);
  auto htbZSCo  = new TH1D("htbZSCo" , "Z t+b (SC) orig." , 240, 0, 10);

  // charge determination done with PLC from DCT track and average thickness correction
  auto hallZ  = new TH1D("hallZ" , "Z t+b+bp" , 240, 0, 10);
  auto htbZ  = new TH1D("htbZ" , "Z t+b" , 240, 0, 10);
  auto hallZSC  = new TH1D("hallZSC" , "Z t+b+bp (SC)" , 240, 0, 10);
  auto hallZSCcorr  = new TH1D("hallZSCcorr" , "Z t+b+bp (SC+fudge)" , 240, 0, 10);
  auto htbZSC  = new TH1D("htbZSC" , "Z t+b (SC)" , 240, 0, 10);

  // charge determination adding thickness correction
  auto hallZ2  = new TH1D("hallZ2" , "Z t+b+bp, t corr." , 240, 0, 10);
  auto htbZ2  = new TH1D("htbZ2" , "Z t+b, t. corr." , 240, 0, 10);
  auto hallZSC2  = new TH1D("hallZSC2" , "Z t+b+bp (SC), t. corr." , 240, 0, 10);
  auto hallZSCcorr2  = new TH1D("hallZSCcorr2" , "Z t+b+bp (SC+fudge), t. corr." , 240, 0, 10);
  auto htbZSC2  = new TH1D("htbZSC2" , "Z t+b (SC), t. corr." , 240, 0, 10);

  // charge determination adding response correction
  auto hallZ3   = new TH1D("hallZ3" , "Z t+b+bp, t corr." , 240, 0, 10);
  auto htbZ3    = new TH1D("htbZ3" , "Z t+b, t. corr." , 240, 0, 10);
  auto hallZSC3 = new TH1D("hallZSC3" , "Z t+b+bp (SC), t. corr." , 240, 0, 10);
  auto hallZSCcorr3  = new TH1D("hallZSCcorr3" , "Z t+b+bp (SC+fudge), t. corr." , 240, 0, 10);
  auto htbZSC3  = new TH1D("htbZSC3" , "Z t+b (SC), t. corr." , 240, 0, 10);
*/
  //charge determination using DCT track for cos(zenith) and thickness
  TH2D* hBvT        = new TH2D("hBvT", "Signal: Bottom vs top;Top signal (MIP);Bottom signal (MIP)",240,0,12,240,0,12);
  TH2D* hZvR        = new TH2D("hZvR", "Signal: Zave vs R (GV);Rigidity (GV);Zave (MIP)",600,0,30,240,0,12);

  //time of flight:  
  auto hTOFdT   = new TH1D("hTOFdT" , "dTOF;ns" , 1000, -20, 20);
  
  //beta:
  auto hbeta   = new TH1D("hbeta" , "beta;(c)" , 1000, -2, 2);
  auto hinvbeta   = new TH1D("hinvbeta" , "1/beta;(1/c)" , 1000, -4, 4);
  TH2D* hZvbeta        = new TH2D("hZvbeta", "Signal: Zave vs beta (c);beta (c);Zave (MIP)",400,0.3,1.3,480,0,12);

  //mass:
  auto hmass   = new TH1D("hmass" , "mass;GV" , 1000, 0, 20);  
  
  //////////////////////////////////
  /////  INPUT EVENT LOOP /////
  //////////////////////////////////
  //loop over events
  auto nevts = toft->GetEntries();
  cout << "Number of events: " << nevts << endl;
  
  //////TEST////
//  nevts = 20;
  /////////////
  for (int ik=0; ik< nevts; ik++)
  {
    toft->GetEntry(ik);
    if (ik%10000 == 0) cout << "Event " << ik << endl;
    
    //get current set of tofcal parameters:
    getSCTOFCalibsForRun(runnum);
    
    
    //////////////////////////////////////////////////////
    //Collect info in shortcut arrays[3] by paddle:
    //////////////////////////////////////////////////////

    ///ORDERING: TOP, BORE, BOTTOM

    //HLX numbering:
    int toppadID = topTOF->fTag.getPaddleID();
    int borpadID = borTOF->fTag.getPaddleID();
    int botpadID = botTOF->fTag.getPaddleID();
    
    //Watch this trick: set up arrays to avoid lots of repeated code.
    int padID[3] = {toppadID, borpadID, botpadID};    
    HLX_Stg2_TOFHit* tofHits[3] = {topTOF, borTOF, botTOF};  //not copies, so use pointer ->
    double tofX[3]    = {tvars.ttofX, tvars.boreX,  tvars.btofX};
    double tofY[3]    = {tvars.ttofY, tvars.boreY,  tvars.btofY};
    double padPLC[3]  = {tvars.topPLC,tvars.borePLC,tvars.botPLC};

    //get paddle center locations in order to normalize DCTY to paddle center
    double padCenter[3] = {0.};
    double padLength[3] = {0.};
    for (int j=0; j<3; j++)
    {
      double west = tofHits[j]->fHitWest.getY();
      double east = tofHits[j]->fHitEast.getY();
      padCenter[j] = (west + east) / 2;
      padLength[j] = (west - east);
    }
    
    //paddle ends:
    HLX_Stg2_TOFPaddleEndHit tE     = topTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit tW     = topTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest);
    HLX_Stg2_TOFPaddleEndHit boreE  = borTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit boreW  = borTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest);
    HLX_Stg2_TOFPaddleEndHit bE     = botTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit bW     = botTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest);
    
    //watch this: an array of paddle ends! Order: Top, bore, bottom
    HLX_Stg2_TOFPaddleEndHit endHit[6] = {tE, tW, boreE, boreW, bE, bW};
    //cool, huh?

    //end signals
    double eSig[6]  = { tE.getSignal(), tW.getSignal(), boreE.getSignal(), boreW.getSignal(), bE.getSignal(), bW.getSignal() };
    double eSigA[6] = { getPadEndArithMeanSignal(tE), getPadEndArithMeanSignal(tW), getPadEndArithMeanSignal(boreE), 
                        getPadEndArithMeanSignal(boreW), getPadEndArithMeanSignal(bE), getPadEndArithMeanSignal(bW) };
    double eSigG[6] = { getPadEndGeomMeanSignal(tE), getPadEndGeomMeanSignal(tW), getPadEndGeomMeanSignal(boreE), 
                        getPadEndGeomMeanSignal(boreW), getPadEndGeomMeanSignal(bE), getPadEndGeomMeanSignal(bW) };
    
    
//  cout << ik << " eN: " << eSig [3] << "\t" << eSig [4] << "\t" << eSig [5] << endl;
//  cout << ik << " eA: " << eSigA[3] << "\t" << eSigA[4] << "\t" << eSigA[5] << endl;
//  cout << ik << " eG: " << eSigG[3] << "\t" << eSigG[4] << "\t" << eSigG[5] << endl;
//  cout << endl;
    
    //fill end-based signal histos
    for (int i=0; i<3; i++) //paddle
    {
      for (int j=0; j<2; j++) //end
      {
        int endID     =  endHit[i*2+j].fTag.getPaddleEndID();
        int globalEndID = padID[i]*2+endID;
        if (eSig[i]> 0) heSig [globalEndID]->Fill(tofY[i]-padCenter[i], log(eSig [i]));
        if (eSigA[i]>0) heSigA[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigA[i]));
        if (eSigG[i]>0) heSigG[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigG[i]));
      }
    }

    ///following distilled down from SC code to correct signal for paddle thicknesses
    //paddle locations:
    double padx[3] = {
      theGeom->fTOF.fTopLocX + theGeom->fTOF.fPaddleMetrology[toppadID]["xloc"],
      theGeom->fTOF.fBorePaddleLocX + theGeom->fTOF.fPaddleMetrology[borpadID]["xloc"],
      theGeom->fTOF.fBottomLocX + theGeom->fTOF.fPaddleMetrology[botpadID]["xloc"]
    };
    double pady[3] = {
      theGeom->fTOF.fTopLocY + theGeom->fTOF.fPaddleMetrology[toppadID]["yloc"],
      theGeom->fTOF.fBorePaddleLocY + theGeom->fTOF.fPaddleMetrology[borpadID]["yloc"],
      theGeom->fTOF.fBottomLocY + theGeom->fTOF.fPaddleMetrology[botpadID]["yloc"]
    };

    // check if track hit each paddle; set thickness
    double thck[3] = {10,10,10}; //mm, default thickness
    double tolx[3] = {100., 303., 100.}; //mm
    double toly[3] = {800., 303., 800.}; //mm
    
    //signals: These were corrected for pathlength (thickness) in TemplateAnal already, and now need to have that uncorrected before we can correct again here
    double padSig[3]  = {topTOF->getSignal(), borTOF->getSignal(), botTOF->getSignal()}; //SPW method
    double padSigA[3] = {getPadMeanArithSignal(topTOF), getPadMeanArithSignal(borTOF), getPadMeanArithSignal(botTOF)}; //arith mean of paddle ends. End fees combined with arith mean.
    double padSigG[3] = {getPadMeanGeomSignal(topTOF),  getPadMeanGeomSignal(borTOF),  getPadMeanGeomSignal(botTOF)}; //geom mean of paddle ends. End fees combined with arith mean.
    
//  cout << "Arith: " << getPadMeanArithSignal(topTOF) << endl;
//  cout << "Geom: " << getPadMeanGeomSignal(topTOF) << endl;
//  
//  cout << ik << " N: " << padSig [0] << "\t" << padSig [1] << "\t" << padSig [2] << endl;
//  cout << ik << " A: " << padSigA[0] << "\t" << padSigA[1] << "\t" << padSigA[2] << endl;
//  cout << ik << " G: " << padSigG[0] << "\t" << padSigG[1] << "\t" << padSigG[2] << endl;
//  cout << endl;
    
    double chg[3]  = {1., 1., 1.};
    double chgA[3] = {1., 1., 1.};
    double chgG[3] = {1., 1., 1.};
    double tcorr[3] = {1., 1., 1.};   // thickness correction factor to adjust the Z value as per impact location thickness
    
    //caclulate correction for paddle Signal by thickness:
    // SC  Average scintillator thickness in mm, which was applied in pulse signal size normalization but which will have to be uncorrected back to raw
    // to apply the correct thickness correction; order here is HELIX ordering (bottom 0-7, top 8-15, bore 16)
    Double_t avgthick[17] = {
      9.79,9.64,10.19,9.56,10.09,9.72,9.84,9.87,  // bottom
      9.56,9.67,9.61,9.88,9.84,9.56,9.97,9.62,  // top
      9.97  // bore
    };

    // SC  average coz(zenith) angle between center of paddle and center of BP; this is done by paddle number, not FEE number ; for BP, use the same coszen as the top paddle hit.
    Double_t coszen[16] = {0.922,0.956,0.977,0.980,0.978,0.961,0.930,0.890, // bottom paddles
    0.922,0.956,0.977,0.980,0.978,0.961,0.930,0.890};   // top paddles average cos(zenith) from toy MC

    //a lamda to check to see if hit is inside paddle:
    auto isInPaddle = [](double x, double xc, double tol) -> bool { 
        return (x <= xc + tol) && (x >= xc - tol); 
    };

    //fill SC hit location histos and get paddle thickness at hit location
    for (int i=0; i<3; i++)
    {
      if (isInPaddle(tofX[i], padx[i], tolx[i]) && isInPaddle(tofY[i], pady[i], toly[i]))
      {
        hXYshift[padID[i]]->Fill(tofX[i]-padx[i],tofY[i]-pady[i]);
        hYXshift[padID[i]]->Fill(tofY[i]-pady[i],tofX[i]-padx[i]); 
        pYXshift[padID[i]]->Fill(tofY[i]-pady[i],tofX[i]-padx[i]);
        double *padt2D = nullptr;
        //from tofcal11:
//        if (i==0) padt2D = &fCal.topt[padID[i]][0][0];
//        if (i==2) padt2D = &fCal.bott[padID[i]-8][0][0];
        //from tofcal12:
        if (i==0 || i==2)
        {
          padt2D = &fCal.padt[padID[i]][0][0];
          thck[i] = thickness(i, padt2D, (tofX[i]-padx[i])/10. + tolx[i]/10., (tofY[i]-pady[i])/10. + toly[i]/10.);  // positions in cm, not mm (/10.), and need to shift by 1/2 paddle dimension (because of thickness map for long paddles)
        } 
        if (i==1) //bore
        {
          padt2D = &fCal.boret[0][0][0]; //bore
          thck[i] = thickness(i, padt2D, (tofX[i]-padx[i])/10. , (tofY[i]-pady[i])/10. );  // positions in cm, not mm (/10.); the BP does NOT need the 1/2 paddle shift
        }
      	hthck[padID[i]] -> Fill(thck[i]);        
      }
//      tcorr[i] = 1./thck[i];   // this tcorr factor will *multiply* the signal size
        tcorr[i] = avgthick[padID[i]] / thck[i];   // this tcorr factor will *multiply* the signal size
    }
    
    /////////////////////////////////////////////////////
    //FILL ALL HISTS AND PROFILES
    /////////////////////////////////////////////////////

    //paddle end fee hists
    //HLX: 0-7 lower, 8-15 upper, 16 is bore  
    for (int i=0; i<6; i++)
    {
      int nFees     = endHit[i].getFEEHitCount();
      int paddleID  = endHit[i].fTag.getPaddleID();
      int endID     =  endHit[i].fTag.getPaddleEndID();
      if (paddleID == 16) continue; //i don't know how to deal with dT for bore paddle
      if ( nFees < 2 ) continue;
      if ( endHit[i].getFEEHit(0).fStatusFlag != HLX_Stg2_TOFFEEHit::kGood ||  endHit[i].getFEEHit(1).fStatusFlag != HLX_Stg2_TOFFEEHit::kGood) continue;
      double dT = 0;
//      cout << "Fee times" << endl;
//      std::streamsize original_precision = std::cout.precision();
//      cout << std::fixed << std::setprecision(3);
      for (int j=0; j<nFees; j++)
      {
        double feeTime = endHit[i].getFEEHit(j).fTime;
        int feeID = endHit[i].getFEEHit(j).fTag.getFEEID();
//        cout << paddleID << " " << endID << " " << feeID << "  " << endHit[i].getFEEHit(j).fTime << endl;
        if (feeID == 0) dT=feeTime; //S or 0
        else      dT = feeTime - dT; //N-S (1-0)
      }
//      std::cout.precision(original_precision); std::cout << std::defaultfloat; 
//      if (dT == 0.) cout << "Pad, times: " << paddleID << "  " << endHit[i].getFEEHit(0).fTime << "  " << endHit[i].getFEEHit(1).fTime << "  " << dT << endl;
      hEnddT[paddleID*2+endID]->Fill(dT); //N-S
//      cout << i << " " << paddleID << " " << endID << " " << dT << endl;
    }
    
    //NOTE: Apply corrections specific to this version of skeletor output:
    for (int i=0; i<3; i++) 
    {
//      padSig[i] *= avgthick[padID[i]];
      padSig [i] *= fCal.mpv[8]/fCal.mpv[0];   //correct for getting the wrong normalization fee.
    } 
    
    //Fill hists by paddle:
    for (int j=0; j<3; j++)
    {
      //do recalcs here if requested to check new parameters:
      double chargeDivPos = tofHits[j]->fChargeDivPos;
      if (recalcChargeDivPos)
      {
        //ChargeDivFrac position from fit to DCTY vs fChargeDivFrac
        double Zdiff = tofHits[j]->fChargeDivFrac;
        chargeDivPos = fCal.ZDm[padID[j]] * Zdiff + fCal.ZDb[padID[j]]; // = DCTY - center
        chargeDivPos += padCenter[j];

      }
      double timeDiffPos = tofHits[j]->fTimeDiffPos;
      if (recalcTimeDiffPos)
      {
        //tofcal12:
        timeDiffPos = fCal.TDm[padID[j]] * tofHits[j]->fTimeDiff + fCal.TDb[padID[j]]; // = DCTY - center
        timeDiffPos += padCenter[j];
      }
      
      //use this if not doing sneaky charge histos after this loop
      //extract signal, correct for scintillator thickness and pathlength
//      padSig [j] *= fCal.mpv[8]/fCal.mpv[0];   //correct for getting the wrong normalization fee.
//      padSigA[j] *= fCal.mpv[8]/fCal.mpv[0];   //correct for getting the wrong normalization fee.
//      padSigG[j] *= fCal.mpv[8]/fCal.mpv[0];   //correct for getting the wrong normalization fee.
//      padSig [j] *= tcorr[j];   //CHANGE AFTER AVG THICKNESS CORRECTION REMOVED FROM TOF TEMPLATEANAL STG2 CODE
//      padSigA[j] *= tcorr[j];   //CHANGE AFTER AVG THICKNESS CORRECTION REMOVED FROM TOF TEMPLATEANAL STG2 CODE
//      padSigG[j] *= tcorr[j];   //CHANGE AFTER AVG THICKNESS CORRECTION REMOVED FROM TOF TEMPLATEANAL STG2 CODE
//      padSig [j] /= padPLC[j];   //correct for pathlength
//      padSigA[j] /= padPLC[j];   //correct for pathlength
//      padSigG[j] /= padPLC[j];   //correct for pathlength
//      padSig [j] /= responseCorr(padID[j],tofY[j]-padCenter[j]);   //correct for response along paddle length
//      padSigA[j] /= responseCorr(padID[j],tofY[j]-padCenter[j]);   //correct for response along paddle length
//      padSigG[j] /= responseCorr(padID[j],tofY[j]-padCenter[j]);   //correct for response along paddle length
//      chg[j]     = sqrt(padSig[j]);
      

      //fill histos:
      hChargeDivPos [padID[j]]->Fill(tofY[j]-padCenter[j],        chargeDivPos-padCenter[j]);  //consistent with calculation of timeDiffPos and chargeDivPos above and in TemplateAnaly
      hTimeDiffPos  [padID[j]]->Fill(tofY[j]-padCenter[j],        timeDiffPos-padCenter[j]);
//      hSignal       [padID[j]]->Fill(tofY[j]-padCenter[j],        padSig[j]);
//      hSignalA      [padID[j]]->Fill(tofY[j]-padCenter[j],        padSigA[j]);
//      hSignalG      [padID[j]]->Fill(tofY[j]-padCenter[j],        padSigG[j]);
//      hChg          [padID[j]]->Fill(tofY[j]-padCenter[j],        chg[j]);
      hPosDiff      [padID[j]]->Fill(tofY[j]-padCenter[j],        timeDiffPos - chargeDivPos);
      hZdT          [padID[j]]->Fill(timeDiffPos,                 chargeDivPos);
      pZdT          [padID[j]]->Fill(timeDiffPos,                 chargeDivPos);
      hZdct         [padID[j]]->Fill(tofHits[j]->fChargeDivFrac,  tofY[j]-padCenter[j]);
      pZdct         [padID[j]]->Fill(tofHits[j]->fChargeDivFrac,  tofY[j]-padCenter[j]);
      hTdct         [padID[j]]->Fill(tofHits[j]->fTimeDiff,       tofY[j]-padCenter[j]);
      pTdct         [padID[j]]->Fill(tofHits[j]->fTimeDiff,       tofY[j]-padCenter[j]);
      //SC hists:
      hXY           [padID[j]]->Fill(tofX[j],                     tofY[j]-padCenter[j]);
      if (DCTChi2<2.)
	      hXYChi      [padID[j]]->Fill(tofX[j],                     tofY[j]-padCenter[j]);
      if (hitCount>60.)
	      hXYHitC     [padID[j]]->Fill(tofX[j],                     tofY[j]-padCenter[j]);
      
      
    //fill end-based hists for each paddle:
      HLX_Stg2_TOFPaddleEndHit ends[2] = {tofHits[j]->getPaddleEndHit(HLX_Geometry_TOF::kEast), tofHits[j]->getPaddleEndHit(HLX_Geometry_TOF::kWest)};
      for (int k = 0; k<2; k++) //ends
      {
        int nFees     = ends[k].getFEEHitCount();
        int paddleID  = ends[k].fTag.getPaddleID();
        int endID     = ends[k].fTag.getPaddleEndID();
        for (int jk=0; jk<nFees; jk++)
        {
          double feeTime = ends[k].getFEEHit(jk).fTime;
          int feeID = ends[k].getFEEHit(jk).fTag.getFEEID();
          auto [laneIdx, chanIdx] = theGeom->fTOF.fGeoMap[{paddleID, endID, feeID}];
          int SCfeeID = laneIdx*8 + chanIdx; //HLX labels. tofcal12
          double feeSignal = ends[k].getFEEHit(jk).fSignalSumRaw;
          double tdiff = tofHits[j]->fTimeDiffPos - (tofY[j]-padCenter[j]);
          hTW[SCfeeID]->Fill(feeSignal, tdiff);
        }
      } //end of end loop
      
    } //end of paddle loop
    
    /*
    //charge determinations:  SC 10/12/2016   Lots of new plots
    //modified 10/23/2026 SLN to make easier to follow and successively apply corrections

    double expo1 = 1.0/1.7;   // nominal power to use in calculating Z
    double expo2 = 0.5;       // SC uses sqrt
    
    ///////
    //I think this section is ripe for deletion
    double ttofZo  = pow(padSig[0] * coszen[padID[0]], expo1);
    double btofZo  = pow(padSig[2] * coszen[padID[2]], expo1);
    double boreZo  = pow(padSig[1] * coszen[padID[0]], expo1);  // use top paddle coszen
    //    cout << "ttofZo, btofZo, boreZo: " << ttofZo << " " << btofZo << " " << boreZo << endl;
    double allZo   = (ttofZo+boreZo+btofZo)/3.;
    double tbZo    = (ttofZo+btofZo)/2.;
    //    cout << "allZo, tbZo: " << allZo << " " << tbZo << endl;

    double ttofZSCo  = padSig[0] * coszen[padID[0]];
    double btofZSCo  = padSig[2] * coszen[padID[2]];
    double boreZSCo  = padSig[1] * coszen[padID[0]];  // use top paddle coszen
    double allZSCo   = pow((ttofZSCo+boreZSCo+btofZSCo)/3.,expo2);   // SC   original
    //    double allZSCo   = sqrt((ttofZSCo+boreZSCo+btofZSCo)/3.);   // SC   original way
    //      double allZSCo   = (sqrt(ttofZSCo)+sqrt(boreZSCo)+sqrt(btofZSCo))/3.;  // SC  this seems to give better Z res.
    double allZSCcorro = (allZSCo - 0.2499) / 0.8385;
    double tbZSCo    = pow((ttofZSCo+btofZSCo)/2.,expo2);   //SC   original
    //    double tbZSCo    = sqrt((ttofZSCo+btofZSCo)/2.);   //SC   original way
    //      double tbZSCo    = (sqrt(ttofZSCo)+sqrt(btofZSCo))/2.;   // SC  this seems to give better Z res.

    if (abs(padSig[0]*coszen[padID[0]]-padSig[2]*coszen[padID[2]])<0.4 ) { // apply top/bottom tight charge agreement
      hallZo -> Fill(allZo);
      htbZo -> Fill(tbZo);

      hallZSCo -> Fill(allZSCo);
      hallZSCcorro -> Fill(allZSCcorro);
      htbZSCo -> Fill(tbZSCo);
    }
    ///////////
*/
    
    //padSig current status: corrected for avg thickness only.
    
    //Zeroth: avg coszen
    double tof0[3], tof0A[3], tof0G[3]; //don't change padSig!
    for (int i=0; i<3; i++) 
    { 
      tof0 [i] = padSig [i] * coszen[padID[i]];
      tof0A[i] = padSigA[i] * coszen[padID[i]];
      tof0G[i] = padSigG[i] * coszen[padID[i]];
    }
    double zall = (tof0 [0] + tof0[1] + tof0[2])/3.;
    double ztb  = (tof0 [0] + tof0[2])/2.;
    if ( zall > 0 && abs(tof0[0] - tof0[2]) < 0.4) hZall[0]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(tof0[0] - tof0[2]) < 0.4) hZtb [0]->Fill(sqrt(ztb));

    zall = (tof0A[0] + tof0A[1] + tof0A[2])/3.;
    ztb  = (tof0A[0] + tof0A[2])/2.;
    if ( zall > 0 && abs(tof0A[0] - tof0A[2]) < 0.4) hZAall[0]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(tof0A[0] - tof0A[2]) < 0.4) hZAtb [0]->Fill(sqrt(ztb));

    zall = (tof0G[0] + tof0G[1] + tof0G[2])/3.;
    ztb  = (tof0G[0] + tof0G[2])/2.;
    if ( zall > 0 && abs(tof0G[0] - tof0G[2]) < 0.4) hZGall[0]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(tof0G[0] - tof0G[2]) < 0.4) hZGtb [0]->Fill(sqrt(ztb));

    //FIRST: apply PLC
    for (int i=0; i<3; i++) 
    {
      padSig [i] /= padPLC[i];
      padSigA[i] /= padPLC[i];
      padSigG[i] /= padPLC[i];
    }

    zall = (padSig [0] + padSig [1] + padSig[2])/3.;
    ztb  = (padSig [0] + padSig [2])/2.;
    if ( zall > 0 && abs(padSig [0] - padSig [2]) < 0.4) hZall[1]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSig [0] - padSig [2]) < 0.4) hZtb [1]->Fill(sqrt(ztb));

    zall = (padSigA[0] + padSigA[1] + padSigA[2])/3.;
    ztb  = (padSigA[0] + padSigA[2])/2.;
    if ( zall > 0 && abs(padSigA[0] - padSigA[2]) < 0.4) hZAall[1]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSigA[0] - padSigA[2]) < 0.4) hZAtb [1]->Fill(sqrt(ztb));

    zall = (padSigG[0] + padSigG[1] + padSigG[2])/3.;
    ztb  = (padSigG[0] + padSigG[2])/2.;
    if ( zall > 0 && abs(padSigG[0] - padSigG[2]) < 0.4) hZGall[1]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSigG[0] - padSigG[2]) < 0.4) hZGtb [1]->Fill(sqrt(ztb));

/*    
    double tofZ1[3] = {pow(padSig[0], expo1), pow(padSig[1], expo1), pow(padSig[2], expo1)};
    double allZ1    = (tofZ1[0]+tofZ1[1]+tofZ1[2])/3.;
    double tbZ1     = (tofZ1[0]+tofZ1[2])/2.;

    double allZ1SC      = pow((tofZ1[0]+tofZ1[1]+tofZ1[2])/3.,expo2);   // SC  original
    double allZ1SCcorr  = (allZ1SC - 0.2499) / 0.8385;
    double tbZ1SC       = pow((tofZ1[0]+tofZ1[2])/2.,expo2);   // SC  original

    if (abs(tofZ1[0]-tofZ1[2])<0.4) {   // apply top/bottom tight charge agreement
      hallZ -> Fill(allZ1);
      htbZ -> Fill(tbZ1);
      
      hallZSC -> Fill(allZ1SC);
      hallZSCcorr -> Fill(allZ1SCcorr);
      htbZSC -> Fill(tbZ1SC);
    }
*/	
    //SECOND: Thickness corrected Z signals
    for (int i=0; i<3; i++)  
    {
      padSig[i] *= tcorr[i];
      padSigA[i] *= tcorr[i];
      padSigG[i] *= tcorr[i];
    }
    
    zall = (padSig [0] + padSig [1] + padSig[2])/3.;
    ztb  = (padSig [0] + padSig [2])/2.;
    if ( zall > 0 && abs(padSig [0] - padSig [2]) < 0.4) hZall[2]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSig [0] - padSig [2]) < 0.4) hZtb [2]->Fill(sqrt(ztb));

    zall = (padSigA[0] + padSigA[1] + padSigA[2])/3.;
    ztb  = (padSigA[0] + padSigA[2])/2.;
    if ( zall > 0 && abs(padSigA[0] - padSigA[2]) < 0.4) hZAall[2]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSigA[0] - padSigA[2]) < 0.4) hZAtb [2]->Fill(sqrt(ztb));

    zall = (padSigG[0] + padSigG[1] + padSigG[2])/3.;
    ztb  = (padSigG[0] + padSigG[2])/2.;
    if ( zall > 0 && abs(padSigG[0] - padSigG[2]) < 0.4) hZGall[2]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSigG[0] - padSigG[2]) < 0.4) hZGtb [2]->Fill(sqrt(ztb));

    /*
    double tofZ2[3] = {pow(padSig[0], expo1), pow(padSig[1], expo1), pow(padSig[2], expo1)};
    double allZ2    = (tofZ2[0]+tofZ2[1]+tofZ2[2])/3.;
    double tbZ2     = (tofZ2[0]+tofZ2[2])/2.;

    double allZSC2      = pow((tofZ2[0]+tofZ2[1]+tofZ2[2])/3.,expo2);   // SC  original
    double allZSCcorr2  = (allZSC2 - 0.2499) / 0.8385;
    double tbZSC2       = pow((tofZ2[0]+tofZ2[2])/2.,expo2);   // SC  original
*/
//fill signal and charge vs dcty hists at this point:
    for (int i=0; i<3; i++)  
    {
      chg[i] = (padSig[i]>0)? sqrt(padSig[i]) : 0;
      hChg[padID[i]]->Fill(tofY[i]-padCenter[i], chg[i]);
      chgA[i] = (padSigA[i]>0)? sqrt(padSigA[i]) : 0;
      hChgA[padID[i]]->Fill(tofY[i]-padCenter[i], chgA[i]);
      chgG[i] = (padSigG[i]>0)? sqrt(padSigG[i]) : 0;
      hChgG[padID[i]]->Fill(tofY[i]-padCenter[i], chgG[i]);
      
      hSignal [padID[i]]->Fill(tofY[i]-padCenter[i], padSig [i]);
      hSignalA[padID[i]]->Fill(tofY[i]-padCenter[i], padSigA[i]);
      hSignalG[padID[i]]->Fill(tofY[i]-padCenter[i], padSigG[i]);
      
    }

    /*
    if (abs(tofZ2[0]-tofZ2[2])<0.4 ) { // apply top/bottom tight charge agreement
      hallZ2 -> Fill(allZ2);
      htbZ2 -> Fill(tbZ2);

      hallZSC2 -> Fill(allZSC2);
      hallZSCcorr2 -> Fill(allZSCcorr2);
      htbZSC2 -> Fill(tbZSC2);
    }
*/

    //THIRD: response corrected version
    //NB: THESE SHOULD BE DIFFERENT PARAMETERS FOR EACH METHOD! RIGHT NOW: SPW method parameters
    for (int i=0; i<3; i++)  
    {
      padSig [i] /= responseCorr(padID[i],tofY[i]-padCenter[i]);   //correct for response along paddle length based on hSignal vs DCTY
      padSigA[i] /= responseCorr(padID[i],tofY[i]-padCenter[i]);   //correct for response along paddle length
      padSigG[i] /= responseCorr(padID[i],tofY[i]-padCenter[i]);   //correct for response along paddle length
    }
    
    zall = (padSig [0] + padSig [1] + padSig[2])/3.;
    ztb  = (padSig [0] + padSig [2])/2.;
    if ( zall > 0 && abs(padSig [0] - padSig [2]) < 0.4) hZall[3]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSig [0] - padSig [2]) < 0.4) hZtb [3]->Fill(sqrt(ztb));

    zall = (padSigA[0] + padSigA[1] + padSigA[2])/3.;
    ztb  = (padSigA[0] + padSigA[2])/2.;
    if ( zall > 0 && abs(padSigA[0] - padSigA[2]) < 0.4) hZAall[3]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSigA[0] - padSigA[2]) < 0.4) hZAtb [3]->Fill(sqrt(ztb));

    zall = (padSigG[0] + padSigG[1] + padSigG[2])/3.;
    ztb  = (padSigG[0] + padSigG[2])/2.;
    if ( zall > 0 && abs(padSigG[0] - padSigG[2]) < 0.4) hZGall[3]->Fill(sqrt(zall));
    if ( ztb  > 0 && abs(padSigG[0] - padSigG[2]) < 0.4) hZGtb [3]->Fill(sqrt(ztb));

/*    
    double tofZ3[3]  = {pow(padSig[0], expo1), pow(padSig[1], expo1), pow(padSig[2], expo1)};
    double allZ3   = (tofZ3[0]+tofZ3[1]+tofZ3[2])/3.;
    double tbZ3    = (tofZ3[0]+tofZ3[2])/2.;

    double allZSC3      = pow((tofZ3[0]+tofZ3[1]+tofZ3[2])/3.,expo2);   // SC  original
    double allZSCcorr3  = (allZSC3 - 0.2499) / 0.8385;
    double tbZSC3       = pow((tofZ3[0]+tofZ3[2])/2.,expo2);   // SC  original
    
    if (abs(tofZ3[0]-tofZ3[2])<0.4 ) { // apply top/bottom tight charge agreement
      hallZ3 -> Fill(allZ3);
      htbZ3 -> Fill(tbZ3);

      hallZSC3 -> Fill(allZSC3);
      hallZSCcorr3 -> Fill(allZSCcorr3);
      htbZSC3 -> Fill(tbZSC3);
    }
*/
    ///TOF Beta and Charge
    double timeDiff = botTOF->fTime - topTOF->fTime;
    hTOFdT->Fill(timeDiff);
    
    // beta distribution
    double beta = tvars.pathLength / 100. / timeDiff / 2.998;    // this will give beta in units of c
    hbeta->Fill(beta);
    if (beta!=0) {
      hinvbeta->Fill(1./beta);
    }

    //this next cut requires the hits to be in all three paddles - top, bore, bot
    bool inPaddle = true;
    for (int i=0; i<3; i++){
      if (!isInPaddle(tofX[i], padx[i], tolx[i]) || !isInPaddle(tofY[i], pady[i], toly[i])) {inPaddle = false;}
    }
  
    //    if (abs(topSig/tvars.topPLC*tcorr[0]-botSig/tvars.botPLC*tcorr[2])<0.2 ) { // apply top/bottom tight charge agreement
    if (abs(chg[0]-chg[2])<0.2 ) { // apply top/bottom tight charge agreement
      if (DCTChi2<2.&&hitCount>60.&&inPaddle)   // DCT cuts, including whether track hits the right paddles
      {
//	      if (beta>0.2&&beta<1.01) {hmass->Fill(allZSC2*rig*sqrt(pow(1./beta,2.)-1.));}   // mass*c^2 in GV units 
	      if (beta>0.2&&beta<1.01) {hmass->Fill(zall*rig*sqrt(pow(1./beta,2.)-1.));}   // mass*c^2 in GV units 
      }
    }
    
//    hZvR->Fill(rig, allZSC2); //rigidity on x; Zave on y
//    hZvbeta->Fill(beta,allZSC2);
    hZvR->Fill(rig, zall); //rigidity on x; Zave on y
    hZvbeta->Fill(beta,zall);
    //    hBvT->Fill(pow(padSig[0],expo2), pow(padSig[2],expo2)); //top signal on x; bottom signal on y
    hBvT->Fill(chgG[0],chgG[2]); //top signal on x; bottom signal on y

  } //end of event loop
  
  cout << "End of event loop" << endl;
  opf->Write();
  opf->Close();


}

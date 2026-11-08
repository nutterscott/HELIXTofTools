//void AnalTOFSkelTree(TTree *toft, const TString &outfile, bool recalcChargeDivPos = true, bool recalcTimeDiffPos = true )  //infile is toft tree output created by skeletor
//input: toft or chain of skeletor toft trees
//outfile: name of output file 

//Driver (at end of this file):
//void AnalTOFSkelList(const char *listName = "filelist.txt", const char *outName  = "combined.AnalTOFSkel.root")
//  std::string dir = "/data/skeletor/";
//  std::string fnameStem = "skeletor.";
//  std::string fnameExt = ".root";


// to run with list of runs:
//void AnalTOFSkelList(const char *listName = "filelist.txt",
//                    const char *outName  = "combined.AnalTOFSkel.root")

//NOW SET TO WORK WITH skeletor.v12.2.NT (new tracking) as input. 
// BUG FIXES TO fSignalSum APPLIED.

//
// SC 10/12/2026   repackage AnalTOFSkel.C into AnalTOFSkelTree.C which will be invoked within AnalTOFSkelChain.C to
// process multiple skeletor files (in a TChain) with a list of input trees read from a text file

////////////////////////////////////////////////

//include PSU tof calib tools:
#include "PSUTofTools13.C"
/////  SET THIS IN PSUTofTools13.C!! ////
//Set up input tofcal ttree:
//string kInputTofcalName = "tofcal13.1.root";

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
  loadSCTOFCalibsv13();
  
  //set up OUTPUT tofcal calibration  
//  createSCTOFCalibs();
  
  //load response correction parameters:
  //from padsigG corrections in Bon.AnalTOFv4.v12.2.NT.2.root
  string kInputResponseCorrFileName = "Bon.AnalTOFv4.v12.2.NT.2.SignalGHiR.PadSigTGraphFits.txt";
  if (outfile.Contains("Boff")) kInputResponseCorrFileName = "Boff.AnalTOFv4.v12.2.NT.2.SignalGHiR.PadSigTGraphFits.txt";
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

  //pos from timing, chargeDivFrac, and difference
  TH2D** hChargeDivPos;                 //position from charge division 
  hChargeDivPos       = new TH2D*[17];
  TH2D** hSignal;                       // Charge (getSignal() = geometric mean of end fees and geometric mean of ends) vs DCTY for all paddles
  hSignal             = new TH2D*[17];
  TH2D** hSignalA;                       // Charge (arithmetic mean of end fees and arith mean of ends) vs DCTY for all paddles
  hSignalA            = new TH2D*[17];
  TH2D** hSignalG;                       // Charge (arithmetic mean of end fees and geometric mean of ends) vs DCTY for all paddles
  hSignalG            = new TH2D*[17];
  TH2D** hSignalHiR;                       // Charge (getSignal() = geometric mean of end fees and geometric mean of ends) vs DCTY for all paddles
  hSignalHiR          = new TH2D*[17];
  TH2D** hSignalAHiR;                       // Charge (arithmetic mean of end fees and arith mean of ends) vs DCTY for all paddles
  hSignalAHiR         = new TH2D*[17];
  TH2D** hSignalGHiR;                       // Charge (arithmetic mean of end fees and geometric mean of ends) vs DCTY for all paddles
  hSignalGHiR         = new TH2D*[17];
  TH2D** hChg;                       // Charge (sqrt(hSignal)) vs DCTY for all paddles
  hChg                = new TH2D*[17];
  TH2D** hChgA;                       // Charge (sqrt(hSignalA)) vs DCTY for all paddles
  hChgA               = new TH2D*[17];
  TH2D** hChgG;                       // Charge (sqrt(hSignalG)) vs DCTY for all paddles
  hChgG               = new TH2D*[17];
  TH2D** hChgHiR;                       // Charge (sqrt(hSignal)) vs DCTY for all paddles
  hChgHiR             = new TH2D*[17];
  TH2D** hChgAHiR;                       // Charge (sqrt(hSignalA)) vs DCTY for all paddles
  hChgAHiR            = new TH2D*[17];
  TH2D** hChgGHiR;                       // Charge (sqrt(hSignalG)) vs DCTY for all paddles
  hChgGHiR            = new TH2D*[17];
  TH2D** hTimeDiffPos;                  //position from end to end time difference 
  hTimeDiffPos        = new TH2D*[17];
  TH2D** hPosDiff;                      //comparison of position from charge division and time diff
  hPosDiff            = new TH2D*[17];
  TH2D** hXY;                           // SC  check (X,Y) for each paddle
  hXY                 = new TH2D*[17];
  TH2D** hXYshift;                      // SC  check (X,Y) for each paddle; (X,Y) shifted to a paddle position (0-20 cm, 0-180 cm)
  hXYshift            = new TH2D*[17];
  TH2D** hYXshift;                      // SC  check (Y,X) for each paddle; (Y,X) shifted to a paddle position (0-20 cm, 0-180 cm);  YX is interesting to look at instead of XY to look for possible paddle rotations
  hYXshift            = new TH2D*[17];
  TProfile** pYXshift;
  pYXshift            = new TProfile*[17];
  TH2D** hXYChi;                        // SC  check (X,Y) for each paddle, after a DCTChi2 cut
  hXYChi              = new TH2D*[17];
  TH2D** hXYHitC;                       // SC  check (X,Y) for each paddle, after a DCT hitCount cut
  hXYHitC             = new TH2D*[17];
  TH1D** hthck;                         // SC  thckness at (X,Y) for each paddle
  hthck               = new TH1D*[17];
  TH2D** hZdT;                          //fChargeDivPos vs fTimeDiffPos timing differences
  hZdT                = new TH2D*[17];
  TProfile** pZdT;
  pZdT                = new TProfile*[17];
  TH2D** hZdct;                         //position differences: dcty vs fChargeDivFrac
  hZdct               = new TH2D*[17];
  TProfile** pZdct;
  pZdct               = new TProfile*[17];
  TH2D** hTdct;                         //position differences: fTimeDiff vs dcty
  hTdct               = new TH2D*[17];
  TProfile** pTdct;
  pTdct               = new TProfile*[17];
  

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
    
    hname1   = TString::Format("%s%d","hSignalHiR",i);
    htitle1  = TString::Format("%s%d%s","fSignal vs DCTY for R>2, paddle",i,";DCTY (mm);Signal (MIPs)");
    if (i<16) hSignalHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hSignalHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hSignalAHiR",i);
    htitle1  = TString::Format("%s%d%s","Arithmetic mean vs DCTY for R>2, paddle",i,";DCTY (mm);Signal (MIPs)");
    if (i<16) hSignalAHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hSignalAHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hSignalGHiR",i);
    htitle1  = TString::Format("%s%d%s","Geometric mean vs DCTY for R>2, paddle",i,";DCTY (mm);Signal (MIPs)");
    if (i<16) hSignalGHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hSignalGHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
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
    
    hname1   = TString::Format("%s%d","hChgHiR",i);
    htitle1  = TString::Format("%s%d%s","Chg vs DCTY for R>2, paddle",i,";DCTY (mm);Chg (Z))");
    if (i<16) hChgHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hChgHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hChgAHiR",i);
    htitle1  = TString::Format("%s%d%s","ChgA vs DCTY for R>2, paddle",i,";DCTY (mm);Chg (Z))");
    if (i<16) hChgAHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hChgAHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
    hname1   = TString::Format("%s%d","hChgGHiR",i);
    htitle1  = TString::Format("%s%d%s","ChgG vs DCTY for R>2, paddle",i,";DCTY (mm);Chg (Z))");
    if (i<16) hChgGHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 800,0,20);
    else      hChgGHiR[i]      = new TH2D(hname1.Data(),htitle1.Data(),50,-300,300, 800,0,20);
    
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
  
  for (int i=0;i<72;i++)
  {
    TString hname   = TString::Format("%s%d","hTW",i);
    TString htitle  = TString::Format("%s%d%s","fTimeDiffPos - DCTy vs FEE signalSumRaw, fee ",i,";signalSumRaw;fTimeDiffPos - DCTY (mm)");
    hTW[i] = new TH2D(hname.Data(),htitle.Data(),400,0,2000, 200,-400,400);
  }

  
  //fee to fee end timing differences
  TH1D** hEnddT;
  hEnddT = new TH1D*[17*2];
  TH2D** heSig;             //ln(end signal) vs dcty, SPW getSignal()
  heSig      = new TH2D*[17*2];  
  TH2D** heSigA;             //ln(end signal) vs dcty, arithmetic mean
  heSigA      = new TH2D*[17*2];  
  TH2D** heSigG;             //ln(end signal) vs dcty, geometric mean
  heSigG      = new TH2D*[17*2];  
  TH2D** heSigHiR;             //ln(end signal) vs dcty, SPW getSignal()
  heSigHiR      = new TH2D*[17*2];  
  TH2D** heSigAHiR;             //ln(end signal) vs dcty, arithmetic mean
  heSigAHiR      = new TH2D*[17*2];  
  TH2D** heSigGHiR;             //ln(end signal) vs dcty, geometric mean
  heSigGHiR      = new TH2D*[17*2];  

//  TProfile** peSig;             //ln(end signal) vs dcty, SPW getSignal()
//  peSig      = new TProfile*[17*2];  
//  TProfile** peSigA;             //ln(end signal) vs dcty, SPW getSignal()
//  peSigA      = new TProfile*[17*2];  
//  TProfile** peSigG;             //ln(end signal) vs dcty, SPW getSignal()
//  peSigG      = new TProfile*[17*2];  
  
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

//    hname1   = TString::Format("%s%d%s%d","peSig",i,"_",j);
//    peSig[i*2+j] = new TProfile(hname1.Data(),htitle1.Data(), 50,-800,800, -2.,8.);
//    peSig[i*2+j]->BuildOptions(-1.5,1.,"");

      hname1   = TString::Format("%s%d%s%d","heSigA",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signalA) vs (DCTY-padCenter),  paddle ", i, " end ",j, ";(DCTY-padCenter) (mm);loge(signalA) (ln(MIPs))");
      heSigA[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 100,-2.,8.);

//    hname1   = TString::Format("%s%d%s%d","peSigA",i,"_",j);
//    peSigA[i*2+j] = new TProfile(hname1.Data(),htitle1.Data(), 50,-800,800, -2.,8.);
//    peSigA[i*2+j]->BuildOptions(-1.5,1.,"");

      hname1   = TString::Format("%s%d%s%d","heSigG",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signalG) vs (DCTY-padCenter),  paddle ", i, " end ",j ,";(DCTY-padCenter) (mm);loge(signalG) (ln(MIPs))");
      heSigG[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 100,-2.,8.);

//    hname1   = TString::Format("%s%d%s%d","peSigG",i,"_",j);
//    peSigG[i*2+j] = new TProfile(hname1.Data(),htitle1.Data(), 50,-800,800, -2.,8.);
//    peSigG[i*2+j]->BuildOptions(-1.5,1.,"");

      hname1   = TString::Format("%s%d%s%d","heSigHiR",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signal) vs (DCTY-padCenter) for R>2, paddle ", i, " end ",j ,";(DCTY-padCenter) (mm);loge(signal) (ln(MIPs))");
      heSigHiR[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(), 50,-800,800, 100,-2.,8.);

      hname1   = TString::Format("%s%d%s%d","heSigAHiR",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signalA) vs (DCTY-padCenter) for R>2,  paddle ", i, " end ",j, ";(DCTY-padCenter) (mm);loge(signalA) (ln(MIPs))");
      heSigAHiR[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 100,-2.,8.);

      hname1   = TString::Format("%s%d%s%d","heSigGHiR",i,"_",j);
      htitle1  = TString::Format("%s%d%s%d%s","loge(signalG) vs (DCTY-padCenter) for R>2,  paddle ", i, " end ",j ,";(DCTY-padCenter) (mm);loge(signalG) (ln(MIPs))");
      heSigGHiR[i*2+j] = new TH2D(hname1.Data(),htitle1.Data(),50,-800,800, 100,-2.,8.);

    }
  }
  
  //charge hists by correction level and signal construction type
  //All use sqrt(padSig)
  // Level 0: avg zenith and cos correction
  // L1: PLC and avg thickness correction
  // L2: thickness correction by position
  // L3: reponse correction
/*  
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
    hZall[i] = new TH1D(hname.Data(),htitle.Data(),400,0.,20.);

    hname   = TString::Format("%s%d","hZAall",i);
    htitle  = TString::Format("%s%d","Arithmetic Charge (t+b+bor), Level ", i);
    hZAall[i] = new TH1D(hname.Data(),htitle.Data(),400,0.,20.);

    hname   = TString::Format("%s%d","hZGall",i);
    htitle  = TString::Format("%s%d","Geometric charge (t+b+bor), Level ", i);
    hZGall[i] = new TH1D(hname.Data(),htitle.Data(),400,0.,20.);

    hname   = TString::Format("%s%d","hZtb",i);
    htitle  = TString::Format("%s%d","SPW Charge (t+b), Level ", i);
    hZtb[i] = new TH1D(hname.Data(),htitle.Data(),400,0.,20.);

    hname   = TString::Format("%s%d","hZAtb",i);
    htitle  = TString::Format("%s%d","Arithmetic Charge (t+b), Level ", i);
    hZAtb[i] = new TH1D(hname.Data(),htitle.Data(),400,0.,20.);

    hname   = TString::Format("%s%d","hZGtb",i);
    htitle  = TString::Format("%s%d","Geometric charge (t+b), Level ", i);
    hZGtb[i] = new TH1D(hname.Data(),htitle.Data(),400,0.,20.);
  }
  */
  auto hZSPWvG  = new TH2D("hZSPWvG","SPW vs G chg;Arith/Geom;SPW Z (G/G)",200,0,10,200,0,10);
  auto hZAvG    = new TH2D("hZAvG","A vs G chg;Arith/Geom;Arith/Arith",200,0,10,200,0,10);
  
  
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
  
  //beta and rig:
  auto hR   = new TH1D("hR" , "Rigidity;(GV)" , 1200, -30, 30);
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
  
  ////// TEST LIMITED SET OF ENTRIES: ////
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

    //Watch this trick: set up arrays to avoid lots of repeated code.
    HLX_Stg2_TOFHit* tofHits[3] = {topTOF, borTOF, botTOF};  //not copies, so use pointer ->
    
    /////////////////
    //Fixes for fuck ups:
    // v12.1.NT (new tracking): screwed up SignalSum. NB: THis takes care of normalization to mpv[8] but not thickness correction
    // Loop over FEEs in each paddle and fix:
    //HLX: 0-7 lower, 8-15 upper, 16 is bore
    // Ensure getFEEHit returns a reference (&), not a copy or a const
    
    
    for (int i=0; i<3; i++) //paddle loop
    {      
      HLX_Stg2_TOFPaddleEndHit* ends[2] = {&(tofHits[i]->fHitEast), &(tofHits[i]->fHitWest)};
      for (int k=0; k<2; k++) //end loop
      {
        int nFees     = ends[k]->getFEEHitCount();
        int paddleID  = ends[k]->fTag.getPaddleID();
        int endID     = ends[k]->fTag.getPaddleEndID();        
        for (int jk=0; jk<nFees; jk++)
        {
          int feeID = ends[k]->getFEEHit(jk).fTag.getFEEID();
          auto& feeHit = ends[k]->fFEEHits[jk];
          auto [laneIdx, chanIdx] = theGeom->fTOF.fGeoMap[{paddleID, endID, feeID}];
          int SCfeeID = laneIdx*8 + chanIdx; //HLX labels. tofcal12
          double feeSignal = feeHit.fSignalSumRaw;
          double chanGain = fCal.mpv[SCfeeID]/fCal.mpv[1];  //normalize to FEE1 (lane 0, ch 1: pad_id=16, end_id=1, fee_id=3) as of tofcal12.2.root
//          double chanGain = fCal.mpv[SCfeeID]/fCal.mpv[8];  //normalize to FEE8 (paddle 0, FEE0)
          double mymippk = 204.6; //constant -- close to the old SC value
          feeHit.fSignalSum = feeSignal/chanGain/mymippk;
//          feeHit.fSignalSum = feeSignal/chanGain/fCal.mippk;
//          feeHit.fSignalSum /= fCal.mippk;  //scale by global mippk to Z=1
//          cout << "FIX: " << feeSignal << "\t" << chanGain << "\t" << fCal.mippk << "\t" << feeHit.fSignalSum << endl;
        }
      }
    }
    
    //spot check:
//    for (int i=0; i<3; i++)
//    {
//      cout << tofHits[i]->getPaddleEndHit(HLX_Geometry_TOF::kEast).getFEEHit(0).fSignalSum << "\t" << tofHits[i]->getPaddleEndHit(HLX_Geometry_TOF::kWest).getFEEHit(0).fSignalSum << endl;
//      cout << tofHits[i]->getPaddleEndHit(HLX_Geometry_TOF::kEast).getFEEHit(1).fSignalSum << "\t" << tofHits[i]->getPaddleEndHit(HLX_Geometry_TOF::kWest).getFEEHit(1).fSignalSum << endl;
//    }
//    cout << botTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast).getFEEHit(0).fSignalSum << "\t" << botTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest).getFEEHit(0).fSignalSum << endl;
    
    // v12.GV
    //NOTE: Apply corrections specific to this version of skeletor output:
//    for (int i=0; i<3; i++) 
//    {
//      padSig[i] *= avgthick[padID[i]];
//      padSig [i] *= fCal.mpv[8]/fCal.mpv[0];   //correct for getting the wrong normalization fee.
//    } 
    
    
    // Fuck-up finishing
    //////////////////
    
    
    
    //HLX numbering:
    int toppadID = topTOF->fTag.getPaddleID();
    int borpadID = borTOF->fTag.getPaddleID();
    int botpadID = botTOF->fTag.getPaddleID();
    
    int padID[3] = {toppadID, borpadID, botpadID};    

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
//    double eSigSC[6]= { getTBPadEndSignalSC(tE), getTBPadEndSignalSC(tW), getBoreEndSignalSC(boreE), 
//                        getBoreEndSignalSC(boreW), getTBPadEndSignalSC(bE), getTBPadEndSignalSC(bW) };
    
  
//  cout << " ****NEW EVENT****" << endl;
//  cout << ik << " eN: " << eSig [0] << "\t" << eSig [1] << "\t" << eSig [2] << "\t" << eSig [3] << "\t" << eSig [4] << "\t" << eSig [5] << endl;
//  cout << ik << " eA: " << eSigA[0] << "\t" << eSigA[1] << "\t" << eSigA[2] << "\t" << eSigA[3] << "\t" << eSigA[4] << "\t" << eSigA[5] << endl;
//  cout << ik << " eG: " << eSigG[0] << "\t" << eSigG[1] << "\t" << eSigG[2] << "\t" << eSigG[3] << "\t" << eSigG[4] << "\t" << eSigG[5] << endl;
//  cout << ik <<" eSC: " <<eSigSC[0] << "\t" <<eSigSC[1] << "\t" <<eSigSC[2] << "\t" <<eSigSC[3] << "\t" <<eSigSC[4] << "\t" <<eSigSC[5] << endl;
//  cout << endl;
    
    //fill end-based signal histos
    for (int i=0; i<3; i++) //paddle
    {
      for (int j=0; j<2; j++) //end
      {
        int endID     =  endHit[i*2+j].fTag.getPaddleEndID();
        int globalEndID = padID[i]*2+endID;
        int k = i*2+j;
        if (eSig [k]>0) {heSig [globalEndID]->Fill(tofY[i]-padCenter[i], log(eSig [k])); /*peSig [globalEndID]->Fill(tofY[i]-padCenter[i], log(eSig [k]));*/}
        if (eSigA[k]>0) {heSigA[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigA[k])); /*peSigA[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigA[k]));*/}
        if (eSigG[k]>0) {heSigG[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigG[k])); /*peSigG[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigG[k]));*/}
        if (rig > 2 || runnum >= 15550 ) //HiR or Boff
        {
          if (eSig [k]>0) {heSigHiR [globalEndID]->Fill(tofY[i]-padCenter[i], log(eSig [k]));}
          if (eSigA[k]>0) {heSigAHiR[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigA[k]));}
          if (eSigG[k]>0) {heSigGHiR[globalEndID]->Fill(tofY[i]-padCenter[i], log(eSigG[k]));}
        }
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
    
    //signals: 
    double padSig[3]  = {topTOF->getSignal(), borTOF->getSignal(), botTOF->getSignal()}; //SPW method
    double padSigA[3] = {getPadMeanArithSignal(topTOF), getPadMeanArithSignal(borTOF), getPadMeanArithSignal(botTOF)}; //arith mean of paddle ends. End fees combined with arith mean.
    double padSigG[3] = {getPadMeanGeomSignal(topTOF),  getPadMeanGeomSignal(borTOF),  getPadMeanGeomSignal(botTOF)}; //geom mean of paddle ends. End fees combined with arith mean.
//    double padSigSC[3]= {getPadSignalSC(topTOF),  getPadSignalSC(borTOF),  getPadSignalSC(botTOF)}; //geom mean of paddle ends. End fees combined with arith mean.
    
//  cout << "Arith: " << getPadMeanArithSignal(topTOF) << endl;
//  cout << "Geom: " << getPadMeanGeomSignal(topTOF) << endl;
//  
//  cout << ik << " N: " << padSig [0] << "\t" << padSig [1] << "\t" << padSig [2] << endl;
//  cout << ik << " A: " << padSigA[0] << "\t" << padSigA[1] << "\t" << padSigA[2] << endl;
//  cout << ik << " G: " << padSigG[0] << "\t" << padSigG[1] << "\t" << padSigG[2] << endl;
////  cout << ik << "SC: " <<padSigSC[0] << "\t" <<padSigSC[1] << "\t" <<padSigSC[2] << endl;
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
          padt2D = &fCal.boret[0][0]; //bore
          thck[i] = thickness(i, padt2D, (tofX[i]-padx[i])/10. , (tofY[i]-pady[i])/10. );  // positions in cm, not mm (/10.); the BP does NOT need the 1/2 paddle shift
        }
      	hthck[padID[i]] -> Fill(thck[i]);        
      }
      tcorr[i] = 10./thck[i];   // this tcorr factor will *multiply* the signal size (thickness in cm)
//        tcorr[i] = avgthick[padID[i]] / thck[i];   // this tcorr factor will *multiply* the signal size
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
      padSig [j] *= tcorr[j];  // v12.1.NT: only accounts for thickness
      padSigA[j] *= tcorr[j];  // v12.1.NT: only accounts for thickness
      padSigG[j] *= tcorr[j];  // v12.1.NT: only accounts for thickness
      
//    cout << "After tcorr, paddle " << j << endl;
//    cout << ik << " N: " << padSig [0] << "\t" << padSig [1] << "\t" << padSig [2] << endl;
//    cout << ik << " A: " << padSigA[0] << "\t" << padSigA[1] << "\t" << padSigA[2] << endl;
//    cout << ik << " G: " << padSigG[0] << "\t" << padSigG[1] << "\t" << padSigG[2] << endl;
//  //  cout << ik << "SC: " <<padSigSC[0] << "\t" <<padSigSC[1] << "\t" <<padSigSC[2] << endl;
//    cout << endl;      
      
      padSig [j] /= padPLC[j];   //correct for pathlength
      padSigA[j] /= padPLC[j];   //correct for pathlength
      padSigG[j] /= padPLC[j];   //correct for pathlength

//    cout << "After PLC, paddle " << j << endl;
//    cout << ik << " N: " << padSig [0] << "\t" << padSig [1] << "\t" << padSig [2] << endl;
//    cout << ik << " A: " << padSigA[0] << "\t" << padSigA[1] << "\t" << padSigA[2] << endl;
//    cout << ik << " G: " << padSigG[0] << "\t" << padSigG[1] << "\t" << padSigG[2] << endl;
//  //  cout << ik << "SC: " <<padSigSC[0] << "\t" <<padSigSC[1] << "\t" <<padSigSC[2] << endl;
//    cout << endl;

      padSig [j] /= applyResponseCorr(padID[j],tofY[j]-padCenter[j]);   //correct for response along paddle length
      padSigA[j] /= applyResponseCorr(padID[j],tofY[j]-padCenter[j]);   //correct for response along paddle length
      padSigG[j] /= applyResponseCorr(padID[j],tofY[j]-padCenter[j]);   //correct for response along paddle length
      chg [j]     = (padSig [j]>0)? sqrt(padSig [j]) : 0;
      chgA[j]     = (padSigA[j]>0)? sqrt(padSigA[j]) : 0;
      chgG[j]     = (padSigG[j]>0)? sqrt(padSigG[j]) : 0;
      

      //fill histos:
      hChargeDivPos [padID[j]]->Fill(tofY[j]-padCenter[j],        chargeDivPos-padCenter[j]);  //consistent with calculation of timeDiffPos and chargeDivPos above and in TemplateAnaly
      hTimeDiffPos  [padID[j]]->Fill(tofY[j]-padCenter[j],        timeDiffPos-padCenter[j]);
      hSignal       [padID[j]]->Fill(tofY[j]-padCenter[j],        padSig[j]);
      hSignalA      [padID[j]]->Fill(tofY[j]-padCenter[j],        padSigA[j]);
      hSignalG      [padID[j]]->Fill(tofY[j]-padCenter[j],        padSigG[j]);
      hChg          [padID[j]]->Fill(tofY[j]-padCenter[j],        chg[j]);
      hChgA         [padID[j]]->Fill(tofY[j]-padCenter[j],        chgA[j]);
      hChgG         [padID[j]]->Fill(tofY[j]-padCenter[j],        chgG[j]);
      if (rig > 2 || runnum >= 15550)
      {
        hSignalHiR       [padID[j]]->Fill(tofY[j]-padCenter[j],        padSig[j]);
        hSignalAHiR      [padID[j]]->Fill(tofY[j]-padCenter[j],        padSigA[j]);
        hSignalGHiR      [padID[j]]->Fill(tofY[j]-padCenter[j],        padSigG[j]);
        hChgHiR          [padID[j]]->Fill(tofY[j]-padCenter[j],        chg[j]);
        hChgAHiR         [padID[j]]->Fill(tofY[j]-padCenter[j],        chgA[j]);
        hChgGHiR         [padID[j]]->Fill(tofY[j]-padCenter[j],        chgG[j]);
      }
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
    
    /////////////////////
    // AnalTOFv3 stuff:
    
/*    
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
    */

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
    
    /*
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
    */

    /*
    double tofZ2[3] = {pow(padSig[0], expo1), pow(padSig[1], expo1), pow(padSig[2], expo1)};
    double allZ2    = (tofZ2[0]+tofZ2[1]+tofZ2[2])/3.;
    double tbZ2     = (tofZ2[0]+tofZ2[2])/2.;

    double allZSC2      = pow((tofZ2[0]+tofZ2[1]+tofZ2[2])/3.,expo2);   // SC  original
    double allZSCcorr2  = (allZSC2 - 0.2499) / 0.8385;
    double tbZSC2       = pow((tofZ2[0]+tofZ2[2])/2.,expo2);   // SC  original
*/
//fill signal and charge vs dcty hists at this point:
    /*
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
*/    

    /*
    if (abs(tofZ2[0]-tofZ2[2])<0.4 ) { // apply top/bottom tight charge agreement
      hallZ2 -> Fill(allZ2);
      htbZ2 -> Fill(tbZ2);

      hallZSC2 -> Fill(allZSC2);
      hallZSCcorr2 -> Fill(allZSCcorr2);
      htbZSC2 -> Fill(tbZSC2);
    }
*/
    /*
    //THIRD: response corrected version
    //NB: THESE SHOULD BE DIFFERENT PARAMETERS FOR EACH METHOD! RIGHT NOW: SPW method parameters
    for (int i=0; i<3; i++)  
    {
      padSig [i] /= applyResponseCorr(padID[i],tofY[i]-padCenter[i]);   //correct for response along paddle length based on hSignal vs DCTY
      padSigA[i] /= applyResponseCorr(padID[i],tofY[i]-padCenter[i]);   //correct for response along paddle length
      padSigG[i] /= applyResponseCorr(padID[i],tofY[i]-padCenter[i]);   //correct for response along paddle length
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
*/
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
    /////////////////////////////////
    
    double zall = (padSigG[0] + padSigA[1] + padSigG[2])/3.; //top, bottom: arith/geometric, bore arith/arith
    double ztb  = (padSigG[0] + padSigG[2])/2.;
    
    
    //fill comparison plots
    hZSPWvG->Fill(chgG[2], chg[2]);
    hZAvG->Fill(chgG[2], chgA[2]);
    
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
  
    double myrig = (std::isfinite(rig))? rig : 0;
    //    if (abs(topSig/tvars.topPLC*tcorr[0]-botSig/tvars.botPLC*tcorr[2])<0.2 ) { // apply top/bottom tight charge agreement
    if (abs(chg[0]-chg[2])<0.2 ) { // apply top/bottom tight charge agreement
      if (DCTChi2<2.&&hitCount>60.&&inPaddle)   // DCT cuts, including whether track hits the right paddles
      {
//	      if (beta>0.2&&beta<1.01) {hmass->Fill(allZSC2*rig*sqrt(pow(1./beta,2.)-1.));}   // mass*c^2 in GV units 
	      if (beta>0.2&&beta<1.01) {hmass->Fill(zall*myrig*sqrt(pow(1./beta,2.)-1.));}   // mass*c^2 in GV units 
      }
    }
    
//    hZvR->Fill(rig, allZSC2); //rigidity on x; Zave on y
//    hZvbeta->Fill(beta,allZSC2);
    hR->Fill(myrig); //rigidity
    hZvR->Fill(myrig, zall); //rigidity on x; Zave on y
    hZvbeta->Fill(beta,zall);
    //    hBvT->Fill(pow(padSig[0],expo2), pow(padSig[2],expo2)); //top signal on x; bottom signal on y
    hBvT->Fill(chgG[0],chgG[2]); //top signal on x; bottom signal on y

  } //end of event loop
  
  cout << "End of event loop" << endl;
  opf->Write();
  opf->Close();


}

//
// SC 10/12/2026   A process to run multiple skeletor files chained together using a variant of AnalTOFSkel.C
// renamed AnalTOFSkelTree.C .
//
// ROOT usage:
//   .L AnalTOFSkelTree.C
//   .L AnalTOFSkelChain.C
//   AnalTOFSkelList("filelist.txt", "all.AnalTOFSkel.root");
//
// where filelist.txt is the list of run numbers, such as 2024_Kiruna_Float_Boff.txt, to chain and process, and all.AnalTOFSkel.root will be
// the output file.

//make the chain global so you can access it from the command line:
TChain chain("toft");

void AnalTOFSkelList(const char *listName = "filelist.txt",
                    const char *outName  = "combined.AnalTOFSkel.root")
{
  std::string dir = "/data/skeletor/";
  std::string fnameStem = "skeletor.";
  std::string fnameExt = ".root";
  
  std::ifstream list(listName);
  if (!list) {
    std::cerr << "Cannot open file list: " << listName << "\n";
    return;
  }

//  TChain chain("toft");
  std::string line;
  int nFiles = 0;

  while (std::getline(list, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::string fname = dir+fnameStem+line.c_str()+fnameExt;
//    if (chain.Add(line.c_str()) <= 0) {
//      std::cerr << "Could not add: " << line << "\n";
    if (chain.Add(fname.c_str()) <= 0) {
      std::cerr << "Could not add: " << fname << "\n";
      continue;
    }
    ++nFiles;
  }

  if (nFiles == 0 || chain.GetEntries() == 0) {
    std::cerr << "No usable files or entries found\n";
    return;
  }

  std::cout << "Added " << nFiles << " files; entries = "
            << chain.GetEntries() << "\n";

  AnalTOFSkelTree(&chain, outName);
}


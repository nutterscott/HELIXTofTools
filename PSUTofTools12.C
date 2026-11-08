//Current status:
// Reads in tofcal12 as created by SC to repair paddle thickness records and establish 

//In this file:
//bool createSCTOFCalibs(string opfName = kOutputTofcalName) //only called once. Creates ttree for writing new tofcals
//void getSCTOFCalibsForRun(int runnum) //Getting fCal for run runnum
//bool loadSCTOFCalibsv12(string ipfName = kInputTofcalName) // Open INPUT SCTOF parameter file v12 and after (root TTree organized by run number). Only called once.
//void fillFitSlopeIntercepts(int runNum) //fill the slopes and intercepts in the new ttree by reading from files. Input is run number.
//void modifyTofcal() //one time pass through all entries and writing out to a new tofcal tree with modified branches and/or values
//void adjustSCTOFCalibsForRun(PSUTOFCalibs oldCal) //do adjustment of bad parameters inside fCal for run

//Deprecated:
//bool loadSCTOFCalibs(string ipfName = kInputTofcalName) // Open INPUT SCTOF parameter file v11 and before (root TTree organized by run number). Only called once.
//void switchfCaltoHLXNumbering()   //one time call to switch SC root tree into HLX numbering scheme for this run.

//can make a new method to read in histos from AnalTOFSkel.C, make fits, and further modify parameters. If that would be something....

//Global shit:
/////  SET THESE FOR READING AND TRANSLATION!! ////
//Set up input and output tofcal ttree
string kInputTofcalName = "tofcal12.2.root";
string kOutputTofcalName = "tofcal12.3.root";

string kMergeInputTofcalName = "tofcal12b.root";


//TREE VARIABLES: Saved in structure for eases of copying
struct PSUTOFCalibs {
  Double_t pup[72][7],epup[72][7],pdn[72][7],epdn[72][7],mippk;
  Int_t runnum,nument;
  Double_t runtime,pedm[72],rms[72],mpv[72],empv[72];
  Double_t etoeshift[17];
  Double_t peakuu[40], widthuu[40], peakud[40], widthud[40], peakdu[40], widthdu[40], peakdd[40], widthdd[40];
  //v11:
  Double_t interT[64],einterT[64],slopeT[64],eslopeT[64],interB[64],einterB[64],slopeB[64],eslopeB[64];
  Double_t interTW[64],einterTW[64],slopeTW[64],eslopeTW[64];
  Int_t padnum[17];
  Double_t topx[8][7],topy[8][17],topt[8][7][17],botx[8][7],boty[8][17],bott[8][7][17]; //,borex[1][11],borey[1][11],boret[1][11][11];
  
  //new branches: tofcal12
  Double_t ZDm[17], ZDb[17];  //paddle slope, intercept, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  Double_t TDm[17], TDb[17];  //paddle slope, intercept, from DCTY vs fTimeDiff:      fTimeDiffPos      = m*fTimeDiff + b
  
  //tofcal12 with HELIX numbering scheme:
  Double_t padx[16][7], pady[16][17], padt[16][7][17]; //combine top and bottom paddles
  Double_t borex[1][11],borey[1][11],boret[1][11][11];
};

TFile*        fInFile;          //INPUT TFile with Calibs
TTree*        fCalTree{};      //INPUT calibration parameters 
TTree*        fMergeCalTree{};      //INPUT merging calibration parameters 
PSUTOFCalibs  fRtCal;         //global struct with INPUT variables that ROOT points to
PSUTOFCalibs  fMergeCal;      //global struct for merging file INPUT calib values
PSUTOFCalibs  fCal;           //global struct for modified OUTPUT calib values

TFile*        fOutFile;       //OUPUT TFile with Calibs
TTree*        fOUTTOFcal{};   //OUTPUT calibration parameters

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//a function we need:
vector<string> getFileListVector(string filename ){
  vector<string> fileListVector;
  std::ifstream file;
  file.open(filename);
  if (!file.good()) {
    cout << "Cannot open file: " << filename << endl;
    exit(-1);
  } else {
    string line;
    int i = 0;
    while (!file.eof()) {
      getline(file, line);
      if (line[0] != '#'  && line != "") {
        fileListVector.push_back(line);
        i++;
      }
    }
  }
  file.close();
  return fileListVector;
}

void fillNewMPVValues(int runNum) //fill the mpv and empv values in the new ttree by reading from txt files. Input is run number.
{
  string baseName="./FEESigs/FEESigs.";

//  cout << "Filling fits for run " << runNum << endl;

  string fname = baseName+std::to_string(runNum)+".txt";
//  cout << "Reading " << fname << endl;
  ifstream input;
  input.open(fname.c_str());
  if (!input.is_open())
  {
    cout << "//////////////////////////////////////////// " << endl;
    cout << "File cannot be found: " << fname << endl;
    cout << "//////////////////////////////////////////// " << endl;
//    exit(-1);
    return;
  }
  
  int feeID;
  double pk, pkerr, chi2;
  
  fCal.runnum = runNum; //needed to insert new values

  for (int i=0; i<72; i++) //fees
  {
    input >> feeID >> pk >> pkerr >> chi2;
    fCal.mpv[feeID] = pk; 
    fCal.empv[feeID] = pkerr; 
  }  
  input.close();  
}

void fillFitSlopeIntercepts(int runNum) //fill the slopes and intercepts in the new ttree by reading from files. Input is run number.
{
//  TString baseName="/data/skeletor/MakeSkeletorPlotsBatch.";
  TString baseName="./PosFits/PosFits.";

//  cout << "Filling fits for run " << runNum << endl;

  int padnum;
  //fChargeDivFracPos:
//  TString fname = baseName+std::to_string(runNum)+".ChargeDivFrac.txt";
  TString fname = baseName+runNum+".ChargeDivFrac.txt";
//  cout << "Reading " << fname << endl;
  ifstream input;
  input.open(fname.Data());
  if (!input.is_open())
  {
    cout << "File cannot be found: " << fname << endl;
//    exit(-1);
    return;
  }
  
  for (int i=0; i<17; i++)
  {
    input >> padnum >> fCal.ZDb[i] >> fCal.ZDm[i]; 
//    cout << fCal.ZDb[i] <<  "  "  << fCal.ZDm[i] << endl; 
  }  
  input.close();
  
  //fTimeDiffPos:
//  fname = baseName+std::to_string(runNum)+".TimeDiff.txt";
  fname = baseName+runNum+".TimeDiff.txt";
//  cout << "Reading " << fname << endl;
  input.open(fname.Data());
  if (!input.is_open())
  {
    cout << "File cannot be found: " << fname << endl;
//    exit(-1);
    return;
  }
  
  for (int i=0; i<17; i++)
  {
    input >> padnum >> fCal.TDb[i] >> fCal.TDm[i]; 
//    cout << fCal.TDb[i] << "  "  << fCal.TDm[i] << endl; 
  }
  input.close();
  
}

bool createSCTOFCalibs(string opfName = kOutputTofcalName) //only called once. Creates ttree for writing new tofcals
{
  ////////  NB: put checks here to make sure file is not already open?  /////////  
  

  fOutFile = new TFile(opfName.c_str(),"RECREATE");
  if (!fOutFile->IsOpen()) {
    cout << "createSCTOFCalibs => could not open file: " << opfName << endl;
    exit(-1);
  }

  cout << "createSCTOFCalibs : creating " << opfName << endl;

  fOUTTOFcal = new TTree("tofcal", opfName.c_str());

  // global ToF parameters, same for any run number
  fOUTTOFcal->Branch("pup",fCal.pup,"pup[72][7]/D");  // for each of 72 FEEs, up ramp parameters (7 of them) for GV function
  //  up ramp:   [0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fOUTTOFcal->Branch("epup",fCal.epup,"epup[72][7]/D");  // for each of 72 FEEs, error on up ramp parameters (7 of them) for GV function
  fOUTTOFcal->Branch("pdn",fCal.pdn,"pdn[72][7]/D");  // for each of 72 FEEs, down ramp parameters (7 of them) for GV function
  //  down ramp:   -[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fOUTTOFcal->Branch("epdn",fCal.epdn,"epdn[72][7]/D");  // for each of 72 FEEs, error on down ramp parameters (7 of them) for GV function
  fOUTTOFcal->Branch("mippk",&fCal.mippk,"mippk/D");  // MIP peak, average most probable Landau fit value, used for determining the charge Z from the pulse integral between samples 11 and 25 (inclusive)
//  fOUTTOFcal->Branch("interTW",&fCal.interTW,"interTW[64]/D");  // intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
//  fOUTTOFcal->Branch("slopeTW",&fCal.slopeTW,"slopeTW[64]/D");  // slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
//  fOUTTOFcal->Branch("einterTW",&fCal.einterTW,"einterTW[64]/D");  // error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
//  fOUTTOFcal->Branch("eslopeTW",&fCal.eslopeTW,"eslopeTW[64]/D");  // error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)	
//  fOUTTOFcal->Branch("padnum",&fCal.padnum,"padnum[17]/I");  // paddle number, HELIX style (0-7 bot, 8-15 top, 16 bore)
  //below: tofcal11
//  fOUTTOFcal->Branch("topx",&fCal.topx,"topx[8][7]/D");  // x values where we have thickness measurements for top paddles
//  fOUTTOFcal->Branch("topy",&fCal.topy,"topy[8][17]/D");  // y values where we have thickness measurements for top paddles
//  fOUTTOFcal->Branch("topt",&fCal.topt,"topt[8][7][17]/D");  // thickness (mm) of top paddle at topx,topy
//  fOUTTOFcal->Branch("botx",&fCal.botx,"botx[8][7]/D");  // x values where we have thickness measurements for bottom paddles
//  fOUTTOFcal->Branch("boty",&fCal.boty,"boty[8][17]/D");  // y values where we have thickness measurements for bottom paddles
//  fOUTTOFcal->Branch("bott",&fCal.bott,"bott[8][7][17]/D");  // thickness (mm) of bottom paddle at topx,topy

  // ToF calibration parameters that change with run number

  fOUTTOFcal->Branch("runnum",&fCal.runnum,"runnum/I");  // run number
  fOUTTOFcal->Branch("nument",&fCal.nument,"nument/I");  // number of entries in that run
  fOUTTOFcal->Branch("runtime",&fCal.runtime,"runtime/D");  // time of the first event in the run, in days using HXA_TimeStamp::DtD()
  fOUTTOFcal->Branch("pedm",&fCal.pedm,"pedm[72]/D");  // for each of 72 FEEs, mean value of the charge pedestal for this run
  fOUTTOFcal->Branch("rms",&fCal.rms,"rms[72]/D");  // for each of 72 FEEs, rms of the charge pedestal peak for this run
  fOUTTOFcal->Branch("mpv",&fCal.mpv,"mpv[72]/D");  // for each of 72 FEEs, most probable charge value of the pedestal-subtracted signal for this run
  fOUTTOFcal->Branch("empv",&fCal.empv,"empv[72]/D");  // for each of 72 FEEs, uncertainty on the mpv value of the pedestal-subtracted charge signal for this run
  fOUTTOFcal->Branch("etoeshift",&fCal.etoeshift,"etoeshift[17]/D");  // for each of 17 paddles, average time shift (in samples) between the two ends for this run (top: 0-7, bore: 8, bottom: 9-16)
  // for the following arrays of size 40: top are 0-15, bore are 16-23, bottom are 24-39;
  // top order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  // bore order is 00E, 01E, 12E, 23E, 00W, 01W, 12W, 23W   where 00E and 00W can be ignored 
  // bottom order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  fOUTTOFcal->Branch("peakuu",  &fCal.peakuu ,"peakuu[40]/D" );  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+up
  fOUTTOFcal->Branch("widthuu", &fCal.widthuu,"widthuu[40]/D" );  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+up
  fOUTTOFcal->Branch("peakud",  &fCal.peakud ,"peakud[40]/D" );  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+dn
  fOUTTOFcal->Branch("widthud", &fCal.widthud,"widthud[40]/D" );  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+dn
  fOUTTOFcal->Branch("peakdu",  &fCal.peakdu ,"peakdu[40]/D" );  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+up
  fOUTTOFcal->Branch("widthdu", &fCal.widthdu,"widthdu[40]/D" );  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+up
  fOUTTOFcal->Branch("peakdd",  &fCal.peakdd ,"peakdd[40]/D" );  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+dn
  fOUTTOFcal->Branch("widthdd", &fCal.widthdd,"widthdd[40]/D" );  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+dn
  
  //new branches from tofcal12:
  //Position from fChargeDivFrac by run:
  fOUTTOFcal->Branch("ZDm",  &fCal.ZDm ,"ZDm[17]/D" );  //paddle slope, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  fOUTTOFcal->Branch("ZDb",  &fCal.ZDb ,"ZDb[17]/D" );  //paddle intercept, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  //Position from fTimeDiff by run:
  fOUTTOFcal->Branch("TDm",  &fCal.TDm ,"TDm[17]/D" );  //paddle slope, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b
  fOUTTOFcal->Branch("TDb",  &fCal.TDb ,"TDb[17]/D" );  //paddle intercept, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b

  //better scheme for top and bottom paddle thickness measurements. NB: EVENTUALLY MOVE THIS TO SINGLE READ IN SEPARATE OBJECT IN tofcal root FILE, NOT IN TREE!
  fOUTTOFcal->Branch("padx",&fCal.padx,"padx[16][7]/D");  // x values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fOUTTOFcal->Branch("pady",&fCal.pady,"pady[16][17]/D");  // y values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fOUTTOFcal->Branch("padt",&fCal.padt,"padt[16][7][17]/D");  // thickness (mm) of top paddle at padx,pady  
  //keep bore stuff:
  fOUTTOFcal->Branch("borex",&fCal.borex,"borex[1][11]/D");  // x values where we have thickness measurements for bore paddle
  fOUTTOFcal->Branch("borey",&fCal.borey,"borey[1][11]/D");  // y values where we have thickness measurements for bore paddle
  fOUTTOFcal->Branch("boret",&fCal.boret,"boret[1][11][11]/D");  // thickness (mm) of bore paddle at topx,topy

  return kTRUE;

}

////////////////////////////////////////////////////////////
void switchfCaltoHLXNumbering()   //one time call to switch SC root tree into HLX numbering scheme for this run.
{
  
  /// insane translation array from lane/chan to SC global FEEID:
  const int fSCFeeIDMap[72] = {
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
  
  int padMap[17] = {9,10,11,12,13,14,15,16,0,1,2,3,4,5,6,7,8}; //padMap[HLX paddle ID] = SC paddle ID
  //HLX: 0-7 lower, 8-15 upper, 16 is bore  
  //SC: 0-7 top, 8=bore, 9-16 bot
  
  //Copy fCal into new structure:
  PSUTOFCalibs SCCal = fCal;
  
  //So painful....
  //loop through each set of variables by paddleID and reset fCal values:
//  auto [laneIdx, chanIdx] = getGeom()->fTOF.fGeoMap[{paddleID, endID,  feeID}];
//  int globFeeID = laneIdx*8+chanIdx;
//  int scFeeID = fSCFeeIDMap[globFeeID];
  
  for (int ipad=0; ipad<17; ipad++) //HLX paddleID
  {
    fCal.etoeshift[ipad]    =SCCal.etoeshift[padMap[ipad]] ;
//    fCal.padnum[ipad]       =SCCal.padnum[padMap[ipad]]    ;
    //NOTE: ZDm, ZDb, and TDm, TDb are already in HELIX numbering scheme since they were added later.
  }
  for (int ifee=0; ifee<72; ifee++) //HLX fee ID
  {
//      auto [laneIdx, chanIdx] = theGeom->fTOF.fGeoMap[{paddleID, endID,  feeID}];
//      int globFeeID = laneIdx*8+chanIdx;
    int scFeeID = fSCFeeIDMap[ifee];
    fCal.pedm[ifee]         =SCCal.pedm[scFeeID]      ;
    fCal.rms[ifee]          =SCCal.rms[scFeeID]       ;
    fCal.mpv[ifee]          =SCCal.mpv[scFeeID]       ;
    fCal.empv[ifee]         =SCCal.empv[scFeeID]      ;
    for (int ip=0; ip<7; ip++)
    {
      fCal.pup[ifee][ip]       =SCCal.pup[scFeeID][ip]    ;
      fCal.epup[ifee][ip]      =SCCal.epup[scFeeID][ip]   ;
      fCal.pdn[ifee][ip]       =SCCal.pdn[scFeeID][ip]    ;
      fCal.epdn[ifee][ip]      =SCCal.epdn[scFeeID][ip]   ;
    }
  }

//    fCal.mippk            =SCCal.mippk         ;
//    fCal.runnum           =SCCal.runnum        ;
//    fCal.nument           =SCCal.nument        ;
//    fCal.runtime          =SCCal.runtime       ;
  
//  cout << "Converting up/down pairs. " << endl;
//  cout << "FEE   scID   HLXID" << endl;
  
  //Copy/translate end-based data in 40 ends in peakuu etc arrays
  // Numbering scheme: 
  // if (paddleID<=15) {scID = paddleID*2 + paddleEndID;}   // so simple now
  // else if (paddleID==16) {scID = 32 + paddleEndID*4 + FEEID;}        // BP, still very simple

  for (int HLXID=0; HLXID<=39; HLXID++)
  {
    int scID = (HLXID <= 15)? HLXID+24 : HLXID-16 ; 
    fCal.peakuu[HLXID]       =SCCal.peakuu[scID]    ;
    fCal.widthuu[HLXID]      =SCCal.widthuu[scID]   ;
    fCal.peakud[HLXID]       =SCCal.peakud[scID]    ;
    fCal.widthud[HLXID]      =SCCal.widthud[scID]   ;
    fCal.peakdu[HLXID]       =SCCal.peakdu[scID]    ;
    fCal.widthdu[HLXID]      =SCCal.widthdu[scID]   ;
    fCal.peakdd[HLXID]       =SCCal.peakdd[scID]    ;
    fCal.widthdd[HLXID]      =SCCal.widthdd[scID]   ;
  }

  //i don't know how to interpret these:
  
//  for (int (i=0; i<64; i++)
//  {
//    fCal.interT[64]       =SCCal.interT[64]    ;
//    fCal.einterT[64]      =SCCal.einterT[64]   ;
//    fCal.slopeT[64]       =SCCal.slopeT[64]    ;
//    fCal.eslopeT[64]      =SCCal.eslopeT[64]   ;
//    fCal.interB[64]       =SCCal.interB[64]    ;
//    fCal.einterB[64]      =SCCal.einterB[64]   ;
//    fCal.slopeB[64]       =SCCal.slopeB[64]    ;
//    fCal.eslopeB[64]      =SCCal.eslopeB[64]   ;
//    fCal.interTW[64]      =SCCal.interTW[64]   ;
//    fCal.einterTW[64]     =SCCal.einterTW[64]  ;
//    fCal.slopeTW[64]      =SCCal.slopeTW[64]   ;
//    fCal.eslopeTW[64]     =SCCal.eslopeTW[64]  ;
//  }

  //switch thickness maps paddle numbering. Use single array for all top and bottom paddles. Make life easier!
  for (int i=0; i<16; i++) //HLX numbers
  {
    for (int ix=0; ix<7; ix++)
    {
      fCal.padx[i][ix]       = (i<8)? SCCal.botx[padMap[i]-9][ix] : SCCal.topx[padMap[i]][ix]    ; //bottom for i<8, note padMap[0] = 9
      for (int iy=0; iy<17; iy++) fCal.padt[i][ix][iy]   = (i<8)? SCCal.bott[padMap[i]-9][ix][iy] : SCCal.topt[padMap[i]][ix][iy];
    }

    for (int iy=0; iy<17; iy++)   fCal.pady[i][iy]       = (i<8)? SCCal.boty[padMap[i]-9][iy] : SCCal.topy[padMap[i]][iy]    ; //bottom for i<8

  }
  
  //bore map does not change? No need to copy since fCal is copy of fRtCal (original read)
// DO THESE CHANGE WITH PADDLE NUMBER CHANGES? I DON'T THINK SO!
//  fCal.borex[1][11]     =SCCal.borex[1][11]  ;
//  fCal.borey[1][11]     =SCCal.borey[1][11]  ;
//  fCal.boret[1][11][11] =SCCal.boret[1][11][11];
    
}

///////////////////////////////////////////////////////////////
//do adjustment of bad parameters inside fCal for run
void adjustSCTOFCalibsForRun(PSUTOFCalibs oldCal) 
{

  ////////  NOTE: IS THE FOLLOWING CHECK STILL NEEDED? /////////////
  for( Int_t j = 0; j < 40; j++) { // check for pathological peak?? corrections, and reset to last good known if needed
    if (fCal.widthuu[j]<0 || fCal.widthuu[j]>0.04) {fCal.peakuu[j] = oldCal.peakuu[j]; }
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

////////////////////////////////////////////////////////
void getSCTOFCalibsForRun(int runnum, bool doAdjustment = true) //Getting fCal for run runnum
{
  static int runcurrent = 0;
//  cout << "getSCTOFCalibsForRun: Getting fCal for run " << runnum << ", runcurrent =  " << runcurrent << endl;
  if (runnum == runcurrent) return;

  //copy previous values
  PSUTOFCalibs oldCal = fCal;

//  cout << "getSCTOFCalibsForRun: fCalTree->GetEntries() " << fCalTree->GetEntries() << endl;

  // Now that we know which run we're dealing with, scroll down fCalTree.root to find the relevant calibration parameters
  for (int ical = 0; ical < fCalTree->GetEntries(); ical++) 
  {	
    fCalTree->GetEntry( ical ); // load up the calibration parameters structure for this run
//    cout << "getSCTOFCalibsForRun: Checking entry: " << ical << ", fRtCal.runnum: " << fRtCal.runnum << ", runnum: " << runnum << endl;
    if (fRtCal.runnum >= runnum) {
//      cout << "getSCTOFCalibsForRun: Setting fCalTree entry ical " << ical << " with fRtCal.runnum = " << fRtCal.runnum << endl;
      runcurrent = fRtCal.runnum;
      fCal = fRtCal;  //copy INPUT to OUTPUT
      break; // now we have the parameters for this run, or at least for the next higher run number in the database
    }
  }
  
//  cout << "getSCTOFCalibsForRun: found entry " << runcurrent << endl;

//adjusting particular values inside fCal for run
  if (doAdjustment) adjustSCTOFCalibsForRun( oldCal ); //adjusting particular values inside fCal for run

}

///////////////////////////////////////////////////////////////
// Open INPUT SCTOF parameter file (root TTree organized by run number)
bool loadSCTOFCalibsv12(string ipfName = kInputTofcalName) //only called once V12 and after!
{
  ////////  NB: put checks here to make sure file is not already open?  /////////  

  if (!std::filesystem::exists(ipfName.c_str())) {
    cout << "loadSCTOFCalibsv12 => file not found: " << ipfName << endl;
    exit(-1);
  }

  fInFile = new TFile(ipfName.c_str());
  if (!fInFile->IsOpen()) {
    cout << "loadSCTOFCalibsv12 => could not open file: " << ipfName << endl;
    exit(-1);
  }

  cout << "loadSCTOFCalibsv12 : load from " << ipfName << endl;

  fCalTree = (TTree*)fInFile->Get("tofcal");

  // global ToF parameters, same for any run number
  fCalTree->SetBranchAddress("pup",&fRtCal.pup);  // for each of 72 FEEs, up ramp parameters (7 of them) for GV function
  //  up ramp:   [0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fCalTree->SetBranchAddress("epup",&fRtCal.epup);  // for each of 72 FEEs, error on up ramp parameters (7 of them) for GV function
  fCalTree->SetBranchAddress("pdn",&fRtCal.pdn);  // for each of 72 FEEs, down ramp parameters (7 of them) for GV function
  //  down ramp:   -[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fCalTree->SetBranchAddress("epdn",&fRtCal.epdn);  // for each of 72 FEEs, error on down ramp parameters (7 of them) for GV function
  fCalTree->SetBranchAddress("mippk",&fRtCal.mippk);  // MIP peak, average most probable Landau fit value, used for determining the charge Z from the pulse integral between samples 11 and 25 (inclusive)
//  fCalTree->SetBranchAddress("interTW",&fRtCal.interTW);  // intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
//  fCalTree->SetBranchAddress("slopeTW",&fRtCal.slopeTW);  // slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
//  fCalTree->SetBranchAddress("einterTW",&fRtCal.einterTW);  // error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
//  fCalTree->SetBranchAddress("eslopeTW",&fRtCal.eslopeTW);  // error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)	
//  fCalTree->SetBranchAddress("padnum",&fRtCal.padnum);  // paddle number, HELIX style (0-7 bot, 8-15 top, 16 bore)
//  fCalTree->SetBranchAddress("topx",&fRtCal.topx);  // x values where we have thickness measurements for top paddles
//  fCalTree->SetBranchAddress("topy",&fRtCal.topy);  // y values where we have thickness measurements for top paddles
//  fCalTree->SetBranchAddress("topt",&fRtCal.topt);  // thickness (mm) of top paddle at topx,topy
//  fCalTree->SetBranchAddress("botx",&fRtCal.botx);  // x values where we have thickness measurements for bottom paddles
//  fCalTree->SetBranchAddress("boty",&fRtCal.boty);  // y values where we have thickness measurements for bottom paddles
//  fCalTree->SetBranchAddress("bott",&fRtCal.bott);  // thickness (mm) of bottom paddle at topx,topy

  // ToF calibration parameters that change with run number

  fCalTree->SetBranchAddress("runnum",&fRtCal.runnum);  // run number
  fCalTree->SetBranchAddress("nument",&fRtCal.nument);  // number of entries in that run
  fCalTree->SetBranchAddress("runtime",&fRtCal.runtime);  // time of the first event in the run, in days since 1st event in run 14753
  fCalTree->SetBranchAddress("pedm",&fRtCal.pedm);  // for each of 72 FEEs, mean value of the charge pedestal for this run
  fCalTree->SetBranchAddress("rms",&fRtCal.rms);  // for each of 72 FEEs, rms of the charge pedestal peak for this run
  fCalTree->SetBranchAddress("mpv",&fRtCal.mpv);  // for each of 72 FEEs, most probable charge value of the pedestal-subtracted signal for this run
  fCalTree->SetBranchAddress("empv",&fRtCal.empv);  // for each of 72 FEEs, uncertainty on the mpv value of the pedestal-subtracted charge signal for this run
  fCalTree->SetBranchAddress("etoeshift",&fRtCal.etoeshift);  // for each of 17 paddles, average time shift (in samples) between the two ends for this run (top: 0-7, bore: 8, bottom: 9-16)
  // for the following arrays of size 40: top are 0-15, bore are 16-23, bottom are 24-39;
  // top order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  // bore order is 00E, 01E, 12E, 23E, 00W, 01W, 12W, 23W   where 00E and 00W can be ignored 
  // bottom order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  fCalTree->SetBranchAddress("peakuu",&fRtCal.peakuu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+up
  fCalTree->SetBranchAddress("widthuu",&fRtCal.widthuu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+up
  fCalTree->SetBranchAddress("peakud",&fRtCal.peakud);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+dn
  fCalTree->SetBranchAddress("widthud",&fRtCal.widthud);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+dn
  fCalTree->SetBranchAddress("peakdu",&fRtCal.peakdu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+up
  fCalTree->SetBranchAddress("widthdu",&fRtCal.widthdu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+up
  fCalTree->SetBranchAddress("peakdd",&fRtCal.peakdd);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+dn
  fCalTree->SetBranchAddress("widthdd",&fRtCal.widthdd);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+dn
    
  //new branches from tofcal12:
  //Position from fChargeDivFrac by run:
  fCalTree->SetBranchAddress("ZDm",  fRtCal.ZDm );  //paddle slope, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  fCalTree->SetBranchAddress("ZDb",  fRtCal.ZDb );  //paddle intercept, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  //Position from fTimeDiff by run:
  fCalTree->SetBranchAddress("TDm",  fRtCal.TDm );  //paddle slope, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b
  fCalTree->SetBranchAddress("TDb",  fRtCal.TDb );  //paddle intercept, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b

  //better scheme for top and bottom paddle thickness measurements. NB: EVENTUALLY MOVE THIS TO SINGLE READ IN SEPARATE OBJECT IN tofcal root FILE, NOT IN TREE!
  fCalTree->SetBranchAddress("padx", fRtCal.padx);  // x values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fCalTree->SetBranchAddress("pady", fRtCal.pady);  // y values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fCalTree->SetBranchAddress("padt", fRtCal.padt);  // thickness (mm) of top paddle at padx,pady  
  //keep bore stuff:
  fCalTree->SetBranchAddress("borex", fRtCal.borex);  // x values where we have thickness measurements for bore paddle
  fCalTree->SetBranchAddress("borey", fRtCal.borey);  // y values where we have thickness measurements for bore paddle
  fCalTree->SetBranchAddress("boret", fRtCal.boret);  // thickness (mm) of bore paddle at topx,topy

  //make copy in case of mods:
  fCalTree->GetEntry( 0 ); // load up the calibration parameters structure for first run
  fCal = fRtCal;

  getSCTOFCalibsForRun( 14765 ); //set default values
  
  return kTRUE;

}

// Open INPUT SCTOF parameter file (root TTree organized by run number)
bool loadSCTOFCalMergeFilev12(string ipfName = kMergeInputTofcalName) //only called once V12 and after!
{
  ////////  NB: put checks here to make sure file is not already open?  /////////  

  if (!std::filesystem::exists(ipfName.c_str())) {
    cout << "loadSCTOFCalMergeFilev12 => file not found: " << ipfName << endl;
    exit(-1);
  }

  fInFile = new TFile(ipfName.c_str());
  if (!fInFile->IsOpen()) {
    cout << "loadSCTOFCalMergeFilev12 => could not open file: " << ipfName << endl;
    exit(-1);
  }

  cout << "loadSCTOFCalMergeFilev12 : load from " << ipfName << endl;

  fMergeCalTree = (TTree*)fInFile->Get("tofcal");

  // global ToF parameters, same for any run number
  fMergeCalTree->SetBranchAddress("pup",&fMergeCal.pup);  // for each of 72 FEEs, up ramp parameters (7 of them) for GV function
  //  up ramp:   [0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fMergeCalTree->SetBranchAddress("epup",&fMergeCal.epup);  // for each of 72 FEEs, error on up ramp parameters (7 of them) for GV function
  fMergeCalTree->SetBranchAddress("pdn",&fMergeCal.pdn);  // for each of 72 FEEs, down ramp parameters (7 of them) for GV function
  //  down ramp:   -[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fMergeCalTree->SetBranchAddress("epdn",&fMergeCal.epdn);  // for each of 72 FEEs, error on down ramp parameters (7 of them) for GV function
  fMergeCalTree->SetBranchAddress("mippk",&fMergeCal.mippk);  // MIP peak, average most probable Landau fit value, used for determining the charge Z from the pulse integral between samples 11 and 25 (inclusive)
//  fMergeCalTree->SetBranchAddress("interTW",&fMergeCal.interTW);  // intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
//  fMergeCalTree->SetBranchAddress("slopeTW",&fMergeCal.slopeTW);  // slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
//  fMergeCalTree->SetBranchAddress("einterTW",&fMergeCal.einterTW);  // error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
//  fMergeCalTree->SetBranchAddress("eslopeTW",&fMergeCal.eslopeTW);  // error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)	
//  fMergeCalTree->SetBranchAddress("padnum",&fMergeCal.padnum);  // paddle number, HELIX style (0-7 bot, 8-15 top, 16 bore)
//  fMergeCalTree->SetBranchAddress("topx",&fMergeCal.topx);  // x values where we have thickness measurements for top paddles
//  fMergeCalTree->SetBranchAddress("topy",&fMergeCal.topy);  // y values where we have thickness measurements for top paddles
//  fMergeCalTree->SetBranchAddress("topt",&fMergeCal.topt);  // thickness (mm) of top paddle at topx,topy
//  fMergeCalTree->SetBranchAddress("botx",&fMergeCal.botx);  // x values where we have thickness measurements for bottom paddles
//  fMergeCalTree->SetBranchAddress("boty",&fMergeCal.boty);  // y values where we have thickness measurements for bottom paddles
//  fMergeCalTree->SetBranchAddress("bott",&fMergeCal.bott);  // thickness (mm) of bottom paddle at topx,topy

  // ToF calibration parameters that change with run number

  fMergeCalTree->SetBranchAddress("runnum",&fMergeCal.runnum);  // run number
  fMergeCalTree->SetBranchAddress("nument",&fMergeCal.nument);  // number of entries in that run
  fMergeCalTree->SetBranchAddress("runtime",&fMergeCal.runtime);  // time of the first event in the run, in days since 1st event in run 14753
  fMergeCalTree->SetBranchAddress("pedm",&fMergeCal.pedm);  // for each of 72 FEEs, mean value of the charge pedestal for this run
  fMergeCalTree->SetBranchAddress("rms",&fMergeCal.rms);  // for each of 72 FEEs, rms of the charge pedestal peak for this run
  fMergeCalTree->SetBranchAddress("mpv",&fMergeCal.mpv);  // for each of 72 FEEs, most probable charge value of the pedestal-subtracted signal for this run
  fMergeCalTree->SetBranchAddress("empv",&fMergeCal.empv);  // for each of 72 FEEs, uncertainty on the mpv value of the pedestal-subtracted charge signal for this run
  fMergeCalTree->SetBranchAddress("etoeshift",&fMergeCal.etoeshift);  // for each of 17 paddles, average time shift (in samples) between the two ends for this run (top: 0-7, bore: 8, bottom: 9-16)
  // for the following arrays of size 40: top are 0-15, bore are 16-23, bottom are 24-39;
  // top order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  // bore order is 00E, 01E, 12E, 23E, 00W, 01W, 12W, 23W   where 00E and 00W can be ignored 
  // bottom order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  fMergeCalTree->SetBranchAddress("peakuu",&fMergeCal.peakuu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+up
  fMergeCalTree->SetBranchAddress("widthuu",&fMergeCal.widthuu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+up
  fMergeCalTree->SetBranchAddress("peakud",&fMergeCal.peakud);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+dn
  fMergeCalTree->SetBranchAddress("widthud",&fMergeCal.widthud);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+dn
  fMergeCalTree->SetBranchAddress("peakdu",&fMergeCal.peakdu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+up
  fMergeCalTree->SetBranchAddress("widthdu",&fMergeCal.widthdu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+up
  fMergeCalTree->SetBranchAddress("peakdd",&fMergeCal.peakdd);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+dn
  fMergeCalTree->SetBranchAddress("widthdd",&fMergeCal.widthdd);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+dn
    
  //new branches from tofcal12:
  //Position from fChargeDivFrac by run:
  fMergeCalTree->SetBranchAddress("ZDm",  fMergeCal.ZDm );  //paddle slope, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  fMergeCalTree->SetBranchAddress("ZDb",  fMergeCal.ZDb );  //paddle intercept, from DCTY vs fChargeDivFrac: fChargeDivFracPos = m*ChargeDivFrac + b
  //Position from fTimeDiff by run:
  fMergeCalTree->SetBranchAddress("TDm",  fMergeCal.TDm );  //paddle slope, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b
  fMergeCalTree->SetBranchAddress("TDb",  fMergeCal.TDb );  //paddle intercept, from DCTY vs fTimeDiff: fTimeDiffPos = m*fTimeDiff + b

  //better scheme for top and bottom paddle thickness measurements. NB: EVENTUALLY MOVE THIS TO SINGLE READ IN SEPARATE OBJECT IN tofcal root FILE, NOT IN TREE!
  fMergeCalTree->SetBranchAddress("padx", fMergeCal.padx);  // x values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fMergeCalTree->SetBranchAddress("pady", fMergeCal.pady);  // y values where we have thickness measurements for top and bot paddles in HELIX numbering scheme
  fMergeCalTree->SetBranchAddress("padt", fMergeCal.padt);  // thickness (mm) of top paddle at padx,pady  
  //keep bore stuff:
  fMergeCalTree->SetBranchAddress("borex", fMergeCal.borex);  // x values where we have thickness measurements for bore paddle
  fMergeCalTree->SetBranchAddress("borey", fMergeCal.borey);  // y values where we have thickness measurements for bore paddle
  fMergeCalTree->SetBranchAddress("boret", fMergeCal.boret);  // thickness (mm) of bore paddle at topx,topy
  
  return kTRUE;

}

///////////////////////////////////////////////////////////////
// Open INPUT SCTOF parameter file (root TTree organized by run number)
bool loadSCTOFCalibsv11(string ipfName = kInputTofcalName) //only called once V11 and before!
{
  ////////  NB: put checks here to make sure file is not already open?  /////////  

  if (!std::filesystem::exists(ipfName.c_str())) {
    cout << "loadSCTOFCalibsv11 => file not found: " << ipfName << endl;
    exit(-1);
  }

  fInFile = new TFile(ipfName.c_str());
  if (!fInFile->IsOpen()) {
    cout << "loadSCTOFCalibs => could not open file: " << ipfName << endl;
    exit(-1);
  }

  cout << "loadSCTOFCalibsv11 : load from " << ipfName << endl;

  fCalTree = (TTree*)fInFile->Get("tofcal");
//  fCalTree = new TChain( "tofcal", "chain" );
//  fCalTree->Add( "tofcal6.root");  // new database with a different TW correction scheme, small bug fixed in ped determination 
//  TTree *fCalTree = new TTree( "tofcal", "tofcal6.root" );

  // global ToF parameters, same for any run number
  fCalTree->SetBranchAddress("pup",&fRtCal.pup);  // for each of 72 FEEs, up ramp parameters (7 of them) for GV function
  //  up ramp:   [0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fCalTree->SetBranchAddress("epup",&fRtCal.epup);  // for each of 72 FEEs, error on up ramp parameters (7 of them) for GV function
  fCalTree->SetBranchAddress("pdn",&fRtCal.pdn);  // for each of 72 FEEs, down ramp parameters (7 of them) for GV function
  //  down ramp:   -[0]+[0]*(1.-[1]*exp(-[2]*(x-[6]))+exp(-[3]*(x-[6]))*(([1]-1.)*cos([4]*(x-[6]))-[5]*sin([4]*(x-[6]))))
  fCalTree->SetBranchAddress("epdn",&fRtCal.epdn);  // for each of 72 FEEs, error on down ramp parameters (7 of them) for GV function
  fCalTree->SetBranchAddress("mippk",&fRtCal.mippk);  // MIP peak, average most probable Landau fit value, used for determining the charge Z from the pulse integral between samples 11 and 25 (inclusive)
  fCalTree->SetBranchAddress("interTW",&fRtCal.interTW);  // intercept for timewalk correction, for 64 pairs of paddles (set to 1. for trivial ones)
  fCalTree->SetBranchAddress("slopeTW",&fRtCal.slopeTW);  // slope for timewalk correction, for 64 pairs of paddles (set to 0. for trivial ones)	
  fCalTree->SetBranchAddress("einterTW",&fRtCal.einterTW);  // error on intercept for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)
  fCalTree->SetBranchAddress("eslopeTW",&fRtCal.eslopeTW);  // error on slope for timewalk correction, for 64 pairs of paddles (set to 0.0001 for trivial ones)	
  fCalTree->SetBranchAddress("padnum",&fRtCal.padnum);  // paddle number, HELIX style (0-7 bot, 8-15 top, 16 bore)
  fCalTree->SetBranchAddress("topx",&fRtCal.topx);  // x values where we have thickness measurements for top paddles
  fCalTree->SetBranchAddress("topy",&fRtCal.topy);  // y values where we have thickness measurements for top paddles
  fCalTree->SetBranchAddress("topt",&fRtCal.topt);  // thickness (mm) of top paddle at topx,topy
  fCalTree->SetBranchAddress("botx",&fRtCal.botx);  // x values where we have thickness measurements for bottom paddles
  fCalTree->SetBranchAddress("boty",&fRtCal.boty);  // y values where we have thickness measurements for bottom paddles
  fCalTree->SetBranchAddress("bott",&fRtCal.bott);  // thickness (mm) of bottom paddle at topx,topy
  fCalTree->SetBranchAddress("borex",&fRtCal.borex);  // x values where we have thickness measurements for bore paddle
  fCalTree->SetBranchAddress("borey",&fRtCal.borey);  // y values where we have thickness measurements for bore paddle
  fCalTree->SetBranchAddress("boret",&fRtCal.boret);  // thickness (mm) of bore paddle at topx,topy

  // ToF calibration parameters that change with run number

  fCalTree->SetBranchAddress("runnum",&fRtCal.runnum);  // run number
  fCalTree->SetBranchAddress("nument",&fRtCal.nument);  // number of entries in that run
  fCalTree->SetBranchAddress("runtime",&fRtCal.runtime);  // time of the first event in the run, in days since 1st event in run 14753
  fCalTree->SetBranchAddress("pedm",&fRtCal.pedm);  // for each of 72 FEEs, mean value of the charge pedestal for this run
  fCalTree->SetBranchAddress("rms",&fRtCal.rms);  // for each of 72 FEEs, rms of the charge pedestal peak for this run
  fCalTree->SetBranchAddress("mpv",&fRtCal.mpv);  // for each of 72 FEEs, most probable charge value of the pedestal-subtracted signal for this run
  fCalTree->SetBranchAddress("empv",&fRtCal.empv);  // for each of 72 FEEs, uncertainty on the mpv value of the pedestal-subtracted charge signal for this run
  fCalTree->SetBranchAddress("etoeshift",&fRtCal.etoeshift);  // for each of 17 paddles, average time shift (in samples) between the two ends for this run (top: 0-7, bore: 8, bottom: 9-16)
  // for the following arrays of size 40: top are 0-15, bore are 16-23, bottom are 24-39;
  // top order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  // bore order is 00E, 01E, 12E, 23E, 00W, 01W, 12W, 23W   where 00E and 00W can be ignored 
  // bottom order is 0E, 0W, 1E, 1W, 2E, 2W, 3E, 3W, 4E, 4W, 5E, 5W, 6E, 6W, 7E, 7W
  fCalTree->SetBranchAddress("peakuu",&fRtCal.peakuu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+up
  fCalTree->SetBranchAddress("widthuu",&fRtCal.widthuu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+up
  fCalTree->SetBranchAddress("peakud",&fRtCal.peakud);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are up+dn
  fCalTree->SetBranchAddress("widthud",&fRtCal.widthud);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are up+dn
  fCalTree->SetBranchAddress("peakdu",&fRtCal.peakdu);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+up
  fCalTree->SetBranchAddress("widthdu",&fRtCal.widthdu);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+up
  fCalTree->SetBranchAddress("peakdd",&fRtCal.peakdd);  // for each pair of side-by-side FEEs, shifted time difference in sample units when ramps are dn+dn
  fCalTree->SetBranchAddress("widthdd",&fRtCal.widthdd);  // for each pair of side-by-side FEEs, sigma of shifted time difference in sample units when ramps are dn+dn
    
  //make copy in case of mods:
  fCalTree->GetEntry( 0 ); // load up the calibration parameters structure for first run
  fCal = fRtCal;

  getSCTOFCalibsForRun( 14765 ); //set default values
  
  return kTRUE;

}

void convertTofcalv11tov12() //one time pass through all *tofcal* entries and writing out to a new tofcal tree
{
  // Load INPUT tofcal calibration data
  loadSCTOFCalibsv11();
//  loadSCTOFCalibsv12();
  
  //set up OUTPUT tofcal calibration  
  createSCTOFCalibs();
  
  int nentries = fCalTree->GetEntries();
  //loop thorugh tofcal entries
  for (int i=0; i<nentries; i++)
  {
    fCalTree->GetEntry(i);
    int runnum = fRtCal.runnum;
//    cout << "Converting entry, run: " << i << ",  " << runnum << endl; 

    //get current set of tofcal parameters:
    //unecessary because fCalTree-GetEntry() fills fCal.
    getSCTOFCalibsForRun(runnum, false);  //don't do adjustment of bad parameters

    //switch numbering scheme for this entry v11->v12:
    switchfCaltoHLXNumbering();

    //fill modified tofcal ttree 
    fOUTTOFcal->Fill();
  }
  
  //finalize output calib file and tree:
  fOUTTOFcal->Write();
  fOutFile->Close();  
}

void modifyTofcal() //one time pass through all *tofcal* entries and writing out to a new tofcal tree
{
  // Load INPUT tofcal calibration data
  loadSCTOFCalibsv12();
  
  //set up OUTPUT tofcal calibration  
  createSCTOFCalibs();
  
  int nentries = fCalTree->GetEntries();
  //loop thorugh tofcal entries
  for (int i=0; i<nentries; i++)
  {
    fCalTree->GetEntry(i);
    int runnum = fRtCal.runnum;
//    cout << "Converting entry, run: " << i << ",  " << runnum << endl; 

    //get current set of tofcal parameters:
    //unecessary because fCalTree-GetEntry() fills fCal.
//    getSCTOFCalibsForRun(runnum, false);  //don't do adjustment of bad parameters

    //add or modify branches:
    fillFitSlopeIntercepts(runnum);
    

    //fill modified tofcal ttree 
    fOUTTOFcal->Fill();
  }
  
  //finalize output calib file and tree:
  fOUTTOFcal->Write();
  fOutFile->Close();  
}

void updateSelectedRuns() //one time pass through all *tofcal* entries and writing out to a new tofcal tree
{
  //list of slected runs:
  vector<int> runList{14788, 14830, 14844, 14851, 14899, 15031, 15164, 15176, 15184};
  
  // Load INPUT tofcal calibration data
  loadSCTOFCalibsv12();
  loadSCTOFCalMergeFilev12();
  
  //set up OUTPUT tofcal calibration  
  createSCTOFCalibs();
  
  int nentries = fCalTree->GetEntries();
  //loop thorugh tofcal entries
  for (int i=0; i<nentries; i++)
  {
    fCalTree->GetEntry(i);
    int runnum = fRtCal.runnum;
//    cout << "Converting entry, run: " << i << ",  " << runnum << endl; 

    //get current set of tofcal parameters:
    getSCTOFCalibsForRun(runnum, false);  //don't do adjustment of bad parameters

    //fill particular branches if special run
    if (std::find(runList.begin(), runList.end(), runnum) != runList.end())
    {
      cout << "Found run " << runnum << endl;
      //scroll through merge tree to find right run
      for (int k=0; k<fMergeCalTree->GetEntries(); k++)
      {
        fMergeCalTree->GetEntry(k);
        if (fMergeCal.runnum == fCal.runnum) break;
      }
      cout << "Loaded merge event with runnum " << fMergeCal.runnum << endl;
      cout << "nument for fCal, fMergeCal: " << fCal.nument << "\t" << fMergeCal.nument << endl;
      cout << "runtime for fCal, fMergeCal: " << fCal.runtime << "\t" << fMergeCal.runtime << endl;
      fCal.runtime = fMergeCal.runtime;
      fCal.nument  = fMergeCal.nument;
      cout << "After: nument for fCal, fMergeCal: " << fCal.nument << "\t" << fMergeCal.nument << endl;
      cout << "After: runtime for fCal, fMergeCal: " << fCal.runtime << "\t" << fMergeCal.runtime << endl;
      for (int ifee=0; ifee<72; ifee++)
      {
        fCal.pedm[ifee]         =fMergeCal.pedm[ifee]      ;
        fCal.rms[ifee]          =fMergeCal.rms[ifee]       ;
        fCal.mpv[ifee]          =fMergeCal.mpv[ifee]       ;
        fCal.empv[ifee]         =fMergeCal.empv[ifee]      ;
        for (int ip=0; ip<7; ip++)
        {
          fCal.pup[ifee][ip]       =fMergeCal.pup[ifee][ip]    ;
          fCal.epup[ifee][ip]      =fMergeCal.epup[ifee][ip]   ;
          fCal.pdn[ifee][ip]       =fMergeCal.pdn[ifee][ip]    ;
          fCal.epdn[ifee][ip]      =fMergeCal.epdn[ifee][ip]   ;
        }        
      }
    }
    
    for (int HLXID=0; HLXID<=39; HLXID++)
    {
      fCal.peakuu[HLXID]       =fMergeCal.peakuu[HLXID]    ;
      fCal.widthuu[HLXID]      =fMergeCal.widthuu[HLXID]   ;
      fCal.peakud[HLXID]       =fMergeCal.peakud[HLXID]    ;
      fCal.widthud[HLXID]      =fMergeCal.widthud[HLXID]   ;
      fCal.peakdu[HLXID]       =fMergeCal.peakdu[HLXID]    ;
      fCal.widthdu[HLXID]      =fMergeCal.widthdu[HLXID]   ;
      fCal.peakdd[HLXID]       =fMergeCal.peakdd[HLXID]    ;
      fCal.widthdd[HLXID]      =fMergeCal.widthdd[HLXID]   ;
    }
    

    //add or modify branches:
//    fillFitSlopeIntercepts(runnum);  

    //fill modified tofcal ttree 
    fOUTTOFcal->Fill();
  }
  
  //finalize output calib file and tree:
  fOUTTOFcal->Write();
  fOutFile->Close();  
}

void MakeNewTofCalTree(string runFileListName) //one time pass through all *runList* entries and writing out to a new tofcal tree
{
  
  //this version makes new tofcal based on Bon and Boff combined list of runs
  //Need to properly fill other values, eg: nument, runtime, pedm
  
  vector<string> runList = getFileListVector(runFileListName);
  
  // Load INPUT tofcal calibration data
  loadSCTOFCalibsv12();
  
  //set up OUTPUT tofcal calibration  
  createSCTOFCalibs();

  //we want a tofcal entry for every file in the Bon and Boff combined run list
  for (int i=0; i<runList.size(); i++)
  {
//    cout << "Analyzing run " << std::stoi(runList[i]) << endl;
    //get current set of tofcal parameters:
    getSCTOFCalibsForRun(std::stoi(runList[i]), false);  //don't do adjustment of bad parameters
//    cout << "Before, after: mpv[5] = " << fCal.mpv[5];
    
    //add or modify branches:
    fillNewMPVValues(std::stoi(runList[i]));
//    cout << "\t" << fCal.mpv[5] << endl;

    //fill modified tofcal ttree 
    fOUTTOFcal->Fill();

  }

  //finalize output calib file and tree:
  fOUTTOFcal->Write();
  fOutFile->Close();  
}


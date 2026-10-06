//void AnalTOFSkel(TString infile) //infile is tof output created by skeletor
//infile: skeletor.15550.root

////////////////////////////////////////////////

//include PSU tof calib tools:
#include "PSUTofTools12.C"
/////  SET THIS IN PSUTofTools12.C!! ////
//Set up input tofcal ttree:
//string kInputTofcalName = "tofcal12.root";

#include "PSUThicknessTools.C"

void AnalTOFSkel(TString infile)
{
  /// Get the geometry singleton - it will autoload the files.
  auto theGeom = HLX_Geometry::getInstance();
  
  /////////////////////////////////////////////////
  //set up input root file for TTree:
  TFile* ipf = TFile::Open(infile.Data());
  if(!ipf || !ipf->IsOpen() || ipf->IsZombie()){
    cout << "Could not open file " << infile << endl;
    exit(0);
  }

  /// INPUT skeletor TOF info
//  vector<HLX_Stg2_TOFHit>* TOFHits = nullptr; //top, bore, bottom
  HLX_Stg2_TOFHit* topTOF = new HLX_Stg2_TOFHit();
  HLX_Stg2_TOFHit* boreTOF = new HLX_Stg2_TOFHit();
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

  TTree *toft = (TTree*)ipf->Get("toft");
  toft = (TTree*)ipf->Get("toft");

//  toft->SetBranchAddress("TOFHits", &TOFHits);
  toft->SetBranchAddress("topTOF", &topTOF);
  toft->SetBranchAddress("boreTOF", &boreTOF);
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
  
  //////////////////////////////////////////
  //book some histograms
  //These 3D histograms don't seem to work as TProfiles, which I prefer to fit over fitting the 2D histo.
//  TH3D* hZdct = new TH3D("hZdct","padID vs DCTY vs fChargeDivFrac;paddle;fChargeDivFrac;dctY (mm)",17,0,17, 50,-1,1, 200,-1000,1000);
//  TH3D* hTdct = new TH3D("hTdct","padID vs DCTY vs fTimeDiff;paddle;fTimeDiff;dctY (mm)",17,0,17, 50,-12,12, 200,-1000,1000);

  //pos from timing, chargeDivFrac, and difference
  TH2D** hChargeDivPos;
  hChargeDivPos = new TH2D*[17];
  TH2D** hTimeDiffPos;
  hTimeDiffPos      = new TH2D*[17];
  TH2D** hPosDiff;
  hPosDiff      = new TH2D*[17];
  TH2D** hXY;   // SC  check (X,Y) for each paddle
  hXY      = new TH2D*[17];
  TH2D** hXYshift;   // SC  check (X,Y) for each paddle; (X,Y) shifted to a paddle position (0-20 cm, 0-180 cm)
  hXYshift      = new TH2D*[17];
  TH2D** hXYChi;   // SC  check (X,Y) for each paddle, after a DCTChi2 cut
  hXYChi      = new TH2D*[17];
  TH2D** hXYHitC;   // SC  check (X,Y) for each paddle, after a DCT hitCount cut
  hXYHitC      = new TH2D*[17];
  TH1D** hthck;   // SC  thckness at (X,Y) for each paddle
  hthck      = new TH1D*[17];

  for (int i=0;i<17;i++)
  {
    TString hname   = TString::Format("%s%d","hChargeDivPos",i);
    TString htitle  = TString::Format("%s%d%s","fChargeDivPos vs DCTY, paddle",i,";DCTY (mm);Pos (mm)");
    if (i<16) hChargeDivPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-1000,1000, 200,-1000,1000);
    else      hChargeDivPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-100,100,  200,-400,400);
    hname.ReplaceAll("hChargeDivPos","hTimeDiffPos");
    htitle.ReplaceAll("ChargeDivPos","TimeDiffPos");
    if (i<16) hTimeDiffPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-1000,1000, 200,-1000,1000);
    else      hTimeDiffPos[i] = new TH2D(hname.Data(),htitle.Data(),200,-400,400,   200,-400,400);
    TString hname2  = TString::Format("%s%d","hPosDiff",i);
    TString htitle2 = TString::Format("%s%d%s","Position diff (fTimeDiffPos-fChargeDivPos) vs DCTY, paddle",i,";DCTY (mm);Pos diff (mm)");
    if (i<16) hPosDiff[i]     = new TH2D(hname2.Data(),htitle2.Data(),200,-500,500, 200,-1000,1000);
    else      hPosDiff[i]     = new TH2D(hname2.Data(),htitle2.Data(),200,-400,400, 200,-400,400);
    // SC create histograms of extrapolated (X,Y) value for each paddle, see if we see the scintillator shapes
    TString hname3  = TString::Format("%s%d","hXY",i);
    TString htitle3 = TString::Format("%s%d%s","DCTY vs DCTX, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXY[i]     = new TH2D(hname3.Data(),htitle3.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXY[i]     = new TH2D(hname3.Data(),htitle3.Data(),200,-400,400, 200,-400,400);
    TString hname6  = TString::Format("%s%d","hXYshift",i);
    TString htitle6 = TString::Format("%s%d%s","shifted DCTY vs DCTX, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXYshift[i]     = new TH2D(hname6.Data(),htitle6.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXYshift[i]     = new TH2D(hname6.Data(),htitle6.Data(),200,-400,400, 200,-400,400);
    TString hname4  = TString::Format("%s%d","hXYChi",i);
    TString htitle4 = TString::Format("%s%d%s","DCTY vs DCTX with DCTChi2 cut, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXYChi[i]     = new TH2D(hname4.Data(),htitle4.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXYChi[i]     = new TH2D(hname4.Data(),htitle4.Data(),200,-400,400, 200,-400,400);
    TString hname5  = TString::Format("%s%d","hXYHitC",i);
    TString htitle5 = TString::Format("%s%d%s","DCTY vs DCTX with DCT hitCount cut, paddle",i,";DCTX (mm);DCTY (mm)");
    if (i<16) hXYHitC[i]     = new TH2D(hname5.Data(),htitle5.Data(),200,-1000,1000, 200,-1000,1000);
    else      hXYHitC[i]     = new TH2D(hname5.Data(),htitle5.Data(),200,-400,400, 200,-400,400);
    TString hname7  = TString::Format("%s%d","hthck",i);
    TString htitle7 = TString::Format("%s%d","Thickness (mm) for paddle",i);
    hthck[i]     = new TH1D(hname7.Data(),htitle7.Data(),200,8.,14.);
  }

  //time walk explorations
  TH2D** hTW;
  hTW      = new TH2D*[72];  

  for (int i=0;i<72;i++)
  {
    TString hname   = TString::Format("%s%d","hTW",i);
    TString htitle  = TString::Format("%s%d%s","fTimeDiffPos - DCTy vs FEE signalSumRaw, fee ",i,";signalSumRaw;fTimeDiffPos - DCTY (mm)");
    hTW[i] = new TH2D(hname.Data(),htitle.Data(),200,0,1000, 200,-400,400);
  }
  
  //fChargeDivPos vs fTimeDiffPos timing differences
  TH2D** hZdT;
  hZdT = new TH2D*[17];
  TProfile** pZdT;
  pZdT = new TProfile*[17];
  
  for (int i=0;i<17;i++)
  {
    TString hname   = TString::Format("%s%d","hZdT",i);
    TString htitle  = TString::Format("%s%d%s","fChargeDivPos vs fTimeDiffPos, paddle",i,";fTimeDiffPos (mm);fChargeDivPos (mm)");
    if (i<16) hZdT[i] = new TH2D(hname.Data(),htitle.Data(),200,-1000,1000, 200,-1000,1000);
    else      hZdT[i] = new TH2D(hname.Data(),htitle.Data(),200,-400,400,   200,-1000,1000);
    hname.ReplaceAll("hZdT","pZdT");
    if (i<16) pZdT[i] = new TProfile(hname.Data(),htitle.Data(),200,-1000,1000, -1000,1000);
    else      pZdT[i] = new TProfile(hname.Data(),htitle.Data(),200,-400,400,   -1000,1000);
  }

  //position differences: dcty vs fChargeDivFrac
  TH2D** hZdct;
  hZdct = new TH2D*[17];
  TProfile** pZdct;
  pZdct = new TProfile*[17];
  
  for (int i=0;i<17;i++)
  {
    TString hname   = TString::Format("%s%d","hZdct",i);
    TString htitle  = TString::Format("%s%d%s","DCTY vs fChargeDivFrac, paddle",i,";fChargeDivFrac;dctY (mm)");
    if (i<16) hZdct[i] = new TH2D(hname.Data(),htitle.Data(),50,-1,1,    200,-1000,1000);
    else      hZdct[i] = new TH2D(hname.Data(),htitle.Data(),50,-0.4,0.4,200,-1000,1000);
    hname.ReplaceAll("hZdct","pZdct");
    if (i<16) pZdct[i] = new TProfile(hname.Data(),htitle.Data(),50,-1,1,    -1000,1000);
    else      pZdct[i] = new TProfile(hname.Data(),htitle.Data(),50,-0.4,0.4,-1000,1000); //bore paddle
  }

  //position differences: fTimeDiff vs dcty
  TH2D** hTdct;
  hTdct = new TH2D*[17];
  TProfile** pTdct;
  pTdct = new TProfile*[17];
  
  for (int i=0;i<17;i++)
  {
    TString hname   = TString::Format("%s%d","hTdct",i);
    TString htitle  = TString::Format("%s%d%s","DCTY vs fTimeDiff, paddle",i,";fTimeDiff;dctY (mm)");
    if (i<16) hTdct[i] = new TH2D(hname.Data(),htitle.Data(),60,-12,12,200,-1000,1000);
    else      hTdct[i] = new TH2D(hname.Data(),htitle.Data(),60,-4,4,  200,-1000,1000);
    hname.ReplaceAll("hTdct","pTdct");
    if (i<16) pTdct[i] = new TProfile(hname.Data(),htitle.Data(),60,-12,12,-1000,1000);
    else      pTdct[i] = new TProfile(hname.Data(),htitle.Data(),60,-4,4,  -1000,1000);
  }

  
  //fee to fee end timing differences
  TH1D** hEnddT;
  hEnddT = new TH1D*[17*2];
  
  for (int i=0;i<17;i++)
  {
    for (int j=0;j<2;j++)
    {
      TString hname = TString::Format("%s%d%s%d","hEnddT",i,"_",j);
      TString htitle = TString::Format("%s%d%s%d%s","Side by side time diff, paddle",i,", end ",j,";dT (ns)");
//      cout << "i*2+j = " << i*2+j << "  " << hname << "   " << htitle << endl;
      hEnddT[i*2+j] = new TH1D(hname.Data(),htitle.Data(),1000,-20,20);
    }
  }
  
  auto hTOFdT   = new TH1D("hTOFdT" , "dTOF;ns" , 1000, -20, 20);
  
  //////////////////////////////////
  /////  INPUT EVENT LOOP /////
  //loop over events
  auto nevts = toft->GetEntries();
  cout << "Number of events: " << nevts << endl;
  
  for (int ik=0; ik< nevts; ik++)
  {
    toft->GetEntry(ik);
    if (ik%10000 == 0) cout << "Event " << ik << endl;
    
    //get current set of tofcal parameters:
    getSCTOFCalibsForRun(runnum);
    
//    cout << ik << "\t" << tvars.ttofX << "\t" << topdircos[0] << "\t" << richZEst << endl;
//    cout << "\t" << topTOF->fTime << "\t" << boreTOF->getSignal() << endl;
    
    HLX_Stg2_TOFPaddleEndHit tE = topTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit tW = topTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest);
    HLX_Stg2_TOFPaddleEndHit bE = botTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit bW = botTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest);
    HLX_Stg2_TOFPaddleEndHit boreE = boreTOF->getPaddleEndHit(HLX_Geometry_TOF::kEast);
    HLX_Stg2_TOFPaddleEndHit boreW = boreTOF->getPaddleEndHit(HLX_Geometry_TOF::kWest);
    
    //watch this: an array of paddle ends!
    HLX_Stg2_TOFPaddleEndHit endHit[6] = {tE, tW, bE, bW, boreE, boreW};
    //cool, huh?
    
    //get fees on paddle ends
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
    
    ///TOF Beta and Charge
    double timeDiff = botTOF->fTime - topTOF->fTime;
    hTOFdT->Fill(timeDiff);

    //Collect results:

    ///ORDERING: TOP, BORE, BOTTOM

    //HLX numbering:
    int tpaddleID=topTOF->fTag.getPaddleID();
    int BpaddleID=boreTOF->fTag.getPaddleID();
    int bpaddleID=botTOF->fTag.getPaddleID();
    
    //Watch this trick: set up arrays to avoid lots of repeated code.
    int padID[3] = {tpaddleID, BpaddleID, bpaddleID};    
    HLX_Stg2_TOFHit* tofHits[3] = {topTOF, boreTOF, botTOF};  //not copies, so use pointer ->
    double tofX[3] = {tvars.ttofX,tvars.boreX,tvars.btofX};
    double tofY[3] = {tvars.ttofY,tvars.boreY,tvars.btofY};
    double padPLC[3] = {tvars.topPLC,tvars.borePLC,tvars.botPLC};

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
    
    ///following distilled down from SC code to correct signal for paddle thicknesses
    //paddle locations:
    double padx[3] = {
      theGeom->fTOF.fTopLocX + theGeom->fTOF.fPaddleMetrology[tpaddleID]["xloc"],
      theGeom->fTOF.fBorePaddleLocX + theGeom->fTOF.fPaddleMetrology[BpaddleID]["xloc"],
      theGeom->fTOF.fBottomLocX + theGeom->fTOF.fPaddleMetrology[bpaddleID]["xloc"]
    };
    double pady[3] = {
      theGeom->fTOF.fTopLocY + theGeom->fTOF.fPaddleMetrology[tpaddleID]["yloc"],
      theGeom->fTOF.fBorePaddleLocY + theGeom->fTOF.fPaddleMetrology[BpaddleID]["yloc"],
      theGeom->fTOF.fBottomLocY + theGeom->fTOF.fPaddleMetrology[bpaddleID]["yloc"]
    };

    //a lamda to check to see if hit is inside paddle:
    auto isInPaddle = [](double x, double xc, double tol) -> bool { 
        return (x <= xc + tol) && (x >= xc - tol); 
    };

    // check if track hit each paddle; set thickness
    double thck[3] = {10,10,10}; //mm, default thickness
    double tolx[3] = {100., 303., 100.}; //mm
    double toly[3] = {800., 303., 800.}; //mm
    
    for (int i=0; i<3; i++)
    {
      if (isInPaddle(tofX[i], padx[i], tolx[i]) && isInPaddle(tofY[i], pady[i], toly[i]))
      {
        hXYshift[padID[i]]->Fill(tofX[i]-padx[i],tofY[i]-pady[i]); 
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
    }
    
    // SC  Average scintillator thickness in mm, which was applied in pulse signal size normalization but which will have to be uncorrected back to raw
    // to apply the correct thickness correction; order here is HELIX ordering (bottom 0-7, top 8-15, bore 16)
    Double_t avgthick[17] = {
      9.79,9.64,10.19,9.56,10.09,9.72,9.84,9.87,  // bottom
      9.56,9.67,9.61,9.88,9.84,9.56,9.97,9.62,  // top
      9.97  // bore
    };

    
    
    //signals:
    double topSig = topTOF->getSignal();
    double botSig = botTOF->getSignal();
    double boreSig = boreTOF->getSignal();

    double ttofZ  = pow(topSig  * tvars.topPLC, 1.0/1.7);
    double boreZ  = pow(boreSig * tvars.borePLC, 1.0/1.7);
    double btofZ  = pow(botSig  * tvars.botPLC, 1.0/1.7);
    double allZ   = (ttofZ+boreZ+btofZ)/3.;
    double tbZ    = (ttofZ+btofZ)/2.;

    double ttofZSC  = topSig  * tvars.topPLC;
    double boreZSC  = boreSig * tvars.borePLC;
    double btofZSC  = botSig  * tvars.botPLC;
    double allZSC   = sqrt((ttofZSC+boreZSC+btofZSC)/3.);
    double allZSCcorr = (allZSC - 0.2499) / 0.8385;
    double tbZSC    = sqrt((ttofZSC+btofZSC)/2.);
    
    
    //fill all paddle-based hists and profiles
    for (int j=0; j<3; j++)
    {
      hChargeDivPos [padID[j]]->Fill(tofY[j]-padCenter[j],        tofHits[j]->fChargeDivPos);
      hTimeDiffPos  [padID[j]]->Fill(tofY[j]-padCenter[j],        tofHits[j]->fTimeDiffPos);
      hPosDiff      [padID[j]]->Fill(tofY[j]-padCenter[j],        tofHits[j]->fTimeDiffPos - tofHits[j]->fChargeDivPos);
      hZdT          [padID[j]]->Fill(tofHits[j]->fTimeDiffPos,    tofHits[j]->fChargeDivPos);
      pZdT          [padID[j]]->Fill(tofHits[j]->fTimeDiffPos,    tofHits[j]->fChargeDivPos);
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
      }
      
    }

  } //end of event loop
  
  cout << "End of event loop" << endl;
    
  ///////////////////////////////////////////
  //output root file for histograms

 //do any fits of profile hists here before writing hists, or in post-anal routine, with screen display.

  TString outfile = infile;
  int insertLoc = infile.Length()-4;
  outfile.Insert(insertLoc,"AnalTOFSkel.");
  TFile* opf = TFile::Open(outfile.Data(), "RECREATE");
  if(!opf || !opf->IsOpen() || opf->IsZombie()){
    cout << "Could not open file " << outfile << endl;
    exit(0);
  }
  
  for (int i=0; i<17; i++)
  {
    hChargeDivPos[i]->Write();
    hTimeDiffPos[i] ->Write();
    hPosDiff[i]     ->Write();    
    hZdT[i]         ->Write();
    pZdT[i]         ->Write();
    hZdct[i]        ->Write();
    pZdct[i]        ->Write();
    hTdct[i]        ->Write();
    pTdct[i]        ->Write();
    //SC:
    hXYshift[i]     ->Write();
    hXY[i]          ->Write();
    hXYChi[i]       ->Write();
    hXYHitC[i]      ->Write();
    hthck[i]        ->Write();
  }

  for (int i=0;i<17;i++)
  {
    for (int j=0;j<2;j++)
    {
      hEnddT[i*2+j]->Write();
    }
  }

  for (int i=0;i<72;i++)
  {
    hTW[i]->Write();
  }
  
  hTOFdT->Write();
 
  opf->Close();


}

// do stuff run by run instead of all runs in list at once with Chain
// 

//include tools to handle reading and writing tofcal
#include "PSUTofTools12.C"

//include PSU thickness tools for paddle thicknesses
#include "PSUThicknessTools.C"

void GetFeeSigHists() //one time pass through all entries and writing out to a new tofcal tree
{
  string datadir = "/data/skeletor/";
  string fstem = "skeletor.";
  string fend = ".root";


  /// Get the geometry singleton - it will autoload the files.
  auto theGeom = HLX_Geometry::getInstance();
  
  ////////////////////////////////////////////
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
  /////////////////////////////////////////////////

  
  // Load INPUT tofcal calibration data
  loadSCTOFCalibsv12();
  
  //set up OUTPUT tofcal calibration  
  //  createSCTOFCalibs();
  
  int nentries = fTOFcal->GetEntries();
  //loop thorugh tofcal entries
  for (int i=0; i<nentries; i++)
  {
    fTOFcal->GetEntry(i);
    int runnum = fRtCal.runnum;
    cout << "Analyzing entry, run: " << i << ",  " << runnum << endl; 
  
  ////////////////////////////  
    //fix missing run in tofcal:
//      runnum = 15630;
//      if (i>=1) continue; 
  //////////////////////////////
    
    //create filename
    string fname = datadir+fstem+std::to_string(runnum)+fend;
    cout << "Opening file " << fname << endl;
    
    //set up input root file for TTree:
    TFile* ipf = TFile::Open(fname.c_str());
    if(!ipf || !ipf->IsOpen() || ipf->IsZombie()){
      cout << "Could not open file " << fname << endl;
      continue; //skip this run, since it is not part of good runs list
    }

    TTree *toft = (TTree*)ipf->Get("toft");
    toft = (TTree*)ipf->Get("toft");

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


    ///////////////////////////////////////////
    //output root file for FEE signal histograms
    string outfile = "./FEESigs/FEESigs."+std::to_string(runnum)+fend;
    TFile *opf = TFile::Open(outfile.c_str(), "RECREATE");

    if(!opf || !opf->IsOpen() || opf->IsZombie())
    {
      cout << "Could not open file " << outfile << endl;
      exit(0);
    }
    
    //create FEE peak histos for this run
    TH1D** hFeeSig;
    hFeeSig      = new TH1D*[72];  
    TH1D** hFeeSigRaw;
    hFeeSigRaw      = new TH1D*[72];  

    for (int i=0;i<72;i++)
    {
      TString hname   = TString::Format("%s%d","hFeeSig",i);
      TString htitle  = TString::Format("%s%d%s","Corrected FEE SignalSum, fee ",i,";Corrected SignalSum");
      hFeeSig[i] = new TH1D(hname.Data(),htitle.Data(),200,0,10);
      hname   = TString::Format("%s%d","hFeeSigRaw",i);
      htitle  = TString::Format("%s%d%s","Raw FEE SignalSum, fee ",i,";SignalSumRaw");
      hFeeSigRaw[i] = new TH1D(hname.Data(),htitle.Data(),500,0,1000);
    }

    ///////////////////////////////////////
    ///  Event loop
    ///////////////////////////////////////
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
//      cout << "Got SCTOFCalibsForRun..." << endl;
    
      //Watch this trick: set up arrays to avoid lots of repeated code.
      HLX_Stg2_TOFHit* tofHits[3] = {topTOF, borTOF, botTOF};  //not copies, so use pointer ->
    
      /////////////////
      //Fixes for fuck ups:
      // v12.1.NT (new tracking): screwed up SignalSum. NB: THis takes care of normalization to mpv[8] but not thickness correction
      // Loop over FEEs in each paddle and fix:
      //HLX: 0-7 lower, 8-15 upper, 16 is bore
      // Ensure getFEEHit returns a reference (&), not a copy or a const
/*    
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
            double chanGain = fCal.mpv[SCfeeID]/fCal.mpv[8];  //normalize to FEE8 (paddle 0, FEE0)
            feeHit.fSignalSum = feeSignal/chanGain/fCal.mippk;
            //          feeHit.fSignalSum /= fCal.mippk;  //scale by global mippk to Z=1
            //          cout << "FIX: " << feeSignal << "\t" << chanGain << "\t" << fCal.mippk << "\t" << feeHit.fSignalSum << endl;
          }
        }
      }
      */
      ///////////////////////

//      cout << "Fixed screwups..." << endl;
      
      //gather information:
      
      int toppadID = topTOF->fTag.getPaddleID();
      int borpadID = borTOF->fTag.getPaddleID();
      int botpadID = botTOF->fTag.getPaddleID();
    
      int padID[3] = {toppadID, borpadID, botpadID};    

      double tofX[3]    = {tvars.ttofX, tvars.boreX,  tvars.btofX};
      double tofY[3]    = {tvars.ttofY, tvars.boreY,  tvars.btofY};
      double padPLC[3]  = {tvars.topPLC,tvars.borePLC,tvars.botPLC};
      
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

      //a lamda to check to see if hit is inside paddle:
      auto isInPaddle = [](double x, double xc, double tol) -> bool { 
        return (x <= xc + tol) && (x >= xc - tol); 
      };
      
      double tcorr[3] = {1., 1., 1.};   // thickness correction factor to adjust the Z value as per impact location thickness
      // check if track hit each paddle; set thickness
      double thck[3] = {10,10,10}; //mm, default thickness
      double tolx[3] = {100., 303., 100.}; //mm
      double toly[3] = {800., 303., 800.}; //mm
      

      //fill SC hit location histos and get paddle thickness at hit location
      for (int i=0; i<3; i++)
      {
        if (isInPaddle(tofX[i], padx[i], tolx[i]) && isInPaddle(tofY[i], pady[i], toly[i]))
        {
          double *padt2D = nullptr;
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
        }
        tcorr[i] = 10./thck[i];   // this tcorr factor will *multiply* the signal size (thickness in cm)
        //        tcorr[i] = avgthick[padID[i]] / thck[i];   // this tcorr factor will *multiply* the signal size
      }
      
//      cout << "Gathered info..." << endl;    
      
      //modify FEE signals with corrections for PLC and thickness and fill hists:
      // borrow code from fuck-up fix
      
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
            if ( ends[k]->getFEEHit(jk).fStatusFlag != HLX_Stg2_TOFFEEHit::kGood ) continue;   //Only fill FEEs flagged as good.
            int feeID = ends[k]->getFEEHit(jk).fTag.getFEEID();
            auto& feeHit = ends[k]->fFEEHits[jk];
            auto [laneIdx, chanIdx] = theGeom->fTOF.fGeoMap[{paddleID, endID, feeID}];
            int SCfeeID = laneIdx*8 + chanIdx; //HLX labels. tofcal12
            double feeSignalRaw = feeHit.fSignalSumRaw;
            //fuck-up fix:
            double chanGain = fCal.mpv[SCfeeID]/fCal.mpv[8];  //normalize to FEE8 (paddle 0, FEE0)
            feeHit.fSignalSum = feeSignalRaw/chanGain/fCal.mippk;
            //back to normal:
            feeHit.fSignalSum *= tcorr[i]/padPLC[i];
            double feeSignal = feeHit.fSignalSum;
            hFeeSig[SCfeeID]->Fill(sqrt(feeSignal));
            feeSignalRaw *= tcorr[i]/padPLC[i];
            hFeeSigRaw[SCfeeID]->Fill(feeSignalRaw);
          }
        }
      }
      
      //fill modified tofcal ttree 
      //    fOUTTOFcal->Fill();
    
    } // event loop
  
    opf->Write(); //write and close FEE hists file
    opf->Close();
    
    ipf->Close(); //close input skeletor file

    //finalize output calib file and tree:
    //  fOUTTOFcal->Write();
    //fOutFile->Close();  
  } //Tofcal tree loop
  
}


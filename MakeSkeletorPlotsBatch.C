// SC  The DrawDownArrowMm thing is from ChatGPT
// Draw a downward arrow at xData. Coordinates are data coordinates, but
// physical dimensions are specified in mm. canvasWmm/Hmm must match the
// intended physical canvas size (especially for printed/PDF output).

void DrawDownArrowMm(double xData,
                     double xAxisY,
                     double canvasWmm = 200.0,
                     double canvasHmm = 200.0,
                     double arrowTopMm = 7.0*10,
                     double arrowTipMm = 2.0*10,
                     double headWidthMm = 2.0*10,
                     double headLengthMm = 1.5*10,
                     Color_t color = kRed)
{
  if (!gPad) return;
  gPad->Update();                         // establish axis/pad transforms

  const double pxPerMmX = gPad->GetCanvas()->GetWw() / canvasWmm;
  const double pxPerMmY = gPad->GetCanvas()->GetWh() / canvasHmm;

  const int xpx = gPad->XtoPixel(xData);
  const int yAxisPx = gPad->YtoPixel(xAxisY);
  const int topPx = yAxisPx - TMath::Nint(arrowTopMm * pxPerMmY);
  const int tipPx = yAxisPx - TMath::Nint(arrowTipMm * pxPerMmY);
  const int basePx = tipPx - TMath::Nint(headLengthMm * pxPerMmY);
  const int halfWidthPx = TMath::Nint(0.5 * headWidthMm * pxPerMmX);

  const double yTop = gPad->PixeltoY(topPx);
  const double yTip = gPad->PixeltoY(tipPx);
  const double yBase = gPad->PixeltoY(basePx);
  const double xLeft = gPad->PixeltoX(xpx - halfWidthPx);
  const double xRight = gPad->PixeltoX(xpx + halfWidthPx);

  auto *shaft = new TLine(xData, yTop, xData, yBase);
  shaft->SetLineColor(color);
  shaft->SetLineWidth(2);
  shaft->Draw();

  double x[4] = {xLeft, xData, xRight, xLeft};
  double y[4] = {yBase, yTip, yBase, yBase};
  auto *head = new TPolyLine(4, x, y, "f");
  head->SetFillColor(color);
  head->SetLineColor(color);
  head->Draw("f");

  gPad->Modified();
  gPad->Update();
}

void MakeSkeletorPlotsBatch(TString fname, TString runNum="", TString TOFtype="", bool writeFits = true)
{
  
  //Examples:
  //fname = "/data/skeletor/skeletor.14802.AnalTOFSkel.root"
  //runNum="14802"
  //TOFtype="GV" or "KS" or blank
  
  
  cout << "Opening file name = " << fname << endl;
  cout << "Run number: " << runNum << endl;
  cout << "TOFtype: " << TOFtype << endl;
  if (TOFtype != "" ) TOFtype += ".";
  
  TFile* f1 = TFile::Open(fname.Data());
  if (!f1 || f1->IsZombie())
  {
    cout << "ERROR ERROR ERROR ERROR" << endl;
    cout << "Unable to open file " << fname << endl;
    exit(0);
  }

//  TString baseName="MakeSkeletorPlotsBatch.";
  TString baseName=fname;
  baseName.ReplaceAll("skeletor.","MakeSkeletorPlotsBatch.");
  baseName.ReplaceAll("AnalTOFSkel.","");
  baseName.ReplaceAll(".root",".pdf");
//  TString fout = baseName+TOFtype+runNum+".pdf";
  TString fout = baseName;
  cout << "Writing " << fout << endl;

  TCanvas** cv;
  cv = new TCanvas*[15];
  int nx = 4, ny = 4;
  int cknt=0;
  int hknt=1;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy();
  gStyle->SetOptFit(1111);
  
  TString cname;

  TH1D** hEnddT;
  hEnddT = new TH1D*[16*2];
  
  for (int i=0;i<16;i++) //bore = 16
  {
    for (int j=0;j<2;j++)
    {
      TString hname = TString::Format("%s%d%s%d","hEnddT",i,"_",j);
      hEnddT[i*2+j] = (TH1D*)f1->Get(hname.Data());
//      cout << "Defining hist " << i << "  " << j << "  " << i*2+j << endl;
    }
  }

  TString toutname = fout;
  toutname.ReplaceAll(".pdf",".EnddT.txt");
  
  ofstream tout;
  
  if (writeFits) tout.open(toutname.Data(),std::ios::out);
  
  for (int i=0;i<16;i++)
  {
    for (int j=0;j<2;j++)
    {
      if( (i*2+j)%(nx*ny) == 0) 
      {
        TString cname = TString::Format("%s%d","c",cknt);
        cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
        cv[cknt]->Divide(nx,ny);
        cout << endl << "NEW CANVAS " << cknt << endl << endl;
        cknt++;
        hknt=1;
      }      
      cv[cknt-1]->cd(hknt);
      //first pass to find centroid:
      hEnddT[i*2+j]->Fit("gaus","","",-8.0,8.0);
      double mean = hEnddT[i*2+j]->GetFunction("gaus")->GetParameter("Mean");
      //second pass:
      hEnddT[i*2+j]->Fit("gaus","","",mean-0.625,mean+0.625);
      double sigma = hEnddT[i*2+j]->GetFunction("gaus")->GetParameter("Sigma");
      mean = hEnddT[i*2+j]->GetFunction("gaus")->GetParameter("Mean");      
      hEnddT[i*2+j]->Draw();
      if (writeFits) tout << i << "\t" << j << "\t" << mean << "\t" << sigma << endl; 
//      cout << endl << "DREW HIST " << hknt << "  " << hEnddT[i*2+j]->GetName() << endl << endl;
      hknt++;
    }
  }
  
  if (writeFits) tout.close();

  //now draw other plots:
  /*
  // time of flight
  gStyle->SetOptLogy(0);
  TH1D* hTOFdT    = (TH1D*)f1->Get("hTOFdT");
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  hTOFdT->Draw();
  cknt++;
  
  // bot vs top charge
  gStyle->SetOptLogy(0);
  TH2D* htvsbTOF    = (TH2D*)f1->Get("htvsbTOF");
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  htvsbTOF->Draw("colz");
  cknt++;

  //TOF hit frequency
  gStyle->SetOptLogy(1);
  TH1D* hTOFHitMap    = (TH1D*)f1->Get("hTOFHitMap");
  TH1D* hSCTOFHitMap  = (TH1D*)f1->Get("hSCTOFHitMap");

  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  cv[cknt]->Divide(1,2);
  cv[cknt]->cd(1);
  hTOFHitMap->Draw();
  cv[cknt]->cd(2);
  hSCTOFHitMap->Draw();
  cknt++;
  
  //signal (charge)
  TH1D* hTTOF   = (TH1D*)f1->Get("hTTOF");
  TH1D* hbTOF   = (TH1D*)f1->Get("hbTOF");
  TH1D* hBTOF   = (TH1D*)f1->Get("hBTOF");
  TH1D* hallTOF = (TH1D*)f1->Get("hallTOF");
  TH1D* htbTOF  = (TH1D*)f1->Get("htbTOF");
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(1);
  cv[cknt]->Divide(3,2);
  cv[cknt]->cd(1);
  hTTOF->Draw();
  cv[cknt]->cd(2);
  hbTOF->Draw();
  cv[cknt]->cd(3);
  hBTOF->Draw();
  cv[cknt]->cd(4);
  hallTOF->Draw();
  cv[cknt]->cd(5);
  htbTOF->Draw();

  cknt++;
*/  

  //pos from timing, chargeDivFrac, and difference
  TH2D** hChargeDivPos;
  hChargeDivPos = new TH2D*[17];

  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries and title only
  cv[cknt]->Divide(nx,ny); //4*4 = 16
  
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hChargeDivPos",i);
    hChargeDivPos[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hChargeDivPos[i]->GetName() << endl << endl;
    hChargeDivPos[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hChargeDivPos[15]->Draw();
    }
  }
  
  cknt++;

  ///////////
  TH2D** hTimeDiffPos;
  hTimeDiffPos      = new TH2D*[17];

  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries and title only
  cv[cknt]->Divide(nx,ny); //4*4 = 16
  
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hTimeDiffPos",i);
    hTimeDiffPos[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hTimeDiffPos[i]->GetName() << endl << endl;
    hTimeDiffPos[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hTimeDiffPos[15]->Draw();
    }
  }
  
  cknt++;
  
  ////////
  TH2D** hPosDiff;
  hPosDiff      = new TH2D*[17];
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries and title only
  cv[cknt]->Divide(nx,ny); //4*4 = 16
  
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hPosDiff",i);
    hPosDiff[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hPosDiff[i]->GetName() << endl << endl;
    hPosDiff[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hPosDiff[15]->Draw();
    }
  }
  
  cknt++;
  
  
  //bore paddle plots
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(2,2);
  
  cv[cknt]->cd(1);
  hChargeDivPos[16]->Draw();

  cv[cknt]->cd(2);
  hTimeDiffPos[16]->Draw();

  cv[cknt]->cd(3);
  hPosDiff[16]->Draw();
  
  cknt++;
  
  
  
  //paddle end to end dT profiles
  TH2D** hZdT;
  hZdT = new TH2D*[17];
  TProfile** pZdT;
  pZdT = new TProfile*[17];
//  TString poutname = baseName+TOFtype+runNum+".PosCompare.txt";  
  TString poutname = baseName;
  poutname.ReplaceAll(".pdf",".PosCompare.txt");  
  ofstream pout;
  if (writeFits) pout.open(poutname.Data(),std::ios::out);
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries and title only
  cv[cknt]->Divide(nx,ny); //4*4 = 16
    
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hZdT",i);
    hZdT[i] = (TH2D*)f1->Get(hname.Data());
    hname.ReplaceAll("hZdT","pZdT");
    pZdT[i] = (TProfile*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
    //fits of profile hists:
    pZdT[i]->Fit("pol1","Q","",-400.,400.);
    TF1* pfit = pZdT[i]->GetFunction("pol1");
    double b = (pfit)? pfit->GetParameter("p0"): -10;
    double m = (pfit)? pfit->GetParameter("p1"): -10;
    if ( isnan(b) || isnan(m) ) b = m = -10;
//    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
    if (writeFits) pout << i << "\t" << b << "\t" << m << endl;
//    cout << endl << "DRAWING " << i << "  " << hZdT[i]->GetName() << endl << endl;
    hZdT[i]->Draw();
    pZdT[i]->Draw("same");
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hZdT[15]->Draw();
      pZdT[15]->Draw("same");
    }    
  }
  
  if (writeFits) pout.close();
  
  cknt++;
  
  
  //DCTY: ChargeDivFrac vs dct Y position
  TH2D** hZdct;
  hZdct = new TH2D*[17];
  TProfile** pZdct;
  pZdct = new TProfile*[17];

//  poutname = baseName+TOFtype+runNum+".ChargeDivFrac.txt";
  poutname = baseName;
  poutname.ReplaceAll(".pdf",".ChargeDivFrac.txt");  
//  ofstream pout;
  if (writeFits) pout.open(poutname.Data(),std::ios::out);
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hZdct",i);
    hZdct[i] = (TH2D*)f1->Get(hname.Data());
    hname.ReplaceAll("hZdct","pZdct");
    pZdct[i] = (TProfile*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
    //fits of profile hists:
    if (i<16) pZdct[i]->Fit("pol1","Q","",-0.4,0.4);
    else      pZdct[i]->Fit("pol1","Q","",-0.1,0.1);
    TF1* pfit = pZdct[i]->GetFunction("pol1");
    double b = (pfit)? pfit->GetParameter("p0"): -10;
    double m = (pfit)? pfit->GetParameter("p1"): -10;
    if ( isnan(b) || isnan(m) ) b = m = -10;
//    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
    if (writeFits) pout << i << "\t" << b << "\t" << m << endl;
//    cout << endl << "DRAWING " << i << "  " << hZdct[i]->GetName() << endl << endl;
    hZdct[i]->Draw();
    pZdct[i]->Draw("same");
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hZdct[15]->Draw();
      pZdct[15]->Draw("same");
    }    
  }

  if (writeFits) pout.close();

  cknt++;

  //DCTY: TimeDiff vs dct Y position
  TH2D** hTdct;
  hTdct = new TH2D*[17];
  TProfile** pTdct;
  pTdct = new TProfile*[17];

//  poutname = baseName+TOFtype+runNum+".TimeDiff.txt";
  poutname = baseName;
  poutname.ReplaceAll(".pdf",".TimeDiff.txt");  
  
  if (writeFits) pout.open(poutname.Data(),std::ios::out);
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hTdct",i);
    hTdct[i] = (TH2D*)f1->Get(hname.Data());
    hname.ReplaceAll("hTdct","pTdct");
    pTdct[i] = (TProfile*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
    //fits of profile hists:
    if (i<16) pTdct[i]->Fit("pol1","Q","",-5.,5.); //regular paddles
    else      pTdct[i]->Fit("pol1","Q","",-2.,2.); //bore paddle
    TF1* pfit = pTdct[i]->GetFunction("pol1");
    double b = (pfit)? pfit->GetParameter("p0"): -10;
    double m = (pfit)? pfit->GetParameter("p1"): -10;
    if ( isnan(b) || isnan(m) ) b = m = -10;
//    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
    if (writeFits) pout << i << "\t" << b << "\t" << m << endl;
//    cout << endl << "DRAWING " << i << "  " << hTdct[i]->GetName() << endl << endl;
    hTdct[i]->Draw();
    pTdct[i]->Draw("same");
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hTdct[15]->Draw();
      pTdct[15]->Draw("same");
    }    
  }

  if (writeFits) pout.close();

  cknt++;
  /*
  
  /// These aren't currently in AnalTOFSkel but can be added if needed ////
  
  //log10(Se/Sw) vs fTimeDiff to get veff
  TH2D** hveff;
  hveff = new TH2D*[17];
  TProfile** pveff;
  pveff = new TProfile*[17];

  poutname = baseName+TOFtype+runNum+".veff.txt";  
  if (writeFits) pout.open(poutname.Data(),std::ios::out);
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hveff",i);
    hveff[i] = (TH2D*)f1->Get(hname.Data());
    hname.ReplaceAll("hveff","pveff");
    pveff[i] = (TProfile*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
    //fits of profile hists:
//    pveff[i]->Fit("pol1","Q","",-450.,450.);
    if (i<16) pveff[i]->Fit("pol1","","",-5.,5.);
    else      pveff[i]->Fit("pol1","","",-2.,2.);
    TF1* pfit = pveff[i]->GetFunction("pol1");
    double b = (pfit)? pfit->GetParameter("p0"): -10;
    double m = (pfit)? pfit->GetParameter("p1"): -10;
    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
    cout << endl << "DRAWING " << i << "  " << hveff[i]->GetName() << endl << endl;
    hveff[i]->Draw();
    pveff[i]->Draw("same");
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hveff[15]->Draw();
      pveff[15]->Draw("same");
    }    
  }

  pout.close();

  cknt++;
  */
  
  //SC plots:
  //XY
  TH2D** hXY;   // SC  check (X,Y) for each paddle
  hXY      = new TH2D*[17];
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  //get all hists but only draw paddles, not bore paddle
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hXY",i);
    hXY[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hXY[i]->GetName() << endl << endl;
    hXY[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hXY[15]->Draw();
    }    
  }

  cknt++;
  


  //XYshift
  TH2D** hXYshift;   // SC  check (X,Y) for each paddle; (X,Y) shifted to a paddle position (0-20 cm, 0-180 cm)
  hXYshift      = new TH2D*[17];
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  //get all hists but only draw paddles, not bore paddle
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hXYshift",i);
    hXYshift[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hXYshift[i]->GetName() << endl << endl;
    hXYshift[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hXYshift[15]->Draw();
    }    
  }

  cknt++;
  
  ////////
  TH2D** hXYChi;   // SC  check (X,Y) for each paddle, after a DCTChi2 cut
  hXYChi      = new TH2D*[17];
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  //get all hists but only draw paddles, not bore paddle
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hXYChi",i);
    hXYChi[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hXYChi[i]->GetName() << endl << endl;
    hXYChi[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hXYChi[15]->Draw();
    }    
  }

  cknt++;
  
  ////////////
  TH2D** hXYHitC;   // SC  check (X,Y) for each paddle, after a DCT hitCount cut
  hXYHitC      = new TH2D*[17];
  
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(nx,ny); //4*4 = 16    

  //get all hists but only draw paddles, not bore paddle
  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hXYHitC",i);
    hXYHitC[i] = (TH2D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hXYHitC[i]->GetName() << endl << endl;
    hXYHitC[i]->Draw();
    if (i==16)  //work around root not obeying "0" goption in Fit()
    {
      hXYHitC[15]->Draw();
    }    
  }

  cknt++;
  
  ////////////
  // Thickness distribution for each paddle
  // SC  Average scintillator thickness in mm
  // order here is HELIX ordering (bottom 0-7, top 8-15, bore 16)
  Double_t avgthick[17] = {
    9.79,9.64,10.19,9.56,10.09,9.72,9.84,9.87,  // bottom
    9.56,9.67,9.61,9.88,9.84,9.56,9.97,9.62,  // top
    9.97  // bore
  };
  
  TH1D** hthck;
  hthck      = new TH1D*[17];

  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries and title only
  cv[cknt]->Divide(nx,ny); //4*4 = 16
  
  for (int i=0;i<16;i++)
  {
    TString hname = TString::Format("%s%d","hthck",i);
    hthck[i] = (TH1D*)f1->Get(hname.Data());
    cv[cknt]->cd(i+1);
//    cout << endl << "DRAWING " << i << "  " << hthck[i]->GetName() << endl << endl;
    hthck[i]->Draw();
    DrawDownArrowMm(avgthick[i], 0.0, 200, 200);
  }
  
  cknt++;

  //////////////////  
  //Now draw bore paddle plots: first canvas
  cname = TString::Format("%s%d","c",cknt);
  cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only
  cv[cknt]->Divide(3,3);
  
  int pad = 1;
  
  cv[cknt]->cd(pad++);
  hZdct[16]->Draw();
  pZdct[16]->Draw("same");

  cv[cknt]->cd(pad++);
  hTdct[16]->Draw();
  pTdct[16]->Draw("same");

//  cv[cknt]->cd(pad++);
//  hveff[16]->Draw();
//  pveff[16]->Draw("same");
  
  cv[cknt]->cd(pad++);
  hZdT[16]->Draw();
  pZdT[16]->Draw("same");
  
  cv[cknt]->cd(pad++);
  hXY[16]->Draw();
  
  cv[cknt]->cd(pad++);
  hXYshift[16]->Draw();
  
  cv[cknt]->cd(pad++);
  hXYChi[16]->Draw();
  
  cv[cknt]->cd(pad++);
  hXYHitC[16]->Draw();
  
  cv[cknt]->cd(pad++);
  TString hname = TString::Format("%s%d","hthck",16);
  hthck[16] = (TH1D*)f1->Get(hname.Data());
  hthck[16]->Draw();
  DrawDownArrowMm(avgthick[16], 0.0, 200, 200);

  cknt++;

  //////////////
  // Write all canvases to pdf

  for (int i=0; i<cknt; i++)
  {
    if (i<cknt-1) 
      cv[i]->Print(fout+"(");
    else
      cv[i]->Print(fout+")");
      
  }
  
}

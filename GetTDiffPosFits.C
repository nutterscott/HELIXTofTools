//Open an AnalTOFSkel file before loading this
//get fits to fTImeDiff vs DCTY plots and draw

void GetTDiffPosFits(bool writeFits = false)
{
  int nx = 4, ny = 4;
  int cknt = -1;
  TCanvas** cv;
  cv = new TCanvas*[10];

  
  //DCTY: TimeDiff vs dct Y position
  TH2D** hTdct;
  hTdct = new TH2D*[17];
  TProfile** pTdct;
  pTdct = new TProfile*[17];

  TString poutname = gFile->GetName();
  poutname.ReplaceAll(".root",".TimeDiff.txt");  
  
  ofstream pout;
  if (writeFits) pout.open(poutname.Data(),std::ios::out);
  

  for (int i=0;i<17;i++)
  {
    if (i%(nx*ny) == 0)
    {
      cknt++;
      TString cname = TString::Format("%s%d","c",cknt);
      cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
      gStyle->SetOptLogy(0);
      gStyle->SetOptStat(11); //# of entries only
      if ( i/(nx*ny) < 1) cv[cknt]->Divide(nx,ny); //4*4 = 16 
    }
    TString hname = TString::Format("%s%d","hTdct",i);
    hTdct[i] = (TH2D*)gFile->Get(hname.Data());
    hname.ReplaceAll("hTdct","pTdct");
    pTdct[i] = (TProfile*)gFile->Get(hname.Data());
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
//    if (i==16)  //work around root not obeying "0" goption in Fit()
//    {
//      hTdct[15]->Draw();
//      pTdct[15]->Draw("same");
//    }    
  }

  if (writeFits) pout.close();

}

void DrawTDm()
{
  
  TTree* tofcal = (TTree*)gFile->Get("tofcal");
  
  int colors[17] = {kOrange,kBlue-6,kGray,kRed,kOrange+6, kOrange-6,  //0-5
                    kCyan-6, kCyan+2, kCyan, kBlue-10, kMagenta,      //6-10
                    kYellow,kSpring,kBlue,kGreen,kBlack, kAzure};     //11-16

  auto cn = new TCanvas("cn","cn",1);
  TH2D* tdm = new TH2D("tdm","Veff vs flight time",100,0,6.5,100,60,120);
  tdm->SetStats(0);
  tdm->Draw();
  tofcal->SetMarkerStyle(20);  // 20 = solid circle (essential so markers are visible)
  tofcal->SetMarkerSize(0.5);  // Adjust marker size  
  for (int i=0; i<17; i++)
  {
    TString dstring = TString::Format("%s%d%s","TDm[",i,"]:runtime");
    tofcal->SetMarkerColorAlpha(colors[i],0.5);
    tofcal->Draw(dstring.Data(),"","same");
  } 
  
  auto cn1 = new TCanvas("cn1","cn1",1);
  cn1->Divide(4,4);
  
  TH1D** htdm;
  htdm = new TH1D*[17];
  double mean[17];
  for (int i=0; i<16; i++)
  {
    cn1->cd(i+1);
    TString dstring = TString::Format("%s%d%s%d%s","TDm[",i,"]>>htdm",i,"(100,60,100)");
    tofcal->Draw(dstring.Data());
    TString hname = TString::Format("%s%d","htdm",i);
    htdm[i] = (TH1D*)gDirectory->Get(hname.Data());
    htdm[i]->Draw();
    mean[i] = htdm[i]->GetMean();
    cout << i << "\t" << mean[i] << endl;
  }

  auto cn2 = new TCanvas("cn2","cn2",1);
  tofcal->Draw("TDm[16]>>htdm16(100,80,120)");
  htdm[16] = (TH1D*)gDirectory->Get("htdm16");
  htdm[16]->Draw();
  mean[16] = htdm[16]->GetMean();
  cout << "16" << "\t" << mean[16] << endl;
   
}

void DrawDCTX()
{
  TH2D** hXY;
  hXY = new TH2D*[17];
  TH1D** hXYx;
  hXYx = new TH1D*[17];
  
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only

  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hXY",i);
    hXY[i] = (TH2D*)gFile->Get(hname.Data());
    hXY[i]->ProjectionX();
    hXYx[i] = (TH1D*)gDirectory->Get((hname+"_px").Data());
    hXYx[i]->SetLineColor(i%8 + 1);
  }
  
  auto cxybot = new TCanvas("cxybot","cxybot",1);

  TH1D* hXYxbot = (TH1D*)hXYx[0]->Clone("hXYxbot");
  hXYxbot->Reset();
  
  //BOT
  for (int i=8;i<16;i++)
  {
    TString opt = "";
    if (i>8) opt = "same";
    hXYx[i]->Draw(opt);
    hXYxbot->Add(hXYx[i]);
  }
  hXYxbot->SetLineColor(9);
  hXYxbot->Draw("same");

  auto cxytop = new TCanvas("cxytop","cxytop",1);

  TH1D* hXYxtop = (TH1D*)hXYx[8]->Clone("hXYxtop");
  hXYxtop->Reset();

  //TOP
  for (int i=0;i<8;i++)
  {
    TString opt = "";
    if (i>0) opt = "same";
    hXYx[i]->Draw(opt);
    hXYxtop->Add(hXYx[i]);
  }
  hXYxtop->SetLineColor(9);
  hXYxtop->Draw("same");
    
}

void DrawDCTXSlice(int lobin = 90, int hibin=110)
{
  TH2D** hXY;
  hXY = new TH2D*[17];
  TH1D** hXYx;
  hXYx = new TH1D*[17];
  
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(11); //# of entries only

  for (int i=0;i<17;i++)
  {
    TString hname = TString::Format("%s%d","hXY",i);
    hXY[i] = (TH2D*)gFile->Get(hname.Data());
    hXY[i]->ProjectionX("_px",lobin,hibin);
    hXYx[i] = (TH1D*)gDirectory->Get((hname+"_px").Data());
    hXYx[i]->SetLineColor(i%8 + 1);
  }
  
  auto cxybot = new TCanvas("cxybot","cxybot",1);

  TH1D* hXYxbot = (TH1D*)hXYx[0]->Clone("hXYxbot");
  hXYxbot->Reset();
  
  //BOT
  for (int i=8;i<16;i++)
  {
    TString opt = "";
    if (i>8) opt = "same";
    hXYx[i]->Draw(opt);
    hXYxbot->Add(hXYx[i]);
  }
  hXYxbot->SetLineColor(9);
  hXYxbot->Draw("same");

  auto cxytop = new TCanvas("cxytop","cxytop",1);

  TH1D* hXYxtop = (TH1D*)hXYx[8]->Clone("hXYxtop");
  hXYxtop->Reset();

  //TOP
  for (int i=0;i<8;i++)
  {
    TString opt = "";
    if (i>0) opt = "same";
    hXYx[i]->Draw(opt);
    hXYxtop->Add(hXYx[i]);
  }
  hXYxtop->SetLineColor(9);
  hXYxtop->Draw("same");
    
}
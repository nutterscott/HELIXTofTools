//analyze all paddles:
//REQUIRED: open ANalTOF output to get hSignalG histograms from

// To use:
// Single paddle: 
// TGraphErrors* SliceAndDice(TH2 *h2, int method = 2, int bincomb=1, double xMin=-800, double xMax=800)

// All paddles: 
// void DoSignalAnal(TString histoName = "hSignal", int method = 2, int bincomb = 2)

// Methods: 
//  0: use bin with max contents
//  1: gaus fit to peak
//  2: landau fit to peak

// bincomb: how many bins in x to combine. Try 2 - 4. Statistics matter!


// int extractPadID(string text) //get paddle number from TString hname
// void DrawSlices(TH1D **hsli, int nsli)  //draw slices
// void DrawEverything(TH1D **hsli, int nsli, TGraphErrors* g1, TH2 *h2) //draw everything:
// TGraphErrors* SliceAndDice(TH2 *h2, int method = 1, double xMin=-800, double xMax=800, int bincomb=1) //Make slices
// void MakeProjectionXPlots(TH2D** hSignal)
// void DoSignalAnal(int method = 2, int bincomb = 2)
// void DrawEndPlots()
// void DrawPadSigPlots()


//GLOBAL hists:
//  for hSignal anal:
TGraphErrors** grSig;
grSig = new TGraphErrors*[17];
// for heSig (ends) anal:
TGraphErrors** geSig;
geSig = new TGraphErrors*[17*2];


int colors[17] = {kOrange,kBlue-6,kGray,kRed,kOrange+6, kOrange-6,  //0-5
                  kCyan-6, kCyan+2, kCyan, kBlue-10, kMagenta,      //6-10
                  kYellow,kSpring,kBlue,kGreen,kBlack, kAzure};     //11-16

//get paddle number from TString hname
std::vector<int> extractPadID(string text) {

    // Match one or more consecutive digits
    std::regex pattern(R"(\d+)"); 
    std::smatch matches;
    
    std::vector<int> numbers;
    
    // Search the string continuously for matches
    std::string::const_iterator search_start(text.cbegin());
    while (std::regex_search(search_start, text.cend(), matches, pattern)) {
        // Convert the string match to an integer and store it
        numbers.push_back(std::stoi(matches[0]));
        
        // Move the search window past the current match
        search_start = matches.suffix().first;
    }

    // Output results
    for (int num : numbers) {
        std::cout << "Paddle ID found: " << num << endl;
    }

    return numbers; //list of numbers found
}

//draw slices:
void DrawSlices(TH1D **hsli, int nsli)
{
  //check whether this is hSignal or heSig
  double xmin = -1.0, xmax = 2.5; //heSig  
  TString hname = hsli[0]->GetName();
  if (hname.Contains("hSignal"))  {    xmin = 0; xmax = 4;      }
  
  int nx = 4, ny = 4;
  const int ncan = static_cast<int>(std::ceil(double(nsli)/nx/ny));
  TCanvas** csli;
  csli = new TCanvas*[ncan];
  int ican = -1;
  int ipad = 1;
  for (int i=0; i<nsli; i++)
  {
    if (i%(nx*ny) == 0)
    {
      ican++;
      cout << "Opening TCanvas " << ican << endl;
      csli[ican] = new TCanvas(Form("%s%d","csli",ican),"Slices",1);
      csli[ican]->Divide(nx,ny);
      ipad = 1;
    }
    csli[ican]->cd(ipad);
    hsli[i]->GetXaxis()->SetRangeUser(xmin, xmax);
    cout << "Drawing pad " << ipad << endl;
    hsli[i]->Draw();
    ++ipad;
  }
}

//draw everything:
void DrawEverything(TH1D **hsli, int nsli, TGraphErrors* g1, TH2 *h2)
{
  gStyle->SetOptStat(11); //# of entries and title only
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas* ce1 = new TCanvas("ce1","Graph1",1);
//  gPad->SetTickx(1);
//  gPad->SetTicky(1);
  g1->Draw("alp");
  TCanvas* ce2 = new TCanvas("ce2","Graph2",1);
  TH1D* h2x = h2->ProjectionX();
  h2x->Draw();
  DrawSlices(hsli, nsli);
}

//Make slices:
TGraphErrors* SliceAndDice(TH2 *h2, int method = 2, int bincomb=1, double xMin=-800, double xMax=800)
{
  //make slices along x from 2D plot
  // in a reduced number of x bins by the factor bincomb, eg bincomb = 1, 2, 3, 5
  // h2 = TH2D 
  // Returns graph of peak vs x

// Methods: 
//  0: use bin with max contents
//  1: gaus fit to peak
//  2: landau fit to peak


  
  int nh = h2->GetXaxis()->GetNbins();
  const int nbins = nh/bincomb;
  
  cout << "# Histogram x bins: " << nh << endl;
  cout << "# Analysis x bins: " << nbins << endl;
  
  TH1D **hsli;
  hsli = new TH1D*[nbins];
  double x[nbins], maxval[nbins], pk[nbins], pkerr[nbins];
  TString hname = h2->GetName();
  
  //check whether this is hSignal or heSig
  double xmin = -1.0, xmax = 2.5; //heSig  
  if (hname.Contains("hSignal"))  {    xmin = 0.1; xmax = 2.0;    }
  
  

  for (int i=0;i<nbins;i++) 
  {
    hsli[i] = h2->ProjectionY(Form("%s_%d",hname.Data(),i+1),i*bincomb+1,(i+1)*bincomb);
    x[i] = (h2->GetXaxis()->GetBinCenter(i*bincomb+1)+h2->GetXaxis()->GetBinCenter((i+1)*bincomb))/2.;
    hsli[i]->SetTitle(Form("%s%3.1f","X: ",x[i]));
    hsli[i]->SetName(Form("%s%s%d",hname.Data(),"_py",i));
  }
  
  //now get pk for each slice
  int knt = 0;
  TF1* fg = new TF1("fg","gaus",xmin,xmax);
  fg->SetParameters(10000.,1.0,0.1);        //initial guess
  TF1* flan = new TF1("flan","landau",xmin,xmax);
  flan->SetParameters(10000.,1.0,0.1);        //initial guess

  for (int i=0;i<nbins;i++) 
  { 
    double bcenter = (h2->GetXaxis()->GetBinCenter(i*bincomb+1)+h2->GetXaxis()->GetBinCenter((i+1)*bincomb))/2.;
    if ( bcenter /*x[i] above*/ < xMin ) continue;
    x[knt] = bcenter;
    hsli[i]->GetXaxis()->SetRangeUser(xmin,xmax);  //limit following search for maxval to this range
    maxval[knt] = hsli[i]->GetBinCenter(hsli[i]->GetMaximumBin());
    double ampl = hsli[i]->GetBinContent(hsli[i]->GetMaximumBin());
    fg->SetParameters(ampl, maxval[knt],0.1);
    flan->SetParameters(ampl, maxval[knt],0.1);
    if (method <= 1) hsli[i]->Fit(fg,"QO","",maxval[knt]-0.4,maxval[knt]+0.3);
    if (method == 2) hsli[i]->Fit(flan,"QO","",maxval[knt]-0.5,maxval[knt]+0.6);
    TF1* fit = hsli[i]->GetFunction("fg");
    if (method == 2) fit = hsli[i]->GetFunction("flan");
    pk[knt] = (fit)? fit->GetParameter(1) : 0;
    pkerr[knt] = (fit)? fit->GetParError(1) : 0;
    if ( x[knt] > xMax ) break;
    cout << x[knt] << "\t" << pk[knt] << endl;
    knt++;
  }
  
  cout << "# Slices: " << knt << endl;
  
  TGraphErrors *g1;
  if (method == 0) 
    g1= new TGraphErrors(knt, x, maxval);
  else if (method >= 1)
    g1 = new TGraphErrors(knt, x, pk,nullptr,pkerr);
  g1->SetMarkerStyle(kFullCircle);
  g1->SetMarkerColor(kGreen);
  std::vector<int> padID = extractPadID( string(hname.Data()) );
  TString gtitle;
  if (padID.size() == 1) //hSignal paddle hists
    gtitle = TString::Format("%s%d%s","Paddle ",padID[0],";DCTY (mm);Signal peak (arb)");
  else  //heSig end hists
    gtitle = TString::Format("%s%d%s%d%s","Paddle ",padID[0]," end ", padID[1], ";DCTY (mm);Signal peak (MIP)");
  g1->SetTitle(gtitle);
  g1->SetName(Form("%s%s",hname.Data(),"_gr"));
/*  g1->Fit("pol3");
  TF1 *f1 = g1->GetFunction("pol3");
  double *pars = f1->GetParameters();
  
  //write out parameters
  cout << endl; //get some space
  cout << charge << "\t";
  for (int i=0; i<4; i++) {cout << pars[i] << "\t";}
  cout << endl;
  cout << endl; //get some space
*/  


  DrawEverything(hsli, knt, g1, h2);
  

  return g1;
}

void MakeProjectionXPlots(TH2D** hSignal)
{
  TH1D** hpx;
  hpx = new TH1D*[17];

  TCanvas *cht, *chb, *chbr;
  
  cht = new TCanvas("cht","Top projX",1);
  for (int i=0; i<16; i++) //top and bot
  {
    if (i==8)
    {
      cht->BuildLegend();
      chb = new TCanvas("chb","Bottom projX",1);
    }
    //get projections
    hpx[i] = (TH1D*)hSignal[i]->ProjectionX();
    hpx[i]->SetLineColor(colors[i]);
    string opt = "same";
    if (i==0 || i==8) opt = "";
    hpx[i]->Draw(opt.c_str());
  }
  chb->BuildLegend();
  
  //bore
  chbr = new TCanvas("chbr","Bore projX",1);
  hpx[16] = (TH1D*)hSignal[16]->ProjectionX();
  hpx[16]->SetLineColor(colors[16]);
  string opt = "same";
  hpx[16]->Draw();

}

void DoSignalAnal(TString histoName = "hSignal", int method = 2, int bincomb = 2)
{
  //histoName options: AnalTOFSkelv4:
  // hSignal, hSignalA, hSignalG
  // hSignalHiR, hSignalAHiR, hSignalGHiR
  
//  int colors[17];
//  for (int i=0; i<17; i++) colors[i] = 1+2*i;

//  int colors[17] = {kOrange,kBlue-6,kGray,kRed,kOrange+6, kOrange-6,  //0-5
//                    kCyan-6, kCyan+2, kCyan, kBlue-10, kMagenta,      //6-10
//                    kYellow,kSpring,kBlue,kGreen,kBlack, kAzure};     //11-16
  
  TH2D** hSignalG;
  hSignalG = new TH2D*[17];

//  TGraphErrors** grSig;
//  grSig = new TGraphErrors*[17];
  
  //get histograms:
  for (int i=0; i<17; i++)
  {
    TString hname = TString::Format("%s%d",histoName.Data(),i);
//    TString hname = TString::Format("%s%d","hSignalGHiR",i);
//    TString hname = TString::Format("%s%d","hSignalG",i);
    hSignalG[i] = (TH2D*)_file0->Get(hname.Data());
    double xmin = -800, xmax = 800; //mm
    if (i==16) { xmin = -300; xmax = 300;} //bore paddle
    grSig[i]   = SliceAndDice(hSignalG[i],method,bincomb,xmin,xmax);
    string grname = string(hname.Data())+"_gr";
    grSig[i]->SetName(grname.c_str());
  }
  
  MakeProjectionXPlots(hSignalG);
  
  //draw them all and save pol2 fits
  ofstream txtout = ofstream("PadSigTGraphFits.txt",std::ios::out);

  TCanvas* cgrtop = new TCanvas("cgrtop","Top graphs",1);
  TCanvas* cgrbot;
//  TLegend* ltop   = new TLegend();
  //TOP + BOT:
  double ymin=0.8, ymax = 1.6;
  for (int i=0; i<16; i++)
  {
    if (i==8)
    {
      cgrtop->BuildLegend();
      cgrbot = new TCanvas("cgrbot","Bottom graphs",1);
    }
    grSig[i]->SetLineColor(colors[i]);
    grSig[i]->SetMarkerColor(colors[i]);
    grSig[i]->Fit("pol2","QO","",-700.,700.);
    TF1* fit = grSig[i]->GetFunction("pol2");
    fit->SetLineColor(colors[i]);
    double* pars = fit->GetParameters();
    cout << "Graph for paddle " << i << " fit params: " << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
    txtout << i << "\t" << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
    string opt = "lpsame";
    if (i==0 || i==8) opt = "alp";
//    cgrtop->cd(i+1);
    if (i==0 || i==8) grSig[i]->GetYaxis()->SetRangeUser(ymin,ymax);
    grSig[i]->Draw(opt.c_str());
  }
  cgrbot->BuildLegend();

/*  
  //BOT:
  TCanvas* cgrbot = new TCanvas("cgrbot","Bottom graphs",1);
  for (int i=8; i<16; i++)
  {
    grSig[i]->SetLineColor(colors[i]);
    grSig[i]->SetMarkerColor(colors[i]);
    grSig[i]->Fit("pol2");
    TF1* fit = grSig[i]->GetFunction("pol2");
    fit->SetLineColor(colors[i]);
    double* pars = fit->GetParameters();
    cout << "Graph for paddle " << i << " fit params: " << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
    txtout << i << "\t" << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
    string opt = "lpsame";
    if (i==8) opt = "alp";
//    cgrbot->cd(i+1);
    if (i==8) grSig[i]->GetYaxis()->SetRangeUser(0.6,1.2);
    grSig[i]->Draw(opt.c_str());
  }
  cgrbot->BuildLegend();
  */
  TCanvas* cgrbore = new TCanvas("cgrbore","Bore graph",1);
  grSig[16]->SetLineColor(colors[16]);
  grSig[16]->SetMarkerColor(colors[16]);
  grSig[16]->Fit("pol2");
  TF1* fit = grSig[16]->GetFunction("pol2");
  fit->SetLineColor(colors[16]);
  double* pars = fit->GetParameters();
  cout << "Graph for paddle 16 fit params: " << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
  txtout << "16\t" << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
  grSig[16]->GetYaxis()->SetRangeUser(ymin,ymax);
  grSig[16]->Draw("alp");
  
  txtout.close();
  
  //the end  
}

void DrawEndPlots(TH2D** heSig)
{
  
  
  //display
  TCanvas* cend[12];
  int cknt = -1;
    
  gStyle->SetOptLogy(0);
  gStyle->SetOptStat(1); //title only
  gStyle->SetOptFit(11); 
  
  //heSig:
  for (int i=0;i<17*2;i++)
  {
    if (i%4 == 0)
    {
      ++cknt;
      TString cname = TString::Format("%s%d","c",cknt);
      cend[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
      cend[cknt]->Divide(2,2); //2*2 = 4    
    }
    cend[cknt]->cd(i%4+1);
    //fits of profile hists:
//  if (i<16*2) peSig[i]->Fit("pol1","Q","",-500.,500.); //regular paddles
//  else        peSig[i]->Fit("pol1","Q","",-250.,250.); //bore paddle
//  TF1* pfit = peSig[i]->GetFunction("pol1");
//  double b = (pfit)? pfit->GetParameter("p0"): -10;
//  double m = (pfit)? pfit->GetParameter("p1"): -10;
//  if ( isnan(b) || isnan(m) ) b = m = -10;
//    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
//    if (writeFits) pout << i << "\t" << b << "\t" << m << endl;
//    cout << endl << "DRAWING " << i << "  " << heSig[i]->GetName() << endl << endl;
//    cout << heSig[i]->GetName() << "\t" << b << "\t" << m << endl;
    heSig[i]->Draw();
//    peSig[i]->Draw("same");
  }
  
//  //heSigA:
//  TCanvas* cend2[12];
//  int cknt2=-1;
//
//  for (int i=0;i<17*2;i++)
//  {
//    if (i%4 == 0)
//    {
//      ++cknt2;
//      TString cname = TString::Format("%s%d","cc",cknt2);
//      cend2[cknt2] = new TCanvas(cname.Data(),cname.Data(),1);
//      cend2[cknt2]->Divide(2,2); //2*2 = 4    
//    }
//    cend2[cknt2]->cd(i%4+1);
//    //fits of profile hists:
//    if (i<16*2) peSigA[i]->Fit("pol1","Q","",-500.,500.); //regular paddles
//    else        peSigA[i]->Fit("pol1","Q","",-250.,250.); //bore paddle
//    TF1* pfit = peSigA[i]->GetFunction("pol1");
//    double b = (pfit)? pfit->GetParameter("p0"): -10;
//    double m = (pfit)? pfit->GetParameter("p1"): -10;
//    if ( isnan(b) || isnan(m) ) b = m = -10;
////    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
////    if (writeFits) pout << i << "\t" << b << "\t" << m << endl;
////    cout << endl << "DRAWING " << i << "  " << heSig[i]->GetName() << endl << endl;
//    cout << heSigA[i]->GetName() << "\t" << b << "\t" << m << endl;
//    heSigA[i]->Draw();
////    peSigA[i]->Draw("same");
//  }
//  
//  //heSigG:
//  TCanvas* cend3[12];
//  int cknt3=-1;
//
//  for (int i=0;i<17*2;i++)
//  {
//    if (i%4 == 0)
//    {
//      ++cknt3;
//      TString cname = TString::Format("%s%d","ccc",cknt3);
//      cend3[cknt3] = new TCanvas(cname.Data(),cname.Data(),1);
//      cend3[cknt3]->Divide(2,2); //2*2 = 4 
//    }
//    cend3[cknt3]->cd(i%4+1);
//    //fits of profile hists:
////    if (i<16*2) peSigG[i]->Fit("pol1","Q","",-500.,500.); //regular paddles
////    else        peSigG[i]->Fit("pol1","Q","",-250.,250.); //bore paddle
////    TF1* pfit = peSigG[i]->GetFunction("pol1");
////    double b = (pfit)? pfit->GetParameter("p0"): -10;
////    double m = (pfit)? pfit->GetParameter("p1"): -10;
////    if ( isnan(b) || isnan(m) ) b = m = -10;
////    if (writeFits) pout << i << "\t" << "{" << b << ", " << m << "}," << endl;
////    if (writeFits) pout << i << "\t" << b << "\t" << m << endl;
////    cout << endl << "DRAWING " << i << "  " << heSig[i]->GetName() << endl << endl;
//    cout << heSigG[i]->GetName() << "\t" << b << "\t" << m << endl;
//    heSigG[i]->Draw();
////    peSigG[i]->Draw("same");
//  }
//  
}

void DoEndSignalAnal(int method = 2, int bincomb = 2)
{
  
//  int colors[17];
//  for (int i=0; i<17; i++) colors[i] = 1+2*i;

//  int colors[17] = {kOrange,kBlue-6,kGray,kRed,kOrange+6, kOrange-6,  //0-5
//                    kCyan-6, kCyan+2, kCyan, kBlue-10, kMagenta,      //6-10
//                    kYellow,kSpring,kBlue,kGreen,kBlack, kAzure};     //11-16
  
  TH2D** heSig;             //ln(end signal) vs dcty, SPW getSignal()
  heSig      = new TH2D*[17*2];  
  TH2D** heSigA;             //ln(end signal) vs dcty, arithmetic mean
  heSigA      = new TH2D*[17*2];  
  TH2D** heSigG;             //ln(end signal) vs dcty, geometric mean
  heSigG      = new TH2D*[17*2];  
  
//  TProfile** peSig;             //ln(end signal) vs dcty, SPW getSignal()
//  peSig      = new TProfile*[17*2];  
//  TProfile** peSigA;             //ln(end signal) vs dcty, arithmetic mean
//  peSigA      = new TProfile*[17*2];  
//  TProfile** peSigG;             //ln(end signal) vs dcty, geometric mean
//  peSigG      = new TProfile*[17*2];  
  
  //get histograms:
  for (int i=0; i<17; i++) //paddles
  {
    for (int j=0; j<2; j++) //ends
    {
      TString hname   = TString::Format("%s%d%s%d","heSigHiR",i,"_",j);
//      TString hname   = TString::Format("%s%d%s%d","heSig",i,"_",j);
      heSig[i*2+j] = (TH2D*)gFile->Get(hname.Data());
      hname   = TString::Format("%s%d%s%d","heSigAHiR",i,"_",j);
//      hname   = TString::Format("%s%d%s%d","heSigA",i,"_",j);
      heSigA[i*2+j] = (TH2D*)gFile->Get(hname.Data());
      hname   = TString::Format("%s%d%s%d","heSigGHiR",i,"_",j);
//      hname   = TString::Format("%s%d%s%d","heSigG",i,"_",j);
      heSigG[i*2+j] = (TH2D*)gFile->Get(hname.Data());

      //create graphs:
      double xmin = -800, xmax = 800; //mm
      if (i==16) { //bore paddle
        xmin = -300; xmax = 300;
        geSig[i*2+j]   = SliceAndDice(heSigA[i*2+j],method,bincomb,xmin,xmax);
      } else {
        xmin = -800, xmax = 800; //mm
        geSig[i*2+j]   = SliceAndDice(heSigG[i*2+j],method,bincomb,xmin,xmax);
      }
      string grname = string(hname.Data())+"_gr";
      geSig[i*2+j]->SetName(grname.c_str());

//    hname   = TString::Format("%s%d%s%d","peSig",i,"_",j);
//    peSig[i*2+j] = (TProfile*)gFile->Get(hname.Data());
//    hname   = TString::Format("%s%d%s%d","peSigA",i,"_",j);
//    peSigA[i*2+j] = (TProfile*)gFile->Get(hname.Data());
//    hname   = TString::Format("%s%d%s%d","peSigG",i,"_",j);
//    peSigG[i*2+j] = (TProfile*)gFile->Get(hname.Data());

    }
  }

  // now GLOBAL:
//  TGraphErrors** geSig;
//  geSig = new TGraphErrors*[17*2];
  
  
  //draw them all and save pol2 fits
  ofstream txtout = ofstream("EndSigTGraphFits.txt",std::ios::out);

  TCanvas* cgetop = new TCanvas("cgetop","Top end graphs",1);
//  cgetop->Divide(1,2);
  TCanvas* cgebot;

  for (int i=0; i<16; i++)
  {
    for (int j=0; j<2; j++)
    {
      if (i==8 && j==0)
      {
        cgetop->BuildLegend();
        cgebot = new TCanvas("cgebot","Bottom end graphs",1);
//        cgebot->Divide(1,2);
      }
//      if (i<8 && j==0) cgetop->cd(1);
//      else if (i<8 && j==0) cgetop->cd(2);
//      else if (i>=8 && j==0) cgebot->cd(1);
//      else if (i>=8 && j==1) cgebot->cd(2);
//      cout << i*2+j << endl;
      geSig[i*2+j]->SetLineColor(colors[i]);
      geSig[i*2+j]->SetMarkerColor(colors[i]);
      geSig[i*2+j]->Fit("pol2","QO","",-700.,700.);
      TF1* fit = geSig[i*2+j]->GetFunction("pol2");
      fit->SetLineColor(colors[i]);
      double* pars = fit->GetParameters();
      cout << "Graph for paddle " << i << ", end " << j << " fit params: " << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
      txtout << i << "\t" << j << "\t" << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
      string opt = "lpsame";
      if ((i==0 || i==8) && j==0) opt = "alp";
  //    cgetop->cd(i+1);
      if ((i==0 || i==8) && j==0) geSig[i*2+j]->GetYaxis()->SetRangeUser(-1.,2.0);
      geSig[i*2+j]->Draw(opt.c_str());
    }
  }
  cgebot->BuildLegend();

  TCanvas* cgebore = new TCanvas("cgebore","Bore ends graph",1);
  for (int j=0; j<2; j++)
  {
    geSig[16*2+j]->SetLineColor(colors[16]);
    geSig[16*2+j]->SetMarkerColor(colors[16]);
    geSig[16*2+j]->Fit("pol2","QO","",-700.,700.);
    TF1* fit = geSig[16*2+j]->GetFunction("pol2");
    fit->SetLineColor(colors[16]);
    double* pars = fit->GetParameters();
    cout << "Graph for paddle 16 end " << j << " fit params: " << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
    txtout << "16\t" << j << "\t" << pars[0] << "\t" << pars[1] << "\t" << pars[2] << endl;
    geSig[16*2+j]->GetYaxis()->SetRangeUser(-1.,2.0);
    if (j==0) geSig[16*2+j]->Draw("alp");
    else geSig[16*2+j]->Draw("same");
  }
  
  txtout.close();
  
  //the end  
}



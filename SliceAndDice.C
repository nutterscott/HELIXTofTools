//analyze all paddles:
//REQUIRED: open ANalTOF output to get hSignal histograms from

//GLOBAL hists:
TGraphErrors** grSig;
grSig = new TGraphErrors*[17];

int colors[17] = {kOrange,kBlue-6,kGray,kRed,kOrange+6, kOrange-6,  //0-5
                  kCyan-6, kCyan+2, kCyan, kBlue-10, kMagenta,      //6-10
                  kYellow,kSpring,kBlue,kGreen,kBlack, kAzure};     //11-16

//get paddle number from TString hname
int extractPadID(string text) {

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

    return numbers[0]; //one number found
}

//draw slices:
void DrawSlices(TH1D **hsli, int nsli)
{
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
    hsli[i]->GetXaxis()->SetRangeUser(0., 4.);
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
TGraphErrors* SliceAndDice(TH2 *h2, int method = 1, double xMin=-800, double xMax=800, int bincomb=1)
{
  //make slices along x from 2D plot
  // in a reduced number of x bins by the factor bincomb, eg bincomb = 1, 2, 3, 5
  // h2 = TH2D 
  // Returns graph of peak vs x
  
  int nh = h2->GetXaxis()->GetNbins();
  const int nbins = nh/bincomb;
  
  cout << "# Histogram x bins: " << nh << endl;
  cout << "# Analysis x bins: " << nbins << endl;
  
  TH1D **hsli;
  hsli = new TH1D*[nbins];
  double x[nbins], maxval[nbins], pk[nbins], pkerr[nbins];
  TString hname = h2->GetName();

  for (int i=0;i<nbins;i++) 
  {
    hsli[i] = h2->ProjectionY(Form("%s_%d",hname.Data(),i+1),i*bincomb+1,(i+1)*bincomb);
    x[i] = (h2->GetXaxis()->GetBinCenter(i*bincomb+1)+h2->GetXaxis()->GetBinCenter((i+1)*bincomb))/2.;
    hsli[i]->SetTitle(Form("%s%3.1f","X: ",x[i]));
    hsli[i]->SetName(Form("%s%s%d",hname.Data(),"_py",i));
  }
  
  //now get pk for each slice
  int knt = 0;
  TF1* fg = new TF1("fg","gaus",0.1,1.5);
  fg->SetParameters(10000.,1.0,0.1);        //initial guess
  TF1* flan = new TF1("flan","landau",0.1,1.6);
  flan->SetParameters(10000.,1.0,0.1);        //initial guess

  for (int i=0;i<nbins;i++) 
  { 
    double bcenter = (h2->GetXaxis()->GetBinCenter(i*bincomb+1)+h2->GetXaxis()->GetBinCenter((i+1)*bincomb))/2.;
    if ( bcenter /*x[i] above*/ < xMin ) continue;
    x[knt] = bcenter;
    hsli[i]->GetXaxis()->SetRangeUser(0.2,1.6);  //limit following search for maxval to this range
    maxval[knt] = hsli[i]->GetBinCenter(hsli[i]->GetMaximumBin());
    double ampl = hsli[i]->GetBinContent(hsli[i]->GetMaximumBin());
    fg->SetParameters(ampl, maxval[knt],0.1);
    flan->SetParameters(ampl, maxval[knt],0.1);
    if (method <= 1) hsli[i]->Fit(fg,"QO","",maxval[knt]-0.3,maxval[knt]+0.2);
    if (method == 2) hsli[i]->Fit(flan,"QO","",maxval[knt]-0.3,maxval[knt]+0.3);
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
  int padID = extractPadID( string(hname.Data()) );
  TString gtitle = TString::Format("%s%d%s","Paddle ",padID,";DCTY (mm);Signal peak (arb)");
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

void DoSignalAnal(int method = 2, int bincomb = 2)
{
  
//  int colors[17];
//  for (int i=0; i<17; i++) colors[i] = 1+2*i;

//  int colors[17] = {kOrange,kBlue-6,kGray,kRed,kOrange+6, kOrange-6,  //0-5
//                    kCyan-6, kCyan+2, kCyan, kBlue-10, kMagenta,      //6-10
//                    kYellow,kSpring,kBlue,kGreen,kBlack, kAzure};     //11-16
  
  TH2D** hSignal;
  hSignal = new TH2D*[17];

//  TGraphErrors** grSig;
//  grSig = new TGraphErrors*[17];
  
  //get histograms:
  for (int i=0; i<17; i++)
  {
    TString hname = TString::Format("%s%d","hSignal",i);
    hSignal[i] = (TH2D*)_file0->Get(hname.Data());
    double xmin = -800, xmax = 800; //mm
    if (i==16) { xmin = -300; xmax = 300;} //bore paddle
    grSig[i]   = SliceAndDice(hSignal[i],method,xmin,xmax,bincomb);
    string grname = string(hname.Data())+"_gr";
    grSig[i]->SetName(grname.c_str());
  }
  
  MakeProjectionXPlots(hSignal);
  
  //draw them all and save pol2 fits
  ofstream txtout = ofstream("PadSigTGraphFits.txt",std::ios::out);

  TCanvas* cgrtop = new TCanvas("cgrtop","Top graphs",1);
  TCanvas* cgrbot;
//  TLegend* ltop   = new TLegend();
  //TOP + BOT:
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
    if (i==0 || i==8) grSig[i]->GetYaxis()->SetRangeUser(0.6,1.2);
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
  grSig[16]->GetYaxis()->SetRangeUser(0.6,1.2);
  grSig[16]->Draw("alp");
  
  txtout.close();
  
  //the end  
}



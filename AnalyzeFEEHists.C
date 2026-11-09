///
//
//  NOTE: TO suppress graphics to the screen, run in batch mode:
// gROOT->SetBatch(kTRUE);
//  OR
// root -b ....
//
//
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


//void AnalyzeFEEHists(  TH1D** hFeeSig, int runnum)
void AnalyzeFEEHists(int runnum)
{
  //get all FEE hists for this run:
  
  //set up input root file for FEEhists:
  string infile = "./FEESigs/FEESigs."+std::to_string(runnum)+".root";
  TFile* ipf = TFile::Open(infile.c_str());
  if(!ipf || !ipf->IsOpen() || ipf->IsZombie()){
    cout << "Could not open file " << infile << endl;
    exit(0);
  }
  
  TString fout = "./FEESigs/FEESigs."+std::to_string(runnum)+".pdf";
//  TString fout = infile.c_str();
//  fout.ReplaceAll(".root",".pdf");

  TString ftxt = fout;
  ftxt.ReplaceAll(".pdf",".txt");

  
  //read in FEE peak histos for this run
  TH1D** hFeeSig;
  hFeeSig      = new TH1D*[72];  

  for (int i=0;i<72;i++)
  {
//    TString hname   = TString::Format("%s%d","hFeeSig",i);
//    TString hname   = TString::Format("%s%d","hFeeSigHiR",i);
    TString hname   = TString::Format("%s%d","hFeeSigRaw",i);
    hFeeSig[i] = (TH1D*)ipf->Get(hname.Data());
  }

  
  //////////////////////
  // Begin analysis here
  //////////////////////

  double xmin = 0.1, xmax = 1.5;
  double parInit[3] = {100.,1.0,0.3};
  double dsmin = 0.3, dsmax=0.3;
  TString hname = hFeeSig[0]->GetName();
  if (hname.Contains("Raw") ) //if analyzing hFeeSigRaw
  {
    xmin = 10, xmax = 800;
    parInit[1]=200; parInit[2]=100;
    dsmin = 100; dsmax=300;
  }
  int method = 2;

  //now get pk for each slice
  int knt = 0;
  TF1* fg = new TF1("fg","gaus",xmin,xmax);
  fg->SetParameters(parInit);        //initial guess
  fg->SetParLimits(2,xmin,xmax);
  
  TF1* flan = new TF1("flan","landau",xmin,xmax);
  flan->SetParameters(parInit);        //initial guess
  flan->SetParLimits(2,xmin,xmax);
  
//  TString opt = "QOB";      //fitting options. B does not seem to work.
  TString opt = "QOB WW";      //fitting options
  
  ofstream txtout = ofstream(ftxt.Data(),std::ios::out);

  for (int i=0;i<72;i++) //FEEs
  { 
    if (hname.Contains("Raw") ) hFeeSig[i]->Rebin(5); //if analyzing hFeeSigRaw
    hFeeSig[i]->GetXaxis()->SetRangeUser(xmin,xmax);  //limit following search for maxval to this range
    double maxval = hFeeSig[i]->GetBinCenter(hFeeSig[i]->GetMaximumBin());
    double ampl = hFeeSig[i]->GetBinContent(hFeeSig[i]->GetMaximumBin());
    fg->SetParameters(ampl, maxval,parInit[2]);
    flan->SetParameters(ampl, maxval,parInit[2]);
    if (method <= 1) hFeeSig[i]->Fit(fg,opt.Data(),"",maxval-dsmin,maxval+dsmax);
    if (method == 2) hFeeSig[i]->Fit(flan,opt.Data(),"",maxval-dsmin,maxval+dsmax);
    TF1* fit = hFeeSig[i]->GetFunction("fg");
    if (method == 2) fit = hFeeSig[i]->GetFunction("flan");
    double* pars = fit->GetParameters();
    // sanity check parameters:
    double chi2 = fit->GetChisquare();  //chi2<100 seems ok for hFeeSigRaw in 15750 but it is not consistently an indicator of bad fits.
    if ( pars[1] < xmin || pars[1] > xmax )// || pars[2]>10.*parInit[2] ) 
//    if ( chi2 > 100 ) 
    {
      //reset fit parameters, try again
      pars[1] = maxval; pars[2]=1.25*parInit[2]; //definition of insanity: doing the same thing again and again and expecting a different result
      cout << "Bad first fit, FEE " << i << endl;
    }
    //fit a second time:
//    pars[2] *= 2.5; //limit to 2.5*width
    //replace xmin,xmax:
    double xmin_2 = max(pars[1]-4.*pars[2], 10.);
    double xmax_2 = min(pars[1]+5.*pars[2], 1000.);
    hFeeSig[i]->GetXaxis()->SetRangeUser(xmin_2,xmax_2);  //limit following search for maxval to this range
    fg->SetParameters(pars); //one of these function's pars will be wrong
    flan->SetParameters(pars);
//    if (method <= 1) hFeeSig[i]->Fit(fg,opt.Data(),"",maxval-dsmin,maxval+dsmax);
//    if (method == 2) hFeeSig[i]->Fit(flan,opt.Data(),"",maxval-dsmin,maxval+dsmax);
    if (method <= 1) hFeeSig[i]->Fit(fg,  opt.Data(),"",xmin,xmax);
    if (method == 2) hFeeSig[i]->Fit(flan,opt.Data(),"",xmin,xmax);
    TF1* fit_2 = hFeeSig[i]->GetFunction("fg");
    if (method == 2) fit_2 = hFeeSig[i]->GetFunction("flan");
    double* pars_2 = fit_2->GetParameters();
        
    double pk = (fit_2)? fit_2->GetParameter(1) : 0;
    double pkerr = (fit_2)? fit_2->GetParError(1) : 0;
    double chi2_2 = fit_2->GetChisquare();
    cout << i << "\t" << pk << "\t" << pkerr << "\t" << chi2 << "\t" << chi2_2 << endl;
    txtout << i << "\t" << pk << "\t" << pkerr << "\t" << chi2_2 << endl;
    
    hFeeSig[i]->GetXaxis()->SetRangeUser(xmin,xmax);      //reset drawing range
    
  }
  
  txtout.close();
    
  //Draw and save to file:
  TCanvas* cv[8];
  int nx = 3, ny=3;
  int ipad=1;
  int cknt = -1;
  TString cname;
    
  gStyle->SetOptStat(11);
  //    gStyle->SetOptLogy();
  gStyle->SetOptFit(1111);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
    
  for (int i=0; i<72; i++)
  {
    if (i%(nx*ny) == 0)
    {
      cknt++;
//      cout << "New TCanvas, i, cknt= " << i << "\t" << cknt << endl;
      cname = TString::Format("%s%d","cv",cknt);
      cv[cknt] = new TCanvas(cname.Data(),cname.Data(),1);
      cv[cknt]->Divide(nx,ny); //3*3=9
      ipad = 1;
    }
    cv[cknt]->cd(ipad);
    hFeeSig[i]->Draw();
    ipad++;
  }
  cknt++;  //last page

  //////////////
  // Write all canvases to pdf

  cout << "Writing " << cknt << " TCanvases to " << fout << endl;

  for (int i=0; i<cknt; i++)
  {
    if (i<cknt-1) 
      cv[i]->Print(fout+"(");
    else
      cv[i]->Print(fout+")");
  }
  
  for (int i=0; i<cknt; i++)
  {
    delete cv[i];
  }
  
  ipf->Close();
  
}

////////////////////////
/////
void RunAnalyzeFEEHists(string runFileListName = "AllGoodRuns.txt")
{

  vector<string> runList = getFileListVector(runFileListName );
  
  for (int i=0; i<runList.size(); i++)
  {
//    if (std::stoi(runList[i]) >= 15630) AnalyzeFEEHists(std::stoi(runList[i]));
    AnalyzeFEEHists(std::stoi(runList[i]));
  }


}


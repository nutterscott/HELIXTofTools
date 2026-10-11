TGraph* SliceAndDice(TH2 *h2, int meth = 0, double xMin=-800, double xMax=800, int bincomb=1)
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
  double x[nbins], maxval[nbins], pk[nbins];
  TString hname = h2->GetName();
  for (int i=0;i<nbins;i++) 
  {
    hsli[i] = h2->ProjectionY(Form("%s_%d",hname.Data(),i+1),i*bincomb+1,(i+1)*bincomb);
    x[i] = (h2->GetXaxis()->GetBinCenter(i*bincomb+1)+h2->GetXaxis()->GetBinCenter((i+1)*bincomb))/2.;
    hsli[i]->SetTitle(Form("%s%3.1f","X: ",x[i]));
  }
  
  //now get pk for each slice
  int knt = 0;
  for (int i=0;i<nbins;i++) 
  { 
    double bcenter = (h2->GetXaxis()->GetBinCenter(i*bincomb+1)+h2->GetXaxis()->GetBinCenter((i+1)*bincomb))/2.;
    if ( bcenter /*x[i] above*/ < xMin ) continue;
    x[knt] = bcenter;
    maxval[knt] = hsli[i]->GetBinCenter(hsli[i]->GetMaximumBin());
    hsli[i]->Fit("gaus","QO","",0.,3.);
    TF1* fit = hsli[i]->GetFunction("gaus");
    pk[knt] = (fit)? fit->GetParameter(1) : 0;
    if ( x[knt] > xMax ) break;
    cout << x[knt] << "\t" << pk[knt] << endl;
    knt++;
  }
  
  TGraph *g1;
  if (meth == 0) 
    g1= new TGraph(knt, x, maxval);
  else 
    g1 = new TGraph(knt, x, pk);
  g1->SetMarkerStyle(kFullCircle);
  g1->SetMarkerColor(kGreen);
  TString gtitle = TString::Format("%s%s",hname.Data()," signal;DCTY (mm);Signal peak (arb)");
  g1->SetTitle(gtitle);
  g1->SetName(Form("%s%s",hname.Data(),"_g"));
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
  return g1;
}
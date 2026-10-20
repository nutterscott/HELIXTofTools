//correct signal for response across length of paddle
// using pol2 fits to sqrt(signal) vs DCTY for each paddle

double p0[17], p1[17], p2[17];

void readResponseCorr(string infilename)
{
    ifstream input;
    input.open(infilename.c_str());
    if (!input.is_open())
    {
      cout << "File cannot be found: " << infilename << endl;
      exit(-1);
      return;
    }
    
    int padnum;
  
    for (int i=0; i<17; i++)
    {
      input >> padnum >> p0[i] >> p1[i] >> p2[i]; 
      cout << padnum << "\t" << p0[i] << "\t" << p1[i] << "\t" << p2[i] << endl; 
    }  
    input.close();
  
}

double responseCorr(int padID, double DCTY)
{
  double corr = 1; //default correction
  
  //position dependent correction to response sqrt(signal) for paddles
  // divide sqrt(signal) by this number
  
  //sanity checks on input:
  if (padID > 16 || padID <0) 
  {
    cout << "ERROR responseCorr: Bad paddle ID " << padID << endl;
    return corr;
  }
  
  double ymin = -800, ymax = 800; //mm
  if (padID == 16) //bore paddle
  {
    ymin = -300; ymax = 300;
  }

  double ypos = DCTY; //copy
  
  //sanity checks on input:
  //default to ymin or ymax if out of range
  if (ypos > ymax ) ypos = ymax;
  if (ypos < ymin ) ypos = ymin;

  corr = p0[padID] + p1[padID]*ypos + p2[padID]*ypos*ypos;
  
  return corr;
}
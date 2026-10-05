/// call this from the command line (while sitting in the build directory):
/// root -l -b -q '../Analysis/Runlists/read_list_example.C("../Analysis/Runlists/2024_Kiruna_Float_Bon.toml")'
void WriteRunNumbers(string runListFilename, string dataDir="./data")
{
  HLX_TOML_Table runListTable(runListFilename);
  vector<int> runNumbers;
  runListTable.getValueAtKeyPath("Runlist", runNumbers);

  cout << "Total number of runs: " << runNumbers.size() << endl;

  //save to a file:
  TString fname = runListFilename;
  fname.ReplaceAll(".toml",".txt");
  ofstream fout(fname.Data());
  for(auto run : runNumbers){
    cout << run << endl;
    fout << run << endl;
  }
  
  fout.close();
  
  

  auto subsetKeys = runListTable.getTableSubKeys("Subset");
  sort(subsetKeys.begin(), subsetKeys.end());
  cout << "Run Subsets: " << subsetKeys.size() << endl;
  if(subsetKeys.size() > 0){
    cout << "Runs per subset: ";
    for(auto key : subsetKeys){
      vector<int> subsetNumbers;
      auto thisKey = format("Subset.{}", key);
      runListTable.getValueAtKeyPath(thisKey, subsetNumbers);
      cout << key<<":"<<subsetNumbers.size() << " ";
    }
    cout << endl;
  }

  exit(0);
}



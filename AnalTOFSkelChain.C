//
// SC 10/12/2026   A process to run multiple skeletor files chained together using a variant of AnalTOFSkel.C
// renamed AnalTOFSkelTree.C .
//
// ROOT usage:
//   .L AnalTOFSkelTree.C
//   .L AnalTOFSkelChain.C
//   AnalTOFSkelList("filelist.txt", "all.AnalTOFSkel.root");
//
// where filelist.txt is the list of run numbers, such as 2024_Kiruna_Float_Boff.txt, to chain and process, and all.AnalTOFSkel.root will be
// the output file.

//make the chain global so you can access it from the command line:
TChain chain("toft");

void AnalTOFSkelList(const char *listName = "filelist.txt",
                    const char *outName  = "combined.AnalTOFSkel.root")
{
  std::string dir = "/data/skeletor/";
  std::string fnameStem = "skeletor.";
  std::string fnameExt = ".root";
  
  std::ifstream list(listName);
  if (!list) {
    std::cerr << "Cannot open file list: " << listName << "\n";
    return;
  }

//  TChain chain("toft");
  std::string line;
  int nFiles = 0;

  while (std::getline(list, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::string fname = dir+fnameStem+line.c_str()+fnameExt;
//    if (chain.Add(line.c_str()) <= 0) {
//      std::cerr << "Could not add: " << line << "\n";
    if (chain.Add(fname.c_str()) <= 0) {
      std::cerr << "Could not add: " << fname << "\n";
      continue;
    }
    ++nFiles;
  }

  if (nFiles == 0 || chain.GetEntries() == 0) {
    std::cerr << "No usable files or entries found\n";
    return;
  }

  std::cout << "Added " << nFiles << " files; entries = "
            << chain.GetEntries() << "\n";

  AnalTOFSkelTree(&chain, outName);
}

// Function to create a TChain from a list of run numbers in a text file
TChain* make_chain(const string datadir, const string runNumList, const string treeName = "TreeName") {

    // Create the TChain object
    TChain* chain = new TChain(treeName.c_str());
    
    // Open the text file containing the list of ROOT files
    std::ifstream infile(runNumList.c_str());
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file list " << runNumList << std::endl;
        delete chain;
        return nullptr;
    }
    
    std::string line;
    int filesAdded = 0;
    
    // Read the file line by line
    while (std::getline(infile, line)) {
        // Strip leading/trailing whitespaces or handle empty lines if necessary
        if (line.empty() || line[0] == '#') continue; // Skip empty lines and comments
        
        //construct filename:
        string fname = datadir+"/"+"run_"+line+".stg1.root";
        cout << fname << endl;
        // Add the file to the TChain
        if (chain->Add(fname.c_str())) {
            filesAdded++;
        } else {
            std::cerr << "Warning: Could not add file " << fname << " to the chain." << std::endl;
        }
    }
    
    infile.close();
    
    std::cout << "Successfully added " << filesAdded << " files to the TChain." << std::endl;
    return chain;
}

void Run(const char* file_path, const char* tree_name) {

    // Open the ROOT file
    TFile* root_file = TFile::Open(file_path);

    // Check if the file is open successfully
    if (!root_file || root_file->IsZombie()) {
        std::cerr << "Failed to open file: " << file_path << std::endl;
        return;
    }

    // Access the TTree
    TTree* tree = dynamic_cast<TTree*>(root_file->Get(tree_name));

    // Check if the tree is found
    if (!tree) {
        std::cerr << "Failed to find tree: " << tree_name << " in file: " << file_path << std::endl;
        return;
    }

    // Print the structure of the TTree
    std::cout << "Structure of TTree '" << tree_name << "':" << std::endl;

    // Loop over branches and print their names and types
    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch) {
            const char* branch_name = branch->GetName();
            const char* branch_type = branch->GetTitle();
            std::cout << "  Branch: " << branch_name << ", Type: " << branch_type << std::endl;
        }
    }

    // Close the ROOT file
    root_file->Close();

}

void PrintTTreeStructure() { 

    Run(
        "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
        ,"TrkAnaExt/trkana" );


    return;
    
}
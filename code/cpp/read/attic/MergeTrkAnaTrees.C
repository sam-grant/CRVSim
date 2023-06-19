void MergeTrkAnaTrees() { 

	// Open file
	TString finName = "../../data/nts.owner.trkana-reco.version.sequencer.2.root";
	TFile *fin = TFile::Open(finName);

	cout<<"---> Opened file "<<finName<<", "<<fin<<endl;

	// Get input trees
	TString tr1Name = "TrkAnaExt/trkana;1";
	TString tr2Name = "TrkAnaExt/trkana;2";

	TTree *tr1 = (TTree*)fin->Get(tr1Name);
	TTree *tr2 = (TTree*)fin->Get(tr2Name);

    //TTree* tr1 = dynamic_cast<TTree*>(fin->Get("TrkAnaExt/trkana;1"));
    //TTree* tr2 = dynamic_cast<TTree*>(fin->Get("TrkAnaExt/trkana;2"));

	cout<<"---> Got trees "<<tr1Name<<" "<<tr1<<" and "<<tr2Name<<" "<<tr2<<endl;

	// Write them both to a second output file
	TString foutName = "../../data/nts.owner.trkana-reco.version.sequencer.merged.root";
	TFile *fout = new TFile(foutName, "RECREATE");

    // Create a new output tree
    TTree* trOut = tr1->CloneTree(0);  // Clone the structure of tr1 without copying the entries
    cout<<"---> Copying entries"<<endl;

    // Merge the entries from tr1 into the merged tree
    trOut->CopyEntries(tr1);

    // Merge the entries from tr2 into the merged tree
    //trOut->CopyEntries(tr2);

    // Write the merged tree to the output file
    //fout->mkdir("TrkAnaExt"); fout->cd("TrkAnaExt");
    trOut->Write();

    // Close the input and output files
    fin->Close();
    fout->Close();

    cout<<"---> Done: written output to "<<foutName<<" "<<fout<<endl;

    // Clean up
    delete fin;
    delete fout;

    return;

}
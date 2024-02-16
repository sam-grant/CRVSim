


void CheckEvent() {

  TString finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"; 
  TFile *fin = TFile::Open(finName); 

  cout << "---> Read " << finName << ", " << fin << endl;

  TTree *tree = (TTree*)fin->Get("TrkAnaExt/trkana"); 

  cout << "---> Get tree " << tree << endl;

  Int_t eventID; 
  
  tree->SetBranchAddress("evtinfo.eventid", &eventID); 

  // Loop over entries in the tree
  for (Long64_t entry = 0; entry < tree->GetEntries(); entry++) {

    tree->GetEntry(entry);

    if (eventID == 476774) { 

      cout << "Event ID: " << eventID << endl;

    }

  }

  fin->Close(); 

  return;

}  
// "evtinfo.eventid" 

//    "crvhit.sectorType" # CRV sector hit
//	    , "crvhit.pos.fCoordinates.fX" # Reconstructed position of the cluster in X
//	    , "crvhit.pos.fCoordinates.fY" # Reconstructed position of the cluster in Y
//	    , "crvhit.pos.fCoordinates.fZ" # Reconstructed position of the cluster in Z
//	    , "crvhit.timeStart" # Earliest time recorded at either end of all bars in the hit
//	    , "crvhit.timeEnd" # Latest time recorded at either end of all bars in the hit
// 	    , "crvhit.time" # average reconstructed hit time of the cluster.
// 	    , "crvhit.PEs" # total number of photoelectrons in this cluser
// 	    , "crvhit.nHits" # Number of individual bar hits combined in this hit
//	    , "crvhit.nLayers" # Number of CRV layers that are part of this cluster
//	    , "crvhit.angle" #

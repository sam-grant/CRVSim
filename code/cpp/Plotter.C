#include <iostream>
#include <vector>

#include "Plotter.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

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

// I can't get this to compile since I'm missing the dictionaries needed. 
// I have tried using rootcling and setting up the links in another header file. 
// I guess I just don't have access to them here. 
// Maybe I need to install some stuff? 

/*#ifdef __CLING__
#pragma link C++ defined_in "LinkDef.h";
#endif*/

// rootcling -f Dictionary.cxx -c Plotter.h LinkDef.h

void Run(TTree *tree, TFile *fout) { 

   	// Book histograms
   	cout<<"---> Booking histograms"<<endl;

   	TH1D *h1_trkTime = new TH1D("h1_trkTime", ";Track time [ns or kct?];Events", 101, 0, 101e3);
  	TH1D *h1_crvTime = new TH1D("h1_crvTime", ";CRV time [ns or kct?];Events", 101, 0, 101e3);

  	TH1F *h1_trkX = new TH1F("h1_trkX", ";Track X-position [mm];Events", 300, -150e3, 150e3);
  	TH1F *h1_crvX = new TH1F("h1_crvX", ";CRV X-position [mm];Events", 300, -150e3, 150e3);

  	TH1F *h1_trkZ = new TH1F("h1_trkZ", ";Track Z-position [mm];Events", 300, -150e3, 150e3); 
  	TH1F *h1_crvZ = new TH1F("h1_crvZ", ";CRV Z-position [mm];Events", 300, -150e3, 150e3);

  	TH1D *h1_deltaTime = new TH1D("h1_deltaTime", ";#Deltat (CRV time #minus track time) [ns or kct?];Events", 300, -50, 250);
  	TH1D *h1_deltaX = new TH1D("h1_deltaX", ";#DeltaX (CRV X-position #minus track X-position) [mm];Events", 300, -150e3, 150e3);
  	TH1D *h1_deltaZ = new TH1D("h1_deltaZ", ";#DeltaZ (CRV Z-position #minus track Z-position) [mm];Events", 300, -150e3, 150e3);

	// Get branches (using header file)
  	InitBranches br(tree);

	for(Long64_t entry = 0; entry < tree->GetEntries(); entry++) {
     
    	tree->GetEntry(entry);

    	// Get variables

       	Float_t trkTime = br.klmcsim_time;   
    	Float_t crvTime = br.crvinfomcplane_time;

       	Float_t trkX = br.klmcsim_pos_fCoordinates_fX;
    	Float_t crvX = br.crvinfomcplane_x;

    	Float_t trkZ = br.klmcsim_pos_fCoordinates_fZ;
    	Float_t crvZ = br.crvinfomcplane_z;

       	Float_t deltaTime = crvTime - trkTime;
    	Float_t deltaX = crvX - trkX;
    	Float_t deltaZ = crvZ - trkZ;

    	// Fill histograms 

    	h1_trkTime->Fill(trkTime);
    	h1_crvTime->Fill(crvTime);

     	h1_trkX->Fill(trkX);
    	h1_crvX->Fill(crvX);

    	h1_trkZ->Fill(trkZ);
    	h1_crvZ->Fill(crvZ);

    	h1_deltaTime->Fill(deltaTime);
    	h1_deltaX->Fill(deltaX);
    	h1_deltaZ->Fill(deltaZ);
  	}

  	fout->mkdir("trkCRVPlots"); fout->cd("trkCRVPlots");

  	h1_trkTime->Write();
  	h1_crvTime->Write();

  	h1_trkX->Write();
  	h1_crvX->Write();

  	h1_trkZ->Write();
  	h1_crvZ->Write();

  	h1_deltaTime->Write();
  	h1_deltaX->Write();
  	h1_deltaZ->Write();

	return;

}

void Plotter() {

		// Open file
	TString finName = "../../data/nts.owner.trkana-reco.version.sequencer.root";
	TFile *fin = TFile::Open(finName);

	// Get tree
	TString treeName = "TrkAnaExt/trkana";
   	TTree *tree = (TTree*)fin->Get(treeName);

   	// Book output
  	TString foutName = "../../plots/trkana-reco.test.root"; 
  	TFile *fout = new TFile(foutName,"RECREATE");

  	// Run
  	Run(tree, fout);

  	cout<<"\nDone. Histogram written to: "<<foutName<<" "<<fout<<endl;

  	// Close files
   	fin->Close();
   	fout->Close();

   	// Clean up pointers
   	delete fin;
   	delete fout;

	return 0;

}
#ifndef Plotter_h
#define Plotter_h

// ROOT includes
#include "TTree.h"

using namespace std;

class InitBranches { 

public: 

   // Declaration of leaf types

   // Tracker variables, I think
   Float_t klmcsim_time;
   Float_t klmcsim_pos_fCoordinates_fX;
   Float_t klmcsim_pos_fCoordinates_fY;
   Float_t klmcsim_pos_fCoordinates_fZ;

   // CRV variables 

   Float_t crvinfomc_time; 
   Float_t crvinfomc_x; 
   Float_t crvinfomc_y; 
   Float_t crvinfomc_z; 

   Float_t crvinfomcplane_time; 
   Float_t crvinfomcplane_x; 
   Float_t crvinfomcplane_y; 
   Float_t crvinfomcplane_z; 

   // Declare constructor
   InitBranches(TTree *tree);

};

// Constructor
InitBranches::InitBranches(TTree* tree) {


   tree->SetBranchAddress("klmcsim.time", &klmcsim_time);
   tree->SetBranchAddress("klmcsim.pos.fCoordinates.fX", &klmcsim_pos_fCoordinates_fX);
   tree->SetBranchAddress("klmcsim.pos.fCoordinates.fY", &klmcsim_pos_fCoordinates_fY);
   tree->SetBranchAddress("klmcsim.pos.fCoordinates.fZ", &klmcsim_pos_fCoordinates_fZ);


   tree->SetBranchAddress("crvinfomc._time", &crvinfomc_time); 
   tree->SetBranchAddress("crvinfomc._x", &crvinfomc_x); 
   tree->SetBranchAddress("crvinfomc._y", &crvinfomc_y); 
   tree->SetBranchAddress("crvinfomc._z", &crvinfomc_z); 

   tree->SetBranchAddress("crvinfomcplane._time", &crvinfomcplane_time); 
   tree->SetBranchAddress("crvinfomcplane._x", &crvinfomcplane_x); 
   tree->SetBranchAddress("crvinfomcplane._y", &crvinfomcplane_y); 
   tree->SetBranchAddress("crvinfomcplane._z", &crvinfomcplane_z); 

}

#endif

